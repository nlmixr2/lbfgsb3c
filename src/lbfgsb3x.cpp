#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <time.h>
#include <Rmath.h>
#include <Rcpp.h>
#include <R_ext/Linpack.h>
#include "lbfgsb_cpp.h"
#define max2( a , b )  ( (a) > (b) ? (a) : (b) )

using namespace Rcpp;

extern "C" void setulb_(int *n, int *m, double *x, double *l, double *u,
                        int *nbd, double *f, double *g, double *factr, double *pgtol,
                        double *wa, int *iwa, int *itask, int *iprint,
                        int *icsave, int *lsave, int *isave, double *dsave);

// Fortran driver; the final solver state goes to `info` (may be NULL)
static void lbfgsb3C_fortran(int n, int lmm, double *x, double *lower,
                          double *upper, int *nbd, double *Fmin, optimfn fn,
                          optimgr gr, int *fail, void *ex, double factr,
                          double pgtol, int *fncount, int *grcount,
                          int maxit, char *msg, int trace, int iprint,
                          double atol, double rtol, double *g,
                          lbfgsb3c_cpp::InfoOut *info) {
  // Optim compatible interface
  fncount[0]=0;
  grcount[0]=0;
  if (!lbfgsb3c_cpp::validLmm(n, lmm)) {
    // the Fortran divides by zero for lmm <= 0
    if (info != NULL) {
      std::memset(info, 0, sizeof(*info));
      info->itask = lbfgsb3c_cpp::kInvalidLmm;
    }
    fail[0] = lbfgsb3c_cpp::kInvalidLmm;
    return;
  }
  int itask= 2;
  // *Fmin=;
  double *lastx = new double[n];
  std::copy(&x[0],&x[0]+n,&lastx[0]);
  int nwa = 2*lmm*n + 11*lmm*lmm + 5*n + 8*lmm;
  double *wa= new double[nwa];
  int niwa = 3*n;
  int *iwa= new int[niwa];
  int icsave = 0;
  int lsave[4] = {0};
  int isave[44] = {0};
  int i=0;
  double dsave[29]= {0};
  // Initial setup
  int doExit=0;
  fncount[0]=0;
  grcount[0]=0;
  int itask2=0;
  while (true){
    if (trace >= 2){
      Rprintf("itask: %d\n", itask);
      Rprintf("computing f and g at prm=\n");
      NumericVector xv(n);
      std::copy(&x[0],&x[0]+n,&xv[0]);
      print(xv);
      // // Calculate f and g
      // Fmin[0] = fn(n, x, ex);
      // fncount[0]++;
      // gr(n, x, g, ex);
      // grcount[0]++;
      Rprintf("\n================================================================================\nBefore call task number %d, or \"%s\"\n", itask, lbfgsb3c_cpp::taskName(itask));
    }
    if (itask==3) doExit=1;
    setulb_(&n, &lmm, x, lower, upper, nbd, Fmin, g, &factr, &pgtol,
            wa, iwa, &itask, &iprint, &icsave, lsave, isave, dsave);

    if (trace > 2) {
      Rprintf("returned from lbfgsb3 \n");
      Rprintf("returned itask is %d or \"%s\"\n",itask,lbfgsb3c_cpp::taskName(itask));
    }
    switch (itask){
    case 4:
    case 20:
    case 21:
      if (trace >= 2) {
        Rprintf("computing f and g at prm=\n");
        NumericVector xv(n);
        std::copy(&x[0],&x[0]+n,&xv[0]);
        print(xv);
      }
      // Calculate f and g
      Fmin[0] = fn(n, x, ex);
      fncount[0]++;
      gr(n, x, g, ex);
      grcount[0]++;
      if (trace > 0) {
        Rprintf("At iteration %d f=%f ", isave[33], *Fmin);
        if (trace > 1) {
          double tmp = 0;
          for (int j = 0; j < n; j++) {
            if (tmp < fabs(g[j])){
              tmp = fabs(g[j]);
            }
          }
          Rprintf("max(abs(g))=%f",tmp);
        }
        Rprintf("\n");
      }
      break;
    case 1:
      // New x;
      if (maxit < fncount[0]){
        itask2=28;
        doExit=1;
        itask=3; // Stop -- gives the right results and restores gradients
        if (trace > 2){
          Rprintf("Exit becuase maximum number of function calls %d met.\n", maxit);
        }
      } else {
        bool converge=fabs(lastx[n-1]-x[n-1]) < fabs(x[n-1])*rtol+atol;
        if (converge){
          for (i=n-1;i--;){
            converge=fabs(lastx[i]-x[i]) < fabs(x[i])*rtol+atol;
            if  (!converge){
              break;
            }
          }
        }
        if (converge){
          itask2=27;
          itask=3; // Stop -- gives the right results and restores gradients
          if (trace > 2){
            Rprintf("CONVERGENCE: Parameters differences below xtol.\n");
          }
          doExit=1;
        }
      }
      std::copy(&x[0],&x[0]+n,&lastx[0]);
      break;
    default:
      doExit=1;
    }
    if (doExit) break;
  }
  if (itask2){
    itask=itask2;
  }
  if (info != NULL) {
    info->itask = itask;
    info->icsave = icsave;
    std::copy(&lsave[0], &lsave[0]+4, info->lsave);
    std::copy(&isave[0], &isave[0]+44, info->isave);
    std::copy(&dsave[0], &dsave[0]+29, info->dsave);
  }
  fail[0]= itask;
  delete[] wa;
  delete[] iwa;
  delete[] lastx;
}

extern "C" void lbfgsb3C_(int n, int lmm, double *x, double *lower,
                          double *upper, int *nbd, double *Fmin, optimfn fn,
                          optimgr gr, int *fail, void *ex, double factr,
                          double pgtol, int *fncount, int *grcount,
                          int maxit, char *msg, int trace, int iprint,
                          double atol, double rtol, double *g) {
  lbfgsb3C_fortran(n, lmm, x, lower, upper, nbd, Fmin, fn, gr, fail, ex,
                   factr, pgtol, fncount, grcount, maxit, msg, trace, iprint,
                   atol, rtol, g, NULL);
}

static void rPrinter(void *data, const char *s) {
  (void)data;
  Rprintf("%s", s);
}

static List infoList(const lbfgsb3c_cpp::InfoOut &info) {
  LogicalVector lsaveR(4);
  NumericVector dsaveR(29);
  IntegerVector isaveR(44);
  std::copy(&info.lsave[0], &info.lsave[0]+4, lsaveR.begin());
  std::copy(&info.dsave[0], &info.dsave[0]+29, dsaveR.begin());
  std::copy(&info.isave[0], &info.isave[0]+44, isaveR.begin());
  return List::create(_["task"] = CharacterVector::create(lbfgsb3c_cpp::taskName(info.itask)),
                      _["itask"]= IntegerVector::create(info.itask),
                      _["lsave"]= lsaveR,
                      _["icsave"]= IntegerVector::create(info.icsave),
                      _["dsave"]= dsaveR,
                      _["isave"] = isaveR);
}

Environment grho;

CharacterVector gnames;

List ev;

double gfn(int n, double *x, void *ex){
  Rcpp::NumericVector par(n);
  std::copy(&x[0], &x[0]+n, &par[0]);
  Function fn = as<Function>(ev["fn"]);
  par.attr("names") = ev["pn"];
  double ret = as<double>(fn(par, grho));
  return ret;
}

void ggr(int n, double *x, double *gr, void *ex){
  Rcpp::NumericVector par(n), ret(n);
  std::copy(&x[0], &x[0]+n, &par[0]);
  Function grad = as<Function>(ev["gr"]);
  par.attr("names") = ev["pn"];
  ret = grad(par, grho);
  std::copy(&ret[0], &ret[0]+n, &gr[0]);
}

//[[Rcpp::export]]
Rcpp::List lbfgsb3cpp(NumericVector par, Function fn, Function gr, NumericVector lower, NumericVector upper, List ctrl, Environment rho){
  Rcpp::List ret;
  ev["fn"] = fn;
  ev["gr"] = gr;
  ev["pn"] = par.attr("names");
  Rcpp::NumericVector g(par.size());
  // CONV in 6, 7, 8; ERROR in 9-19; WARN in 23-26
  IntegerVector traceI = as<IntegerVector>(ctrl["trace"]);
  if (traceI.size() != 1) stop("trace has to have one element in it.");
  int trace = traceI[0];
  NumericVector factrN = as<NumericVector>(ctrl["factr"]);
  if (factrN.size() != 1) stop("factr has to have one element in it.");
  double factr = factrN[0];
  NumericVector pgtolN = as<NumericVector>(ctrl["pgtol"]);
  if (pgtolN.size() != 1) stop("pgtol has to have one element in it.");
  double pgtol = pgtolN[0];
  NumericVector atolN = as<NumericVector>(ctrl["abstol"]);
  if (atolN.size() != 1) stop("abstol has to have one element in it.");
  double atol = atolN[0];
  NumericVector rtolN = as<NumericVector>(ctrl["reltol"]);
  if (atolN.size() != 1) stop("reltol has to have one element in it.");
  double rtol = rtolN[0];
  LogicalVector infoN = as<LogicalVector>(ctrl["info"]);
  if (infoN.size() != 1) stop("info has to have one element in it.");
  bool addInfo = infoN[0];
  IntegerVector lmmN = as<IntegerVector>(ctrl["lmm"]);
  if (lmmN.size() != 1) stop("lmm has to have one element in it.");
  int lmm = lmmN[0];//lmmN.size();
  int n = par.size();
  if (lmmN[0] == NA_INTEGER || !lbfgsb3c_cpp::validLmm(n, lmm))
    stop("lmm must be a whole number >= 1 small enough for the workspace to fit in memory (got %d).", lmm);
  IntegerVector maxitN = as<IntegerVector>(ctrl["maxit"]);
  if (maxitN.size() != 1) stop("maxit has to have one element in it.");
  int maxit = maxitN[0];
  IntegerVector iprintN = as<IntegerVector>(ctrl["iprint"]);
  if (iprintN.size() != 1) stop("iprint has to have one element in it.");
  int iprint = iprintN[0];
  // 0 = Fortran (default), 1 = thread-safe C++ port
  int engine = 0;
  if (ctrl.containsElementNamed("engine")) {
    IntegerVector engineN = as<IntegerVector>(ctrl["engine"]);
    if (engineN.size() != 1) stop("engine has to have one element in it.");
    engine = engineN[0];
  }
  // double *g = new double[par.size()];
  double *low = new double[par.size()];
  if (lower.size() == 1){
    std::fill_n(&low[0],par.size(),lower[0]);
  } else if (lower.size() == par.size()){
    std::copy(lower.begin(),lower.end(),&low[0]);
  } else {
    delete [] low;
    stop("Lower bound must match the size of par or only have one element.");
  }
  double *up = new double[par.size()];
  if (upper.size() == 1){
    std::fill_n(&up[0],par.size(),upper[0]);
  } else if (upper.size() == par.size()){
    std::copy(upper.begin(),upper.end(),&up[0]);
  } else {
    delete [] low;
    delete [] up;
    stop("Upper bound must match the size of par or only have one element.");
  }
  double *x = new double[par.size()];
  std::copy(par.begin(),par.end(),&x[0]);
  int *nbd = new int[par.size()];
  int i;
  for (i = par.size();i--;){
    /*
      nbd(i)=0 if x(i) is unbounded,
      1 if x(i) has only a lower bound,
      2 if x(i) has both lower and upper bounds,
      3 if x(i) has only an upper bound.
    */
    nbd[i] = 0;
    if (R_FINITE(low[i])) nbd[i] = 1;
    if (R_FINITE(up[i]))  nbd[i] = 3 - nbd[i];
  }
  double fmin=std::numeric_limits<double>::max();
  int fail = 0, fncount=0, grcount=0;
  grho=rho;
  void *ex =NULL;
  char msg[120];
  lbfgsb3c_cpp::InfoOut info = {};
  if (engine == 1) {
    lbfgsb3c_cpp::Printer pr = {rPrinter, NULL};
    lbfgsb3c_cpp::lbfgsb3Cts_core(n, lmm, x, low, up, nbd, &fmin, gfn, ggr,
                                  &fail, ex, factr, pgtol, &fncount,
                                  &grcount, maxit, trace, iprint, atol, rtol,
                                  &g[0], &pr, &info);
  } else {
    lbfgsb3C_fortran(n, lmm, x, low, up, nbd, &fmin, gfn, ggr,
                     &fail, ex, factr, pgtol, &fncount,
                     &grcount, maxit, msg, trace, iprint , atol, rtol, &g[0],
                     &info);
  }
  NumericVector parf(par.size());
  std::copy(&x[0],&x[0]+par.size(),parf.begin());
  parf.attr("names")=ev["pn"];
  g.attr("names")=ev["pn"];
  ret["par"]=parf;
  ret["grad"]=g;
  ret["value"] = fmin;
  IntegerVector cnt = IntegerVector::create(fncount,grcount);
  ret["counts"] = cnt;
  switch (fail){
  case 6:
  case 7:
  case 8:
  case 27:
    ret["convergence"]=IntegerVector::create(0);
    break;
  case 28:
    ret["convergence"]=IntegerVector::create(1);
    break;
  case 23:
  case 24:
  case 25:
  case 26:
    ret["convergence"] = IntegerVector::create(51);
    break;
  case 9:
  case 10:
  case 11:
  case 12:
  case 13:
  case 14:
  case 15:
  case 16:
  case 17:
  case 18:
  case 19:
  case 29:
    ret["convergence"] = IntegerVector::create(52);
    break;
  default:
    ret["convergence"] = IntegerVector::create(NA_INTEGER);
  }
  ret["message"]= CharacterVector::create(lbfgsb3c_cpp::taskName(fail));
  if (addInfo) ret["info"] = infoList(info);
  delete [] x;
  delete [] low;
  delete [] up;
  delete [] nbd;
  return ret;
}
