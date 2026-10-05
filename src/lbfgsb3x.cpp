#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <vector>
#include <time.h>
#include <Rmath.h>
#include <Rcpp.h>
#include <R_ext/Linpack.h>
#include "lbfgsb_cpp.h"

using namespace Rcpp;

extern "C" void setulb_(int *n, int *m, double *x, double *l, double *u,
                        int *nbd, double *f, double *g, double *factr, double *pgtol,
                        double *wa, int *iwa, int *itask, int *iprint,
                        int *icsave, int *lsave, int *isave, double *dsave);

// lbfgsb_cpp's driver step for the Fortran setulb_ (which does its own
// printing through R, so `pr` is unused)
static void fortranStep(int n, int m, double *x, const double *l,
                        const double *u, const int *nbd, double &f,
                        double *g, double factr, double pgtol, double *wa,
                        int *iwa, int &itask, int iprint, int &icsave,
                        int *lsave, int *isave, double *dsave,
                        const lbfgsb3c_cpp::Printer *pr) {
  (void)pr;
  setulb_(&n, &m, x, const_cast<double*>(l), const_cast<double*>(u),
          const_cast<int*>(nbd), &f, g, &factr, &pgtol, wa, iwa, &itask,
          &iprint, &icsave, lsave, isave, dsave);
}

static void rPrinter(void *data, const char *s) {
  (void)data;
  Rprintf("%s", s);
}

extern "C" void lbfgsb3C_(int n, int lmm, double *x, double *lower,
                          double *upper, int *nbd, double *Fmin, optimfn fn,
                          optimgr gr, int *fail, void *ex, double factr,
                          double pgtol, int *fncount, int *grcount,
                          int maxit, char *msg, int trace, int iprint,
                          double atol, double rtol, double *g) {
  (void)msg;
  lbfgsb3c_cpp::Printer pr = {rPrinter, NULL};
  lbfgsb3c_cpp::lbfgsbDriver(fortranStep, n, lmm, x, lower, upper, nbd, Fmin,
                             fn, gr, fail, ex, factr, pgtol, fncount,
                             grcount, maxit, trace, iprint, atol, rtol, g,
                             &pr, NULL);
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

// R callbacks for one lbfgsb3cpp() call, passed to the solver through
// `ex` so that nested optimizations do not share them
struct RCallbacks {
  Function fn;
  Function gr;
  RObject pn;
  Environment rho;
};

double gfn(int n, double *x, void *ex){
  RCallbacks *cb = static_cast<RCallbacks*>(ex);
  Rcpp::NumericVector par(n);
  std::copy(&x[0], &x[0]+n, &par[0]);
  par.attr("names") = cb->pn;
  double ret = as<double>(cb->fn(par, cb->rho));
  return ret;
}

void ggr(int n, double *x, double *gr, void *ex){
  RCallbacks *cb = static_cast<RCallbacks*>(ex);
  Rcpp::NumericVector par(n), ret(n);
  std::copy(&x[0], &x[0]+n, &par[0]);
  par.attr("names") = cb->pn;
  ret = cb->gr(par, cb->rho);
  if (ret.size() != n) stop("gradient must have %d elements (got %d).", n, (int)ret.size());
  std::copy(&ret[0], &ret[0]+n, &gr[0]);
}

//[[Rcpp::export]]
Rcpp::List lbfgsb3cpp(NumericVector par, Function fn, Function gr, NumericVector lower, NumericVector upper, List ctrl, Environment rho){
  Rcpp::List ret;
  RCallbacks cb = {fn, gr, par.attr("names"), rho};
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
  if (rtolN.size() != 1) stop("reltol has to have one element in it.");
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
  // vectors so an R error (bad bounds, or in fn/gr) does not leak them
  std::vector<double> lowV(par.size());
  double *low = lowV.data();
  if (lower.size() == 1){
    std::fill_n(&low[0],par.size(),lower[0]);
  } else if (lower.size() == par.size()){
    std::copy(lower.begin(),lower.end(),&low[0]);
  } else {
    stop("Lower bound must match the size of par or only have one element.");
  }
  std::vector<double> upV(par.size());
  double *up = upV.data();
  if (upper.size() == 1){
    std::fill_n(&up[0],par.size(),upper[0]);
  } else if (upper.size() == par.size()){
    std::copy(upper.begin(),upper.end(),&up[0]);
  } else {
    stop("Upper bound must match the size of par or only have one element.");
  }
  std::vector<double> xV(par.size());
  double *x = xV.data();
  std::copy(par.begin(),par.end(),&x[0]);
  std::vector<int> nbdV(par.size());
  int *nbd = nbdV.data();
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
  void *ex = &cb;
  lbfgsb3c_cpp::InfoOut info = {};
  lbfgsb3c_cpp::Printer pr = {rPrinter, NULL};
  lbfgsb3c_cpp::lbfgsbDriver(engine == 1 ? lbfgsb3c_cpp::setulb : fortranStep,
                             n, lmm, x, low, up, nbd, &fmin, gfn, ggr, &fail,
                             ex, factr, pgtol, &fncount, &grcount, maxit,
                             trace, iprint, atol, rtol, &g[0], &pr, &info);
  NumericVector parf(par.size());
  std::copy(&x[0],&x[0]+par.size(),parf.begin());
  parf.attr("names")=cb.pn;
  g.attr("names")=cb.pn;
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
  return ret;
}
