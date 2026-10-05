// Test helper: solve many independent problems with the thread-safe C++
// port inside an OpenMP loop and serially, so the results can be compared.
#include <Rcpp.h>
#include <vector>
#include "lbfgsb_cpp.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

namespace {

// Generalized Rosenbrock with a per-problem scale
struct RosenProb {
  double p;
};

double rosenFn(int n, double *x, void *ex) {
  double p = static_cast<RosenProb *>(ex)->p;
  double f = 1.0;
  for (int i = 1; i < n; ++i) {
    double a = x[i] - x[i - 1] * x[i - 1];
    double b = 1.0 - x[i];
    f += p * a * a + b * b;
  }
  return f;
}

void rosenGr(int n, double *x, double *g, void *ex) {
  double p = static_cast<RosenProb *>(ex)->p;
  for (int i = 0; i < n; ++i) g[i] = 0.0;
  for (int i = 1; i < n; ++i) {
    double a = x[i] - x[i - 1] * x[i - 1];
    g[i - 1] += -4.0 * p * a * x[i - 1];
    g[i] += 2.0 * p * a - 2.0 * (1.0 - x[i]);
  }
}

// Solve problem k; odd problems are bounded to exercise the bound logic.
// Output row: par (n), value, fncount, grcount, fail
void solveOne(int k, int n, double *out) {
  RosenProb prob = {10.0 + 5.0 * k};
  std::vector<double> x(n), lower(n), upper(n), g(n, 0.0);
  std::vector<int> nbd(n);
  for (int i = 0; i < n; ++i) {
    x[i] = (i % 2 == 0) ? -1.2 + 0.01 * k : 1.0 + 0.01 * k;
    if (k % 2 == 1) {
      lower[i] = -1.5;
      upper[i] = 0.9;
      nbd[i] = 2;
    } else {
      lower[i] = R_NegInf;
      upper[i] = R_PosInf;
      nbd[i] = 0;
    }
  }
  double fmin = 0.0;
  int fail = 0, fncount = 0, grcount = 0;
  lbfgsb3Cts_(n, 5, x.data(), lower.data(), upper.data(), nbd.data(), &fmin,
              rosenFn, rosenGr, &fail, &prob, 1e7, 0.0, &fncount, &grcount,
              1000, NULL, 0, -1, 0.0, 1e-8, g.data());
  for (int i = 0; i < n; ++i) out[i] = x[i];
  out[n] = fmin;
  out[n + 1] = fncount;
  out[n + 2] = grcount;
  out[n + 3] = fail;
}

} // namespace

//[[Rcpp::export(name=".lbfgsb3cThreadTest", rng=false)]]
List lbfgsb3cThreadTest_(int nprob, int nthreads) {
  if (nprob < 1 || nprob == NA_INTEGER) stop("nprob must be a positive integer");
  if (nthreads < 1 || nthreads == NA_INTEGER) stop("nthreads must be a positive integer");
  const int n = 10;
  const int ncol = n + 4;
  std::vector<double> par((size_t)nprob * ncol), ser((size_t)nprob * ncol);
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(dynamic)
#endif
  for (int k = 0; k < nprob; ++k) solveOne(k, n, &par[(size_t)k * ncol]);
  for (int k = 0; k < nprob; ++k) solveOne(k, n, &ser[(size_t)k * ncol]);
  NumericMatrix parM(ncol, nprob, par.begin()), serM(ncol, nprob, ser.begin());
  bool openmp = false;
#ifdef _OPENMP
  openmp = true;
#else
  (void)nthreads;
#endif
  return List::create(_["parallel"] = transpose(parM),
                      _["serial"] = transpose(serM),
                      _["openmp"] = openmp);
}

// Test helper: call lbfgsb3Cts_ directly with the given (n, lmm) pairs;
// n is at most 3
//[[Rcpp::export(name=".lbfgsb3cLmmTest", rng=false)]]
List lbfgsb3cLmmTest_(IntegerVector lmm, IntegerVector n) {
  int k = lmm.size();
  if (n.size() != k) stop("n and lmm must have the same length");
  IntegerVector fail(k), fncount(k);
  for (int j = 0; j < k; ++j) {
    RosenProb prob = {100.0};
    double x[3] = {-1.2, 1.0, -1.2}, l[3] = {0, 0, 0}, u[3] = {0, 0, 0},
           g[3] = {0, 0, 0}, f = 0.0;
    int nbd[3] = {0, 0, 0}, fl = 0, fc = -1, gc = -1;
    if (n[j] > 3) stop("n must be at most 3");
    lbfgsb3Cts_(n[j], lmm[j], x, l, u, nbd, &f, rosenFn, rosenGr, &fl, &prob,
                1e7, 0.0, &fc, &gc, 1000, NULL, 0, -1, 0.0, 1e-8, g);
    fail[j] = fl;
    fncount[j] = fc;
  }
  return List::create(_["fail"] = fail, _["fncount"] = fncount);
}
