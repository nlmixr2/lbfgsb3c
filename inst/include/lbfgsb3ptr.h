#ifndef __LBFGS3PTR_H__
#define __LBFGS3PTR_H__
#include <stdlib.h> // for NULL
#include <Rmath.h>
#include <R_ext/Rdynload.h>

#if defined(__cplusplus)
extern "C" {
#endif

  typedef double optimfn(int n, double *par, void *ex);

  typedef void optimgr(int n, double *par, double *gr, void *ex);

  typedef void (*lbfgsb3_fn)(int n, int lmm, double *x, double *lower,
                             double *upper, int *nbd, double *Fmin, optimfn fn,
                             optimgr gr, int *fail, void *ex, double factr,
                             double pgtol, int *fncount, int *grcount,
                             int maxit, char *msg, int trace, int nREPORT, double atol,
                             double rtol, double *g);
  // Fortran L-BFGS-B driver.  Not thread safe.
  extern lbfgsb3_fn lbfgsb3C;

  // Thread-safe C++ port with the same signature as lbfgsb3C.  It keeps
  // no shared state and never calls the R API, so it may be
  // called concurrently (for example inside an OpenMP loop) as long as
  // `fn` and `gr` are themselves thread safe.  `trace`, `nREPORT` and `msg`
  // are ignored; nothing is printed.  Its workspace is allocated once per
  // calling thread and reused, so memory grows with the number of threads,
  // not the number of optimizations.  An invalid `lmm` (< 1, or so large
  // the workspace size overflows an int) returns *fail = 29
  // ("ERROR: INVALID LMM") without calling fn; lbfgsb3C does the same.
  // It stays NULL when the installed
  // lbfgsb3c is too old to provide it.
  extern lbfgsb3_fn lbfgsb3Cts;

  static inline SEXP iniLbfgsb3ptr0(SEXP p) {
    if (lbfgsb3C == NULL) {
      lbfgsb3C = (lbfgsb3_fn) R_ExternalPtrAddrFn(VECTOR_ELT(p, 0));
    }
    if (lbfgsb3Cts == NULL && Rf_length(p) > 1) {
      lbfgsb3Cts = (lbfgsb3_fn) R_ExternalPtrAddrFn(VECTOR_ELT(p, 1));
    }
    return R_NilValue;
  }

#define iniLbfgsb3                             \
  lbfgsb3_fn lbfgsb3C = NULL;                   \
  lbfgsb3_fn lbfgsb3Cts = NULL;                 \
  SEXP iniLbfgsb3ptr(SEXP ptr) {                \
    return iniLbfgsb3ptr0(ptr);                 \
  }

#if defined(__cplusplus)
}
#endif

#endif
