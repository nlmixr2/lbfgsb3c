// Thread-safe C++ port of L-BFGS-B 3.0 (lbfgsb.f, blas.f, and R's
// LINPACK dpofa/dtrsl).  The Fortran is kept alongside; this port has
// no static or global mutable state and never calls the R API, so
// independent problems may be solved concurrently (e.g. OpenMP).
#ifndef LBFGSB3C_LBFGSB_CPP_H
#define LBFGSB3C_LBFGSB_CPP_H

#if defined(__cplusplus)
extern "C" {
#endif

typedef double optimfn(int n, double *par, void *ex);
typedef void optimgr(int n, double *par, double *gr, void *ex);

// Thread-safe entry point with the same signature as lbfgsb3C_.
// `trace`, `iprint` and `msg` are ignored: nothing is printed and msg is
// not written.  fn and gr must themselves be thread safe.  An invalid lmm
// (< 1, or too large for the workspace) returns *fail = 29 without
// evaluating fn.
void lbfgsb3Cts_(int n, int lmm, double *x, double *lower,
                 double *upper, int *nbd, double *Fmin, optimfn fn,
                 optimgr gr, int *fail, void *ex, double factr,
                 double pgtol, int *fncount, int *grcount,
                 int maxit, char *msg, int trace, int iprint,
                 double atol, double rtol, double *g);

#if defined(__cplusplus)
}

namespace lbfgsb3c_cpp {

// Optional sink for diagnostic output; print == nullptr means silent.
struct Printer {
  void (*print)(void *data, const char *s);
  void *data;
};

// Solver state after a run (mirrors the Fortran working arrays).
struct InfoOut {
  int itask;
  int icsave;
  int lsave[4];
  int isave[44];
  double dsave[29];
};

// itask code returned when lmm is invalid (see validLmm)
const int kInvalidLmm = 29;

// Name of an itask code (1..29), or "UNKNOWN".
const char *taskName(int itask);

// True when lmm >= 1 and the solver workspace for (n, lmm) fits in an int
// (the Fortran indexes it with default integers).
bool validLmm(int n, int lmm);

// Reverse-communication setulb (port of lbfgsb.f).
void setulb(int n, int m, double *x, const double *l, const double *u,
            const int *nbd, double &f, double *g, double factr,
            double pgtol, double *wa, int *iwa, int &itask, int iprint,
            int &icsave, int *lsave, int *isave, double *dsave,
            const Printer *pr);

// Driver with the lbfgsb3C_ loop; `pr` (may be null) receives trace and
// iprint output, `info` (may be null) receives the final solver state.
void lbfgsb3Cts_core(int n, int lmm, double *x, double *lower,
                     double *upper, int *nbd, double *Fmin, optimfn fn,
                     optimgr gr, int *fail, void *ex, double factr,
                     double pgtol, int *fncount, int *grcount,
                     int maxit, int trace, int iprint,
                     double atol, double rtol, double *g,
                     const Printer *pr, InfoOut *info);

} // namespace lbfgsb3c_cpp
#endif

#endif
