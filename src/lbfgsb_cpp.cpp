// Thread-safe C++ port of L-BFGS-B 3.0 (Nocedal & Morales, 2011) as
// modified for lbfgsb3c (src/lbfgsb.f, src/blas.f) plus R's LINPACK
// dpofa/dtrsl.  The translation is routine-for-routine and keeps the
// Fortran arithmetic order; each routine notes its Fortran source line.
//
// Thread safety: no shared mutable state, no R API.  All solver state
// lives in the wa/iwa/isave/dsave/lsave arrays exactly as in the Fortran
// reverse-communication design; wa/iwa come from a per-thread workspace
// that is allocated once (growing as needed) and zeroed for each run.
// Diagnostic output goes to an optional Printer callback.
//
// L-BFGS-B is released under the "New BSD License" (see
// inst/License-lbfgsb-orig.txt).
#include "lbfgsb_cpp.h"
#include <algorithm>
#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <vector>

namespace lbfgsb3c_cpp {

namespace {

// 1-based vector view, v(i) == p[i-1]
template <class T> struct V1 {
  T *p;
  explicit V1(T *p_) : p(p_) {}
  T &operator()(int i) const { return p[i - 1]; }
};

// 1-based column-major matrix view, a(i,j) == p[(j-1)*ld + i-1]
template <class T> struct V2 {
  T *p;
  int ld;
  V2(T *p_, int ld_) : p(p_), ld(ld_) {}
  T &operator()(int i, int j) const {
    return p[(size_t)(j - 1) * (size_t)ld + (size_t)(i - 1)];
  }
};

inline double dmax(double a, double b) { return (a > b) ? a : b; }
inline double dmin(double a, double b) { return (a < b) ? a : b; }
inline double dmax3(double a, double b, double c) { return dmax(dmax(a, b), c); }

void say(const Printer *pr, const char *fmt, ...) {
  if (pr == nullptr || pr->print == nullptr) return;
  char buf[512];
  va_list ap;
  va_start(ap, fmt);
  vsnprintf(buf, sizeof(buf), fmt, ap);
  va_end(ap);
  pr->print(pr->data, buf);
}

void sayVec(const Printer *pr, const char *label, const double *x, int n) {
  if (pr == nullptr || pr->print == nullptr) return;
  say(pr, "%s", label);
  for (int i = 0; i < n; ++i) say(pr, " %g", x[i]);
  say(pr, "\n");
}

void sayRestart(const Printer *pr, const char *what) {
  say(pr, "%s\n   refresh the lbfgs memory and restart the iteration.\n", what);
}

// ---------------------------------------------------------------- BLAS
// (blas.f; only unit strides are used by L-BFGS-B)

// f: blas.f:166
double ddot(int n, const double *dx, const double *dy) {
  double dtemp = 0.0;
  for (int i = 0; i < n; ++i) dtemp = dtemp + dx[i] * dy[i];
  return dtemp;
}

// f: blas.f:64
void daxpy(int n, double da, const double *dx, double *dy) {
  if (n <= 0) return;
  if (da == 0.0) return;
  for (int i = 0; i < n; ++i) dy[i] = dy[i] + da * dx[i];
}

// f: blas.f:114 (forward element-by-element copy, as in the Fortran)
void dcopy(int n, const double *dx, double *dy) {
  for (int i = 0; i < n; ++i) dy[i] = dx[i];
}

// f: blas.f:217
void dscal(int n, double da, double *dx) {
  if (n <= 0) return;
  for (int i = 0; i < n; ++i) dx[i] = da * dx[i];
}

// ------------------------------------------------------------- LINPACK
// (R's src/appl/dpofa.f with its positive-definiteness tolerance, and
// src/appl/dtrsl.f)

void dpofa(double *a_, int lda, int n, int &info) {
  V2<double> a(a_, lda);
  const double eps = 1.0e-14;
  for (int j = 1; j <= n; ++j) {
    info = j;
    double s = 0.0;
    int jm1 = j - 1;
    for (int k = 1; k <= jm1; ++k) {
      double t = a(k, j) - ddot(k - 1, &a(1, k), &a(1, j));
      t = t / a(k, k);
      a(k, j) = t;
      s = s + t * t;
    }
    s = a(j, j) - s;
    if (s <= eps * fabs(a(j, j))) return;
    a(j, j) = sqrt(s);
  }
  info = 0;
}

void dtrsl(double *t_, int ldt, int n, double *b_, int job, int &info) {
  V2<double> t(t_, ldt);
  V1<double> b(b_);
  for (info = 1; info <= n; ++info) {
    if (t(info, info) == 0.0) return;
  }
  info = 0;
  int kase = 1;
  if (job % 10 != 0) kase = 2;
  if ((job % 100) / 10 != 0) kase = kase + 2;
  int j, jj;
  double temp;
  switch (kase) {
  case 1: // job = 00, T lower, T*x = b
    b(1) = b(1) / t(1, 1);
    for (j = 2; j <= n; ++j) {
      temp = -b(j - 1);
      daxpy(n - j + 1, temp, &t(j, j - 1), &b(j));
      b(j) = b(j) / t(j, j);
    }
    break;
  case 2: // job = 01, T upper, T*x = b
    b(n) = b(n) / t(n, n);
    for (jj = 2; jj <= n; ++jj) {
      j = n - jj + 1;
      temp = -b(j + 1);
      daxpy(j, temp, &t(1, j + 1), &b(1));
      b(j) = b(j) / t(j, j);
    }
    break;
  case 3: // job = 10, T lower, trans(T)*x = b
    b(n) = b(n) / t(n, n);
    for (jj = 2; jj <= n; ++jj) {
      j = n - jj + 1;
      b(j) = b(j) - ddot(jj - 1, &t(j + 1, j), &b(j + 1));
      b(j) = b(j) / t(j, j);
    }
    break;
  case 4: // job = 11, T upper, trans(T)*x = b
    b(1) = b(1) / t(1, 1);
    for (j = 2; j <= n; ++j) {
      b(j) = b(j) - ddot(j - 1, &t(1, j), &b(1));
      b(j) = b(j) / t(j, j);
    }
    break;
  }
}

// ------------------------------------------------------------- L-BFGS-B

// Project x onto [l, u] for one variable (the first loop of active).
// Returns 1 when x is at (or was moved to) a bound.
int activeProject(int nbd, double l, double u, double &x, int &prjctd) {
  if (nbd <= 0) return 0;
  if (nbd <= 2 && x <= l) {
    if (x < l) {
      prjctd = 1;
      x = l;
    }
    return 1;
  }
  if (nbd >= 2 && x >= u) {
    if (x > u) {
      prjctd = 1;
      x = u;
    }
    return 1;
  }
  return 0;
}

// iwhere for one variable (the second loop of active)
int activeWhere(int nbd, double l, double u) {
  if (nbd == 0) return -1;
  if (nbd == 2 && u - l <= 0.0) return 3;
  return 0;
}

// f: lbfgsb.f:1091
void active(int n, const double *l_, const double *u_, const int *nbd_,
            double *x_, int *iwhere_, int iprint, int &prjctd, int &cnstnd,
            int &boxed, const Printer *pr) {
  V1<const double> l(l_), u(u_);
  V1<const int> nbd(nbd_);
  V1<double> x(x_);
  V1<int> iwhere(iwhere_);
  int nbdd = 0;
  prjctd = 0;
  cnstnd = 0;
  boxed = 1;
  for (int i = 1; i <= n; ++i)
    nbdd = nbdd + activeProject(nbd(i), l(i), u(i), x(i), prjctd);
  for (int i = 1; i <= n; ++i) {
    if (nbd(i) != 2) boxed = 0;
    if (nbd(i) != 0) cnstnd = 1;
    iwhere(i) = activeWhere(nbd(i), l(i), u(i));
  }
  if (iprint >= 0) {
    if (prjctd == 1) say(pr, "initial X infeasible. Restart with projection.\n");
    if (cnstnd == 0) say(pr, "This problem is unconstrained.\n");
  }
  if (iprint > 0) say(pr, " Variables exactly at bounds for X0  %d\n", nbdd);
}

// f: lbfgsb.f:1199
void bmv(int m, double *sy_, double *wt_, int col, double *v_, double *p_,
         int &info) {
  V2<double> sy(sy_, m);
  V1<double> v(v_), p(p_);
  if (col == 0) return;
  p(col + 1) = v(col + 1);
  for (int i = 2; i <= col; ++i) {
    int i2 = col + i;
    double sum = 0.0;
    for (int k = 1; k <= i - 1; ++k) sum = sum + sy(i, k) * v(k) / sy(k, k);
    p(i2) = v(i2) + sum;
  }
  dtrsl(wt_, m, col, &p(col + 1), 11, info);
  if (info != 0) return;
  for (int i = 1; i <= col; ++i) p(i) = v(i) / sqrt(sy(i, i));
  dtrsl(wt_, m, col, &p(col + 1), 1, info);
  if (info != 0) return;
  for (int i = 1; i <= col; ++i) p(i) = -p(i) / sqrt(sy(i, i));
  for (int i = 1; i <= col; ++i) {
    double sum = 0.0;
    for (int k = i + 1; k <= col; ++k)
      sum = sum + sy(k, i) * p(col + k) / sy(i, i);
    p(i) = p(i) + sum;
  }
}

// f: lbfgsb.f:2506
void hpsolb(int n, double *t_, int *iorder_, int iheap) {
  V1<double> t(t_);
  V1<int> iorder(iorder_);
  int i, j, indxin, indxou;
  double ddum, out;
  if (iheap == 0) {
    for (int k = 2; k <= n; ++k) {
      ddum = t(k);
      indxin = iorder(k);
      i = k;
      while (i > 1) {
        j = i / 2;
        if (ddum < t(j)) {
          t(i) = t(j);
          iorder(i) = iorder(j);
          i = j;
        } else {
          break;
        }
      }
      t(i) = ddum;
      iorder(i) = indxin;
    }
  }
  if (n > 1) {
    i = 1;
    out = t(1);
    indxou = iorder(1);
    ddum = t(n);
    indxin = iorder(n);
    for (;;) {
      j = i + i;
      if (j <= n - 1) {
        if (t(j + 1) < t(j)) j = j + 1;
        if (t(j) < ddum) {
          t(i) = t(j);
          iorder(i) = iorder(j);
          i = j;
          continue;
        }
      }
      break;
    }
    t(i) = ddum;
    iorder(i) = indxin;
    t(n) = out;
    iorder(n) = indxou;
  }
}

// State of cauchy (lbfgsb.f:1315), whose steps are split into methods
// to keep each one small.  The operations and their order are unchanged.
struct Cauchy {
  int n, m, col, head, col2, iprint;
  double theta, epsmch;
  V1<double> x, g, t, d, xcp, p, wbp;
  V1<const double> l, u;
  V1<const int> nbd;
  V1<int> iorder, iwhere;
  V2<double> wy, ws;
  double *sy_, *wt_, *p_, *c_, *wbp_, *v_;
  const Printer *pr;
  bool bnded = true;
  int nfree = 0, nbreak = 0, ibkmin = 0;
  double bkmin = 0.0, f1 = 0.0, f2 = 0.0, f2_org = 0.0, dtm = 0.0,
         tsum = 0.0, tl = 0.0, tu = 0.0;

  Cauchy(int n_, double *x_, const double *l_, const double *u_,
         const int *nbd_, double *g_, int *iorder_, int *iwhere_,
         double *t_, double *d_, double *xcp_, int m_, double *wy_,
         double *ws_, double *sy, double *wt, double theta_, int col_,
         int head_, double *pw, double *cw, double *wbpw, double *vw,
         int iprint_, double epsmch_, const Printer *pr_)
      : n(n_), m(m_), col(col_), head(head_), col2(2 * col_),
        iprint(iprint_), theta(theta_), epsmch(epsmch_), x(x_), g(g_),
        t(t_), d(d_), xcp(xcp_), p(pw), wbp(wbpw), l(l_), u(u_), nbd(nbd_),
        iorder(iorder_), iwhere(iwhere_), wy(wy_, n_), ws(ws_, n_),
        sy_(sy), wt_(wt), p_(pw), c_(cw), wbp_(wbpw), v_(vw), pr(pr_) {}

  // Reset iwhere(i) from the bounds and the sign of -g(i)
  void classify(int i, double neggi) {
    if (iwhere(i) == 3 || iwhere(i) == -1) return;
    if (nbd(i) <= 2) tl = x(i) - l(i);
    if (nbd(i) >= 2) tu = u(i) - x(i);
    bool xlower = nbd(i) <= 2 && tl <= 0.0;
    bool xupper = nbd(i) >= 2 && tu <= 0.0;
    iwhere(i) = 0;
    if (xlower) {
      if (neggi <= 0.0) iwhere(i) = 1;
    } else if (xupper) {
      if (neggi >= 0.0) iwhere(i) = 2;
    } else {
      if (fabs(neggi) <= 0.0) iwhere(i) = -3;
    }
  }

  // Record breakpoint tbreak for variable i
  void addBreak(int i, double tbreak) {
    nbreak = nbreak + 1;
    iorder(nbreak) = i;
    t(nbreak) = tbreak;
    if (nbreak == 1 || t(nbreak) < bkmin) {
      bkmin = t(nbreak);
      ibkmin = nbreak;
    }
  }

  // Direction, p = W'd and the breakpoint of a free variable i
  void collect(int i, double neggi) {
    int pointr = head;
    d(i) = neggi;
    f1 = f1 - neggi * neggi;
    for (int j = 1; j <= col; ++j) {
      p(j) = p(j) + wy(i, pointr) * neggi;
      p(col + j) = p(col + j) + ws(i, pointr) * neggi;
      pointr = pointr % m + 1;
    }
    if (nbd(i) <= 2 && nbd(i) != 0 && neggi < 0.0) {
      addBreak(i, tl / (-neggi));
    } else if (nbd(i) >= 2 && neggi > 0.0) {
      addBreak(i, tu / neggi);
    } else {
      nfree = nfree - 1;
      iorder(nfree) = i;
      if (fabs(neggi) > 0.0) bnded = false;
    }
  }

  void scan() {
    for (int i = 1; i <= n; ++i) {
      double neggi = -g(i);
      classify(i, neggi);
      if (iwhere(i) != 0 && iwhere(i) != -1) {
        d(i) = 0.0;
      } else {
        collect(i, neggi);
      }
    }
  }

  // Scan the variables and set xcp = x.  False when there are no
  // breakpoints and no free variables (xcp = x is the answer).
  bool start() {
    nfree = n + 1;
    if (iprint >= 99) say(pr, "--- CAUCHY entered---\n");
    for (int i = 1; i <= col2; ++i) p(i) = 0.0;
    scan();
    if (theta != 1.0) dscal(col, theta, &p(col + 1));
    dcopy(n, &x(1), &xcp(1));
    if (nbreak == 0 && nfree == n + 1) {
      printXcp();
      return false;
    }
    return true;
  }

  // First and second derivative of the piecewise quadratic at t = 0 and
  // its stationary point dtm.  False on a bmv failure.
  bool initSlopes(int &info) {
    for (int j = 0; j < col2; ++j) c_[j] = 0.0;
    f2 = -theta * f1;
    f2_org = f2;
    if (col > 0) {
      bmv(m, sy_, wt_, col, p_, v_, info);
      if (info != 0) return false;
      f2 = f2 - ddot(col2, v_, p_);
    }
    dtm = -f1 / f2;
    tsum = 0.0;
    return true;
  }

  void printXcp() const {
    if (iprint > 100) sayVec(pr, "Cauchy X[1:5] = ", &xcp(1), n > 5 ? 5 : n);
  }

  // Next breakpoint (in increasing order) for pass `iter`; sets tj
  int nextBreak(int iter, int nleft, double &tj) {
    if (iter == 1) {
      tj = bkmin;
      return iorder(ibkmin);
    }
    if (iter == 2 && ibkmin != nbreak) {
      t(ibkmin) = t(nbreak);
      iorder(ibkmin) = iorder(nbreak);
    }
    hpsolb(nleft, &t(1), &iorder(1), iter - 2);
    tj = t(nleft);
    return iorder(nleft);
  }

  // Fix variable ibp at the bound it reaches; returns zibp
  double fixAtBound(int ibp, double dibp) {
    double zibp;
    if (dibp > 0.0) {
      zibp = u(ibp) - x(ibp);
      xcp(ibp) = u(ibp);
      iwhere(ibp) = 2;
    } else {
      zibp = l(ibp) - x(ibp);
      xcp(ibp) = l(ibp);
      iwhere(ibp) = 1;
    }
    return zibp;
  }

  // Limited-memory part of the f1, f2 update after fixing ibp
  void updateCurvature(int ibp, double dt, double dibp, double dibp2,
                       int &info) {
    daxpy(col2, dt, p_, c_);
    int pointr = head;
    for (int j = 1; j <= col; ++j) {
      wbp(j) = wy(ibp, pointr);
      wbp(col + j) = theta * ws(ibp, pointr);
      pointr = pointr % m + 1;
    }
    bmv(m, sy_, wt_, col, wbp_, v_, info);
    if (info != 0) return;
    double wmc = ddot(col2, c_, v_);
    double wmp = ddot(col2, p_, v_);
    double wmw = ddot(col2, wbp_, v_);
    daxpy(col2, -dibp, wbp_, p_);
    f1 = f1 + dibp * wmc;
    f2 = f2 + 2.0 * dibp * wmp - dibp2 * wmw;
  }

  void printPiece(int nseg, double dt) const {
    say(pr, "Piece %d\n", nseg);
    say(pr, "f1 at start point = %g\n", f1);
    say(pr, "f2 at start point = %g\n", f1);
    say(pr, "Distance to the next break point =  %g\n", dt);
    say(pr, "Distance to the stationary point =  %g\n", dtm);
  }

  // Step through the breakpoints.  Returns 1 when every variable ended
  // up at a bound (Fortran label 999), 0 when the minimizer lies within
  // a segment (label 888), and -1 on a bmv failure.
  int breakpoints(int &nseg, int &info) {
    int nleft = nbreak, iter = 1;
    double tj = 0.0;
    for (;;) {
      double tj0 = tj;
      int ibp = nextBreak(iter, nleft, tj);
      double dt = tj - tj0;
      if (dt != 0.0 && iprint >= 100) printPiece(nseg, dt);
      if (dtm < dt) return 0;
      tsum = tsum + dt;
      nleft = nleft - 1;
      iter = iter + 1;
      double dibp = d(ibp);
      d(ibp) = 0.0;
      double zibp = fixAtBound(ibp, dibp);
      if (iprint >= 100) say(pr, "Variable fixed, index  %d\n", ibp);
      if (nleft == 0 && nbreak == n) {
        dtm = dt;
        return 1;
      }
      nseg = nseg + 1;
      double dibp2 = dibp * dibp;
      f1 = f1 + dt * f2 + dibp2 - theta * dibp * zibp;
      f2 = f2 - theta * dibp2;
      if (col > 0) {
        updateCurvature(ibp, dt, dibp, dibp2, info);
        if (info != 0) return -1;
      }
      f2 = dmax(epsmch * f2_org, f2);
      if (nleft == 0) break;
      dtm = -f1 / f2;
    }
    if (bnded) {
      f1 = 0.0;
      f2 = 0.0;
      dtm = 0.0;
    } else {
      dtm = -f1 / f2;
    }
    return 0;
  }

  // Move xcp to the minimizer along the last segment (label 888)
  void finishSegment(int nseg) {
    if (iprint >= 99) {
      say(pr, "Piece %d\n", nseg);
      say(pr, "f1 at start point = %g\n", f1);
      say(pr, "f2 at start point = %g\n", f1);
      say(pr, "Distance to the stationary point =  %g\n", dtm);
    }
    if (dtm <= 0.0) dtm = 0.0;
    tsum = tsum + dtm;
    daxpy(n, tsum, &d(1), &xcp(1));
  }
};

// f: lbfgsb.f:1315
void cauchy(int n, double *x_, const double *l_, const double *u_,
            const int *nbd_, double *g_, int *iorder_, int *iwhere_,
            double *t_, double *d_, double *xcp_, int m, double *wy_,
            double *ws_, double *sy_, double *wt_, double theta, int col,
            int head, double *p_, double *c_, double *wbp_, double *v_,
            int &nseg, int iprint, double sbgnrm, int &info, double epsmch,
            const Printer *pr) {
  if (sbgnrm <= 0.0) {
    if (iprint >= 0) say(pr, "Subgnorm =0, GCP = X.\n");
    dcopy(n, x_, xcp_);
    return;
  }
  Cauchy cp(n, x_, l_, u_, nbd_, g_, iorder_, iwhere_, t_, d_, xcp_, m, wy_,
            ws_, sy_, wt_, theta, col, head, p_, c_, wbp_, v_, iprint,
            epsmch, pr);
  if (!cp.start()) return;
  if (!cp.initSlopes(info)) return;
  nseg = 1;
  if (iprint >= 99) say(pr, "no. of breakpoints = %d\n", cp.nbreak);

  int atBounds = (cp.nbreak == 0) ? 0 : cp.breakpoints(nseg, info);
  if (atBounds < 0) return;
  if (atBounds == 0) cp.finishSegment(nseg);

  if (col > 0) daxpy(cp.col2, cp.dtm, p_, c_);
  cp.printXcp();
  if (iprint >= 99) say(pr, "--- exit CAUCHY---\n");
}

// f: lbfgsb.f:1848
void cmprlb(int n, int m, double *x_, double *g_, double *ws_, double *wy_,
            double *sy_, double *wt_, double *z_, double *r_, double *wa_,
            int *index_, double theta, int col, int head, int nfree,
            int cnstnd, int &info) {
  V1<double> x(x_), g(g_), z(z_), r(r_), wa(wa_);
  V1<int> index(index_);
  V2<double> ws(ws_, n), wy(wy_, n);
  if (cnstnd == 0 && col > 0) {
    for (int i = 1; i <= n; ++i) r(i) = -g(i);
  } else {
    for (int i = 1; i <= nfree; ++i) {
      int k = index(i);
      r(i) = -theta * (z(k) - x(k)) - g(k);
    }
    bmv(m, sy_, wt_, col, &wa(2 * m + 1), &wa(1), info);
    if (info != 0) {
      info = -8;
      return;
    }
    int pointr = head;
    for (int j = 1; j <= col; ++j) {
      double a1 = wa(j);
      double a2 = theta * wa(col + j);
      for (int i = 1; i <= nfree; ++i) {
        int k = index(i);
        r(i) = r(i) + wy(k, pointr) * a1 + ws(k, pointr) * a2;
      }
      pointr = pointr % m + 1;
    }
  }
}

// f: lbfgsb.f:1916
void errclb(int n, int m, double factr, const double *l_, const double *u_,
            const int *nbd_, int &itask, int &info, int &k,
            const Printer *pr) {
  V1<const double> l(l_), u(u_);
  V1<const int> nbd(nbd_);
  if (n <= 0) {
    itask = 13;
    say(pr, "  ERROR: N .LE. 0\n");
    return;
  }
  if (m <= 0) {
    say(pr, "  ERROR: M .LE. 0\n");
    return;
  }
  if (factr <= 0.0) {
    say(pr, "  ERROR: FACTR .LT. 0\n");
    return;
  }
  for (int i = 1; i <= n; ++i) {
    if (nbd(i) < 0 || nbd(i) > 3) {
      itask = 12;
      info = -6;
      k = i;
    }
    if (nbd(i) == 2) {
      if (l(i) > u(i)) {
        itask = 14;
        info = -7;
        k = i;
      }
    }
  }
}

// Pieces of formk (lbfgsb.f:1991), which forms the LEL^T factorization
// of the indefinite middle matrix K.  wn1 holds the inner products that
// are updated in place between iterations.
struct FormkMats {
  int n, m, head;
  V2<double> wn1, ws, wy;
  FormkMats(int n_, int m_, int head_, double *wn1_, double *ws_,
            double *wy_)
      : n(n_), m(m_), head(head_), wn1(wn1_, 2 * m_), ws(ws_, n_),
        wy(wy_, n_) {}
};

// Pointer to the newest column of the circular storage
int newestPointer(int head, int col, int m) {
  int ptr = head + col - 1;
  if (ptr > m) ptr = ptr - m;
  return ptr;
}

// Shift the old part of wn1 once the memory is full
void formkShift(FormkMats &k) {
  int m = k.m;
  for (int jy = 1; jy <= m - 1; ++jy) {
    int js = m + jy;
    dcopy(m - jy, &k.wn1(jy + 1, jy + 1), &k.wn1(jy, jy));
    dcopy(m - jy, &k.wn1(js + 1, js + 1), &k.wn1(js, js));
    dcopy(m - 1, &k.wn1(k.m + 2, jy + 1), &k.wn1(k.m + 1, jy));
  }
}

// New row (and the last column of the lower-left block) of wn1 from the
// free variables ind(1:nsub) and active variables ind(nsub+1:n)
void formkNewRow(FormkMats &k, const V1<int> &ind, int nsub, int col) {
  int m = k.m, iy = col, is = m + col;
  int ipntr = newestPointer(k.head, col, m);
  int jpntr = k.head;
  for (int jy = 1; jy <= col; ++jy) {
    int js = m + jy;
    double temp1 = 0.0, temp2 = 0.0, temp3 = 0.0;
    for (int kk = 1; kk <= nsub; ++kk) {
      int k1 = ind(kk);
      temp1 = temp1 + k.wy(k1, ipntr) * k.wy(k1, jpntr);
    }
    for (int kk = nsub + 1; kk <= k.n; ++kk) {
      int k1 = ind(kk);
      temp2 = temp2 + k.ws(k1, ipntr) * k.ws(k1, jpntr);
      temp3 = temp3 + k.ws(k1, ipntr) * k.wy(k1, jpntr);
    }
    k.wn1(iy, jy) = temp1;
    k.wn1(is, js) = temp2;
    k.wn1(is, jy) = temp3;
    jpntr = jpntr % m + 1;
  }
  int jy = col;
  jpntr = newestPointer(k.head, col, m);
  ipntr = k.head;
  for (int i = 1; i <= col; ++i) {
    is = m + i;
    double temp3 = 0.0;
    for (int kk = 1; kk <= nsub; ++kk) {
      int k1 = ind(kk);
      temp3 = temp3 + k.ws(k1, ipntr) * k.wy(k1, jpntr);
    }
    ipntr = ipntr % m + 1;
    k.wn1(is, jy) = temp3;
  }
}

// Update the diagonal blocks of wn1 for variables entering
// (indx2(1:nenter)) and leaving (indx2(ileave:n)) the free set
void formkDiagonal(FormkMats &k, const V1<int> &indx2, int nenter,
                   int ileave, int upcl) {
  int m = k.m, ipntr = k.head;
  for (int iy = 1; iy <= upcl; ++iy) {
    int is = m + iy, jpntr = k.head;
    for (int jy = 1; jy <= iy; ++jy) {
      int js = m + jy;
      double temp1 = 0.0, temp2 = 0.0, temp3 = 0.0, temp4 = 0.0;
      for (int kk = 1; kk <= nenter; ++kk) {
        int k1 = indx2(kk);
        temp1 = temp1 + k.wy(k1, ipntr) * k.wy(k1, jpntr);
        temp2 = temp2 + k.ws(k1, ipntr) * k.ws(k1, jpntr);
      }
      for (int kk = ileave; kk <= k.n; ++kk) {
        int k1 = indx2(kk);
        temp3 = temp3 + k.wy(k1, ipntr) * k.wy(k1, jpntr);
        temp4 = temp4 + k.ws(k1, ipntr) * k.ws(k1, jpntr);
      }
      k.wn1(iy, jy) = k.wn1(iy, jy) + temp1 - temp3;
      k.wn1(is, js) = k.wn1(is, js) - temp2 + temp4;
      jpntr = jpntr % m + 1;
    }
    ipntr = ipntr % m + 1;
  }
}

// Update the off-diagonal block of wn1 for entering/leaving variables
void formkOffDiagonal(FormkMats &k, const V1<int> &indx2, int nenter,
                      int ileave, int upcl) {
  int m = k.m, ipntr = k.head;
  for (int is = m + 1; is <= m + upcl; ++is) {
    int jpntr = k.head;
    for (int jy = 1; jy <= upcl; ++jy) {
      double temp1 = 0.0, temp3 = 0.0;
      for (int kk = 1; kk <= nenter; ++kk) {
        int k1 = indx2(kk);
        temp1 = temp1 + k.ws(k1, ipntr) * k.wy(k1, jpntr);
      }
      for (int kk = ileave; kk <= k.n; ++kk) {
        int k1 = indx2(kk);
        temp3 = temp3 + k.ws(k1, ipntr) * k.wy(k1, jpntr);
      }
      if (is <= jy + m) {
        k.wn1(is, jy) = k.wn1(is, jy) + temp1 - temp3;
      } else {
        k.wn1(is, jy) = k.wn1(is, jy) - temp1 + temp3;
      }
      jpntr = jpntr % m + 1;
    }
    ipntr = ipntr % m + 1;
  }
}

// Form the upper triangle of wn from wn1
void formkAssemble(V2<double> &wn, const V2<double> &wn1, const V2<double> &sy,
                   int m, int col, double theta) {
  for (int iy = 1; iy <= col; ++iy) {
    int is = col + iy, is1 = m + iy;
    for (int jy = 1; jy <= iy; ++jy) {
      int js = col + jy, js1 = m + jy;
      wn(jy, iy) = wn1(iy, jy) / theta;
      wn(js, is) = wn1(is1, js1) * theta;
    }
    for (int jy = 1; jy <= iy - 1; ++jy) wn(jy, is) = -wn1(is1, jy);
    for (int jy = iy; jy <= col; ++jy) wn(jy, is) = wn1(is1, jy);
    wn(iy, iy) = wn(iy, iy) + sy(iy, iy);
  }
}

// Factor wn in place: the (1,1) block, then the (2,2) block after the
// Schur complement update
void formkFactor(double *wn_, int m, int col, int &info) {
  V2<double> wn(wn_, 2 * m);
  int m2 = 2 * m, col2 = 2 * col;
  dpofa(wn_, m2, col, info);
  if (info != 0) {
    info = -1;
    return;
  }
  for (int js = col + 1; js <= col2; ++js)
    dtrsl(wn_, m2, col, &wn(1, js), 11, info);
  for (int is = col + 1; is <= col2; ++is)
    for (int js = is; js <= col2; ++js)
      wn(is, js) = wn(is, js) + ddot(col, &wn(1, is), &wn(1, js));
  dpofa(&wn(col + 1, col + 1), m2, col, info);
  if (info != 0) info = -2;
}

// f: lbfgsb.f:1991
void formk(int n, int nsub, int *ind_, int nenter, int ileave, int *indx2_,
           int iupdat, int updatd, double *wn_, double *wn1_, int m,
           double *ws_, double *wy_, double *sy_, double theta, int col,
           int head, int &info) {
  V1<int> ind(ind_), indx2(indx2_);
  V2<double> wn(wn_, 2 * m), sy(sy_, m);
  FormkMats k(n, m, head, wn1_, ws_, wy_);
  int upcl = col;
  if (updatd == 1) {
    if (iupdat > m) formkShift(k);
    formkNewRow(k, ind, nsub, col);
    upcl = col - 1;
  }
  formkDiagonal(k, indx2, nenter, ileave, upcl);
  formkOffDiagonal(k, indx2, nenter, ileave, upcl);
  formkAssemble(wn, k.wn1, sy, m, col, theta);
  formkFactor(wn_, m, col, info);
}

// f: lbfgsb.f:2316
void formt(int m, double *wt_, double *sy_, double *ss_, int col,
           double theta, int &info) {
  V2<double> wt(wt_, m), sy(sy_, m), ss(ss_, m);
  for (int j = 1; j <= col; ++j) wt(1, j) = theta * ss(1, j);
  for (int i = 2; i <= col; ++i) {
    for (int j = i; j <= col; ++j) {
      int k1 = (i < j ? i : j) - 1;
      double ddum = 0.0;
      for (int k = 1; k <= k1; ++k) ddum = ddum + sy(i, k) * sy(j, k) / sy(k, k);
      wt(i, j) = ddum + theta * ss(i, j);
    }
  }
  dpofa(wt_, m, col, info);
  if (info != 0) info = -3;
}

// Variables leaving (indx2(ileave:n)) and entering (indx2(1:nenter))
// the free set since the last iteration (first part of freev)
void freevChanges(int n, int nfree, const V1<int> &index, int &nenter,
                  int &ileave, V1<int> &indx2, const V1<int> &iwhere,
                  int iprint, const Printer *pr) {
  for (int i = 1; i <= nfree; ++i) {
    int k = index(i);
    if (iwhere(k) > 0) {
      ileave = ileave - 1;
      indx2(ileave) = k;
      if (iprint >= 100)
        say(pr, "Variable k leaves the set of free variables for k = %d\n", k);
    }
  }
  for (int i = 1 + nfree; i <= n; ++i) {
    int k = index(i);
    if (iwhere(k) <= 0) {
      nenter = nenter + 1;
      indx2(nenter) = k;
      if (iprint >= 100) say(pr, "Var entering free vars is k= %d\n", k);
    }
  }
  if (iprint >= 99) {
    say(pr, " no. variables leaving  = %d\n", n + 1 - ileave);
    say(pr, " no. variables entering = %d\n", nenter);
  }
}

// f: lbfgsb.f:2384
void freev(int n, int &nfree, int *index_, int &nenter, int &ileave,
           int *indx2_, int *iwhere_, int &wrk, int updatd, int cnstnd,
           int iprint, int iter, const Printer *pr) {
  V1<int> index(index_), indx2(indx2_), iwhere(iwhere_);
  nenter = 0;
  ileave = n + 1;
  if (iter > 0 && cnstnd == 1)
    freevChanges(n, nfree, index, nenter, ileave, indx2, iwhere, iprint, pr);
  wrk = ((ileave < n + 1) || (nenter > 0) || updatd == 1) ? 1 : 0;
  // free variables first, then the active ones from the end
  nfree = 0;
  int iact = n + 1;
  for (int i = 1; i <= n; ++i) {
    if (iwhere(i) <= 0) {
      nfree = nfree + 1;
      index(nfree) = i;
    } else {
      iact = iact - 1;
      index(iact) = i;
    }
  }
  if (iprint >= 99) {
    say(pr, " no. variables free = %d\n", nfree);
    say(pr, " at GCP  %d\n", iter + 1);
  }
}

// The four cases of dcstep (lbfgsb.f:3622); each returns the new trial
// step stpf.

// Case 1: a higher function value; the minimum is bracketed
double dcstepHigher(double stx, double fx, double dx, double stp, double fp,
                    double dp) {
  const double two = 2.0, three = 3.0;
  double theta = three * (fx - fp) / (stp - stx) + dx + dp;
  double s = dmax3(fabs(theta), fabs(dx), fabs(dp));
  double a = theta / s;
  double gamma = s * sqrt(a * a - (dx / s) * (dp / s));
  if (stp < stx) gamma = -gamma;
  double p = (gamma - dx) + theta;
  double q = ((gamma - dx) + gamma) + dp;
  double r = p / q;
  double stpc = stx + r * (stp - stx);
  double stpq =
      stx + ((dx / ((fx - fp) / (stp - stx) + dx)) / two) * (stp - stx);
  if (fabs(stpc - stx) < fabs(stpq - stx)) return stpc;
  return stpc + (stpq - stpc) / two;
}

// Case 2: lower function value, derivatives of opposite sign
double dcstepOpposite(double stx, double fx, double dx, double stp,
                      double fp, double dp) {
  const double three = 3.0;
  double theta = three * (fx - fp) / (stp - stx) + dx + dp;
  double s = dmax3(fabs(theta), fabs(dx), fabs(dp));
  double a = theta / s;
  double gamma = s * sqrt(a * a - (dx / s) * (dp / s));
  if (stp > stx) gamma = -gamma;
  double p = (gamma - dp) + theta;
  double q = ((gamma - dp) + gamma) + dx;
  double r = p / q;
  double stpc = stp + r * (stx - stp);
  double stpq = stp + (dp / (dp - dx)) * (stx - stp);
  if (fabs(stpc - stp) > fabs(stpq - stp)) return stpc;
  return stpq;
}

// Case 3: lower function value, same-sign derivatives, and the
// magnitude of the derivative decreases
double dcstepDecreasing(double stx, double fx, double dx, double sty,
                        double stp, double fp, double dp, bool brackt,
                        double stpmin, double stpmax) {
  const double p66 = 0.66, three = 3.0;
  double theta = three * (fx - fp) / (stp - stx) + dx + dp;
  double s = dmax3(fabs(theta), fabs(dx), fabs(dp));
  double a = theta / s;
  double b = a * a - (dx / s) * (dp / s);
  double gamma = s * sqrt(dmax(0.0, b));
  if (stp > stx) gamma = -gamma;
  double p = (gamma - dp) + theta;
  double q = (gamma + (dx - dp)) + gamma;
  double r = p / q;
  double stpc;
  if (r < 0.0 && gamma != 0.0) {
    stpc = stp + r * (stx - stp);
  } else if (stp > stx) {
    stpc = stpmax;
  } else {
    stpc = stpmin;
  }
  double stpq = stp + (dp / (dp - dx)) * (stx - stp);
  double stpf;
  if (brackt) {
    stpf = (fabs(stpc - stp) < fabs(stpq - stp)) ? stpc : stpq;
    if (stp > stx) return dmin(stp + p66 * (sty - stp), stpf);
    return dmax(stp + p66 * (sty - stp), stpf);
  }
  stpf = (fabs(stpc - stp) > fabs(stpq - stp)) ? stpc : stpq;
  stpf = dmin(stpmax, stpf);
  return dmax(stpmin, stpf);
}

// Case 4: lower function value, same-sign derivatives, and the
// magnitude of the derivative does not decrease
double dcstepNotDecreasing(double stx, double sty, double fy, double dy,
                           double stp, double fp, double dp, bool brackt,
                           double stpmin, double stpmax) {
  const double three = 3.0;
  if (brackt) {
    double theta = three * (fp - fy) / (sty - stp) + dy + dp;
    double s = dmax3(fabs(theta), fabs(dy), fabs(dp));
    double a = theta / s;
    double gamma = s * sqrt(a * a - (dy / s) * (dp / s));
    if (stp > sty) gamma = -gamma;
    double p = (gamma - dp) + theta;
    double q = ((gamma - dp) + gamma) + dy;
    double r = p / q;
    return stp + r * (sty - stp);
  }
  if (stp > stx) return stpmax;
  return stpmin;
}

// f: lbfgsb.f:3622
void dcstep(double &stx, double &fx, double &dx, double &sty, double &fy,
            double &dy, double &stp, double fp, double dp, bool &brackt,
            double stpmin, double stpmax) {
  double sgnd = dp * (dx / fabs(dx));
  double stpf;
  if (fp > fx) {
    stpf = dcstepHigher(stx, fx, dx, stp, fp, dp);
    brackt = true;
  } else if (sgnd < 0.0) {
    stpf = dcstepOpposite(stx, fx, dx, stp, fp, dp);
    brackt = true;
  } else if (fabs(dp) < fabs(dx)) {
    stpf = dcstepDecreasing(stx, fx, dx, sty, stp, fp, dp, brackt, stpmin,
                            stpmax);
  } else {
    stpf = dcstepNotDecreasing(stx, sty, fy, dy, stp, fp, dp, brackt,
                               stpmin, stpmax);
  }

  // update the interval that contains a minimizer
  if (fp > fx) {
    sty = stp;
    fy = fp;
    dy = dp;
  } else {
    if (sgnd < 0.0) {
      sty = stx;
      fy = fx;
      dy = dx;
    }
    stx = stp;
    fx = fp;
    dx = dp;
  }
  stp = stpf;
}

// Line-search state kept between dcsrch (lbfgsb.f:3255) calls in
// isave(1:2) and dsave(1:13)
struct DcsrchState {
  bool brackt = false;
  int stage = 0;
  double ginit = 0, gtest = 0, gx = 0, gy = 0, finit = 0, fx = 0, fy = 0,
         stx = 0, sty = 0, stmin = 0, stmax = 0, width = 0, width1 = 0;

  void load(const V1<int> &isave, const V1<double> &dsave) {
    brackt = (isave(1) == 1);
    stage = isave(2);
    ginit = dsave(1);
    gtest = dsave(2);
    gx = dsave(3);
    gy = dsave(4);
    finit = dsave(5);
    fx = dsave(6);
    fy = dsave(7);
    stx = dsave(8);
    sty = dsave(9);
    stmin = dsave(10);
    stmax = dsave(11);
    width = dsave(12);
    width1 = dsave(13);
  }

  void save(V1<int> &isave, V1<double> &dsave) const {
    isave(1) = brackt ? 1 : 0;
    isave(2) = stage;
    dsave(1) = ginit;
    dsave(2) = gtest;
    dsave(3) = gx;
    dsave(4) = gy;
    dsave(5) = finit;
    dsave(6) = fx;
    dsave(7) = fy;
    dsave(8) = stx;
    dsave(9) = sty;
    dsave(10) = stmin;
    dsave(11) = stmax;
    dsave(12) = width;
    dsave(13) = width1;
  }
};

// Error code for invalid dcsrch arguments, or itask unchanged
int dcsrchArgCheck(int itask, double g, double stp, double ftol,
                   double gtol, double xtol, double stpmin, double stpmax) {
  if (stp < stpmin) itask = 16;
  if (stp > stpmax) itask = 15;
  if (g >= 0.0) itask = 11;
  if (ftol < 0.0) itask = 9;
  if (gtol < 0.0) itask = 10;
  if (xtol < 0.0) itask = 19;
  if (stpmin < 0.0) itask = 18;
  if (stpmax < stpmin) itask = 17;
  return itask;
}

void dcsrchStart(DcsrchState &st, double f, double g, double stp,
                 double ftol, double stpmin, double stpmax) {
  const double p5 = 0.5, xtrapu = 4.0;
  st.brackt = false;
  st.stage = 1;
  st.finit = f;
  st.ginit = g;
  st.gtest = ftol * st.ginit;
  st.width = stpmax - stpmin;
  st.width1 = st.width / p5;
  st.stx = 0.0;
  st.fx = st.finit;
  st.gx = st.ginit;
  st.sty = 0.0;
  st.fy = st.finit;
  st.gy = st.ginit;
  st.stmin = 0.0;
  st.stmax = stp + xtrapu * stp;
}

// Convergence and warning tests; returns the (possibly updated) itask
int dcsrchTests(const DcsrchState &st, int itask, double f, double g,
                double stp, double ftest, double gtol, double xtol,
                double stpmin, double stpmax) {
  if (st.brackt && (stp <= st.stmin || stp >= st.stmax)) itask = 23;
  if (st.brackt && st.stmax - st.stmin <= xtol * st.stmax) itask = 26;
  if (stp == stpmax && f <= ftest && g <= st.gtest) itask = 24;
  if (stp == stpmin && (f > ftest || g >= st.gtest)) itask = 25;
  if (f <= ftest && fabs(g) <= gtol * (-st.ginit)) itask = 6;
  return itask;
}

// Use dcstep to compute a new step (with the modified function in the
// first stage)
void dcsrchStep(DcsrchState &st, double f, double g, double &stp,
                double ftest) {
  if (st.stage == 1 && f <= st.fx && f > ftest) {
    double fm = f - stp * st.gtest;
    double fxm = st.fx - st.stx * st.gtest;
    double fym = st.fy - st.sty * st.gtest;
    double gm = g - st.gtest;
    double gxm = st.gx - st.gtest;
    double gym = st.gy - st.gtest;
    dcstep(st.stx, fxm, gxm, st.sty, fym, gym, stp, fm, gm, st.brackt,
           st.stmin, st.stmax);
    st.fx = fxm + st.stx * st.gtest;
    st.fy = fym + st.sty * st.gtest;
    st.gx = gxm + st.gtest;
    st.gy = gym + st.gtest;
  } else {
    dcstep(st.stx, st.fx, st.gx, st.sty, st.fy, st.gy, stp, f, g,
           st.brackt, st.stmin, st.stmax);
  }
}

// Safeguard the new step and update the interval of uncertainty
void dcsrchSafeguard(DcsrchState &st, double &stp, double xtol,
                     double stpmin, double stpmax) {
  const double p5 = 0.5, p66 = 0.66, xtrapl = 1.1, xtrapu = 4.0;
  if (st.brackt) {
    if (fabs(st.sty - st.stx) >= p66 * st.width1)
      stp = st.stx + p5 * (st.sty - st.stx);
    st.width1 = st.width;
    st.width = fabs(st.sty - st.stx);
  }
  if (st.brackt) {
    st.stmin = dmin(st.stx, st.sty);
    st.stmax = dmax(st.stx, st.sty);
  } else {
    st.stmin = stp + xtrapl * (stp - st.stx);
    st.stmax = stp + xtrapu * (stp - st.stx);
  }
  stp = dmax(stp, stpmin);
  stp = dmin(stp, stpmax);
  if ((st.brackt && (stp <= st.stmin || stp >= st.stmax)) ||
      (st.brackt && st.stmax - st.stmin <= xtol * st.stmax))
    stp = st.stx;
}

// f: lbfgsb.f:3255
void dcsrch(double f, double g, double &stp, double ftol, double gtol,
            double xtol, double stpmin, double stpmax, int &itask,
            int *isave_, double *dsave_) {
  V1<int> isave(isave_);
  V1<double> dsave(dsave_);
  DcsrchState st;

  if (itask == 2) {
    itask = dcsrchArgCheck(itask, g, stp, ftol, gtol, xtol, stpmin, stpmax);
    if ((itask >= 9) && (itask <= 19)) return;
    dcsrchStart(st, f, g, stp, ftol, stpmin, stpmax);
    itask = 4;
    st.save(isave, dsave);
    return;
  }
  st.load(isave, dsave);

  double ftest = st.finit + stp * st.gtest;
  if (st.stage == 1 && f <= ftest && g >= 0.0) st.stage = 2;

  itask = dcsrchTests(st, itask, f, g, stp, ftest, gtol, xtol, stpmin,
                      stpmax);
  if (!((itask >= 23) || ((itask <= 8) && (itask >= 6)))) {
    dcsrchStep(st, f, g, stp, ftest);
    dcsrchSafeguard(st, stp, xtol, stpmin, stpmax);
    itask = 4;
  }
  st.save(isave, dsave);
}

// Largest step to the boundary of the feasible region along d (part of
// lnsrlb)
double lnsrlbMaxStep(int n, const V1<const double> &l,
                     const V1<const double> &u, const V1<const int> &nbd,
                     const V1<double> &x, const V1<double> &d) {
  double stpmx = 1.0e10;
  for (int i = 1; i <= n; ++i) {
    double a1 = d(i);
    if (nbd(i) == 0) continue;
    if (a1 < 0.0 && nbd(i) <= 2) {
      double a2 = l(i) - x(i);
      if (a2 >= 0.0) {
        stpmx = 0.0;
      } else if (a1 * stpmx < a2) {
        stpmx = a2 / a1;
      }
    } else if (a1 > 0.0 && nbd(i) >= 2) {
      double a2 = u(i) - x(i);
      if (a2 <= 0.0) {
        stpmx = 0.0;
      } else if (a1 * stpmx > a2) {
        stpmx = a2 / a1;
      }
    }
  }
  return stpmx;
}

// f: lbfgsb.f:2618
void lnsrlb(int n, const double *l_, const double *u_, const int *nbd_,
            double *x_, double &f, double &fold, double &gd, double &gdold,
            double *g_, double *d_, double *r_, double *t_, double *z_,
            double &stp, double &dnorm, double &dtd, double &xstep,
            double &stpmx, int iter, int &ifun, int &iback, int &nfgv,
            int &info, int &itask, int boxed, int cnstnd, int &icsave,
            int *isave, double *dsave, const Printer *pr) {
  V1<const double> l(l_), u(u_);
  V1<const int> nbd(nbd_);
  V1<double> x(x_), d(d_), t(t_);
  const double one = 1.0, big = 1.0e10;
  const double ftol = 1.0e-3, gtol = 0.9, xtol = 0.1;

  if (itask != 20) {
    // start a new line search
    dtd = ddot(n, d_, d_);
    dnorm = sqrt(dtd);
    stpmx = big;
    if (cnstnd == 1) stpmx = (iter == 0) ? one : lnsrlbMaxStep(n, l, u, nbd, x, d);
    stp = (iter == 0 && boxed == 0) ? dmin(one / dnorm, stpmx) : one;
    dcopy(n, x_, t_);
    dcopy(n, g_, r_);
    fold = f;
    ifun = 0;
    iback = 0;
    icsave = 2;
  }

  gd = ddot(n, g_, d_);
  if (ifun == 0) {
    gdold = gd;
    if (gd >= 0.0) {
      say(pr, " ascent direction in projection gd =  %g\n", gd);
      info = -4;
      return;
    }
  }

  dcsrch(f, gd, stp, ftol, gtol, xtol, 0.0, stpmx, icsave, isave, dsave);

  xstep = stp * dnorm;
  if ((icsave < 6) || ((icsave > 8) && (icsave < 23))) {
    itask = 20;
    ifun = ifun + 1;
    nfgv = nfgv + 1;
    iback = ifun - 1;
    if (stp == one) {
      dcopy(n, z_, x_);
    } else {
      for (int i = 1; i <= n; ++i) x(i) = stp * d(i) + t(i);
    }
  } else {
    itask = 1;
  }
}

// f: lbfgsb.f:2756
void matupd(int n, int m, double *ws_, double *wy_, double *sy_,
            double *ss_, double *d_, double *r_, int &itail, int iupdat,
            int &col, int &head, double &theta, double rr, double dr,
            double stp, double dtd) {
  V2<double> ws(ws_, n), wy(wy_, n), sy(sy_, m), ss(ss_, m);
  if (iupdat <= m) {
    col = iupdat;
    itail = (head + iupdat - 2) % m + 1;
  } else {
    itail = itail % m + 1;
    head = head % m + 1;
  }
  dcopy(n, d_, &ws(1, itail));
  dcopy(n, r_, &wy(1, itail));
  theta = rr / dr;
  if (iupdat > m) {
    for (int j = 1; j <= col - 1; ++j) {
      dcopy(j, &ss(2, j + 1), &ss(1, j));
      dcopy(col - j, &sy(j + 1, j + 1), &sy(j, j));
    }
  }
  int pointr = head;
  for (int j = 1; j <= col - 1; ++j) {
    sy(col, j) = ddot(n, d_, &wy(1, pointr));
    ss(j, col) = ddot(n, &ws(1, pointr), d_);
    pointr = pointr % m + 1;
  }
  if (stp == 1.0) {
    ss(col, col) = dtd;
  } else {
    ss(col, col) = stp * stp * dtd;
  }
  sy(col, col) = dr;
}

// f: lbfgsb.f:2844
void projgr(int n, const double *l_, const double *u_, const int *nbd_,
            const double *x_, const double *g_, double &sbgnrm) {
  V1<const double> l(l_), u(u_), x(x_), g(g_);
  V1<const int> nbd(nbd_);
  sbgnrm = 0.0;
  for (int i = 1; i <= n; ++i) {
    double gi = g(i);
    if (nbd(i) != 0) {
      if (gi < 0.0) {
        if (nbd(i) >= 2) gi = dmax((x(i) - u(i)), gi);
      } else {
        if (nbd(i) <= 2) gi = dmin((x(i) - l(i)), gi);
      }
    }
    sbgnrm = dmax(sbgnrm, fabs(gi));
  }
}

// Subspace Newton direction d = (1/theta) r + (1/theta^2) Z'W K^-1 W'Z r
// (first part of subsm)
void subsmDirection(int n, int m, int nsub, const V1<int> &ind,
                    double *d_, double *ws_, double *wy_, double theta,
                    int col, int head, double *wv_, double *wn_,
                    int &info) {
  V1<double> d(d_), wv(wv_);
  V2<double> ws(ws_, n), wy(wy_, n);
  int pointr = head;
  for (int i = 1; i <= col; ++i) {
    double temp1 = 0.0, temp2 = 0.0;
    for (int j = 1; j <= nsub; ++j) {
      int k = ind(j);
      temp1 = temp1 + wy(k, pointr) * d(j);
      temp2 = temp2 + ws(k, pointr) * d(j);
    }
    wv(i) = temp1;
    wv(col + i) = theta * temp2;
    pointr = pointr % m + 1;
  }
  int m2 = 2 * m, col2 = 2 * col;
  dtrsl(wn_, m2, col2, wv_, 11, info);
  if (info != 0) return;
  for (int i = 1; i <= col; ++i) wv(i) = -wv(i);
  dtrsl(wn_, m2, col2, wv_, 1, info);
  if (info != 0) return;
  pointr = head;
  for (int jy = 1; jy <= col; ++jy) {
    int js = col + jy;
    for (int i = 1; i <= nsub; ++i) {
      int k = ind(i);
      d(i) = d(i) + wy(k, pointr) * wv(jy) / theta + ws(k, pointr) * wv(js);
    }
    pointr = pointr % m + 1;
  }
  dscal(nsub, 1.0 / theta, d_);
}

// x(k) + dk projected onto its bounds; sets iword = 1 when it lands on
// a bound
double subsmProject(int nbd, double l, double u, double xk, double dk,
                    double xcur, int &iword) {
  if (nbd == 0) return xk + dk;
  if (nbd == 1) {
    double xn = dmax(l, xk + dk);
    if (xn == l) iword = 1;
    return xn;
  }
  if (nbd == 2) {
    xk = dmax(l, xk + dk);
    double xn = dmin(u, xk);
    if (xn == l || xn == u) iword = 1;
    return xn;
  }
  if (nbd == 3) {
    double xn = dmin(u, xk + dk);
    if (xn == u) iword = 1;
    return xn;
  }
  return xcur;
}

// Step-length limit for one variable in the backtracking step; temp1
// carries over between variables as in the Fortran loop
void subsmLimit(int nbd, double l, double u, double xk, double dk,
                double alpha, double &temp1) {
  if (dk < 0.0 && nbd <= 2) {
    double temp2 = l - xk;
    if (temp2 >= 0.0) {
      temp1 = 0.0;
    } else if (dk * alpha < temp2) {
      temp1 = temp2 / dk;
    }
  } else if (dk > 0.0 && nbd >= 2) {
    double temp2 = u - xk;
    if (temp2 <= 0.0) {
      temp1 = 0.0;
    } else if (dk * alpha > temp2) {
      temp1 = temp2 / dk;
    }
  }
}

// Backtrack to the largest feasible step along d (end of subsm)
void subsmBacktrack(int nsub, const V1<int> &ind, const V1<const double> &l,
                    const V1<const double> &u, const V1<const int> &nbd,
                    V1<double> &x, V1<double> &d) {
  double alpha = 1.0, temp1 = alpha;
  int ibd = 0;
  for (int i = 1; i <= nsub; ++i) {
    int k = ind(i);
    if (nbd(k) == 0) continue;
    subsmLimit(nbd(k), l(k), u(k), x(k), d(i), alpha, temp1);
    if (temp1 < alpha) {
      alpha = temp1;
      ibd = i;
    }
  }
  if (alpha < 1.0) {
    double dk = d(ibd);
    int k = ind(ibd);
    if (dk > 0.0) {
      x(k) = u(k);
      d(ibd) = 0.0;
    } else if (dk < 0.0) {
      x(k) = l(k);
      d(ibd) = 0.0;
    }
  }
  for (int i = 1; i <= nsub; ++i) {
    int k = ind(i);
    x(k) = x(k) + alpha * d(i);
  }
}

// f: lbfgsb.f:2893
void subsm(int n, int m, int nsub, int *ind_, const double *l_,
           const double *u_, const int *nbd_, double *x_, double *d_,
           double *xp_, double *ws_, double *wy_, double theta,
           double *xx_, double *gg_, int col, int head, int &iword,
           double *wv_, double *wn_, int iprint, int &info,
           const Printer *pr) {
  V1<int> ind(ind_);
  V1<const double> l(l_), u(u_);
  V1<const int> nbd(nbd_);
  V1<double> x(x_), d(d_), xx(xx_), gg(gg_);

  if (nsub <= 0) return;
  if (iprint >= 99) say(pr, " ----- SUBSM entered -----\n");

  subsmDirection(n, m, nsub, ind, d_, ws_, wy_, theta, col, head, wv_, wn_,
                 info);
  if (info != 0) return;

  // project the subspace minimizer onto the bounds
  iword = 0;
  dcopy(n, x_, xp_);
  for (int i = 1; i <= nsub; ++i) {
    int k = ind(i);
    x(k) = subsmProject(nbd(k), l(k), u(k), x(k), d(i), x(k), iword);
  }

  if (iword != 0) {
    double dd_p = 0.0;
    for (int i = 1; i <= n; ++i) dd_p = dd_p + (x(i) - xx(i)) * gg(i);
    if (dd_p > 0.0) {
      dcopy(n, xp_, x_);
      say(pr, " Positive dir derivative in projection \n");
      say(pr, " Using the backtracking step \n");
      subsmBacktrack(nsub, ind, l, u, nbd, x, d);
    }
  }

  if (iprint >= 99) say(pr, " exit SUBSM \n");
}

// mainlb (lbfgsb.f:305) as a small state machine.  The Fortran labels
// become steps: 111 = firstEval, 222 = iterate, 666 = lineSearch,
// 777 = newX, 999/1000 = save and return to the driver.  Local variables
// that survive between calls are saved in lsave/isave/dsave exactly as
// in the Fortran.  The CPU timers are disabled (always zero) as in the
// R version of the Fortran, so their bookkeeping is reduced to the saved
// zeros.
struct Mainlb {
  enum Step { kStart, kFirstEval, kIterate, kLineSearch, kNewX, kSave };

  int n, m;
  double *x;
  const double *l, *u;
  const int *nbd;
  double &f;
  double *g;
  double factr, pgtol;
  double *ws, *wy, *sy, *ss, *wt, *wn, *snd, *z, *r, *d, *t, *xp, *wa;
  int *index, *iwhere, *indx2;
  int &itask;
  int iprint;
  int &icsave;
  V1<int> lsave, isave;
  V1<double> dsave;
  const Printer *pr;

  int prjctd = 0, cnstnd = 0, boxed = 0, updatd = 0, wrk = 0;
  int nintol = 0, iback = 0, nskip = 0, head = 0, col = 0, iter = 0,
      itail = 0, iupdat = 0, nseg = 0, nfgv = 0, info = 0, ifun = 0,
      iword = 0, nfree = 0, nact = 0, ileave = 0, nenter = 0;
  double theta = 0, fold = 0, tol = 0, xstep = 0, sbgnrm = 0, dnorm = 0,
         dtd = 0, epsmch = 0, cpu1 = 0, sbtime = 0, time1 = 0, gd = 0,
         gdold = 0, stp = 0, stpmx = 0;

  Mainlb(int n_, int m_, double *x_, const double *l_, const double *u_,
         const int *nbd_, double &f_, double *g_, double factr_,
         double pgtol_, double *ws_, double *wy_, double *sy_, double *ss_,
         double *wt_, double *wn_, double *snd_, double *z_, double *r_,
         double *d_, double *t_, double *xp_, double *wa_, int *index_,
         int *iwhere_, int *indx2_, int &itask_, int iprint_,
         int &icsave_, int *lsave_, int *isave_, double *dsave_,
         const Printer *pr_)
      : n(n_), m(m_), x(x_), l(l_), u(u_), nbd(nbd_), f(f_), g(g_),
        factr(factr_), pgtol(pgtol_), ws(ws_), wy(wy_), sy(sy_), ss(ss_),
        wt(wt_), wn(wn_), snd(snd_), z(z_), r(r_), d(d_), t(t_), xp(xp_),
        wa(wa_), index(index_), iwhere(iwhere_), indx2(indx2_),
        itask(itask_), iprint(iprint_), icsave(icsave_), lsave(lsave_),
        isave(isave_), dsave(dsave_), pr(pr_) {}

  // itask == 2 (START); false when errclb found an input error
  bool init() {
    epsmch = DBL_EPSILON;
    head = 1;
    theta = 1.0;
    nfree = n;
    tol = factr * epsmch;
    int k = 0;
    errclb(n, m, factr, l, u, nbd, itask, info, k, pr);
    if ((itask >= 9) && (itask <= 19)) return false;
    active(n, l, u, nbd, x, iwhere, iprint, prjctd, cnstnd, boxed, pr);
    return true;
  }

  void load() {
    prjctd = lsave(1);
    cnstnd = lsave(2);
    boxed = lsave(3);
    updatd = lsave(4);
    nintol = isave(1);
    iback = isave(4);
    nskip = isave(5);
    head = isave(6);
    col = isave(7);
    itail = isave(8);
    iter = isave(9);
    iupdat = isave(10);
    nseg = isave(12);
    nfgv = isave(13);
    info = isave(14);
    ifun = isave(15);
    iword = isave(16);
    nfree = isave(17);
    nact = isave(18);
    ileave = isave(19);
    nenter = isave(20);
    theta = dsave(1);
    fold = dsave(2);
    tol = dsave(3);
    dnorm = dsave(4);
    epsmch = dsave(5);
    cpu1 = dsave(6);
    sbtime = dsave(8);
    time1 = dsave(10);
    gd = dsave(11);
    stpmx = dsave(12);
    sbgnrm = dsave(13);
    stp = dsave(14);
    gdold = dsave(15);
    dtd = dsave(16);
  }

  void save() {
    lsave(1) = prjctd;
    lsave(2) = cnstnd;
    lsave(3) = boxed;
    lsave(4) = updatd;
    isave(1) = nintol;
    isave(4) = iback;
    isave(5) = nskip;
    isave(6) = head;
    isave(7) = col;
    isave(8) = itail;
    isave(9) = iter;
    isave(10) = iupdat;
    isave(12) = nseg;
    isave(13) = nfgv;
    isave(14) = info;
    isave(15) = ifun;
    isave(16) = iword;
    isave(17) = nfree;
    isave(18) = nact;
    isave(19) = ileave;
    isave(20) = nenter;
    dsave(1) = theta;
    dsave(2) = fold;
    dsave(3) = tol;
    dsave(4) = dnorm;
    dsave(5) = epsmch;
    dsave(6) = cpu1;
    dsave(8) = sbtime;
    dsave(10) = time1;
    dsave(11) = gd;
    dsave(12) = stpmx;
    dsave(13) = sbgnrm;
    dsave(14) = stp;
    dsave(15) = gdold;
    dsave(16) = dtd;
  }

  // Discard the limited-memory matrices and restart from B = I
  void resetMemory() {
    info = 0;
    col = 0;
    head = 1;
    theta = 1.0;
    iupdat = 0;
    updatd = 0;
  }

  // Where to continue after returning from the driver with `itask`
  Step resumeAt() const {
    if (itask == 20) return kLineSearch;
    if (itask == 1) return kNewX;
    if (itask == 21) return kFirstEval;
    return kStart;
  }

  // Label 111: f and g at the starting point are available
  Step firstEval() {
    nfgv = 1;
    projgr(n, l, u, nbd, x, g, sbgnrm);
    if (iprint >= 1)
      say(pr, "At iterate %d f= %g |proj g|=  %g\n", iter, f, sbgnrm);
    if (sbgnrm <= pgtol) {
      itask = 7;
      return kSave;
    }
    return kIterate;
  }

  // Formk, cmprlb and subsm; false when the memory had to be reset
  bool subspaceMin() {
    if (wrk == 1)
      formk(n, nfree, index, nenter, ileave, indx2, iupdat, updatd, wn, snd,
            m, ws, wy, sy, theta, col, head, info);
    if (info != 0) {
      if (iprint >= 1)
        sayRestart(pr, " Nonpositive definiteness in Cholesky factorization in formk;");
      resetMemory();
      return false;
    }
    cmprlb(n, m, x, g, ws, wy, sy, wt, z, r, wa, index, theta, col, head,
           nfree, cnstnd, info);
    if (info == 0)
      subsm(n, m, nfree, index, l, u, nbd, z, r, xp, ws, wy, theta, x, g,
            col, head, iword, wa, wn, iprint, info, pr);
    if (info != 0) {
      sayRestart(pr, " Singular triangular system detected;");
      resetMemory();
      return false;
    }
    return true;
  }

  // Label 222: generalized Cauchy point, subspace minimization and the
  // search direction d = z - x
  Step iterate() {
    if (iprint >= 99) say(pr, "ITERATION  %d\n", iter + 1);
    iword = -1;
    if (cnstnd == 0 && col > 0) {
      dcopy(n, x, z);
      wrk = updatd;
      nseg = 0;
    } else {
      cauchy(n, x, l, u, nbd, g, indx2, iwhere, t, d, z, m, wy, ws, sy, wt,
             theta, col, head, &wa[0], &wa[2 * m], &wa[4 * m], &wa[6 * m],
             nseg, iprint, sbgnrm, info, epsmch, pr);
      if (info != 0) {
        sayRestart(pr, " Singular triangular system detected;");
        resetMemory();
        return kIterate;
      }
      nintol = nintol + nseg;
      freev(n, nfree, index, nenter, ileave, indx2, iwhere, wrk, updatd,
            cnstnd, iprint, iter, pr);
      nact = n - nfree;
    }
    if (nfree != 0 && col != 0 && !subspaceMin()) return kIterate;
    for (int i = 0; i < n; ++i) d[i] = z[i] - x[i];
    return kLineSearch;
  }

  // The line search failed: restore the previous iterate
  Step lineSearchFailed() {
    dcopy(n, t, x);
    dcopy(n, r, g);
    f = fold;
    if (col == 0) {
      if (info == 0) {
        info = -9;
        nfgv = nfgv - 1;
        ifun = ifun - 1;
        iback = iback - 1;
      }
      itask = 5;
      iter = iter + 1;
      return kSave;
    }
    if (iprint >= 1) sayRestart(pr, " Bad direction in the line search;");
    if (info == 0) nfgv = nfgv - 1;
    resetMemory();
    itask = 22;
    return kIterate;
  }

  // Label 666
  Step lineSearch() {
    lnsrlb(n, l, u, nbd, x, f, fold, gd, gdold, g, d, r, t, z, stp, dnorm,
           dtd, xstep, stpmx, iter, ifun, iback, nfgv, info, itask, boxed,
           cnstnd, icsave, &isave(22), &dsave(17), pr);
    if (info != 0 || iback >= 20) return lineSearchFailed();
    if (itask != 20) {
      iter = iter + 1;
      projgr(n, l, u, nbd, x, g, sbgnrm);
    }
    return kSave;
  }

  // BFGS update of the limited-memory matrices (end of label 777)
  void updateMemory() {
    V1<double> rv(r), gv(g);
    for (int i = 1; i <= n; ++i) rv(i) = gv(i) - rv(i);
    double rr = ddot(n, r, r), dr, ddum;
    if (stp == 1.0) {
      dr = gd - gdold;
      ddum = -gdold;
    } else {
      dr = (gd - gdold) * stp;
      dscal(n, stp, d);
      ddum = -gdold * stp;
    }
    if (dr <= epsmch * ddum) {
      nskip = nskip + 1;
      updatd = 0;
      if (iprint >= 1) {
        say(pr, " ys = %g\n", dr);
        say(pr, " BFGS update skipped for +gs= %g\n", ddum);
      }
      return;
    }
    updatd = 1;
    iupdat = iupdat + 1;
    matupd(n, m, ws, wy, sy, ss, d, r, itail, iupdat, col, head, theta, rr,
           dr, stp, dtd);
    formt(m, wt, sy, ss, col, theta, info);
    if (info != 0) {
      if (iprint >= 1)
        sayRestart(pr, " Nonpositive definiteness in Cholesky factorization in formt;");
      resetMemory();
    }
  }

  // Label 777: convergence tests at a new iterate, then the update
  Step newX() {
    if (sbgnrm <= pgtol) {
      itask = 7;
      return kSave;
    }
    double ddum = dmax3(fabs(fold), fabs(f), 1.0);
    if ((fold - f) <= tol * ddum) {
      itask = 8;
      if (iback >= 10) info = -5;
      return kSave;
    }
    updateMemory();
    return kIterate;
  }

  void run() {
    Step next;
    if (itask == 2) {
      if (!init()) return;
      next = kStart;
    } else {
      load();
      next = resumeAt();
    }
    for (;;) {
      switch (next) {
      case kStart:
        itask = 21;  // ask the driver for f0 and g0
        next = kSave;
        break;
      case kFirstEval:
        next = firstEval();
        break;
      case kIterate:
        next = iterate();
        break;
      case kLineSearch:
        next = lineSearch();
        break;
      case kNewX:
        next = newX();
        break;
      case kSave:
        save();
        return;
      }
    }
  }
};

// f: lbfgsb.f:305
void mainlb(int n, int m, double *x, const double *l, const double *u,
            const int *nbd, double &f, double *g, double factr, double pgtol,
            double *ws, double *wy, double *sy, double *ss, double *wt,
            double *wn, double *snd, double *z_, double *r_, double *d_,
            double *t, double *xp, double *wa, int *index, int *iwhere,
            int *indx2, int &itask, int iprint, int &icsave, int *lsave_,
            int *isave_, double *dsave_, const Printer *pr) {
  Mainlb lb(n, m, x, l, u, nbd, f, g, factr, pgtol, ws, wy, sy, ss, wt, wn,
            snd, z_, r_, d_, t, xp, wa, index, iwhere, indx2, itask, iprint,
            icsave, lsave_, isave_, dsave_, pr);
  lb.run();
}

const char *const kTaskNames[29] = {
    "NEW_X",
    "START",
    "STOP",
    "FG",
    "ABNORMAL_TERMINATION_IN_LNSRCH",
    "CONVERGENCE",
    "CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL",
    "CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH",
    "ERROR: FTOL .LT. ZERO",
    "ERROR: GTOL .LT. ZERO",
    "ERROR: INITIAL G .GE. ZERO",
    "ERROR: INVALID NBD",
    "ERROR: N .LE. 0",
    "ERROR: NO FEASIBLE SOLUTION",
    "ERROR: STP .GT. STPMAX",
    "ERROR: STP .LT. STPMIN",
    "ERROR: STPMAX .LT. STPMIN",
    "ERROR: STPMIN .LT. ZERO",
    "ERROR: XTOL .LT. ZERO",
    "FG_LNSRCH",
    "FG_START",
    "RESTART_FROM_LNSRCH",
    "WARNING: ROUNDING ERRORS PREVENT PROGRESS",
    "WARNING: STP .eq. STPMAX",
    "WARNING: STP .eq. STPMIN",
    "WARNING: XTOL TEST SATISFIED",
    "CONVERGENCE: Parameters differences below xtol",
    "Maximum number of iterations reached",
    "ERROR: INVALID LMM"};

} // namespace

const char *taskName(int itask) {
  if (itask < 1 || itask > 29) return "UNKNOWN";
  return kTaskNames[itask - 1];
}

bool validLmm(int n, int lmm) {
  if (lmm < 1) return false;
  // 11*lmm^2 alone must fit in an int; this bound also keeps the 64-bit
  // arithmetic below from overflowing for any int n
  if (lmm > 13972) return false;
  long long ln = n < 0 ? 0 : n, lm = lmm;
  long long nwa = 2 * lm * ln + 11 * lm * lm + 5 * ln + 8 * lm;
  return nwa <= (long long)INT_MAX;
}

// f: lbfgsb.f:49
void setulb(int n, int m, double *x, const double *l, const double *u,
            const int *nbd, double &f, double *g, double factr,
            double pgtol, double *wa, int *iwa, int &itask, int iprint,
            int &icsave, int *lsave, int *isave_, double *dsave,
            const Printer *pr) {
  V1<int> isave(isave_);
  if ((itask < 1) || (itask > 26)) {
    say(pr, "TASK NOT IN VALID RANGE\n");
    itask = -999;
    return;
  }
  if (itask == 2) {
    isave(1) = m * n;
    isave(2) = m * m;
    isave(3) = 4 * m * m;
    isave(4) = 1;                     // ws      m*n
    isave(5) = isave(4) + isave(1);   // wy      m*n
    isave(6) = isave(5) + isave(1);   // wsy     m**2
    isave(7) = isave(6) + isave(2);   // wss     m**2
    isave(8) = isave(7) + isave(2);   // wt      m**2
    isave(9) = isave(8) + isave(2);   // wn      4*m**2
    isave(10) = isave(9) + isave(3);  // wsnd    4*m**2
    isave(11) = isave(10) + isave(3); // wz      n
    isave(12) = isave(11) + n;        // wr      n
    isave(13) = isave(12) + n;        // wd      n
    isave(14) = isave(13) + n;        // wt      n
    isave(15) = isave(14) + n;        // wxp     n
    isave(16) = isave(15) + n;        // wa      8*m
  }
  int lws = isave(4), lwy = isave(5), lsy = isave(6), lss = isave(7),
      lwt = isave(8), lwn = isave(9), lsnd = isave(10), lz = isave(11),
      lr = isave(12), ld = isave(13), lt = isave(14), lxp = isave(15),
      lwa = isave(16);
  mainlb(n, m, x, l, u, nbd, f, g, factr, pgtol, &wa[lws - 1], &wa[lwy - 1],
         &wa[lsy - 1], &wa[lss - 1], &wa[lwt - 1], &wa[lwn - 1],
         &wa[lsnd - 1], &wa[lz - 1], &wa[lr - 1], &wa[ld - 1], &wa[lt - 1],
         &wa[lxp - 1], &wa[lwa - 1], &iwa[0], &iwa[n], &iwa[2 * n], itask,
         iprint, icsave, lsave, &isave(22), dsave, pr);
}

namespace {

// Solver workspace.  Each thread keeps one that only grows, so memory
// scales with the number of threads running the optimizer rather than
// the number of optimizations.
struct Workspace {
  std::vector<double> wa;
  std::vector<int> iwa;
  std::vector<double> lastx;
  bool inUse = false;

  // Size for (n, lmm) and zero the part a run uses
  void prepare(int n, int lmm) {
    size_t nwa = (size_t)(2 * lmm * n + 11 * lmm * lmm + 5 * n + 8 * lmm);
    size_t niwa = (size_t)(3 * n);
    if (nwa < 1) nwa = 1;
    if (niwa < 1) niwa = 1;
    if (wa.size() < nwa) wa.resize(nwa);
    if (iwa.size() < niwa) iwa.resize(niwa);
    if (lastx.size() < (size_t)n) lastx.resize(n);
    std::fill(wa.begin(), wa.begin() + nwa, 0.0);
    std::fill(iwa.begin(), iwa.begin() + niwa, 0);
  }
};

thread_local Workspace tlWorkspace;

// Claims this thread's workspace, or a private one when the thread's is
// already in use (an objective function that itself runs lbfgsb3Cts).
class WorkspaceLease {
public:
  WorkspaceLease() : ws_(tlWorkspace.inUse ? &own_ : &tlWorkspace) {
    ws_->inUse = true;
  }
  ~WorkspaceLease() { ws_->inUse = false; }
  WorkspaceLease(const WorkspaceLease &) = delete;
  WorkspaceLease &operator=(const WorkspaceLease &) = delete;
  Workspace &get() { return *ws_; }

private:
  Workspace own_;
  Workspace *ws_;
};

} // namespace

namespace {

// Set the outputs for input rejected before running; false if accepted
bool rejectInput(int n, int lmm, int *fail, const Printer *pr,
                 InfoOut *info) {
  int code = 0;
  if (n <= 0) {
    // errclb's check, done before the workspace is sized from n
    say(pr, "  ERROR: N .LE. 0\n");
    code = 13;
  } else if (!validLmm(n, lmm)) {
    // lmm <= 0 would divide by zero in matupd; a huge lmm overflows the
    // workspace size
    say(pr, "  ERROR: INVALID LMM (%d)\n", lmm);
    code = kInvalidLmm;
  }
  if (code == 0) return false;
  if (info != nullptr) {
    std::memset(info, 0, sizeof(*info));
    info->itask = code;
  }
  fail[0] = code;
  return true;
}

// True when every |lastx - x| < |x|*rtol + atol.  Checks x[n-1] and then
// x[n-2..0], as the original lbfgsb3C_ loop did.
bool xConverged(const double *lastx, const double *x, int n, double rtol,
                double atol) {
  bool converge =
      fabs(lastx[n - 1] - x[n - 1]) < fabs(x[n - 1]) * rtol + atol;
  for (int i = n - 1; converge && i--;)
    converge = fabs(lastx[i] - x[i]) < fabs(x[i]) * rtol + atol;
  return converge;
}

void traceBeforeStep(const Printer *pr, int itask, const double *x, int n) {
  say(pr, "itask: %d\n", itask);
  sayVec(pr, "computing f and g at prm=\n", x, n);
  say(pr, "\n================================================================================\nBefore call task number %d, or \"%s\"\n",
      itask, taskName(itask));
}

void traceEval(const Printer *pr, int trace, int iter, double f,
               const double *g, int n) {
  say(pr, "At iteration %d f=%f ", iter, f);
  if (trace > 1) {
    double tmp = 0;
    for (int j = 0; j < n; j++) {
      if (tmp < fabs(g[j])) tmp = fabs(g[j]);
    }
    say(pr, "max(abs(g))=%f", tmp);
  }
  say(pr, "\n");
}

// Driver-side stopping rules at a new iterate (itask == 1): maxit on the
// number of function evaluations and the x tolerance.  Returns the final
// itask code (28 or 27), or 0 to continue.
int newXStop(int maxit, int fncount, const double *lastx, const double *x,
             int n, double rtol, double atol, int trace, const Printer *pr) {
  if (maxit < fncount) {
    if (trace > 2)
      say(pr, "Exit becuase maximum number of function calls %d met.\n", maxit);
    return 28;
  }
  if (xConverged(lastx, x, n, rtol, atol)) {
    if (trace > 2) say(pr, "CONVERGENCE: Parameters differences below xtol.\n");
    return 27;
  }
  return 0;
}

} // namespace

namespace {

// Evaluate f and g at x for an FG request from setulb
void evalFG(int n, double *x, double *Fmin, double *g, optimfn fn,
            optimgr gr, void *ex, int *fncount, int *grcount, int trace,
            int iter, const Printer *pr) {
  if (trace >= 2) sayVec(pr, "computing f and g at prm=\n", x, n);
  Fmin[0] = fn(n, x, ex);
  fncount[0]++;
  gr(n, x, g, ex);
  grcount[0]++;
  if (trace > 0) traceEval(pr, trace, iter, *Fmin, g, n);
}

void fillInfo(InfoOut *info, int itask, int icsave, const int *lsave,
              const int *isave, const double *dsave) {
  if (info == nullptr) return;
  info->itask = itask;
  info->icsave = icsave;
  std::memcpy(info->lsave, lsave, sizeof(info->lsave));
  std::memcpy(info->isave, isave, sizeof(info->isave));
  std::memcpy(info->dsave, dsave, sizeof(info->dsave));
}

} // namespace

// The lbfgsb3C_ driver loop for either implementation of setulb, with
// output routed to `pr` and the final state to `info`.
void lbfgsbDriver(SetulbStep step, int n, int lmm, double *x, double *lower,
                  double *upper, int *nbd, double *Fmin, optimfn fn,
                  optimgr gr, int *fail, void *ex, double factr,
                  double pgtol, int *fncount, int *grcount, int maxit,
                  int trace, int iprint, double atol, double rtol,
                  double *g, const Printer *pr, InfoOut *info) {
  fncount[0] = 0;
  grcount[0] = 0;
  if (rejectInput(n, lmm, fail, pr, info)) return;
  WorkspaceLease lease;
  Workspace &ws = lease.get();
  ws.prepare(n, lmm);
  std::vector<double> &lastx = ws.lastx;
  std::copy(x, x + n, lastx.begin());
  int itask = 2, icsave = 0, itask2 = 0;
  int lsave[4] = {0};
  int isave[44] = {0};
  double dsave[29] = {0};
  bool doExit = false;
  while (!doExit) {
    if (trace >= 2) traceBeforeStep(pr, itask, x, n);
    // a STOP request gets one more call to restore the best point
    if (itask == 3) doExit = true;
    step(n, lmm, x, lower, upper, nbd, *Fmin, g, factr, pgtol, ws.wa.data(),
         ws.iwa.data(), itask, iprint, icsave, lsave, isave, dsave, pr);
    if (trace > 2) {
      say(pr, "returned from lbfgsb3 \n");
      say(pr, "returned itask is %d or \"%s\"\n", itask, taskName(itask));
    }
    if (itask == 4 || itask == 20 || itask == 21) {
      // FG, FG_LNSRCH, FG_START
      evalFG(n, x, Fmin, g, fn, gr, ex, fncount, grcount, trace, isave[33],
             pr);
    } else if (itask == 1) {
      // NEW_X
      itask2 = newXStop(maxit, fncount[0], lastx.data(), x, n, rtol, atol,
                        trace, pr);
      if (itask2 != 0) {
        itask = 3;
        doExit = true;
      }
      std::copy(x, x + n, lastx.begin());
    } else {
      doExit = true;
    }
  }
  if (itask2) itask = itask2;
  fillInfo(info, itask, icsave, lsave, isave, dsave);
  fail[0] = itask;
}

void lbfgsb3Cts_core(int n, int lmm, double *x, double *lower,
                     double *upper, int *nbd, double *Fmin, optimfn fn,
                     optimgr gr, int *fail, void *ex, double factr,
                     double pgtol, int *fncount, int *grcount,
                     int maxit, int trace, int iprint,
                     double atol, double rtol, double *g,
                     const Printer *pr, InfoOut *info) {
  lbfgsbDriver(setulb, n, lmm, x, lower, upper, nbd, Fmin, fn, gr, fail, ex,
               factr, pgtol, fncount, grcount, maxit, trace, iprint, atol,
               rtol, g, pr, info);
}

} // namespace lbfgsb3c_cpp

extern "C" void lbfgsb3Cts_(int n, int lmm, double *x, double *lower,
                            double *upper, int *nbd, double *Fmin,
                            optimfn fn, optimgr gr, int *fail, void *ex,
                            double factr, double pgtol, int *fncount,
                            int *grcount, int maxit, char *msg, int trace,
                            int iprint, double atol, double rtol,
                            double *g) {
  (void)msg;
  (void)trace;
  (void)iprint;
  lbfgsb3c_cpp::lbfgsb3Cts_core(n, lmm, x, lower, upper, nbd, Fmin, fn, gr,
                                fail, ex, factr, pgtol, fncount, grcount,
                                maxit, 0, -1, atol, rtol, g, nullptr,
                                nullptr);
}
