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
  for (int i = 1; i <= n; ++i) {
    if (nbd(i) > 0) {
      if (nbd(i) <= 2 && x(i) <= l(i)) {
        if (x(i) < l(i)) {
          prjctd = 1;
          x(i) = l(i);
        }
        nbdd = nbdd + 1;
      } else if (nbd(i) >= 2 && x(i) >= u(i)) {
        if (x(i) > u(i)) {
          prjctd = 1;
          x(i) = u(i);
        }
        nbdd = nbdd + 1;
      }
    }
  }
  for (int i = 1; i <= n; ++i) {
    if (nbd(i) != 2) boxed = 0;
    if (nbd(i) == 0) {
      iwhere(i) = -1;
    } else {
      cnstnd = 1;
      if (nbd(i) == 2 && u(i) - l(i) <= 0.0) {
        iwhere(i) = 3;
      } else {
        iwhere(i) = 0;
      }
    }
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

// f: lbfgsb.f:1315
void cauchy(int n, double *x_, const double *l_, const double *u_,
            const int *nbd_, double *g_, int *iorder_, int *iwhere_,
            double *t_, double *d_, double *xcp_, int m, double *wy_,
            double *ws_, double *sy_, double *wt_, double theta, int col,
            int head, double *p_, double *c_, double *wbp_, double *v_,
            int &nseg, int iprint, double sbgnrm, int &info, double epsmch,
            const Printer *pr) {
  V1<double> x(x_), g(g_), t(t_), d(d_), xcp(xcp_), p(p_), wbp(wbp_);
  V1<const double> l(l_), u(u_);
  V1<const int> nbd(nbd_);
  V1<int> iorder(iorder_), iwhere(iwhere_);
  V2<double> wy(wy_, n), ws(ws_, n);
  bool xlower, xupper, bnded;
  int i, j, col2, nfree, nbreak, pointr, ibp = 0, nleft, ibkmin, iter;
  double f1, f2, dt, dtm, tsum, dibp, zibp, dibp2, bkmin, tu = 0.0, tl = 0.0,
      wmc, wmp, wmw, tj, tj0, neggi, f2_org;

  if (sbgnrm <= 0.0) {
    if (iprint >= 0) say(pr, "Subgnorm =0, GCP = X.\n");
    dcopy(n, x_, xcp_);
    return;
  }
  bnded = true;
  nfree = n + 1;
  nbreak = 0;
  ibkmin = 0;
  bkmin = 0.0;
  col2 = 2 * col;
  f1 = 0.0;
  if (iprint >= 99) say(pr, "--- CAUCHY entered---\n");

  for (i = 1; i <= col2; ++i) p(i) = 0.0;

  for (i = 1; i <= n; ++i) {
    neggi = -g(i);
    if (iwhere(i) != 3 && iwhere(i) != -1) {
      if (nbd(i) <= 2) tl = x(i) - l(i);
      if (nbd(i) >= 2) tu = u(i) - x(i);
      xlower = nbd(i) <= 2 && tl <= 0.0;
      xupper = nbd(i) >= 2 && tu <= 0.0;
      iwhere(i) = 0;
      if (xlower) {
        if (neggi <= 0.0) iwhere(i) = 1;
      } else if (xupper) {
        if (neggi >= 0.0) iwhere(i) = 2;
      } else {
        if (fabs(neggi) <= 0.0) iwhere(i) = -3;
      }
    }
    pointr = head;
    if (iwhere(i) != 0 && iwhere(i) != -1) {
      d(i) = 0.0;
    } else {
      d(i) = neggi;
      f1 = f1 - neggi * neggi;
      for (j = 1; j <= col; ++j) {
        p(j) = p(j) + wy(i, pointr) * neggi;
        p(col + j) = p(col + j) + ws(i, pointr) * neggi;
        pointr = pointr % m + 1;
      }
      if (nbd(i) <= 2 && nbd(i) != 0 && neggi < 0.0) {
        nbreak = nbreak + 1;
        iorder(nbreak) = i;
        t(nbreak) = tl / (-neggi);
        if (nbreak == 1 || t(nbreak) < bkmin) {
          bkmin = t(nbreak);
          ibkmin = nbreak;
        }
      } else if (nbd(i) >= 2 && neggi > 0.0) {
        nbreak = nbreak + 1;
        iorder(nbreak) = i;
        t(nbreak) = tu / neggi;
        if (nbreak == 1 || t(nbreak) < bkmin) {
          bkmin = t(nbreak);
          ibkmin = nbreak;
        }
      } else {
        nfree = nfree - 1;
        iorder(nfree) = i;
        if (fabs(neggi) > 0.0) bnded = false;
      }
    }
  }

  if (theta != 1.0) dscal(col, theta, &p(col + 1));

  dcopy(n, x_, xcp_);

  if (nbreak == 0 && nfree == n + 1) {
    if (iprint > 100) sayVec(pr, "Cauchy X[1:5] = ", xcp_, n > 5 ? 5 : n);
    return;
  }

  for (j = 1; j <= col2; ++j) c_[j - 1] = 0.0;

  f2 = -theta * f1;
  f2_org = f2;
  if (col > 0) {
    bmv(m, sy_, wt_, col, p_, v_, info);
    if (info != 0) return;
    f2 = f2 - ddot(col2, v_, p_);
  }
  dtm = -f1 / f2;
  tsum = 0.0;
  nseg = 1;
  if (iprint >= 99) say(pr, "no. of breakpoints = %d\n", nbreak);

  if (nbreak == 0) goto L888;

  nleft = nbreak;
  iter = 1;
  tj = 0.0;

L777:
  tj0 = tj;
  if (iter == 1) {
    tj = bkmin;
    ibp = iorder(ibkmin);
  } else {
    if (iter == 2) {
      if (ibkmin != nbreak) {
        t(ibkmin) = t(nbreak);
        iorder(ibkmin) = iorder(nbreak);
      }
    }
    hpsolb(nleft, t_, iorder_, iter - 2);
    tj = t(nleft);
    ibp = iorder(nleft);
  }

  dt = tj - tj0;

  if (dt != 0.0 && iprint >= 100) {
    say(pr, "Piece %d\n", nseg);
    say(pr, "f1 at start point = %g\n", f1);
    say(pr, "f2 at start point = %g\n", f1);
    say(pr, "Distance to the next break point =  %g\n", dt);
    say(pr, "Distance to the stationary point =  %g\n", dtm);
  }

  if (dtm < dt) goto L888;

  tsum = tsum + dt;
  nleft = nleft - 1;
  iter = iter + 1;
  dibp = d(ibp);
  d(ibp) = 0.0;
  if (dibp > 0.0) {
    zibp = u(ibp) - x(ibp);
    xcp(ibp) = u(ibp);
    iwhere(ibp) = 2;
  } else {
    zibp = l(ibp) - x(ibp);
    xcp(ibp) = l(ibp);
    iwhere(ibp) = 1;
  }
  if (iprint >= 100) say(pr, "Variable fixed, index  %d\n", ibp);
  if (nleft == 0 && nbreak == n) {
    dtm = dt;
    goto L999;
  }

  nseg = nseg + 1;
  dibp2 = dibp * dibp;

  f1 = f1 + dt * f2 + dibp2 - theta * dibp * zibp;
  f2 = f2 - theta * dibp2;

  if (col > 0) {
    daxpy(col2, dt, p_, c_);
    pointr = head;
    for (j = 1; j <= col; ++j) {
      wbp(j) = wy(ibp, pointr);
      wbp(col + j) = theta * ws(ibp, pointr);
      pointr = pointr % m + 1;
    }
    bmv(m, sy_, wt_, col, wbp_, v_, info);
    if (info != 0) return;
    wmc = ddot(col2, c_, v_);
    wmp = ddot(col2, p_, v_);
    wmw = ddot(col2, wbp_, v_);
    daxpy(col2, -dibp, wbp_, p_);
    f1 = f1 + dibp * wmc;
    f2 = f2 + 2.0 * dibp * wmp - dibp2 * wmw;
  }

  f2 = dmax(epsmch * f2_org, f2);
  if (nleft > 0) {
    dtm = -f1 / f2;
    goto L777;
  } else if (bnded) {
    f1 = 0.0;
    f2 = 0.0;
    dtm = 0.0;
  } else {
    dtm = -f1 / f2;
  }

L888:
  if (iprint >= 99) {
    say(pr, "Piece %d\n", nseg);
    say(pr, "f1 at start point = %g\n", f1);
    say(pr, "f2 at start point = %g\n", f1);
    say(pr, "Distance to the stationary point =  %g\n", dtm);
  }
  if (dtm <= 0.0) dtm = 0.0;
  tsum = tsum + dtm;
  daxpy(n, tsum, d_, xcp_);

L999:
  if (col > 0) daxpy(col2, dtm, p_, c_);
  if (iprint > 100) sayVec(pr, "Cauchy X[1:5] = ", xcp_, n > 5 ? 5 : n);
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

// f: lbfgsb.f:1991
void formk(int n, int nsub, int *ind_, int nenter, int ileave, int *indx2_,
           int iupdat, int updatd, double *wn_, double *wn1_, int m,
           double *ws_, double *wy_, double *sy_, double theta, int col,
           int head, int &info) {
  V1<int> ind(ind_), indx2(indx2_);
  V2<double> wn(wn_, 2 * m), wn1(wn1_, 2 * m), ws(ws_, n), wy(wy_, n),
      sy(sy_, m);
  int m2, ipntr, jpntr, iy, is, jy, js, is1, js1, k1, i, k, col2, pbegin,
      pend, dbegin, dend, upcl;
  double temp1, temp2, temp3, temp4;

  if (updatd == 1) {
    if (iupdat > m) {
      for (jy = 1; jy <= m - 1; ++jy) {
        js = m + jy;
        dcopy(m - jy, &wn1(jy + 1, jy + 1), &wn1(jy, jy));
        dcopy(m - jy, &wn1(js + 1, js + 1), &wn1(js, js));
        dcopy(m - 1, &wn1(m + 2, jy + 1), &wn1(m + 1, jy));
      }
    }
    pbegin = 1;
    pend = nsub;
    dbegin = nsub + 1;
    dend = n;
    iy = col;
    is = m + col;
    ipntr = head + col - 1;
    if (ipntr > m) ipntr = ipntr - m;
    jpntr = head;
    for (jy = 1; jy <= col; ++jy) {
      js = m + jy;
      temp1 = 0.0;
      temp2 = 0.0;
      temp3 = 0.0;
      for (k = pbegin; k <= pend; ++k) {
        k1 = ind(k);
        temp1 = temp1 + wy(k1, ipntr) * wy(k1, jpntr);
      }
      for (k = dbegin; k <= dend; ++k) {
        k1 = ind(k);
        temp2 = temp2 + ws(k1, ipntr) * ws(k1, jpntr);
        temp3 = temp3 + ws(k1, ipntr) * wy(k1, jpntr);
      }
      wn1(iy, jy) = temp1;
      wn1(is, js) = temp2;
      wn1(is, jy) = temp3;
      jpntr = jpntr % m + 1;
    }
    jy = col;
    jpntr = head + col - 1;
    if (jpntr > m) jpntr = jpntr - m;
    ipntr = head;
    for (i = 1; i <= col; ++i) {
      is = m + i;
      temp3 = 0.0;
      for (k = pbegin; k <= pend; ++k) {
        k1 = ind(k);
        temp3 = temp3 + ws(k1, ipntr) * wy(k1, jpntr);
      }
      ipntr = ipntr % m + 1;
      wn1(is, jy) = temp3;
    }
    upcl = col - 1;
  } else {
    upcl = col;
  }

  ipntr = head;
  for (iy = 1; iy <= upcl; ++iy) {
    is = m + iy;
    jpntr = head;
    for (jy = 1; jy <= iy; ++jy) {
      js = m + jy;
      temp1 = 0.0;
      temp2 = 0.0;
      temp3 = 0.0;
      temp4 = 0.0;
      for (k = 1; k <= nenter; ++k) {
        k1 = indx2(k);
        temp1 = temp1 + wy(k1, ipntr) * wy(k1, jpntr);
        temp2 = temp2 + ws(k1, ipntr) * ws(k1, jpntr);
      }
      for (k = ileave; k <= n; ++k) {
        k1 = indx2(k);
        temp3 = temp3 + wy(k1, ipntr) * wy(k1, jpntr);
        temp4 = temp4 + ws(k1, ipntr) * ws(k1, jpntr);
      }
      wn1(iy, jy) = wn1(iy, jy) + temp1 - temp3;
      wn1(is, js) = wn1(is, js) - temp2 + temp4;
      jpntr = jpntr % m + 1;
    }
    ipntr = ipntr % m + 1;
  }

  ipntr = head;
  for (is = m + 1; is <= m + upcl; ++is) {
    jpntr = head;
    for (jy = 1; jy <= upcl; ++jy) {
      temp1 = 0.0;
      temp3 = 0.0;
      for (k = 1; k <= nenter; ++k) {
        k1 = indx2(k);
        temp1 = temp1 + ws(k1, ipntr) * wy(k1, jpntr);
      }
      for (k = ileave; k <= n; ++k) {
        k1 = indx2(k);
        temp3 = temp3 + ws(k1, ipntr) * wy(k1, jpntr);
      }
      if (is <= jy + m) {
        wn1(is, jy) = wn1(is, jy) + temp1 - temp3;
      } else {
        wn1(is, jy) = wn1(is, jy) - temp1 + temp3;
      }
      jpntr = jpntr % m + 1;
    }
    ipntr = ipntr % m + 1;
  }

  m2 = 2 * m;
  for (iy = 1; iy <= col; ++iy) {
    is = col + iy;
    is1 = m + iy;
    for (jy = 1; jy <= iy; ++jy) {
      js = col + jy;
      js1 = m + jy;
      wn(jy, iy) = wn1(iy, jy) / theta;
      wn(js, is) = wn1(is1, js1) * theta;
    }
    for (jy = 1; jy <= iy - 1; ++jy) wn(jy, is) = -wn1(is1, jy);
    for (jy = iy; jy <= col; ++jy) wn(jy, is) = wn1(is1, jy);
    wn(iy, iy) = wn(iy, iy) + sy(iy, iy);
  }

  dpofa(wn_, m2, col, info);
  if (info != 0) {
    info = -1;
    return;
  }
  col2 = 2 * col;
  for (js = col + 1; js <= col2; ++js) dtrsl(wn_, m2, col, &wn(1, js), 11, info);

  for (is = col + 1; is <= col2; ++is)
    for (js = is; js <= col2; ++js)
      wn(is, js) = wn(is, js) + ddot(col, &wn(1, is), &wn(1, js));

  dpofa(&wn(col + 1, col + 1), m2, col, info);
  if (info != 0) {
    info = -2;
    return;
  }
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

// f: lbfgsb.f:2384
void freev(int n, int &nfree, int *index_, int &nenter, int &ileave,
           int *indx2_, int *iwhere_, int &wrk, int updatd, int cnstnd,
           int iprint, int iter, const Printer *pr) {
  V1<int> index(index_), indx2(indx2_), iwhere(iwhere_);
  int iact, i, k;
  nenter = 0;
  ileave = n + 1;
  if (iter > 0 && cnstnd == 1) {
    for (i = 1; i <= nfree; ++i) {
      k = index(i);
      if (iwhere(k) > 0) {
        ileave = ileave - 1;
        indx2(ileave) = k;
        if (iprint >= 100)
          say(pr, "Variable k leaves the set of free variables for k = %d\n", k);
      }
    }
    for (i = 1 + nfree; i <= n; ++i) {
      k = index(i);
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
  if ((ileave < n + 1) || (nenter > 0) || updatd == 1) {
    wrk = 1;
  } else {
    wrk = 0;
  }
  nfree = 0;
  iact = n + 1;
  for (i = 1; i <= n; ++i) {
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

// f: lbfgsb.f:3622
void dcstep(double &stx, double &fx, double &dx, double &sty, double &fy,
            double &dy, double &stp, double fp, double dp, bool &brackt,
            double stpmin, double stpmax) {
  const double p66 = 0.66, two = 2.0, three = 3.0;
  double gamma, p, q, r, s, sgnd, stpc, stpf, stpq, theta, a, b;

  sgnd = dp * (dx / fabs(dx));

  if (fp > fx) {
    theta = three * (fx - fp) / (stp - stx) + dx + dp;
    s = dmax3(fabs(theta), fabs(dx), fabs(dp));
    a = theta / s;
    gamma = s * sqrt(a * a - (dx / s) * (dp / s));
    if (stp < stx) gamma = -gamma;
    p = (gamma - dx) + theta;
    q = ((gamma - dx) + gamma) + dp;
    r = p / q;
    stpc = stx + r * (stp - stx);
    stpq = stx + ((dx / ((fx - fp) / (stp - stx) + dx)) / two) * (stp - stx);
    if (fabs(stpc - stx) < fabs(stpq - stx)) {
      stpf = stpc;
    } else {
      stpf = stpc + (stpq - stpc) / two;
    }
    brackt = true;
  } else if (sgnd < 0.0) {
    theta = three * (fx - fp) / (stp - stx) + dx + dp;
    s = dmax3(fabs(theta), fabs(dx), fabs(dp));
    a = theta / s;
    gamma = s * sqrt(a * a - (dx / s) * (dp / s));
    if (stp > stx) gamma = -gamma;
    p = (gamma - dp) + theta;
    q = ((gamma - dp) + gamma) + dx;
    r = p / q;
    stpc = stp + r * (stx - stp);
    stpq = stp + (dp / (dp - dx)) * (stx - stp);
    if (fabs(stpc - stp) > fabs(stpq - stp)) {
      stpf = stpc;
    } else {
      stpf = stpq;
    }
    brackt = true;
  } else if (fabs(dp) < fabs(dx)) {
    theta = three * (fx - fp) / (stp - stx) + dx + dp;
    s = dmax3(fabs(theta), fabs(dx), fabs(dp));
    a = theta / s;
    b = a * a - (dx / s) * (dp / s);
    gamma = s * sqrt(dmax(0.0, b));
    if (stp > stx) gamma = -gamma;
    p = (gamma - dp) + theta;
    q = (gamma + (dx - dp)) + gamma;
    r = p / q;
    if (r < 0.0 && gamma != 0.0) {
      stpc = stp + r * (stx - stp);
    } else if (stp > stx) {
      stpc = stpmax;
    } else {
      stpc = stpmin;
    }
    stpq = stp + (dp / (dp - dx)) * (stx - stp);
    if (brackt) {
      if (fabs(stpc - stp) < fabs(stpq - stp)) {
        stpf = stpc;
      } else {
        stpf = stpq;
      }
      if (stp > stx) {
        stpf = dmin(stp + p66 * (sty - stp), stpf);
      } else {
        stpf = dmax(stp + p66 * (sty - stp), stpf);
      }
    } else {
      if (fabs(stpc - stp) > fabs(stpq - stp)) {
        stpf = stpc;
      } else {
        stpf = stpq;
      }
      stpf = dmin(stpmax, stpf);
      stpf = dmax(stpmin, stpf);
    }
  } else {
    if (brackt) {
      theta = three * (fp - fy) / (sty - stp) + dy + dp;
      s = dmax3(fabs(theta), fabs(dy), fabs(dp));
      a = theta / s;
      gamma = s * sqrt(a * a - (dy / s) * (dp / s));
      if (stp > sty) gamma = -gamma;
      p = (gamma - dp) + theta;
      q = ((gamma - dp) + gamma) + dy;
      r = p / q;
      stpc = stp + r * (sty - stp);
      stpf = stpc;
    } else if (stp > stx) {
      stpf = stpmax;
    } else {
      stpf = stpmin;
    }
  }

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

// f: lbfgsb.f:3255
void dcsrch(double f, double g, double &stp, double ftol, double gtol,
            double xtol, double stpmin, double stpmax, int &itask,
            int *isave_, double *dsave_) {
  V1<int> isave(isave_);
  V1<double> dsave(dsave_);
  const double p5 = 0.5, p66 = 0.66, xtrapl = 1.1, xtrapu = 4.0;
  bool brackt;
  int stage;
  double finit, ftest, fm, fx, fxm, fy, fym, ginit, gtest, gm, gx, gxm, gy,
      gym, stx, sty, stmin, stmax, width, width1;

  if (itask == 2) {
    if (stp < stpmin) itask = 16;
    if (stp > stpmax) itask = 15;
    if (g >= 0.0) itask = 11;
    if (ftol < 0.0) itask = 9;
    if (gtol < 0.0) itask = 10;
    if (xtol < 0.0) itask = 19;
    if (stpmin < 0.0) itask = 18;
    if (stpmax < stpmin) itask = 17;
    if ((itask >= 9) && (itask <= 19)) return;

    brackt = false;
    stage = 1;
    finit = f;
    ginit = g;
    gtest = ftol * ginit;
    width = stpmax - stpmin;
    width1 = width / p5;
    stx = 0.0;
    fx = finit;
    gx = ginit;
    sty = 0.0;
    fy = finit;
    gy = ginit;
    stmin = 0.0;
    stmax = stp + xtrapu * stp;
    itask = 4;
    goto L1000;
  } else {
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

  ftest = finit + stp * gtest;
  if (stage == 1 && f <= ftest && g >= 0.0) stage = 2;

  if (brackt && (stp <= stmin || stp >= stmax)) itask = 23;
  if (brackt && stmax - stmin <= xtol * stmax) itask = 26;
  if (stp == stpmax && f <= ftest && g <= gtest) itask = 24;
  if (stp == stpmin && (f > ftest || g >= gtest)) itask = 25;
  if (f <= ftest && fabs(g) <= gtol * (-ginit)) itask = 6;

  if ((itask >= 23) || ((itask <= 8) && (itask >= 6))) goto L1000;

  if (stage == 1 && f <= fx && f > ftest) {
    fm = f - stp * gtest;
    fxm = fx - stx * gtest;
    fym = fy - sty * gtest;
    gm = g - gtest;
    gxm = gx - gtest;
    gym = gy - gtest;
    dcstep(stx, fxm, gxm, sty, fym, gym, stp, fm, gm, brackt, stmin, stmax);
    fx = fxm + stx * gtest;
    fy = fym + sty * gtest;
    gx = gxm + gtest;
    gy = gym + gtest;
  } else {
    dcstep(stx, fx, gx, sty, fy, gy, stp, f, g, brackt, stmin, stmax);
  }

  if (brackt) {
    if (fabs(sty - stx) >= p66 * width1) stp = stx + p5 * (sty - stx);
    width1 = width;
    width = fabs(sty - stx);
  }

  if (brackt) {
    stmin = dmin(stx, sty);
    stmax = dmax(stx, sty);
  } else {
    stmin = stp + xtrapl * (stp - stx);
    stmax = stp + xtrapu * (stp - stx);
  }

  stp = dmax(stp, stpmin);
  stp = dmin(stp, stpmax);

  if ((brackt && (stp <= stmin || stp >= stmax)) ||
      (brackt && stmax - stmin <= xtol * stmax))
    stp = stx;

  itask = 4;

L1000:
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
  double a1, a2;

  if (itask == 20) goto L556;

  dtd = ddot(n, d_, d_);
  dnorm = sqrt(dtd);

  stpmx = big;
  if (cnstnd == 1) {
    if (iter == 0) {
      stpmx = one;
    } else {
      for (int i = 1; i <= n; ++i) {
        a1 = d(i);
        if (nbd(i) != 0) {
          if (a1 < 0.0 && nbd(i) <= 2) {
            a2 = l(i) - x(i);
            if (a2 >= 0.0) {
              stpmx = 0.0;
            } else if (a1 * stpmx < a2) {
              stpmx = a2 / a1;
            }
          } else if (a1 > 0.0 && nbd(i) >= 2) {
            a2 = u(i) - x(i);
            if (a2 <= 0.0) {
              stpmx = 0.0;
            } else if (a1 * stpmx > a2) {
              stpmx = a2 / a1;
            }
          }
        }
      }
    }
  }

  if (iter == 0 && boxed == 0) {
    stp = dmin(one / dnorm, stpmx);
  } else {
    stp = one;
  }

  dcopy(n, x_, t_);
  dcopy(n, g_, r_);
  fold = f;
  ifun = 0;
  iback = 0;
  icsave = 2;

L556:
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
  V1<double> x(x_), d(d_), xx(xx_), gg(gg_), wv(wv_);
  V2<double> ws(ws_, n), wy(wy_, n);
  int pointr, m2, col2, ibd, js, i, k;
  double alpha, xk, dk, temp1, temp2, dd_p;

  if (nsub <= 0) return;
  if (iprint >= 99) say(pr, " ----- SUBSM entered -----\n");

  pointr = head;
  for (i = 1; i <= col; ++i) {
    temp1 = 0.0;
    temp2 = 0.0;
    for (int j = 1; j <= nsub; ++j) {
      k = ind(j);
      temp1 = temp1 + wy(k, pointr) * d(j);
      temp2 = temp2 + ws(k, pointr) * d(j);
    }
    wv(i) = temp1;
    wv(col + i) = theta * temp2;
    pointr = pointr % m + 1;
  }

  m2 = 2 * m;
  col2 = 2 * col;
  dtrsl(wn_, m2, col2, wv_, 11, info);
  if (info != 0) return;
  for (i = 1; i <= col; ++i) wv(i) = -wv(i);
  dtrsl(wn_, m2, col2, wv_, 1, info);
  if (info != 0) return;

  pointr = head;
  for (int jy = 1; jy <= col; ++jy) {
    js = col + jy;
    for (i = 1; i <= nsub; ++i) {
      k = ind(i);
      d(i) = d(i) + wy(k, pointr) * wv(jy) / theta + ws(k, pointr) * wv(js);
    }
    pointr = pointr % m + 1;
  }

  dscal(nsub, 1.0 / theta, d_);

  iword = 0;
  dcopy(n, x_, xp_);

  for (i = 1; i <= nsub; ++i) {
    k = ind(i);
    dk = d(i);
    xk = x(k);
    if (nbd(k) != 0) {
      if (nbd(k) == 1) {
        x(k) = dmax(l(k), xk + dk);
        if (x(k) == l(k)) iword = 1;
      } else {
        if (nbd(k) == 2) {
          xk = dmax(l(k), xk + dk);
          x(k) = dmin(u(k), xk);
          if (x(k) == l(k) || x(k) == u(k)) iword = 1;
        } else {
          if (nbd(k) == 3) {
            x(k) = dmin(u(k), xk + dk);
            if (x(k) == u(k)) iword = 1;
          }
        }
      }
    } else {
      x(k) = xk + dk;
    }
  }

  if (iword == 0) goto L911;

  dd_p = 0.0;
  for (i = 1; i <= n; ++i) dd_p = dd_p + (x(i) - xx(i)) * gg(i);
  if (dd_p > 0.0) {
    dcopy(n, xp_, x_);
    say(pr, " Positive dir derivative in projection \n");
    say(pr, " Using the backtracking step \n");
  } else {
    goto L911;
  }

  alpha = 1.0;
  temp1 = alpha;
  ibd = 0;
  for (i = 1; i <= nsub; ++i) {
    k = ind(i);
    dk = d(i);
    if (nbd(k) != 0) {
      if (dk < 0.0 && nbd(k) <= 2) {
        temp2 = l(k) - x(k);
        if (temp2 >= 0.0) {
          temp1 = 0.0;
        } else if (dk * alpha < temp2) {
          temp1 = temp2 / dk;
        }
      } else if (dk > 0.0 && nbd(k) >= 2) {
        temp2 = u(k) - x(k);
        if (temp2 <= 0.0) {
          temp1 = 0.0;
        } else if (dk * alpha > temp2) {
          temp1 = temp2 / dk;
        }
      }
      if (temp1 < alpha) {
        alpha = temp1;
        ibd = i;
      }
    }
  }

  if (alpha < 1.0) {
    dk = d(ibd);
    k = ind(ibd);
    if (dk > 0.0) {
      x(k) = u(k);
      d(ibd) = 0.0;
    } else if (dk < 0.0) {
      x(k) = l(k);
      d(ibd) = 0.0;
    }
  }
  for (i = 1; i <= nsub; ++i) {
    k = ind(i);
    x(k) = x(k) + alpha * d(i);
  }

L911:
  if (iprint >= 99) say(pr, " exit SUBSM \n");
}

// f: lbfgsb.f:305
void mainlb(int n, int m, double *x, const double *l, const double *u,
            const int *nbd, double &f, double *g, double factr, double pgtol,
            double *ws, double *wy, double *sy, double *ss, double *wt,
            double *wn, double *snd, double *z_, double *r_, double *d_,
            double *t, double *xp, double *wa, int *index, int *iwhere,
            int *indx2, int &itask, int iprint, int &icsave, int *lsave_,
            int *isave_, double *dsave_, const Printer *pr) {
  V1<int> lsave(lsave_), isave(isave_);
  V1<double> dsave(dsave_), z(z_), r(r_), d(d_), gv(g), xv(x);
  const double one = 1.0, zero = 0.0;
  int prjctd = 0, cnstnd = 0, boxed = 0, updatd = 0, wrk = 0;
  int k = 0, nintol = 0, iback = 0, nskip = 0, head = 0, col = 0, iter = 0,
      itail = 0, iupdat = 0, nseg = 0, nfgv = 0, info = 0, ifun = 0,
      iword = 0, nfree = 0, nact = 0, ileave = 0, nenter = 0;
  double theta = 0, fold = 0, dr = 0, rr = 0, tol = 0, xstep = 0,
         sbgnrm = 0, ddum = 0, dnorm = 0, dtd = 0, epsmch = 0, cpu1 = 0,
         cpu2 = 0, sbtime = 0, time1 = 0, gd = 0, gdold = 0, stp = 0,
         stpmx = 0;

  if (itask == 2) {
    epsmch = DBL_EPSILON;
    time1 = 0.0;
    col = 0;
    head = 1;
    theta = one;
    iupdat = 0;
    updatd = 0;
    iback = 0;
    itail = 0;
    iword = 0;
    nact = 0;
    ileave = 0;
    nenter = 0;
    fold = zero;
    dnorm = zero;
    cpu1 = zero;
    gd = zero;
    stpmx = zero;
    sbgnrm = zero;
    stp = zero;
    gdold = zero;
    dtd = zero;
    iter = 0;
    nfgv = 0;
    nseg = 0;
    nintol = 0;
    nskip = 0;
    nfree = n;
    ifun = 0;
    tol = factr * epsmch;
    sbtime = 0;
    info = 0;

    errclb(n, m, factr, l, u, nbd, itask, info, k, pr);
    if ((itask >= 9) && (itask <= 19)) return;

    active(n, l, u, nbd, x, iwhere, iprint, prjctd, cnstnd, boxed, pr);
  } else {
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

    if (itask == 20) goto L666;
    if (itask == 1) goto L777;
    if (itask == 21) goto L111;
  }

  // Compute f0 and g0.
  itask = 21;
  goto L1000;

L111:
  nfgv = 1;
  projgr(n, l, u, nbd, x, g, sbgnrm);
  if (iprint >= 1) say(pr, "At iterate %d f= %g |proj g|=  %g\n", iter, f, sbgnrm);
  if (sbgnrm <= pgtol) {
    itask = 7;
    goto L999;
  }

L222:
  if (iprint >= 99) say(pr, "ITERATION  %d\n", iter + 1);
  iword = -1;
  if (cnstnd == 0 && col > 0) {
    dcopy(n, x, z_);
    wrk = updatd;
    nseg = 0;
    goto L333;
  }

  cpu1 = 0.0;
  cauchy(n, x, l, u, nbd, g, indx2, iwhere, t, d_, z_, m, wy, ws, sy, wt,
         theta, col, head, &wa[0], &wa[2 * m], &wa[4 * m], &wa[6 * m], nseg,
         iprint, sbgnrm, info, epsmch, pr);
  if (info != 0) {
    sayRestart(pr, " Singular triangular system detected;");
    info = 0;
    col = 0;
    head = 1;
    theta = one;
    iupdat = 0;
    updatd = 0;
    cpu2 = 0.0;
    goto L222;
  }
  cpu2 = 0.0;
  nintol = nintol + nseg;

  freev(n, nfree, index, nenter, ileave, indx2, iwhere, wrk, updatd, cnstnd,
        iprint, iter, pr);
  nact = n - nfree;

L333:
  if (nfree == 0 || col == 0) goto L555;

  cpu1 = 0.0;
  if (wrk == 1)
    formk(n, nfree, index, nenter, ileave, indx2, iupdat, updatd, wn, snd, m,
          ws, wy, sy, theta, col, head, info);
  if (info != 0) {
    if (iprint >= 1)
      sayRestart(pr, " Nonpositive definiteness in Cholesky factorization in formk;");
    info = 0;
    col = 0;
    head = 1;
    theta = one;
    iupdat = 0;
    updatd = 0;
    cpu2 = 0.0;
    sbtime = sbtime + cpu2 - cpu1;
    goto L222;
  }

  cmprlb(n, m, x, g, ws, wy, sy, wt, z_, r_, wa, index, theta, col, head,
         nfree, cnstnd, info);
  if (info != 0) goto L444;

  subsm(n, m, nfree, index, l, u, nbd, z_, r_, xp, ws, wy, theta, x, g, col,
        head, iword, wa, wn, iprint, info, pr);

L444:
  if (info != 0) {
    sayRestart(pr, " Singular triangular system detected;");
    info = 0;
    col = 0;
    head = 1;
    theta = one;
    iupdat = 0;
    updatd = 0;
    cpu2 = 0.0;
    sbtime = sbtime + cpu2 - cpu1;
    goto L222;
  }
  cpu2 = 0.0;
  sbtime = sbtime + cpu2 - cpu1;

L555:
  for (int i = 1; i <= n; ++i) d(i) = z(i) - xv(i);
  cpu1 = 0.0;

L666:
  lnsrlb(n, l, u, nbd, x, f, fold, gd, gdold, g, d_, r_, t, z_, stp, dnorm,
         dtd, xstep, stpmx, iter, ifun, iback, nfgv, info, itask, boxed,
         cnstnd, icsave, &isave(22), &dsave(17), pr);
  if (info != 0 || iback >= 20) {
    dcopy(n, t, x);
    dcopy(n, r_, g);
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
      goto L999;
    } else {
      if (iprint >= 1) sayRestart(pr, " Bad direction in the line search;");
      if (info == 0) nfgv = nfgv - 1;
      info = 0;
      col = 0;
      head = 1;
      theta = one;
      iupdat = 0;
      updatd = 0;
      itask = 22;
      cpu2 = 0.0;
      goto L222;
    }
  } else if (itask == 20) {
    goto L1000;
  } else {
    cpu2 = 0.0;
    iter = iter + 1;
    projgr(n, l, u, nbd, x, g, sbgnrm);
    goto L1000;
  }

L777:
  if (sbgnrm <= pgtol) {
    itask = 7;
    goto L999;
  }
  ddum = dmax3(fabs(fold), fabs(f), one);
  if ((fold - f) <= tol * ddum) {
    itask = 8;
    if (iback >= 10) info = -5;
    goto L999;
  }

  for (int i = 1; i <= n; ++i) r(i) = gv(i) - r(i);
  rr = ddot(n, r_, r_);
  if (stp == one) {
    dr = gd - gdold;
    ddum = -gdold;
  } else {
    dr = (gd - gdold) * stp;
    dscal(n, stp, d_);
    ddum = -gdold * stp;
  }

  if (dr <= epsmch * ddum) {
    nskip = nskip + 1;
    updatd = 0;
    if (iprint >= 1) {
      say(pr, " ys = %g\n", dr);
      say(pr, " BFGS update skipped for +gs= %g\n", ddum);
    }
    goto L888;
  }

  updatd = 1;
  iupdat = iupdat + 1;
  matupd(n, m, ws, wy, sy, ss, d_, r_, itail, iupdat, col, head, theta, rr,
         dr, stp, dtd);
  formt(m, wt, sy, ss, col, theta, info);
  if (info != 0) {
    if (iprint >= 1)
      sayRestart(pr, " Nonpositive definiteness in Cholesky factorization in formt;");
    info = 0;
    col = 0;
    head = 1;
    theta = one;
    iupdat = 0;
    updatd = 0;
    goto L222;
  }

L888:
  goto L222;

L999:
  // time2 = 0; time = time2 - time1 (timing disabled, unused)

L1000:
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
  (void)xstep;
  (void)k;
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

// Port of the lbfgsb3C_ driver loop in lbfgsb3x.cpp, with output routed
// to `pr` and the final state to `info` instead of R objects.
void lbfgsb3Cts_core(int n, int lmm, double *x, double *lower,
                     double *upper, int *nbd, double *Fmin, optimfn fn,
                     optimgr gr, int *fail, void *ex, double factr,
                     double pgtol, int *fncount, int *grcount,
                     int maxit, int trace, int iprint,
                     double atol, double rtol, double *g,
                     const Printer *pr, InfoOut *info) {
  fncount[0] = 0;
  grcount[0] = 0;
  if (!validLmm(n, lmm)) {
    // lmm <= 0 would divide by zero in matupd; a huge lmm overflows the
    // workspace size.  Fail before allocating or evaluating anything.
    say(pr, "  ERROR: INVALID LMM (%d)\n", lmm);
    if (info != nullptr) {
      std::memset(info, 0, sizeof(*info));
      info->itask = kInvalidLmm;
    }
    fail[0] = kInvalidLmm;
    return;
  }
  int itask = 2;
  WorkspaceLease lease;
  Workspace &ws = lease.get();
  ws.prepare(n, lmm);
  std::vector<double> &wa = ws.wa;
  std::vector<int> &iwa = ws.iwa;
  std::vector<double> &lastx = ws.lastx;
  std::copy(x, x + n, lastx.begin());
  int icsave = 0;
  int lsave[4] = {0};
  int isave[44] = {0};
  double dsave[29] = {0};
  int doExit = 0;
  int itask2 = 0;
  fncount[0] = 0;
  grcount[0] = 0;
  while (true) {
    if (trace >= 2) {
      say(pr, "itask: %d\n", itask);
      sayVec(pr, "computing f and g at prm=\n", x, n);
      say(pr, "\n================================================================================\nBefore call task number %d, or \"%s\"\n",
          itask, taskName(itask));
    }
    if (itask == 3) doExit = 1;
    setulb(n, lmm, x, lower, upper, nbd, *Fmin, g, factr, pgtol, wa.data(),
           iwa.data(), itask, iprint, icsave, lsave, isave, dsave, pr);
    if (trace > 2) {
      say(pr, "returned from lbfgsb3 \n");
      say(pr, "returned itask is %d or \"%s\"\n", itask, taskName(itask));
    }
    switch (itask) {
    case 4:
    case 20:
    case 21:
      if (trace >= 2) sayVec(pr, "computing f and g at prm=\n", x, n);
      Fmin[0] = fn(n, x, ex);
      fncount[0]++;
      gr(n, x, g, ex);
      grcount[0]++;
      if (trace > 0) {
        say(pr, "At iteration %d f=%f ", isave[33], *Fmin);
        if (trace > 1) {
          double tmp = 0;
          for (int j = 0; j < n; j++) {
            if (tmp < fabs(g[j])) tmp = fabs(g[j]);
          }
          say(pr, "max(abs(g))=%f", tmp);
        }
        say(pr, "\n");
      }
      break;
    case 1:
      if (maxit < fncount[0]) {
        itask2 = 28;
        doExit = 1;
        itask = 3;
        if (trace > 2)
          say(pr, "Exit becuase maximum number of function calls %d met.\n", maxit);
      } else {
        bool converge = fabs(lastx[n - 1] - x[n - 1]) < fabs(x[n - 1]) * rtol + atol;
        if (converge) {
          // keeps lbfgsb3C_'s loop (visits n-2..0) for identical results
          for (int i = n - 1; i--;) {
            converge = fabs(lastx[i] - x[i]) < fabs(x[i]) * rtol + atol;
            if (!converge) break;
          }
        }
        if (converge) {
          itask2 = 27;
          itask = 3;
          if (trace > 2) say(pr, "CONVERGENCE: Parameters differences below xtol.\n");
          doExit = 1;
        }
      }
      std::copy(x, x + n, lastx.begin());
      break;
    default:
      doExit = 1;
    }
    if (doExit) break;
  }
  if (itask2) itask = itask2;
  if (info != nullptr) {
    info->itask = itask;
    info->icsave = icsave;
    std::memcpy(info->lsave, lsave, sizeof(lsave));
    std::memcpy(info->isave, isave, sizeof(isave));
    std::memcpy(info->dsave, dsave, sizeof(dsave));
  }
  fail[0] = itask;
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
