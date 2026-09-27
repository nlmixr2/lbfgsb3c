#include "../inst/include/lbfgsb3c.h"
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

void lbfgsb3C_(int n, int lmm, double *x, double *lower,
	       double *upper, int *nbd, double *Fmin, optimfn fn,
	       optimgr gr, int *fail, void *ex, double factr,
	       double pgtol, int *fncount, int *grcount,
	       int maxit, char *msg, int trace, int iprint,
	       double atol, double rtol, double *g);

void lbfgsb3Cts_(int n, int lmm, double *x, double *lower,
		 double *upper, int *nbd, double *Fmin, optimfn fn,
		 optimgr gr, int *fail, void *ex, double factr,
		 double pgtol, int *fncount, int *grcount,
		 int maxit, char *msg, int trace, int iprint,
		 double atol, double rtol, double *g);

SEXP _lbfgsb3c_lbfgsb3cpp(SEXP, SEXP, SEXP, SEXP, SEXP,
			  SEXP, SEXP);
SEXP _lbfgsb3c_lbfgsb3cThreadTest_(SEXP, SEXP);
SEXP _lbfgsb3c_lbfgsb3cLmmTest_(SEXP);

SEXP _lbfgsb3c_ptr(void) {
  int pro = 0;  // Counter for the number of PROTECT calls

  // Create an external pointer
  SEXP lbfgsb3Cptr = PROTECT(R_MakeExternalPtrFn((DL_FUNC)&lbfgsb3C_, R_NilValue, R_NilValue)); pro++;
  SEXP lbfgsb3Ctsptr = PROTECT(R_MakeExternalPtrFn((DL_FUNC)&lbfgsb3Cts_, R_NilValue, R_NilValue)); pro++;

  // Element order is part of the interface (see inst/include/lbfgsb3ptr.h);
  // only append new entries.
#define nVec 2

  SEXP ret = PROTECT(Rf_allocVector(VECSXP, nVec)); pro++;
  SEXP retN = PROTECT(Rf_allocVector(STRSXP, nVec)); pro++;

  SET_VECTOR_ELT(ret, 0, lbfgsb3Cptr);
  SET_STRING_ELT(retN, 0, Rf_mkChar("lbfgsb3C"));
  SET_VECTOR_ELT(ret, 1, lbfgsb3Ctsptr);
  SET_STRING_ELT(retN, 1, Rf_mkChar("lbfgsb3Cts"));

#undef nVec

  // Set the names attribute of the list
  Rf_setAttrib(ret, R_NamesSymbol, retN);

  // Unprotect all protected objects
  UNPROTECT(pro);

  // Return the list of external pointers
  return ret;

}

void R_init_lbfgsb3c(DllInfo *info){
  R_CallMethodDef callMethods[]  = {
    {"_lbfgsb3c_lbfgsb3cpp", (DL_FUNC) &_lbfgsb3c_lbfgsb3cpp, 7},
    {"_lbfgsb3c_ptr", (DL_FUNC) &_lbfgsb3c_ptr, 0},
    {"_lbfgsb3c_lbfgsb3cThreadTest_", (DL_FUNC) &_lbfgsb3c_lbfgsb3cThreadTest_, 2},
    {"_lbfgsb3c_lbfgsb3cLmmTest_", (DL_FUNC) &_lbfgsb3c_lbfgsb3cLmmTest_, 1},
    {NULL, NULL, 0}
  };
  // C callable to assign environments.
  R_RegisterCCallable("lbfgsb3c", "lbfgsb3C_", (DL_FUNC) lbfgsb3C_);
  R_RegisterCCallable("lbfgsb3c", "lbfgsb3Cts_", (DL_FUNC) lbfgsb3Cts_);
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}
