# lbfgsb3c 2024-3.6 changes

* Added a thread-safe C++ port of the L-BFGS-B 3.0 Fortran code (the
  Fortran is kept and remains the default).  It keeps no shared
  mutable state and does not call the R API, so independent problems can
  be solved concurrently (for example as an inner optimizer in an OpenMP
  loop).  The port reproduces the Fortran results exactly.  Its
  workspace is allocated once per thread and reused (zeroed per run), so
  memory scales with the number of threads, not the number of
  optimizations.

* The C++ port is available from C/C++ as `lbfgsb3Cts` (the second
  element of `.lbfgsb3cPtr()`, set up by `iniLbfgsb3` in
  `lbfgsb3ptr.h`) and from R with `lbfgsb3c(..., control =
  list(engine = "cpp"))`.

* `lbfgsb3C_` no longer stores its `$info` result in a global R list.

* `lmm` is now validated.  `lbfgsb3c()` errors unless it is a whole
  number >= 1; the C entry points return `fail = 29` ("ERROR: INVALID
  LMM") instead of dividing by zero (`lmm <= 0`) or overflowing the
  workspace size (very large `lmm`).  `n <= 0` returns `fail = 13`
  before any work arrays are allocated.

* Nested `lbfgsb3c()` calls (an objective that runs its own
  optimization) now work: the R callbacks are kept per call instead of
  in globals that the inner call overwrote.  A gradient of the wrong
  length is now an error instead of a read past its end.

* Added `src/Makevars` (OpenMP for the internal thread-safety test
  helper, and explicit `$(FLIBS)`).

# Version: 2024-3.5 changes

* Added function pointer interface (instead of only low level abi interface)

# lbfgsb3c 2024-3.4 changes

* LTO fixes and remove unused code for Fortran fixes

# lbfgsb3c 2024-3.3 changes

* Fixed `lmm`.  In prior version with the R interface `lmm` was not
  being passed through correctly (though it was passed through in C
  correctly)

* Fixed some bugs in printout

* Reverted the code to fix some of the issues in the FORTRAN code

# lbfgsb3c 2020-3.3 changes

* Now allow `rho=NULL` to work the same as if `rho` was not supplied

* Be More careful about `$convergence` by adding a default value of `NA_INTEGER`

* `$convergence` is now an integer instead of a real number

* Fix too many arguments for format as requested by CRAN

* Added a `NEWS.md` file to track changes to the package.

# lbfgsb3c 2020-3.2 prior changes and information

To do

Add test using a plain C function for optimization. lbfgsb3c is
supposed to handle this.

--------------------------------------------------------------
2019-03-19
    o Packages lbfgsb3 and lbfgsb3c merged into latter. Vignette added.
    o Suppressed printout when trace>2 and starting (f not defined)

2015-01-20
    o Fixup line longer than 72 chars in lbfgsb.f. Undeclared
      integer itask in errclb subroutine. Thanks to Berend Hasselman.

New package lbfgsb3 2014.7.31
