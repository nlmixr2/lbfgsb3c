# Changelog

## lbfgsb3c 2024-3.6 changes

- Added a thread-safe C++ port of the L-BFGS-B 3.0 Fortran code (the
  Fortran is kept and remains the default). It keeps no shared mutable
  state and does not call the R API, so independent problems can be
  solved concurrently (for example as an inner optimizer in an OpenMP
  loop). The port reproduces the Fortran results exactly. Its workspace
  is allocated once per thread and reused (zeroed per run), so memory
  scales with the number of threads, not the number of optimizations.

- The C++ port is available from C/C++ as `lbfgsb3Cts` (the second
  element of
  [`.lbfgsb3cPtr()`](https://nlmixr2.github.io/lbfgsb3c/reference/dot-lbfgsb3cPtr.md),
  set up by `iniLbfgsb3` in `lbfgsb3ptr.h`) and from R with
  `lbfgsb3c(..., control = list(engine = "cpp"))`.

- `lbfgsb3C_` no longer stores its `$info` result in a global R list.

- `lmm` is now validated.
  [`lbfgsb3c()`](https://nlmixr2.github.io/lbfgsb3c/reference/lbfgsb3c.md)
  errors unless it is a whole number \>= 1; the C entry points return
  `fail = 29` (“ERROR: INVALID LMM”) instead of dividing by zero
  (`lmm <= 0`) or overflowing the workspace size (very large `lmm`).
  `n <= 0` returns `fail = 13` before any work arrays are allocated.

- Nested
  [`lbfgsb3c()`](https://nlmixr2.github.io/lbfgsb3c/reference/lbfgsb3c.md)
  calls (an objective that runs its own optimization) now work: the R
  callbacks are kept per call instead of in globals that the inner call
  overwrote. A gradient of the wrong length is now an error instead of a
  read past its end.

- Added `src/Makevars` (OpenMP for the internal thread-safety test
  helper, and explicit `$(FLIBS)`).

- `control$reltol` is now checked for length one (the check tested
  `abstol` instead, so an empty `reltol` was read out of bounds).

## lbfgsb3c 2024-3.4 changes

CRAN release: 2024-04-04

- LTO fixes and remove unused code for Fortran fixes

## lbfgsb3c 2024-3.3 changes

CRAN release: 2024-04-01

- Fixed `lmm`. In prior version with the R interface `lmm` was not being
  passed through correctly (though it was passed through in C correctly)

- Fixed some bugs in printout

- Reverted the code to fix some of the issues in the FORTRAN code

## lbfgsb3c 2020-3.3 changes

CRAN release: 2023-11-28

- Now allow `rho=NULL` to work the same as if `rho` was not supplied

- Be More careful about `$convergence` by adding a default value of
  `NA_INTEGER`

- `$convergence` is now an integer instead of a real number

- Fix too many arguments for format as requested by CRAN

- Added a `NEWS.md` file to track changes to the package.

## lbfgsb3c 2020-3.2 prior changes and information

CRAN release: 2020-03-03

To do

Add test using a plain C function for optimization. lbfgsb3c is supposed
to handle this.

------------------------------------------------------------------------

2019-03-19 o Packages lbfgsb3 and lbfgsb3c merged into latter. Vignette
added. o Suppressed printout when trace\>2 and starting (f not defined)

2015-01-20 o Fixup line longer than 72 chars in lbfgsb.f. Undeclared
integer itask in errclb subroutine. Thanks to Berend Hasselman.

New package lbfgsb3 2014.7.31
