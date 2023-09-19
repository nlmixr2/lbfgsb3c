<!-- badges: start -->
[![R-CMD-check](https://github.com/nlmixr2/lbfgsb3c/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/nlmixr2/lbfgsb3c/actions/workflows/R-CMD-check.yaml)
[![CRAN total downloads](https://cranlogs.r-pkg.org/badges/grand-total/lbfgsb3c)](https://cran.r-project.org/package=lbfgsb3c)
[![CRAN total downloads](https://cranlogs.r-pkg.org/badges/lbfgsb3c)](https://cran.r-project.org/package=lbfgsb3c)
[![Codecov test coverage](https://codecov.io/gh/nlmixr2/lbfgsb3c/branch/main/graph/badge.svg)](https://app.codecov.io/gh/nlmixr2/lbfgsb3c?branch=main)
<!-- badges: end -->

# libfgsb3c interface from C
This is the fork of the libfgsb3 from cran with the following differences:
- The return type has changed is is very similar to what `optim` returns
- Allows a direct C/C++ interface through a R registered function,
  similar to C interface to `optim` with 2 additional arguments.
- Allows adjustment of tolerances for minimization success.
- Added `xtolAtol` and `xtolRtol` minimization success criterion.
- Added `maxit` termination
