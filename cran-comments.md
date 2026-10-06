## Release summary

- Adds a thread-safe C++ port of the L-BFGS-B 3.0 code (opt-in; the
  Fortran remains the default), validates `lmm`, and fixes several
  memory-safety issues (see NEWS.md).

## Test environments

- Local: Ubuntu, R 4.6.1, `R CMD check --as-cran`
- R-hub v2: linux, windows, macos, macos-arm64, m1-san, clang-asan,
  clang-ubsan, clang16-23, gcc13-16, gcc-asan, atlas, c23, intel, mkl,
  lto, nold, noremap, donttest, ubuntu-* (R-devel and release)
- Local valgrind (memcheck) of tests and examples: no errors

## R CMD check results

0 errors | 0 warnings | 0 notes on R-hub.

Locally there is 1 NOTE about `-mno-omit-leaf-frame-pointer`, which comes
from the system R build flags, not the package.

## Reverse dependencies

We checked 14 reverse dependencies (revdepcheck), comparing R CMD check
results across CRAN and dev versions of this package.

- We saw 0 new problems
- We failed to check 0 packages
