# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when
working with code in this repository.

## Overview

**lbfgsb3c** is an R package wrapping the Nocedal-Morales L-BFGS-B 3.0
limited-memory BFGS minimizer with bounds on parameters. It is a fork of
CRAN's `lbfgsb3` that returns `optim()`-like results, adds extra stopping
criteria (`abstol`/`reltol`, `maxit`), and exposes an `optim()`-style C
interface so other packages (nlmixr2 ecosystem) can call the optimizer
from compiled code. The package mixes R, Rcpp/C++, C and Fortran.

Version scheme is date-like (`2026-3.5`), not semantic versioning.

## Build and Development Commands

```r
devtools::install()          # or from shell: R CMD INSTALL .
devtools::document()         # roxygen2 -> man/, NAMESPACE
Rcpp::compileAttributes(".") # regenerate src/RcppExports.cpp, R/RcppExports.R
```

### Run Tests
```sh
rm -f src/*.o src/*.so; Rscript -e "devtools::test()"
Rscript -e "devtools::test(filter='cpp-engine')"   # single file
```

Exact (bit-for-bit) Fortran vs C++ comparison, only meaningful with R's
reference BLAS (see `.github/workflows/exact-engine.yaml`):
```sh
LBFGSB3C_EXACT_ENGINE_TEST=true Rscript -e "devtools::test(filter='cpp-engine')"
```

### R CMD Check
```sh
rm -f src/*.o src/*.so; R CMD build . && R CMD check --as-cran lbfgsb3c_*.tar.gz
```

## Architecture

### R layer (`R/lbfgsb3.R`)
- `lbfgsb3c(par, fn, gr, lower, upper, control, ..., rho)` is the only real
  entry point; `lbfgsb3`, `lbfgsb3f`, `lbfgsb3x` are aliases.  Calling it as
  `lbfgsb3`/`lbfgsb3f` turns on `control$info` (detected from the call stack).
- Validates `control` names and `lmm` (whole number >= 1), maps
  `control$engine` (`"fortran"` default, `"cpp"`) to 0/1, builds a
  `numDeriv::grad` gradient when `gr` is NULL, then calls the Rcpp
  `lbfgsb3cpp()` (`src/lbfgsb3x.cpp`).
- `.lbfgsb3cPtr()` returns the exported function-pointer list.

### Compiled code (`src/`)
- `lbfgsb.f`, `blas.f` -- original L-BFGS-B 3.0 Fortran (`setulb_`),
  modified to print through R. Default engine; NOT thread safe.
- `lbfgsb_cpp.cpp`/`.h` -- thread-safe C++ port (namespace `lbfgsb3c_cpp`,
  entry `lbfgsb3Cts_`). No static/global mutable state and no R API calls;
  workspace allocated once per thread and reused. Must reproduce the
  Fortran results exactly.
- `lbfgsb3x.cpp` -- `lbfgsb3C_` (Fortran-backed C entry) and the Rcpp
  `lbfgsb3cpp()` R entry; both drive `lbfgsb3c_cpp::lbfgsbDriver` with
  either `fortranStep` or `lbfgsb3c_cpp::setulb`.
- `threadtest.cpp` -- internal test helpers `.lbfgsb3cThreadTest` (OpenMP
  vs serial runs) and `.lbfgsb3cLmmTest`.
- `init.c` -- MANUAL `.Call` registration (`R_init_lbfgsb3c`) with
  hardcoded arities, `_lbfgsb3c_ptr` pointer table, and
  `R_RegisterCCallable` for `lbfgsb3C_`/`lbfgsb3Cts_`.
- `Makevars`/`Makevars.win` -- OpenMP flags and `$(FLIBS)`.

Error codes: invalid `lmm` returns `fail = 29` ("ERROR: INVALID LMM");
`n <= 0` returns `fail = 13` before allocating work arrays.

### Public C interface (`inst/include/`)
- `lbfgsb3ptr.h` -- preferred: `iniLbfgsb3` macro + `iniLbfgsb3ptr0()`
  fill `lbfgsb3C` (slot 0) and `lbfgsb3Cts` (slot 1) from `.lbfgsb3cPtr()`.
- `lbfgsb3c.h`, `lbfgsb3x.h` -- older `R_GetCCallable` wrappers.

> [!IMPORTANT]
> The `_lbfgsb3c_ptr` list in `src/init.c` is APPEND-ONLY: never reorder,
> rename or remove slots (`lbfgsb3C`, `lbfgsb3Cts`) or change the
> `lbfgsb3_fn` signature. Reverse dependencies read slots by position and
> released versions cannot be patched. Keep the `R_RegisterCCallable` names
> too. Add new entry points by bumping `nVec` and appending.

### Adding or changing an Rcpp export
Run `Rcpp::compileAttributes(".")` AND hand-edit `src/init.c` (prototype
plus the arg count in `callMethods[]`). A stale arity compiles fine but
fails at runtime with "Incorrect number of arguments".

### Key Files

| File | Purpose |
|------|---------|
| `R/lbfgsb3.R` | `lbfgsb3c()` and aliases, control defaults, roxygen docs |
| `src/lbfgsb3x.cpp` | Rcpp entry `lbfgsb3cpp`, `lbfgsb3C_`, engine dispatch |
| `src/lbfgsb_cpp.cpp` | Thread-safe C++ port of L-BFGS-B |
| `src/lbfgsb.f` | Original Fortran L-BFGS-B 3.0 |
| `src/init.c` | Routine registration and pointer table |
| `inst/include/lbfgsb3ptr.h` | Header for downstream C/C++ users |
| `vignettes/lbfgsb3c.Rmd` | Vignette (Fortran-objective chunk is `eval=FALSE`) |

### Generated Files (do not edit by hand)
- `NAMESPACE`, `man/*.Rd` -- roxygen2
- `R/RcppExports.R`, `src/RcppExports.cpp` -- `Rcpp::compileAttributes()`

## Testing

- testthat edition 3, files `tests/testthat/test-*.R`.
- `test-cpp-engine.R` compares engines (`.cmpEngines()`): identical
  convergence always; full bit-for-bit only with
  `LBFGSB3C_EXACT_ENGINE_TEST=true`. Optimized BLAS/FMA cause rounding
  differences, so do not make default tests require exact equality.
- Thread-safety test uses `.lbfgsb3cThreadTest()`; OpenMP-specific
  expectation is skipped when OpenMP is unavailable.
- `test-240319.R` uses `adagio` only if installed (Suggests).
- Inside tests, call internal functions by bare name (no `:::`).

## CI (`.github/workflows/`)
`R-CMD-check.yaml`, `exact-engine.yaml` (reference BLAS, exact engine
comparison), `test-coverage.yaml`, `pkgdown.yaml`, `rhub.yaml`.

## Code Style and CRAN Cautions

- ASCII only in all files (CRAN): use `--`, `->`, straight quotes.
- Internal R helpers and locals use a leading dot (`.cmpEngines`,
  `.engines`); exported functions are camelCase.
- Do not add global/static mutable state or R API calls to the C++ port;
  it must stay callable concurrently from OpenMP threads.
- Fortran must not use `WRITE`/`PRINT` or non-portable intrinsics (the
  unused `timer.f` was removed because `cpu_time` broke flang 21+).
- The vignette must build outside the source tree (no
  `pkgbuild::compile_dll()`).
- `NEWS.md`: one `# lbfgsb3c <version> changes` section per release,
  past-tense bullets.
- Release files: `cran-comments.md`, `CRAN-SUBMISSION` (both
  Rbuildignored). Reverse dependencies are checked with revdepcheck
  (`revdep/`, git-ignored).
- Do not commit build artifacts (`*.o`, `*.so`, `*.dll`).

## CRAN tests: never assert on timing

- Never put a wall-clock timing expectation in a test that can run on CRAN:
  no `expect_lt()`/`expect_true()` on `Sys.time()` differences,
  `system.time()`, `proc.time()`, `bench`/`microbenchmark` results, or any
  other elapsed-time threshold. CRAN check machines are shared and heavily
  loaded, so elapsed time there is random and a timing assertion fails at
  random. Example: nlmixr2plot `test-model-diagram.R` asserted a
  `modelGraph()` call took < 30 s; on CRAN it took 48.4 s and failed.
- Test the result, not the speed. If a speed regression guard is really
  needed, put it only in a test that starts with `skip_on_cran()` (or the
  package's own CRAN skip helper), and never in examples or vignettes.
- Tests that take more than a few seconds should also skip on CRAN.
