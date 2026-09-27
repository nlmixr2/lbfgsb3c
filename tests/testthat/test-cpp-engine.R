# The C++ port (engine="cpp") must reproduce the Fortran (engine="fortran")
.cmpEngines <- function(par, fn, gr, lower = -Inf, upper = Inf,
                        control = list(), ...) {
  control$info <- TRUE
  .f <- lbfgsb3c(par, fn, gr, lower = lower, upper = upper,
                 control = c(control, list(engine = "fortran")), ...)
  .c <- lbfgsb3c(par, fn, gr, lower = lower, upper = upper,
                 control = c(control, list(engine = "cpp")), ...)
  expect_equal(.c$par, .f$par, tolerance = 1e-12)
  expect_equal(.c$value, .f$value, tolerance = 1e-12)
  expect_equal(.c$grad, .f$grad, tolerance = 1e-10)
  expect_identical(.c$counts, .f$counts)
  expect_identical(.c$convergence, .f$convergence)
  expect_identical(.c$message, .f$message)
  expect_identical(.c$info$itask, .f$info$itask)
  expect_identical(.c$info$isave, .f$info$isave)
  expect_identical(.c$info$lsave, .f$info$lsave)
  expect_equal(.c$info$dsave, .f$info$dsave, tolerance = 1e-10)
  invisible(list(fortran = .f, cpp = .c))
}

test_that("engine control is validated", {
  expect_error(lbfgsb3c(c(1, 2), function(x) sum(x^2), function(x) 2 * x,
                        control = list(engine = "julia")))
})

.rosen.f <- function(x) {
  n <- length(x)
  1 + sum(100 * (x[2:n] - x[1:(n - 1)]^2)^2 + (1 - x[2:n])^2)
}
.rosen.g <- function(x) {
  n <- length(x)
  g <- double(n)
  tn <- 2:n
  z1 <- x[tn] - x[tn - 1]^2
  g[tn] <- 2 * (100 * z1 - (1 - x[tn]))
  g[tn - 1] <- g[tn - 1] - 400 * x[tn - 1] * z1
  g
}

test_that("cpp engine matches fortran: unconstrained Rosenbrock", {
  .cmpEngines(c(a = -1.2, b = 1, c = -1.2, d = 1), .rosen.f, .rosen.g)
  .cmpEngines(rep(3, 25), .rosen.f, .rosen.g)
  .cmpEngines(c(1.02, 1.02, 1.02), .rosen.f, .rosen.g)
})

test_that("cpp engine matches fortran: two-sided, one-sided and fixed bounds", {
  .cmpEngines(rep(3, 25), .rosen.f, .rosen.g, lower = 2, upper = 4)
  .cmpEngines(rep(3, 10), .rosen.f, .rosen.g, lower = 2)
  .cmpEngines(rep(-3, 10), .rosen.f, .rosen.g, upper = -0.5)
  .cmpEngines(rep(0.5, 6), .rosen.f, .rosen.g,
              lower = c(-Inf, 0, 0.5, -1, -Inf, 0.2),
              upper = c(Inf, Inf, 0.5, 2, 0.8, Inf))
  # infeasible start is projected onto the box
  .cmpEngines(rep(10, 5), .rosen.f, .rosen.g, lower = -1, upper = 2)
})

test_that("cpp engine matches fortran: bounds test", {
  bt.f <- function(x) sum(x * x)
  bt.g <- function(x) 2.0 * x
  n <- 4
  lower <- (seq_len(n) - 1) * (n - 1) / n
  upper <- seq_len(n) * (n + 1) / n
  .r <- .cmpEngines(0.5 * (lower + upper), bt.f, bt.g, lower = lower,
                    upper = upper)
  expect_equal(.r$cpp$par, c(0, 0.75, 1.5, 2.25))
})

test_that("cpp engine matches fortran: Chebyquad", {
  cyq.res <- function(x) {
    n <- length(x)
    res <- rep(0, n)
    for (i in 1:n) {
      rr <- 0
      for (k in 1:n) {
        z7 <- 1
        z2 <- 2 * x[k] - 1
        z8 <- z2
        j <- 1
        while (j < i) {
          z6 <- z7
          z7 <- z8
          z8 <- 2 * z2 * z7 - z6
          j <- j + 1
        }
        rr <- rr + z8
      }
      rr <- rr / n
      if (2 * trunc(i / 2) == i) rr <- rr + 1 / (i * i - 1)
      res[i] <- rr
    }
    res
  }
  cyq.f <- function(x) sum(cyq.res(x)^2)
  cyq.g <- function(x) numDeriv::grad(cyq.f, x)
  for (n in c(2, 3, 5, 8, 10)) {
    .cmpEngines(seq_len(n) / (n + 1), cyq.f, cyq.g, lower = -10, upper = 10)
  }
})

test_that("cpp engine matches fortran: control settings", {
  # maximum number of iterations
  .r <- .cmpEngines(rep(3, 25), .rosen.f, .rosen.g, control = list(maxit = 5))
  expect_identical(.r$cpp$convergence, 1L)
  # x-tolerance stop, pgtol, factr and lmm
  .cmpEngines(rep(3, 25), .rosen.f, .rosen.g,
              control = list(abstol = 1e-4, reltol = 0))
  .cmpEngines(rep(3, 25), .rosen.f, .rosen.g, control = list(pgtol = 1e-3))
  .cmpEngines(rep(3, 25), .rosen.f, .rosen.g, control = list(factr = 1e12))
  .cmpEngines(rep(3, 25), .rosen.f, .rosen.g, control = list(lmm = 1))
  .cmpEngines(rep(3, 25), .rosen.f, .rosen.g, control = list(lmm = 17))
  # extra arguments pass through
  .cmpEngines(rep(3, 10), function(x, s) .rosen.f(x) * s,
              function(x, s) .rosen.g(x) * s, s = 2)
})

test_that("cpp engine matches fortran: invalid bounds", {
  .r <- suppressWarnings(
    .cmpEngines(c(1, 1), function(x) sum(x^2), function(x) 2 * x,
                lower = c(0, 2), upper = c(1, 1)))
  expect_identical(.r$cpp$message, "ERROR: NO FEASIBLE SOLUTION")
})

test_that("cpp engine prints trace output", {
  expect_output(
    lbfgsb3c(c(-1.2, 1), .rosen.f, .rosen.g,
             control = list(engine = "cpp", trace = 1)),
    "At iteration")
  expect_output(
    lbfgsb3c(c(-1.2, 1), .rosen.f, .rosen.g,
             control = list(engine = "cpp", iprint = 1)),
    "At iterate")
})

test_that("lbfgsb3Cts is thread safe", {
  .r <- lbfgsb3c:::.lbfgsb3cThreadTest(64L, 4L)
  expect_identical(.r$parallel, .r$serial)
  # every problem ran to a convergence code
  expect_true(all(.r$serial[, 14] %in% c(6, 7, 8, 27)))
  skip_if_not(.r$openmp, "OpenMP not available")
  expect_true(.r$openmp)
})

test_that(".lbfgsb3cPtr exposes the thread-safe pointer", {
  .p <- .lbfgsb3cPtr()
  expect_identical(names(.p), c("lbfgsb3C", "lbfgsb3Cts"))
  expect_true(all(vapply(.p, typeof, character(1)) == "externalptr"))
})

test_that("nested optimizations do not share workspace or R callbacks", {
  # The inner problem is solved inside the outer objective, while the
  # outer run holds this thread's workspace and R callbacks.  The inner
  # minimum of s * rosen is s (at x = 1), so the outer objective is
  # (x1 - 1)^2 + (x2 - 1)^2 + 1 + x1^2, minimized at (0.5, 1).
  .inner <- function(s, engine) {
    lbfgsb3c(rep(3, 4), function(x) .rosen.f(x) * s,
             function(x) .rosen.g(x) * s,
             control = list(engine = engine, factr = 10))$value
  }
  for (.engine in c("fortran", "cpp")) {
    .of <- function(x) sum((x - 1)^2) + .inner(1 + x[1]^2, .engine)
    .og <- function(x) c(2 * (x[1] - 1) + 2 * x[1], 2 * (x[2] - 1))
    .r <- lbfgsb3c(c(a = 2, b = 2), .of, .og,
                   control = list(engine = .engine))
    expect_equal(unname(.r$par), c(0.5, 1), tolerance = 1e-4,
                 info = .engine)
    expect_identical(names(.r$par), c("a", "b"))
    expect_identical(.r$convergence, 0L)
  }
})

test_that("nonsensical lmm values error", {
  .f <- function(x) sum(x^2)
  .g <- function(x) 2 * x
  for (.engine in c("fortran", "cpp")) {
    for (.lmm in list(0, -1, 2.5, NA, NA_integer_, Inf, c(3, 4), "5")) {
      expect_error(lbfgsb3c(c(1, 2), .f, .g,
                            control = list(lmm = .lmm, engine = .engine)),
                   "lmm", info = paste(.engine, format(.lmm)))
    }
    # passes the R check but the workspace size overflows an int
    expect_error(lbfgsb3c(c(1, 2), .f, .g,
                          control = list(lmm = 20000, engine = .engine)),
                 "lmm")
    expect_equal(lbfgsb3c(c(1, 2), .f, .g,
                          control = list(lmm = 1, engine = .engine))$par,
                 c(0, 0))
  }
})

test_that("C entry points reject invalid lmm and n without evaluating", {
  # lmm = 1e9 would overflow even 64-bit workspace-size arithmetic
  .r <- lbfgsb3c:::.lbfgsb3cLmmTest(c(0L, -3L, 20000L, 1000000000L, 5L, 5L, 5L),
                                    c(3L, 3L, 3L, 3L, 0L, -1L, 3L))
  expect_identical(.r$fail[1:6], c(rep(29L, 4), 13L, 13L))
  expect_identical(.r$fncount[1:6], rep(0L, 6))
  expect_true(.r$fail[7] %in% c(6L, 7L, 8L, 27L))
})

test_that("a gradient of the wrong length is an error", {
  for (.engine in c("fortran", "cpp")) {
    expect_error(lbfgsb3c(c(1, 2, 3), function(x) sum(x^2),
                          function(x) 2 * x[1:2],
                          control = list(engine = .engine)),
                 "gradient must have 3 elements")
  }
})
