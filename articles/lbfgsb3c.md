# lbfgsb3c: Using the 2011 version of L-BFGS-B.

## Abstract

In 2011 the authors of the L-BFGSB program published a correction and
update to their 1995 code. The latter is the basis of the L-BFGS-B
method of the [`optim()`](https://rdrr.io/r/stats/optim.html) function
in base-R. The package `lbfgsb3` wrapped the updated code using a
`.Fortran` call after removing a very large number of Fortran output
statements. Matthew Fidler used this Fortran code and an `Rcpp`
interface to produce package `lbfgsb3c` where the function
[`lbfgsb3c()`](https://nlmixr2.github.io/lbfgsb3c/reference/lbfgsb3c.md)
returns an object similar to that of base-R
[`optim()`](https://rdrr.io/r/stats/optim.html) and that of
[`optimx::optimr()`](https://rdrr.io/pkg/optimx/man/optimr.html).
Subsequently, in a fine example of the collaborations that have made
**R** so useful, we have merged the functionality of package `lbfgsb3`
into `lbfgsb3c`, as explained in this vignette. Note that this document
is intended primarily to document our efforts to check the differences
in variants of the code rather than be expository.

### Provenance of the R optim::L-BFGS-B and related solvers

The base-R code lbfgsb.c (at writing in R-3.5.2/src/appl/) is commented:

    /* l-bfgs-b.f -- translated by f2c (version 19991025).

      From ?optim:
      The code for method ‘"L-BFGS-B"’ is based on Fortran code by Zhu,
      Byrd, Lu-Chen and Nocedal obtained from Netlib (file 'opt/lbfgs_bcm.shar')

      The Fortran files contained no copyright information.

      Byrd, R. H., Lu, P., Nocedal, J. and Zhu, C.  (1995) A limited
      memory algorithm for bound constrained optimization.
      \emph{SIAM J. Scientific Computing}, \bold{16}, 1190--1208.
    */

The paper Byrd, Lu, Nocedal, and C. Y. Zhu (1995) builds on Lu et al.
(1994). There have been a number of other workers who have followed-up
on this work, but **R** code and packages seem to have largely stayed
with codes derived from these original papers. Though the date of the
paper is 1995, the ideas it embodies were around for a decade and a half
at least, in particular in Nocedal80 and LiuN89. The definitive Fortran
code was published as Zhu et al. (1997). This is available as
`toms/778.zip` on www.netlib.org. A side-by-side comparison of the main
subroutines in the two downloads from Netlib unfortunately shows a lot
of differences. I have not tried to determine if these affect
performance or are simply cosmetic.

More seriously perhaps, there were some deficiencies in the code(s), and
in 2011 Nocedal’s team published a Fortran code with some corrections
(Morales and Nocedal (2011)). Since the **R** code predates this, I
prepared package `lbfgsb3` (Nash et al. (2015)) to wrap the Fortran
code. However, I did not discover any test cases where the
`optim::L-BFGS-B` and `lbfgsb3` were different, though I confess to only
running some limited tests. There are, in fact, more in this vignette.

In 2016, I was at a Fields Institute optimization conference in Toronto
for the 70th birthday of Andy Conn. By sheer serendipity, Nocedal did
not attend the conference, but sat down next to me at the conference
dinner. When I asked him about the key changes, he said that the most
important one was to fix the computation of the machine precision, which
was not always correct in the 1995 code. Since **R** gets this number as
`.Machine$double.eps`, the offending code is irrelevant.

Within Morales and Nocedal (2011), there is also reported an improvement
in the subspace minimization that is applied in cases of bounds
constraints. Since few of the tests I have applied impose such
constraints, it is reasonable that I will not have observed performance
differences between the base-R `optim` code and my `lbfsgb3` package.
More appropriate tests are welcome, and on my agenda.

Besides the ACM TOMS code, there are two related codes from the
Northwestern team on NETLIB: <https://netlib.org/opt/lbfgs_um.shar> is
for unconstrained minimization, while
<https://netlib.org/opt/lbfgs_bcm.shar> handles bounds constrained
problems. To these are attached references Liu and Nocedal (1989) and
Byrd, Lu, Nocedal, and C. Zhu (1995) respectively, most likely
reflecting the effort required to implement the constraints.

The unconstrained code has been converted to **C** under the leadership
of Naoaki Okazaki (see <http://www.chokkan.org/software/liblbfgs/>, or
the fork at <https://github.com/MIRTK/LBFGS>). This has been wrapped for
**R** as Coppola et al. (2014) as the `lbfgs` package. This can be
called from
[`optimx::optimr()`](https://rdrr.io/pkg/optimx/man/optimr.html).

Using Rcpp (see Eddelbuettel and François (2011)) and the Fortran code
in package `lbfgs3`, Matthew Fidler developed package `lbfgsb3c` (Fidler
et al. (2018)). As this provides a more standard call and return than
`lbfgsb3` Fidler and I are unified the two packages and released them
both under the same name `lbfgsb3c`.

### Functions in package `lbfgsb3c`

There is really only one optimizer function in the package, but it may
be called by four (4) names:

- [`lbfgsb3c()`](https://nlmixr2.github.io/lbfgsb3c/reference/lbfgsb3c.md)
  uses Rcpp (Eddelbuettel (2013), Eddelbuettel and François (2011),
  Eddelbuettel and Balamuta (2017)) to streamline the call to the
  underlying Fortran. This is the base function used.
- [`lbfgsb3x()`](https://nlmixr2.github.io/lbfgsb3c/reference/lbfgsb3c.md)
  is an alias of
  [`lbfgsb3c()`](https://nlmixr2.github.io/lbfgsb3c/reference/lbfgsb3c.md).
  We were using this name for a while, and have kept the alias to avoid
  having to edit test scripts.
- `lbfgsb3`, which imitates a `.Fortran` call of the compiled 2011
  Fortran code. The object returned by this routine is NOT equivalent to
  the object returned by base-R
  [`optim()`](https://rdrr.io/r/stats/optim.html) or by
  [`optimx::optimr()`](https://rdrr.io/pkg/optimx/man/optimr.html).
  Instead, it includes a structure `info` which contains the detailed
  diagnostic information of the Fortran code. For most users, this is
  not of interest, and I only recommend use of this function for those
  needing to examine how the optimization has been carried out.
- [`lbfgsb3f()`](https://nlmixr2.github.io/lbfgsb3c/reference/lbfgsb3c.md)
  is an alias of `lbfsgb3()`.

We recommend using the `lbfsgb3c()` call for most uses.

#### Candlestick function

``` r

# candlestick function
# J C Nash 2011-2-3
cstick.f<-function(x,alpha=100){
  x<-as.vector(x)
  r2<-crossprod(x)
  f<-as.double(r2+alpha/r2)
  return(f)
}

cstick.g<-function(x,alpha=100){
  x<-as.vector(x)
  r2<-as.numeric(crossprod(x))
  g1<-2*x
  g2 <- (-alpha)*2*x/(r2*r2)
  g<-as.double(g1+g2)
  return(g)
}
library(lbfgsb3c)
nn <- 2
x0 <- c(10,10)
lo <- c(1, 1)
up <- c(10,10)
print(x0)
```

    ## [1] 10 10

``` r

## c2o <- opm(x0, cstick.f, cstick.g, lower=lo, upper=up, method=meths, control=list(trace=0))
## print(summary(c2o, order=value))
c2l1 <- lbfgsb3c(x0, cstick.f, cstick.g, lower=lo, upper=up)
c2l1
```

    ## $par
    ## [1] 2.236068 2.236068
    ## 
    ## $grad
    ## [1] -1.121272e-06 -1.121272e-06
    ## 
    ## $value
    ## [1] 20
    ## 
    ## $counts
    ## [1] 14 14
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $message
    ## [1] "CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH"

``` r

## meths <- c("L-BFGS-B", "lbfgsb3c", "Rvmmin", "Rcgmin", "Rtnmin")
## require(optimx)

## cstick2a <- opm(x0, cstick.f, cstick.g, method=meths, upper=up, lower=lo, control=list(kkt=FALSE))
## print(summary(cstick2a, par.select=1:2, order=value))
lo <- c(4, 4)
## c2ob <- opm(x0, cstick.f, cstick.g, lower=lo, upper=up, method=meths, control=list(trace=0))
## print(summary(c2ob, order=value))
c2l2 <- lbfgsb3c(x0, cstick.f, cstick.g, lower=lo, upper=up)
c2l2
```

    ## $par
    ## [1] 4 4
    ## 
    ## $grad
    ## [1] 7.21875 7.21875
    ## 
    ## $value
    ## [1] 35.125
    ## 
    ## $counts
    ## [1] 2 2
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $message
    ## [1] "CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL"

``` r

## cstick2b <- opm(x0, cstick.f, cstick.g, method=meths, upper=up, lower=lo, control=list(kkt=FALSE))
## print(summary(cstick2b, par.select=1:2, order=value))

## nn <- 100
## x0 <- rep(10, nn)
## up <- rep(10, nn)
## lo <- rep(1e-4, nn)
## cco <- opm(x0, cstick.f, cstick.g, lower=lo, upper=up, method=meths, control=list(trace=0, kkt=FALSE))
## print(summary(cco, par.select=1:4, order=value))
```

#### Extended Rosenbrock function (from funconstrain)

``` r

# require(funconstrain) ## not in CRAN, so explicit inclusion of this function
# exrosen <- ex_rosen()
# exrosenf <- exrosen$fn
exrosenf <- function (par) {
    n <- length(par)
    if (n%%2 != 0) {
        stop("Extended Rosenbrock: n must be even")
    }
    fsum <- 0
    for (i in 1:(n/2)) {
        p2 <- 2 * i
        p1 <- p2 - 1
        f_p1 <- 10 * (par[p2] - par[p1]^2)
        f_p2 <- 1 - par[p1]
        fsum <- fsum + f_p1 * f_p1 + f_p2 * f_p2
    }
    fsum
}
# exroseng <- exrosen$gr
exroseng <- function (par) {
    n <- length(par)
    if (n%%2 != 0) {
        stop("Extended Rosenbrock: n must be even")
    }
    grad <- rep(0, n)
    for (i in 1:(n/2)) {
        p2 <- 2 * i
        p1 <- p2 - 1
        xx <- par[p1] * par[p1]
        yx <- par[p2] - xx
        f_p1 <- 10 * yx
        f_p2 <- 1 - par[p1]
        grad[p1] <- grad[p1] - 400 * par[p1] * yx - 2 * f_p2
        grad[p2] <- grad[p2] + 200 * yx
    }
    grad
}

exrosenx0 <- function (n = 20) {
    if (n%%2 != 0) {
        stop("Extended Rosenbrock: n must be even")
    }
    rep(c(-1.2, 1), n/2)
}


require(lbfgsb3c)
## require(optimx)

## require(optimx)
for (n in seq(2,12, by=2)) {
  cat("ex_rosen try for n=",n,"\n")
  x0 <- exrosenx0(n)
  lo <- rep(-1.5, n)
  up <- rep(3, n)
  print(x0)
  cat("optim L-BFGS-B\n")
  eo <- optim(x0, exrosenf, exroseng, lower=lo, upper=up, method="L-BFGS-B", control=list(trace=0))
  print(eo)
  cat("lbfgsb3c\n")
  el <- lbfgsb3c(x0, exrosenf, exroseng, lower=lo, upper=up, control=list(trace=0))
  print(el)
##    erfg <- opm(x0, exrosenf, exroseng, method=meths, lower=lo, upper=up)
##    print(summary(erfg, par.select=1:2, order=value))
}
```

    ## ex_rosen try for n= 2 
    ## [1] -1.2  1.0
    ## optim L-BFGS-B
    ## $par
    ## [1] 1 1
    ## 
    ## $value
    ## [1] 3.84442e-14
    ## 
    ## $counts
    ## function gradient 
    ##       51       51 
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $message
    ## [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
    ## 
    ## lbfgsb3c
    ## $par
    ## [1] 1 1
    ## 
    ## $grad
    ## [1] -8.65458e-07  6.26284e-07
    ## 
    ## $value
    ## [1] 3.844419e-14
    ## 
    ## $counts
    ## [1] 51 51
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $message
    ## [1] "CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH"
    ## 
    ## ex_rosen try for n= 4 
    ## [1] -1.2  1.0 -1.2  1.0
    ## optim L-BFGS-B
    ## $par
    ## [1] 1 1 1 1
    ## 
    ## $value
    ## [1] 7.688832e-14
    ## 
    ## $counts
    ## function gradient 
    ##       51       51 
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $message
    ## [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
    ## 
    ## lbfgsb3c
    ## $par
    ## [1] 1 1 1 1
    ## 
    ## $grad
    ## [1] -8.65458e-07  6.26284e-07 -8.65458e-07  6.26284e-07
    ## 
    ## $value
    ## [1] 7.688837e-14
    ## 
    ## $counts
    ## [1] 51 51
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $message
    ## [1] "CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH"
    ## 
    ## ex_rosen try for n= 6 
    ## [1] -1.2  1.0 -1.2  1.0 -1.2  1.0
    ## optim L-BFGS-B
    ## $par
    ## [1] 1 1 1 1 1 1
    ## 
    ## $value
    ## [1] 1.153326e-13
    ## 
    ## $counts
    ## function gradient 
    ##       51       51 
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $message
    ## [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
    ## 
    ## lbfgsb3c
    ## $par
    ## [1] 1 1 1 1 1 1
    ## 
    ## $grad
    ## [1] -8.654581e-07  6.262840e-07 -8.654581e-07  6.262840e-07 -8.654581e-07
    ## [6]  6.262840e-07
    ## 
    ## $value
    ## [1] 1.153324e-13
    ## 
    ## $counts
    ## [1] 51 51
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $message
    ## [1] "CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH"
    ## 
    ## ex_rosen try for n= 8 
    ## [1] -1.2  1.0 -1.2  1.0 -1.2  1.0 -1.2  1.0
    ## optim L-BFGS-B
    ## $par
    ## [1] 1 1 1 1 1 1 1 1
    ## 
    ## $value
    ## [1] 1.537766e-13
    ## 
    ## $counts
    ## function gradient 
    ##       51       51 
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $message
    ## [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
    ## 
    ## lbfgsb3c
    ## $par
    ## [1] 1 1 1 1 1 1 1 1
    ## 
    ## $grad
    ## [1] -8.65458e-07  6.26284e-07 -8.65458e-07  6.26284e-07 -8.65458e-07
    ## [6]  6.26284e-07 -8.65458e-07  6.26284e-07
    ## 
    ## $value
    ## [1] 1.537766e-13
    ## 
    ## $counts
    ## [1] 51 51
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $message
    ## [1] "CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH"
    ## 
    ## ex_rosen try for n= 10 
    ##  [1] -1.2  1.0 -1.2  1.0 -1.2  1.0 -1.2  1.0 -1.2  1.0
    ## optim L-BFGS-B
    ## $par
    ##  [1] 1 1 1 1 1 1 1 1 1 1
    ## 
    ## $value
    ## [1] 1.922208e-13
    ## 
    ## $counts
    ## function gradient 
    ##       51       51 
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $message
    ## [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
    ## 
    ## lbfgsb3c
    ## $par
    ##  [1] 1 1 1 1 1 1 1 1 1 1
    ## 
    ## $grad
    ##  [1] -8.654580e-07  6.262841e-07 -8.654580e-07  6.262841e-07 -8.654580e-07
    ##  [6]  6.262841e-07 -8.654580e-07  6.262841e-07 -8.654580e-07  6.262841e-07
    ## 
    ## $value
    ## [1] 1.92221e-13
    ## 
    ## $counts
    ## [1] 51 51
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $message
    ## [1] "CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH"
    ## 
    ## ex_rosen try for n= 12 
    ##  [1] -1.2  1.0 -1.2  1.0 -1.2  1.0 -1.2  1.0 -1.2  1.0 -1.2  1.0
    ## optim L-BFGS-B
    ## $par
    ##  [1] 1 1 1 1 1 1 1 1 1 1 1 1
    ## 
    ## $value
    ## [1] 2.306651e-13
    ## 
    ## $counts
    ## function gradient 
    ##       51       51 
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $message
    ## [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
    ## 
    ## lbfgsb3c
    ## $par
    ##  [1] 1 1 1 1 1 1 1 1 1 1 1 1
    ## 
    ## $grad
    ##  [1] -8.65458e-07  6.26284e-07 -8.65458e-07  6.26284e-07 -8.65458e-07
    ##  [6]  6.26284e-07 -8.65458e-07  6.26284e-07 -8.65458e-07  6.26284e-07
    ## [11] -8.65458e-07  6.26284e-07
    ## 
    ## $value
    ## [1] 2.306649e-13
    ## 
    ## $counts
    ## [1] 51 51
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $message
    ## [1] "CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH"

### Using compiled function code

While you may use the same interface as described in the writing R
extensions to interface compiled code with this function, see
[L-BFGS-B](https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Optimization),
it is sometimes more convenient to use your own compiled code.

The following example shows how this is done using the file `jrosen.f`.
We have unfortunately found that compilation is not always portable
across systems, so this example is presented without execution.

           subroutine rosen(n, x, fval)
           double precision x(n), fval, dx
           integer n, i
           fval = 0.0D0
           do 10 i=1,(n-1)
              dx = x(i + 1) - x(i) * x(i)
              fval = fval + 100.0 * dx * dx
              dx = 1.0 - x(i)
              fval = fval + dx * dx
     10    continue
           return
           end

Here is the example script. Note that we must have the file `jrosen.f`
available. Because the executable files on different systems use
different conventions and structures, we have turned evaluation off here
so this vignette can be built on multiple platforms. However, we wished
to provide examples of how compiled code could be used.

``` r

system("R CMD SHLIB jrosen.f")
dyn.load("jrosen.so")
is.loaded("rosen")
x0 <- as.double(c(-1.2,1))
fv <- as.double(-999)
n <- as.double(2)
testf <- .Fortran("rosen", n=as.integer(n), x=as.double(x0), fval=as.double(fv))
testf

rrosen <- function(x) {
  fval <- 0.0
  for (i in 1:(n-1)) {
    dx <- x[i + 1] - x[i] * x[i]
    fval <- fval + 100.0 * dx * dx
    dx <- 1.0 - x[i]
    fval <- fval + dx * dx
  }
  fval
}

(rrosen(x0))

frosen <- function(x){
  nn <- length(x)
  if (nn > 100) { stop("max number of parameters is 100")}
  fv <- -999.0
  val <- .Fortran("rosen", n=as.integer(nn), x=as.double(x), fval=as.double(fv))
  val$fval # NOTE--need ONLY function value returned
}
# Test the funcion
tval <- frosen(x0)
str(tval)

cat("Run with Nelder-Mead using R function\n")
mynm <- optim(x0, rrosen, control=list(trace=0))
print(mynm)
cat("\n\n Run with Nelder-Mead using Fortran function")
mynmf <- optim(x0, frosen, control=list(trace=0))
print(mynmf)


library(lbfgsb3c)
library(microbenchmark)
cat("try lbfgsb3c, no Gradient \n")
cat("R function\n")
tlR<-microbenchmark(myopR <- lbfgsb3c(x0, rrosen, gr=NULL, control=list(trace=0)))
print(tlR)
print(myopR)
cat("Fortran function\n")
tlF<-microbenchmark(myop <- lbfgsb3c(x0, frosen, gr=NULL, control=list(trace=0)))
print(tlF)
print(myop)
```

In this example, Fortran execution was actually SLOWER than plain R on
the system where it was run.

### References

Byrd, Richard H., Peihuang Lu, Jorge Nocedal, and Ci You Zhu. 1995. “A
Limited Memory Algorithm for Bound Constrained Optimization.” *SIAM
Journal on Scientific Computing* 16 (5): 1190–208.

Byrd, Richard H., Peihuang Lu, Jorge Nocedal, and Ciyou Zhu. 1995. “A
Limited Memory Algorithm for Bound Constrained Optimization.” *SIAM J.
Sci. Comput.* (Philadelphia, PA, USA) 16 (5): 1190–208.
<https://doi.org/10.1137/0916069>.

Coppola, Antonio, Brandon Stewart, and Naoaki Okazaki. 2014. *Lbfgs:
Limited-Memory BFGS Optimization*.
<https://CRAN.R-project.org/package=lbfgs>.

Eddelbuettel, Dirk. 2013. *Seamless R and C++ Integration with Rcpp*.
Springer. <https://doi.org/10.1007/978-1-4614-6868-4>.

Eddelbuettel, Dirk, and James Joseph Balamuta. 2017. “Extending extitR
with extitC++: A Brief Introduction to extitRcpp.” *PeerJ Preprints* 5
(August): e3188v1. <https://doi.org/10.7287/peerj.preprints.3188v1>.

Eddelbuettel, Dirk, and Romain François. 2011. “Rcpp: Seamless R and C++
Integration.” *Journal of Statistical Software* 40 (8): 1–18.
<https://doi.org/10.18637/jss.v040.i08>.

Fidler, Matthew L, John C Nash, Ciyou Zhu, Richard Byrd, Jorge Nocedal,
and Jose Luis Morales. 2018. *lbfgsb3c: Limited Memory BFGS Minimizer
with Bounds on Parameters with Optim() ’c’ Interface*.
<https://CRAN.R-project.org/package=lbfgsb3c>.

Liu, Dong C., and Jorge Nocedal. 1989. “On the Limited Memory BFGS
Method for Large Scale Optimization.” *Math. Program.* 45 (1-3): 503–28.
<https://doi.org/10.1007/BF01589116>.

Lu, Peihuang, Jorge Nocedal, Ciyou Zhu, and Richard H. Byrd. 1994. “A
Limited-Memory Algorithm for Bound Constrained Optimization.” *SIAM
Journal on Scientific Computing* 16: 1190–208.

Morales, José Luis, and Jorge Nocedal. 2011. “Remark on Algorithm 778:
L-BFGS-B: Fortran subroutines for large-scale bound constrained
optimization.” *ACM Trans. Math. Softw.* (New York, NY, USA) 38 (1):
7:1–4.

Nash, John C, Ciyou Zhu, Richard Byrd, Jorge Nocedal, and Jose Luis
Morales. 2015. *lbfgsb3: Limited Memory BFGS Minimizer with Bounds on
Parameters*. <https://CRAN.R-project.org/package=lbfgsb3>.

Zhu, Ciyou, Richard H. Byrd, Peihuang Lu, and Jorge Nocedal. 1997.
“Algorithm 778: L-BFGS-B: Fortran Subroutines for Large-Scale
Bound-Constrained Optimization.” *ACM Trans. Math. Softw.* (New York,
NY, USA) 23 (4): 550–60. <https://doi.org/10.1145/279232.279236>.
