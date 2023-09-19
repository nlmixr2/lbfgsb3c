# lbfgsb3c development changes

* Now allow `rho=NULL` to work the same as if `rho` was not supplied

* Be More careful about `$convergence` by adding a default value of `NA_INTEGER`

* `$convergence` is now an integer instead of a real number

* Added a `NEWS.md` file to track changes to the package.

# News before lbfgsb3c

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
