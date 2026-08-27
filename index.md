# libfgsb3c interface from C

This is the fork of the libfgsb3 from cran with the following
differences: - The return type has changed is is very similar to what
`optim` returns - Allows a direct C/C++ interface through a R registered
function, similar to C interface to `optim` with 2 additional
arguments. - Allows adjustment of tolerances for minimization success. -
Added `xtolAtol` and `xtolRtol` minimization success criterion. - Added
`maxit` termination
