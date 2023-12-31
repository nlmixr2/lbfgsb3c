test_that("bounds test", {
    bt.f<-function(x){
        sum(x*x)
    }

    bt.g<-function(x){
        gg<-2.0*x
    }

    n<-4
    lower<-rep(0,n)
    upper<-lower ## to get arrays set
    ## bdmsk<-rep(1,n)
    ## bdmsk[(trunc(n/2)+1)]<-0
    for (i in 1:n) {
        lower[i]<-1.0*(i-1)*(n-1)/n
        upper[i]<-1.0*i*(n+1)/n
    }
    xx<-0.5*(lower+upper)

    abt <- lbfgsb3c(xx, bt.f, bt.g, lower=lower, upper=upper,
                    control=list(trace=0, iprint=0))
    expect_equal(abt$par, c(0, 0.75, 1.5, 2.25))
})
