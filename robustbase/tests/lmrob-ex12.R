
library(robustbase)

## EX 1
data(coleman)
(m1 <- lmrob(Y ~ ., data = coleman))
summary(m1)
## Values will change once we use R's random number generator !
##>> stopifnot(
all.equal(unname(coef(m1)),
          c(29.5123314250361, -1.66091845261954, 0.083453132305779,
            0.665713483136765, 1.17852545547180, -4.01036910730544))
##>> )

## EX 2
gen <- function(n,p, n0, y0, x0, beta = rep(1, p))
{
    stopifnot(n >= 1, p >= 1, n0 >= 0, length(beta) == p)
    x <- matrix(rnorm(n*p),n,p) # iid x's
    y <- x %*% beta + rnorm(n)
    xc <- matrix(0,n0,p)
    xc[,1] <- x0
    xc <- xc + 0.1*matrix(rnorm(n0*p),n0,p)
    x[1:n0,] <- xc
    y[1:n0] <- y0 + .001*rnorm(n0)
    list(x=x, y=y)
}

## generate --a sample of  n  observations with  p  variables
## and 10% of outliers near (x1,y) = (10,10)
n <- 500 ; n0 <- n %/% 10
p <- 7 ## p = 20 is more impressive but too slow for "standard test"
set.seed(17)
a <- gen(n=n, p=p, n0= n0, y0=10, x0=10)
plot(a$x[,1], a$y, col = c(rep(2, n0), rep(1, n-n0)))
system.time( tmp <- lmrob(y~x, data = a) )
plot(tmp) #- currently 5 plots; MM: don't like #3 (Response vs fitted)

## don't compute robust distances --> much faster (factor of  ~ 8 !)
system.time(t2 <- lmrob(y~x, data = a,
                        control = lmrob.control(compute.rd = FALSE)))
## ==> most of the CPU time is spent in cov.rob()!
(st2 <- summary(t2))
l1 <- lm(y~x, data = a)
cbind(robust = coef(st2)[,1:2],
      lm = coef(summary(l1))[,1:2])
## rm(a,tmp, t2)



cat('Time elapsed: ', proc.time(),'\n') # "stats"
