library(robustbase)

## minimal testing only
data(ruspini, package = "cluster")

### NOTA BENE:  scale.tau2() is *not* consistent {by constant factor}
rub1 <- covOGK(ruspini, 1, scaleTau2, covGK, hard.rejection)
rub2 <- covOGK(ruspini, 2, scaleTau2, covGK, hard.rejection)

AE <- function(x,y) all.equal(x,y, tol = 2e-15)
## The following test is already fulfilled by Kjell Konis'  original code:
stopifnot(AE(c(rub1$wcov)[c(1,3:4)],
             c(917.99893333333, 94.9232, 2340.319288888888)),
          all.equal(rub1$wcov, rub2$wcov, tol=0)
          ,
          AE(c(rub1$cov)[c(1,3:4)],
             c(923.5774514441657, 91.5385216376565, 2342.4556232436971))
          ,
          AE(c(rub2$cov)[c(1,3:4)],
             c(927.2465953711782, 91.8009184487779, 2346.5790105548940))
          )

q()

## More tests: ------------ not yet / too slow --> do other place ------------

set.seed(101)
library(MASS)# mvrnorm()
N <- 1000
p <- 50

mu <- rep(0, p)
## p x p  cov.matrix with const. sigma=10
S0 <- toeplitz(10* ARMAacf(0.8, lag.max = p - 1))


nSim <- 20
r <- vector("list", nSim)
cpu <- numeric(nSim)
for(i in 1:nSim) {
    X <- mvrnorm(N, mu = mu, Sigma = S0)
    cpu[i] <-
        system.time(cc <- covOGK(X, n.it = 2,
                                 sigmamu = scaleTau2, weight.fn = hard.rejection))
    r[[i]] <- cc
}


cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
