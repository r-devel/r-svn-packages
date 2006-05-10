
library(robustbase)

## EX 1
data(coleman)
## "empty model" (not really a lot of sense)
(m0 <- lmrob(Y ~ 0, data = coleman))
summary(m0)
## "Intercept" only: robust mean
(m1 <- lmrob(Y ~ 1, data = coleman))
summary(m1)


(mC <- lmrob(Y ~ ., data = coleman))
summary(mC)
## Values will change once we use R's random number generator !
stopifnot(
all.equal(unname(coef(mC)),
          c(29.45, -1.650, 0.0830, 0.6656, 1.177, -3.995),
          ## tol = 0.00243 for 64-bit, 0.00242 for 32-bit
          ## --
          ## was c(29.52, -1.662, 0.0835, 0.6657, 1.179, -4.012),
          ##     tol =0.000146 for 64-bit, 0.000596 for 32-bit
          tol = 0.004)
)
dput(signif(unname(coef(mC)), 7))
## 2006-04-28 :
## 64b (0.1-6): c(29.53885, -1.664495, 0.08360545, 0.6657679, 1.179005, -4.016079)
## 32b (0.1-6): c(29.35569, -1.635266, 0.08242489, 0.6653669, 1.174706, -3.973863)

##-> mean       c(29.44727, -1.64988,  0.08301517, 0.6655674, 1.176856, -3.994971)

## by repeating the above, we get *different* results !!!
##  [even though the seed is the same ?????   :
## 64b also had c(29.4531,  -1.651344, 0.08305924, 0.6655888, 1.177105, -3.996584)
##              c(29.49953, -1.658865, 0.08336762, 0.6656872, 1.178221, -4.007395
## 32b also had c(29.39156, -1.64123, 0.08265623, 0.6654513, 1.175597, -3.982232)

## 2006-04-24 :
## 32-bit (0.1-6): c(29.411951606792, -1.64459850945731, 0.0827889726545072,
##                   0.665497831767958, 1.17609979960314, -3.98699004062665)
## 64-bit (0.1-6): c(29.5388508044274, -1.66449462058860, 0.0836054450771018,
##                   0.665767853354177, 1.17900459411589, -4.01607879779799)


## 2006-04-20 :
## 64-bit (0.1-6): c(29.5388508044274, -1.66449462058860, 0.0836054450771018,
##                   0.665767853354177, 1.17900459411589, -4.01607879779799)
## 32-bit (0.1-6): c(29.5044520305338, -1.65965482653590, 0.0834004535344066,
##                   0.665697308050525, 1.17833839348011, -4.00853793688203)

str(mC)


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
system.time( m1 <- lmrob(y~x, data = a) )
plot(m1, ask=FALSE)
##-> currently 5 plots; MM:I  don't like #3 (Response vs fitted)

## don't compute robust distances --> faster by factor of two:
system.time(m2 <- lmrob(y~x, data = a,
                        control = lmrob.control(compute.rd = FALSE)))
## ==> half of the CPU time is spent in covMcd()!
(sm2 <- summary(m2))
l1 <- lm(y~x, data = a)
cbind(robust = coef(sm2)[,1:2],
      lm = coef(summary(l1))[,1:2])

##--- Now use n > 2000 --> so we use C internal fast_s_large_n(...)
n <- 2500 ; n0 <- n %/% 10
a2 <- gen(n=n, p = 3, n0= n0, y0=10, x0=10)
plot(a2$x[,1], a2$y, col = c(rep(2, n0), rep(1, n-n0)))
system.time( m3 <- lmrob(y~x, data = a2) )
m3
system.time( m4 <- lmrob(y~x, data = a2, compute.rd = FALSE))
(sm4 <- summary(m4))

dput(signif(cf <- unname(coef(m3)), 7))
## 2006-05-12(64-b): c(0.03801423, 0.9965265, 1.005553, 0.9998042)
## 2006-05-12(64-b): c(0.03842982, 0.996453,  1.005504, 1.000424) #< R CMD check

## 2006-05-12(32-b):
## 2006-04-28(64-b): c(0.03799075, 0.9965306, 1.005554, 0.9997679)
## 2006-04-28(32-b):
## 2006-04-28(64-b): c(0.03842982, 0.996453,  1.005504, 1.000424)
## 2006-04-28(32-b): c(0.03844032, 0.9964511, 1.005502, 1.00044)
dput(signif(sd <- unname(coef(sm4)[, "Std. Error"]), 7))
## 2006-05-12(64-b): c(0.02519034, 0.002868783, 0.02602468, 0.02628035)
## 2006-05-12(64-b): c(0.02546277, 0.002883934, 0.02635756, 0.02667132)#<.. check
## 2006-05-12(32-b):
## 2006-04-28(64-b): c(0.02517215, 0.002867516, 0.02600593, 0.0262572)
## 2006-04-28(32-b): c(0.02547035, 0.002884445, 0.02636594, 0.0266819)
stopifnot(identical(coef(m3), coef(m4)),
	  all.equal(cf, c(0.0382, 0.9965, 1.00555, 0.9998), tol= 5e-4),
	  all.equal(sd, c(0.025321, 0.002876, 0.026186, 0.02647), tol=0.01)
	  )


## rm(a,m1, m2, m3, m4, sm2, l1)



cat('Time elapsed: ', proc.time(),'\n') # "stats"
