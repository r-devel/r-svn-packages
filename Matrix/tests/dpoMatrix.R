### Testing positive definite matrices

library(Matrix)

h9 <- Hilbert(9)
stopifnot(c(0,0) == dim(Hilbert(0)),
          c(9,9) == dim(h9))
str(h9)
all.equal(c(determinant(h9)$modulus), -96.7369456, tol= 2e-8)
stopifnot(0 == length(h9@factors))# nothing yet
round(ch9 <- chol(h9), 3) ## round() preserves 'triangular' !
str(f9 <- as(chol(h9), "dtrMatrix"))
## h9 now has factorization
stopifnot(names(h9@factors) == "Cholesky",
          all.equal(rcond(h9), 9.0938e-13),
          all.equal(rcond(f9), 9.1272e-7, tol = 1e-6))# more precision fails
str(h9)# has 'factors'
options(digits=4)
(cf9 <- crossprod(f9))# looks the same as  h9 :
stopifnot(all.equal(as.matrix(h9),
                    as.matrix(cf9), tol= 1e-15))

h9. <- round(h9, 2)# actually loses pos.def. "slightly"
h9.p <- as(h9., "dppMatrix")
h4  <- h9.[1:4, 1:4] # this and the next
h9.[1,1] <- 10       # had failed in 0.995-14
h9.p. <- as(h9., "dppMatrix")
h9.p[1,1] <- 10 # failed in 0.995-14

stopifnot(is(h9., "symmetricMatrix"),
          is(h9.p, "symmetricMatrix"),
          is(h4,   "symmetricMatrix"))

h9.p[1,2] <- 99 #-> becomes "dgeMatrix"

str(hp9 <- as(h9, "dppMatrix"))# packed
stopifnot(is(thp9 <- t(hp9), "dppMatrix"))

hs <- as(hp9, "dspMatrix")
hs@x <- 1/hp9@x # is not pos.def. anymore
validObject(hs)
stopifnot(diag(hs) == seq(1, by = 2, length = 9))

s9 <- solve(hp9, seq(nrow(hp9)))
signif(t(s9)/10000, 4)# only rounded numbers are platform-independent
(I9 <- hp9 %*% s9)
m9 <- matrix(1:9, dimnames = list(NULL,NULL))
stopifnot(all.equal(m9, as.matrix(I9), tol = 2e-9))

### Testing nearPD() --- this is partly in  ../man/nearPD.Rd :
pr <- Matrix(c(1,     0.477, 0.644, 0.478, 0.651, 0.826,
               0.477, 1,     0.516, 0.233, 0.682, 0.75,
               0.644, 0.516, 1,     0.599, 0.581, 0.742,
               0.478, 0.233, 0.599, 1,     0.741, 0.8,
               0.651, 0.682, 0.581, 0.741, 1,     0.798,
               0.826, 0.75,  0.742, 0.8,   0.798, 1),
             nrow = 6, ncol = 6)

nL <-
    list(r   = nearPD(pr, conv.tol = 1e-7), # default
         r.1 = nearPD(pr, conv.tol = 1e-7, corr = TRUE),
         rH  = nearPD(pr, conv.tol = 1e-15),
         rH.1= nearPD(pr, conv.tol = 1e-15, corr = TRUE))

sapply(nL, `[`, c("iterations", "normF"))

allnorms <- function(d) sapply(c("1","I","F","M"), function(typ) norm(d, typ))

## "F" and "M" distances are larger for the (corr=TRUE) constrained:
100 * sapply(nL, function(rr) allnorms((pr - rr  $ mat)))

## But indeed, the 'corr = TRUE' constraint yield a better solution,
## if you need the constraint :  cov2cor() does not just fix it up :
100 * (nn <- sapply(nL, function(rr) allnorms((pr - cov2cor(rr  $ mat)))))

stopifnot(
all.equal(nn["1",],
          c(r = 0.099944428698, r.1 =0.087461417994,
            rH= 0.099944428698, rH.1=0.087461430806), tol=1e-9))

nr <- nL $rH.1 $mat
stopifnot(
    all.equal(nr[lower.tri(nr)],
	      c(0.48796803265083, 0.64265188295401, 0.49063868812228, 0.64409905497094,
		0.80871120142824, 0.51411473401472, 0.25066882763262, 0.67235131534931,
		0.72583206922437, 0.59682778611131, 0.58219178154582, 0.7449631866236,
		0.72988206459063, 0.77215024062758, 0.81319175546212), tol = 1e-10))

set.seed(27)
m9 <- h9 + rnorm(9^2)/1000 ; m9 <- (m9 + t(m9))/2
nm9 <- nearPD(m9)

CPU <- 0
for(M in c(5, 12))
    for(i in 1:50) {
	m <- matrix(round(rnorm(M^2),2), M, M)
	m <- m + t(m)
	diag(m) <- pmax(0, diag(m)) + 1
	m <- cov2cor(m)
	CPU <- CPU + system.time(n.m <- nearPD(m))[1]
	X <- as(n.m$mat, "matrix")
	stopifnot(all.equal(X, (X + t(X))/2, tol = 8*.Machine$double.eps),
		  all.equal(eigen(n.m$mat, only.values=TRUE)$values,
			    n.m$eigenvalues, tol = 4e-8))
    }
cat('Time elapsed for nearPD(): ', CPU,'\n')

## cov2cor()
m <- diag(6:1) %*% as(pr,"matrix") %*% diag(6:1) # so we can "vector-index"
m[upper.tri(m)] <- 0
ltm <- which(lower.tri(m))
ne <- length(ltm)
set.seed(17)
m[ltm[sample(ne, 3/4*ne)]] <- 0
m <- (m + t(m))/2 # now is a covariance matrix with many 0 entries
(spr <- Matrix(m))
cspr <- cov2cor(spr)
ev <- eigen(cspr, only.v = TRUE)$values
stopifnot(is(spr, "dsCMatrix"),
          is(cspr,"dsCMatrix"),
          all.equal(ev, c(1.5901626099,  1.1902658504, 1, 1,
                          0.80973414959, 0.40983739006), tol=1e-10))


cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''

