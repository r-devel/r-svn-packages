
#### Testing  medcouple  and related functions

library(robustbase)
source(system.file("mcnaive.R", package = "robustbase"))# mcNaive()

system.time(stopifnot(0 == sapply(1:100, function(n) mc(seq_len(n)))))
system.time(stopifnot(0 == sapply(1:100, function(n) mc(seq_len(n), doRefl=FALSE))))


## This is somewhat interesting {diff the output !}
## particularly since *most* give  the  'not found' diagnostic
set.seed(17)
for(n in 1:100) {
    cat(sprintf("n =%3d:\n------\n", n))
    mcval <- mc(rlnorm(n), trace=TRUE, doRefl=FALSE)
    cat(sprintf(" --> mc(rlnorm(%d)) = %g\n", n, mcval))
}

system.time(
stopifnot(0 == sapply(1:100, function(n) mcNaive(seq_len(n), "simple")))
)
system.time(
stopifnot(0 == sapply(1:100, function(n) mcNaive(seq_len(n), "h.use")))
)
## 0.26 (on nb-mm; just slightly slower than "simple"

x1 <- c(1, 2, 7, 9, 10)
mcNaive(x1) # = -1/3
stopifnot(all.equal(-1/3, mcNaive(x1), tol=1e-12),
          all.equal(-1/3, mcNaive(x1, "h.use"), tol=1e-12),
          all.equal(-1/3, mc(x1), tol = 1e-12))

x2 <- c(-1, 0, 0, 0, 1, 2)
mcNaive(x2, meth="simple") # = 0 - which is wrong
mcNaive(x2, meth="h.use")  # = 1/6 = 0.16666
stopifnot(all.equal(1/6, mc(x2),               tol=1e-12),
          all.equal(1/6, mcNaive(x2, "h.use"), tol=1e-12))

x3 <- c(-2, rep(-1,4), rep(0,6), 2, 2, 2:4)
mcNaive(x3,"h.use") # 1/3
mcNaive(x3,"simple")#  0

try( mc(x3, doRefl = FALSE, maxit = 15, trace = 3)) ## "non-convergence"
str(robustbase:::mcComp(-x3, doRefl = FALSE, maxit = 15, trace = 3))

set.seed(17)
for(n in 3:50) {
    cat(" ")
    x <- rlnorm(n)
    mc1 <- mc(x)
    mc2 <- mcNaive(x)
    mc3 <- mcNaive(x, method = "h.use")
    stopifnot(all.equal(mc1, mc3, tol = 1e-12),
              mc2 == mc3)
## BUG: Fails (as 'x3' above)
##     x <- round(2 * rnorm(n)) # many ties, often at median -- not converging
##     mc1 <- mc(x)
##     mc2 <- mcNaive(x)
##     mc3 <- mcNaive(x, method = "h.use")
##     stopifnot(mc1 == mc3)
##     if(mc2 != mc3) cat("D!")
    cat(".")
};  cat("\n")
