
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
## 0.26 (on nb-mm; just slightly slower than "simple")
## 0.18 (on lynne[2007])

allEQ <- function(x,y) all.equal(x,y, tol = 1e-12)

x1 <- c(1, 2, 7, 9, 10)
mcNaive(x1) # = -1/3
stopifnot(allEQ(-1/3, mcNaive(x1)),
          allEQ(-1/3, mcNaive(x1, "h.use")),
          allEQ(-1/3, mc(x1)))

x2 <- c(-1, 0, 0, 0, 1, 2)
mcNaive(x2, meth="simple") # = 0 - which is wrong
mcNaive(x2, meth="h.use")  # = 1/6 = 0.16666
stopifnot(allEQ(1/6, mc(x2)),
          allEQ(1/6, mcNaive(x2, "h.use")))

x3 <- c(-2, rep(-1,4), rep(0,6), 2, 2, 2:4)
mcNaive(x3,"h.use") # 1/3
mcNaive(x3,"simple")#  0

try( mc(x3, doRefl = FALSE, maxit = 15, trace = 3)) ## "non-convergence" (32-bit)
str(robustbase:::mcComp(-x3, doRefl = FALSE, maxit = 15, trace = 4))

### And here is the "real" problem of the whole 'eps' idea:
x4 <- c(1:5,7,10,15,25, 1e15)
mcNaive(x4,"h.use") # 0.5833333
mcNaive(x4,"simple")# == " == 7/12
try( mc(x4) )# not converged  !!
str(robustbase:::mcComp( x4, doRefl= FALSE, maxit = 15, trace= 3))## = 0: conv.quickly
str(robustbase:::mcComp(-x4, doRefl= FALSE, maxit = 15, trace= 3)) # *not* conv!

## a much more extreme eps seems the cure:
str(robustbase:::mcComp( x4, doRefl= FALSE, eps1=.Machine$double.xmin))
str(robustbase:::mcComp(-x4, doRefl= FALSE, eps1=.Machine$double.xmin))

stopifnot(0 == sapply(1:100, function(n)
          mc(seq_len(n), doRefl=FALSE, eps1=.Machine$double.xmin)))

### Examples "like x3" (non-convergence on 32-bit)
xClist <- list(## length 5 :
               c(0,0, 1, 3,3),
               c(0,0, 1, 3:4),
               ##
               ## length 6 :
               c(0,0, 2, 4:6),    c(0,0, 2, 3, 4, 6), c(0,0, 4, 5, 7, 8),
               c(0, 1,1, 2, 6,6), c(0, 3,3, 4, 6,6),
               c(0,0, 1, 3, 5,5),
               c(0,0, 1, 4,4, 6), c(0,0, 1, 4,4, 7),  c(0,0, 1, 5,5, 6)
               )

rlis <- lapply(xClist, function(x)
               try(mc(x, maxit=9, eps1=.Machine$double.xmin), silent=TRUE))
table(sapply(rlis, class))
if(R.version$arch == "x86_64") {
    print(unlist(rlis))
    rl2 <- lapply(xClist, mc, maxit=9, eps1= 1e-10)
    stopifnot(allEQ(rlis, rl2),
              allEQ(unlist(rlis), sapply(xClist, mcNaive)))
}


set.seed(17)
for(n in 3:50) {
    cat(" ")
    x <- rlnorm(n)
    mc1 <- mc(x)
    mc2 <- mcNaive(x)
    mc3 <- mcNaive(x, method = "h.use")
    stopifnot(allEQ(mc1, mc3),
              mc2 == mc3)
    x <- round(2 * rnorm(n)) # many ties, often at median -- not converging
    if(R.version$arch == "x86_64") {
        ## non-convergence BUG  rarely and only on 32-bit
        mc1 <- mc(x, eps1 = .Machine$double.xmin)
        mc2 <- mcNaive(x, method = "simple")
        mc3 <- mcNaive(x, method = "h.use")
        stopifnot(allEQ(mc1, mc3))
        if(mc2 != mc3) {
            cat("d"); if(!isTRUE(allEQ(mc2, mc3))) cat("!!")
        }
    }
    cat(".")
};  cat("\n")

quit('no')
##  ------

## Find short example of non-convergence (32-bit)
n <- 9
for(ii in 1:100) {
    x <- round(2 * rnorm(n)) # many ties, often at median -- not converging
    mc1 <- mc(x, eps1 = .Machine$double.xmin)
}
##
x5 <- c(-3, -3, -2, -1, -1, 0, 0, 1, 2, 2, 3, 4)
x6 <- c(-5, -2, -1, -1, -1, 0, 0, 0, 2, 2, 2, 4)

x7 <- c(-2, -2, -2, -1, -1, 1, 1, 1, 3)
x8 <- c(-3, -1, -1,  0,  1, 2, 2, 2, 2)
