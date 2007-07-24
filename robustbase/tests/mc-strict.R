
#### Testing  medcouple	 and related functions

### here, we do "strict tests" -- hence no *.Rout.save
### hence, can also produce non-reproducible output such as timing

library(robustbase)
source(system.file("mcnaive.R", package = "robustbase"))# mcNaive()

allEQ <- function(x,y) all.equal(x,y, tol = 1e-12)
DO <- function(...) system.time(stopifnot(...))

DO(0 == sapply(1:100, function(n) mc(seq_len(n))))
DO(0 == sapply(1:100, function(n) mc(seq_len(n), doRefl=FALSE)))


DO(0 == sapply(1:100, function(n) mcNaive(seq_len(n), "simple")))
DO(0 == sapply(1:100, function(n) mcNaive(seq_len(n), "h.use" )))


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

x4 <- c(1:5,7,10,15,25, 1e15) ## - bombed in orignal algo
mcNaive(x4,"h.use") # 0.5833333
stopifnot(allEQ( 7/12, mcNaive(x4, "h.use")),
	  allEQ( 7/12, mc( x4, doRefl= FALSE, eps1=.Machine$double.xmin)),
	  allEQ(-7/12, mc(-x4, doRefl= FALSE, eps1=.Machine$double.xmin)))


set.seed(17)
for(n in 3:50) {
    cat(" ")
    for(k in 1:5) {
	x <- rlnorm(n)
	mc1 <- mc(x)
	mc2 <- mcNaive(x, method = "simple")
	mc3 <- mcNaive(x, method = "h.use" )
	stopifnot(all.equal(mc1, mc3, tol = 1e-10),# 1e-12 not quite ok
		  mc2 == mc3)
	cat(".")
    }
};  cat("\n")


DO(0 == sapply(1:100, function(n)
   mc(seq_len(n), doRefl=FALSE, eps1=.Machine$double.xmin)))

###----  Strict tests of adjOutlyingness():


set.seed(1); system.time(a1 <- adjOutlyingness(longley))
set.seed(2); system.time(a2 <- adjOutlyingness(hbk))
set.seed(3); system.time(a3 <- adjOutlyingness(hbk[, 1:3]))# the 'X' space
set.seed(4); system.time(a4 <- adjOutlyingness(milk))
set.seed(5); system.time(a5 <- adjOutlyingness(wood))
set.seed(5); system.time(a5 <- adjOutlyingness(wood[, 1:5]))# the 'X' space

## Not yet strict checking;  the results are currently  *very* doubtful!


### Some platform info :
local({ nms <- names(Si <- Sys.info())
        dropNms <- c("nodename", "machine", "login")
        structure(Si[c("nodename", nms[is.na(match(nms, dropNms))])],
                  class="simple.list") })

if(identical(1L, grep("linux", R.version[["os"]]))) { ##----- Linux - only ----
    ##
    Sys.procinfo <- function(procfile)
    {
        l2 <- strsplit(readLines(procfile),"[ \t]*:[ \t]*")
        r <- sapply(l2[sapply(l2, length) == 2],
                    function(c2)structure(c2[2], names= c2[1]))
        attr(r,"Name") <- procfile
        class(r) <- "simple.list"
        r
    }
    ##
    Scpu <- Sys.procinfo("/proc/cpuinfo")
    Smem <- Sys.procinfo("/proc/meminfo")
    print(Scpu[c("model name", "cpu MHz", "cache size", "bogomips")])
    print(Smem[c("MemTotal", "SwapTotal")])
}
