if(FALSE)
    library(cluster, lib.loc="/u/maechler/R/Pkgs/cluster.Rcheck")
library(cluster)

x <- cbind(c(0, -4, -22, -14, 0, NA, -28, 1, 10, -1,
             100 + c(13, 0, 2, 4, 7, 8, 1)),
           c(-5, -14, NA, -35, -30, NA, 7, 2, -18, 13,
             47, 64, 48, NA, NA, 44, 65))
x
(d <- dist(x,'manhattan'))
summary(d, na.rm = TRUE) # max = 270
## First call with "trace" (seg.fault typically later ...):
try( clara(x, k=2, metric="manhattan", samples=2, sampsize=10, trace = 3) )
## Originally:already shows the problem:  nbest[] = c(0,0,...,0) must be WRONG!!
## Now: gives the proper error message.

## S-plus 6.1.2 (rel.2 for Linux, 2002) gives
##> cc <- clara(x, k=2, metric="manhattan", samples=2, sampsize=10)
## Problem in .Fortran("clara",: Internal error: data for decrementing
## ref. count didn't point to a valid arena (0x0), while calling subroutine clara
