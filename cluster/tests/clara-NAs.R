if(FALSE)
    library(cluster, lib.loc="/u/maechler/R/Pkgs/cluster.Rcheck")
library(cluster)

set.seed(47)
x <- round(rbind(cbind(rnorm(17, 0,3), rnorm(17, 0,3)),
                 cbind(rnorm( 8,50,3), rnorm( 8,10,3))), 2)
x[c(7, 13, 25, 35, 38, 46, 47, 50)] <- NA
x
## call this repeatedly -> seg.fault (at "17" on lynne)
for(N in 1:20){ cat(N,":\n"); cc <- clara(x, 2, sampsize=10)}

