if(FALSE)
    library(cluster, lib.loc="/u/maechler/R/Pkgs/cluster.Rcheck")
library(cluster)

x <- cbind(c(0.05, -0.76, -4.4, -2.77, 0.12,  NA,  -5.48,  0.27, 2.01,
             -0.24, 52.66, 50, 47.92, 50.78, 51.52, 51.69, 50.18),
           c(-0.98, -2.67, NA, -6.97, -5.9, NA, 1.45, 0.43, -3.6,
             2.66, 9.54, 12.85, 8.54,  NA,  NA, 8.93, 13.1))
x
## call this repeatedly -> seg.fault (at "17" on lynne)
for(N in 1:20){ cat(N,":\n"); cc <- clara(x, k=2, samples=2, sampsize=10)}

