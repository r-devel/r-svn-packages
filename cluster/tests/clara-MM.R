### Clara with really LARGE data set --- CPU testing
library(cluster)


## generate 15*N objects, divided into 2 clusters.
N <- 10000
## For back-compatibility:
if(R.version$major != "1" || as.numeric(R.version$minor) >= 7) RNGversion("1.6")
set.seed(521)
x <- rbind(cbind(rnorm(7*N, 0,8), rnorm(7*N, 0,8)),
           cbind(rnorm(8*N,50,8), rnorm(8*N,10,8)))
.proctime00 <- proc.time()
for(nn in 1:3)
    print(clara2 <- clara(x[sample(nrow(x)),], 2, sampsize = 128, samples= 500))
cat('Time elapsed: ', proc.time() - .proctime00,'\n')
## nb-mm, with cluster 1.6-4 :  85.2 (elapsed 91.7)

###-- Larger example: 20*N objects, divided into 5 clusters.
N <- 20000 # 100'000 would give swapping like crazy for 256 MB RAM
x5 <- rbind(cbind(rnorm(4*N, 0,4), rnorm(4*N, 0,4)),
            cbind(rnorm(4*N,10,8), rnorm(4*N,40,6)),
            cbind(rnorm(4*N,30,4), rnorm(4*N, 0,4)),
            cbind(rnorm(4*N,40,4), rnorm(4*N,20,2)),
            cbind(rnorm(4*N,50,4), rnorm(4*N,50,4)))
## plus 1 random dimension
x5 <- cbind(x5, rnorm(nrow(x5)))

.proctime00 <- proc.time()
for(nn in 1:3)
    print(clara(x5[sample(nrow(x5)),], 5, samples= 100))
## Last Line:
cat("Time elapsed: ", proc.time() - .proctime00,"\n")
## nb-mm, with cluster 1.6-4 : 74.4 (elapsed 88.2)
## ~~~~~  1.7-* quite faster : 67.6 (elapsed 68.7)
