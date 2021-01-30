library(cluster)
##      =========>> This will *NOT* be in final tarball, i.e, only tested by maintainer
                                        #               vvvvvvvvvvvvvvvvvvvvv
sessionInfo()# {for the package maintainer; do here, as don't use *.Rout.save}

## generate 1500 objects, divided into 2 clusters.
set.seed(264)
x <- rbind(cbind(rnorm(700, 0,8), rnorm(700, 0,8)),
           cbind(rnorm(800,50,8), rnorm(800,10,8)))

.proctime00 <- proc.time()

## full size sample {should be = pam()}:
n0 <- length(iSml <- c(1:70, 701:720))
summary(clara0 <- clara(x[iSml,], k = 2, sampsize = n0))
          pam0 <- pam  (x[iSml,], k = 2)
stopifnot(identical(clara0$clustering, pam0$clustering),
          all.equal(clara0$objective, unname(pam0$objective[2]),
                    tol = 100*.Machine$double.eps))

summary(clara2 <- clara(x, 2))

clInd <- c("objective", "medoids", "clusinfo")
clIn4 <- c(clInd, "sample")
## clara() {as original code} always draws the *same* random samples !!!!
clara(x, 2, samples = 50)[clInd]

for(i in 1:20)
    print(clara(x[sample(nrow(x)),], 2, samples = 50)[clInd])

clara(x, 2, samples = 101)[clInd]
clara(x, 2, samples = 149)[clInd]
clara(x, 2, samples = 200)[clInd]
## Note that this last one is practically identical to the slower  pam() one

x[print(sample(length(x), 20))] <- NA
clara(x, 2, samples = 50)[clInd]

###-- Larger example: 2000 objects, divided into 5 clusters.
x5 <- rbind(cbind(rnorm(400, 0,4), rnorm(400, 0,4)),
            cbind(rnorm(400,10,8), rnorm(400,40,6)),
            cbind(rnorm(400,30,4), rnorm(400, 0,4)),
            cbind(rnorm(400,40,4), rnorm(400,20,2)),
            cbind(rnorm(400,50,4), rnorm(400,50,4)))
## plus 1 random dimension
x5 <- cbind(x5, rnorm(nrow(x5)))

(cl5def <- clara(x5, 5))
summary(cl5.50 <- clara(x5, 5, samples = 50)) # instead of default 'samples = 5'
## 3 times many samples of default sampsize = 50
cx5many <- list(clara(x5, 5, samples  = 999)
              , clara(x5, 5, samples = 1000)
              , clara(x5, 5, samples = 1001))

cm <- sapply(cx5many, `[[`, "clustering") # 2000 x 3
stopifnot(cm == cm[,1]) # all columns are the same
## not perfect:
iClus <- cx5many[[1]]$clustering
clTrue.p <- rep(c(1:2,4:3,5), each=400) # <<-- permuted (only 3 <--> 4 ) true clustering
sunflowerplot(iClus ~ clTrue.p)
abline(h=0, v=0, lty=3) ; abline(0,1, col=adjustcolor("blue", .5))
(i.off <- which(iClus != clTrue.p))
cbind(iClus, clTrue.p)[i.off, ]
stopifnot(length(i.off) == 2,
          i.off == c(489, 797))

system.time(cl5F <- clara(x5, 5, sampsize = 2000, samples = 1))# full sample
##=> should be equivalent to pam()
system.time(pam5 <- pam(x5, 5))
stopifnot(identical(cl5F$clustering,       pam5$clustering),
          all.equal(cl5F$objective, unname(pam5$objective[2]),
                    tol = 100*.Machine$double.eps))


## Last Line:
cat('Time elapsed: ', proc.time() - .proctime00,'\n')
## Lynne (P IV, 1.6 GHz):  7.5
## nb-mm (P III,700 MHz): 29.2
