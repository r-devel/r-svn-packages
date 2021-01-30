library(cluster)
##      =========>> This will *NOT* be in final tarball, i.e, only tested by maintainer
## generate 1500 objects, divided into 2 clusters.
suppressWarnings(RNGversion("3.5.0")) # << as long as we don't have R >= 3.6.0
set.seed(264)
x <- rbind(cbind(rnorm(700, 0,8), rnorm(700, 0,8)),
           cbind(rnorm(800,50,8), rnorm(800,10,8)))

isEq <- function(x,y, epsF = 100)
    is.logical(r <- all.equal(x,y, tol = epsF * .Machine$double.eps)) && r


.proctime00 <- proc.time()

## full size sample {should be = pam()}:
n0 <- length(iSml <- c(1:70, 701:720))
summary(clara0 <- clara(x[iSml,], k = 2, sampsize = n0))
          pam0 <- pam  (x[iSml,], k = 2)
stopifnot(identical(clara0$clustering, pam0$clustering)
        , isEq(clara0$objective, unname(pam0$objective[2]))
          )

summary(clara2 <- clara(x, 2))

clInd <- c("objective", "i.med", "medoids", "clusinfo")
clInS <- c(clInd, "sample")
## clara() {as original code} always draws the *same* random samples !!!!
clara(x, 2, samples = 50)[clInd]

## The for() code below *did* hang (in some 64-bit versions of R / clara / ...):
## but this has been fixed *long* ago
xx <- x[sample(nrow(x)),]
if(FALSE) { ## only run manually
    ff <- "/u/maechler/R/MM/Pkg-ex/cluster/clara2-ex.rds"
    if(file.exists(ff))
        xx <- readRDS(ff)
    else
        saveRDS(xx, file=ff)
}
print(clara(xx, 2, samples= 5, trace = 3)[clInd])

print(clara(xx, 2, samples=50, trace = 2)[clInd])
##                             ^^^^^^^^^ a bit less output ...
## for sample 45
## From that output I gather the sample indices nsel[] and all k{ran} are
ii <- c(194, 1411,
        1430,398,570,72,27,62,1301,368,1390,991,296,431,1019,186,558,258,413,
        647,585,1352,1073,873,377,711,1498,865,1436,1335,189,622,760,226,146,
        145,1349,382,1368,934,204,303,856,1489)
## Trying
pp <- pam(x[ii,], k =2) # gives no problem

for(i in 1:19)
    print(clara(x[sample(nrow(x)),], 2, samples = 50)[clInd])

clara(x, 2, samples = 101)[clInd]
clara(x, 2, samples = 149)[clInd]
clara(x, 2, samples = 200)[clInd]
## Note that this last one is practically identical to the slower  pam() one

x[print(sample(length(x), 20))] <- NA
clara(x, 2, samples = 50)[clInd]

###--- Start a version of  example(clara) -------

## xclara : artificial data with 3 clusters of 1000 bivariate objects each.
data(xclara)
(clx3 <- clara(xclara, 3))
## Plot similar to Figure 5 in Struyf et al (1996)
plot(clx3)

## The  rngR = TRUE case is currently in the non-strict tests
## ./clara-ex.R
## ~~~~~~~~~~~~

###--- End version of example(clara) -------

## Last Line:
cat('Time elapsed: ', proc.time() - .proctime00,'\n')
