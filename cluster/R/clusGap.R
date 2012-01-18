#### Originally from orphaned package SLmisc
#### (Version: 1.4.1, 2007-04-12, Maintainer: Matthias Kohl <kohl@sirs-lab.com>)
#### License: GPL (version 2 or later)
####
#### which said
####  "function corresponds to function gap in package SAGx"

## MM: SAGx is now in Bioconductor --- 1.10.1{devel} or 1.11.1{release}
##     had gap() *corrected* to re-cluster using FUNcluster
##
## MM: Package 'lga' -- has gap() and lga and robust lga [-> UBC]
##    - it uses  boot() nicely  [2012-01: ORPHANED because  Justin Harrington is amiss]
## MM: renamed arguments, and changed almost everything

clusGap <- function (x, FUNcluster, K.max, B=500, verbose = 1, ...)
{
    stopifnot(is.function(FUNcluster), length(dim(x)) == 2, K.max >= 2,
              (n <- nrow(x)) >= 1, (p <- ncol(x)) >= 1)
    if(B != (B. <- as.integer(B)) || (B <- B.) <= 0)
        stop("'B' has to be a positive integer")

    if(is.data.frame(x))
        x <- as.matrix(x)

    W.k <- function(X, kk) {
        clus <- if(kk > 1) FUNcluster(X, kk, ...)$cluster else rep.int(1L, nrow(X))
        ##                 ---------- =  =       -------- kmeans() has 'cluster'; pam() 'clustering'
        ## FIXME: instead of by(), use tapply(1:n, factor(clustering), function(x) ..)
        ## i.e.   vapply(split(x, clustering), function(x) ..,  0.)
        0.5* sum(by(X, factor(clus),
                    function(x) ##D { cat(sprintf(" inside W.k: (n, p) = (%d, %d)\n", nrow(x),ncol(x)))
                                 sum(dist(x)/nrow(x)) ## MM{FIXED}: had ncol(x) !
                             ##D}
                    ))
    }

    logW <- E.logW <- SE.sim <- numeric(K.max)
    if(verbose) cat("Clustering k = 1,2,..., K.max (= ",K.max,"): .. ", sep='')
    for(k in 1:K.max)
        logW[k] <- log(W.k(x, k))
    if(verbose) cat("done\n")

    xs <- scale(x, center=TRUE, scale=FALSE)
    m.x <- rep(attr(xs,"scaled:center"), each = n)# for back transforming
    V.sx <- svd(xs)$v
    rng.x1 <- apply(xs %*% V.sx, # = transformed(x)
                    2, range)
    logWks <- matrix(0., B, K.max)
    if(verbose) cat("Bootstrapping, b = 1,2,..., B (= ", B,
                    ")  [one \".\" per sample]:\n", sep="")
    for (b in 1:B) {
        z1 <- apply(rng.x1, 2,
                    function(M, nn) runif(nn, min=M[1], max=M[2]),
                    nn=n)
        z <- tcrossprod(z1, V.sx) + m.x # back transformed
        for(k in 1:K.max) {
            logWks[b,k] <- log(W.k(z, k))
        }
        if(verbose) cat(".", if(b %% 50 == 0) paste(b,"\n"))
    }
    if(verbose && (B %% 50 != 0)) cat("",B,"\n")
    E.logW <- colMeans(logWks)
    SE.sim <- sqrt((1 + 1/B) * apply(logWks, 2, var))
    structure(class = "clusGap",
              list(Tab = cbind(logW, E.logW, gap = E.logW - logW, SE.sim),
                   ## K.max == nrow(T)
                   n = n, B = B, FUNcluster=FUNcluster))
}

if(FALSE)##-- SAGx/R/gap.r (May 31, 2007 SAGX_1.10-1.tar.gz {development version})
# Tibshirani, Walther and Hastie (2000) #
# Check added 02JUL04;correction to uniform sampling 15MAY06    #
# Tibshirani, Walther and Hastie (2000) #
# Check added 02JUL04;correction to uniform sampling 15MAY06    #
# se added 10SEP06, bug correction 14NOV06
# bug correction 22MAY07
gap <- function(data=swiss,class=g,B=500, cluster.func = myclus){
# data = swiss; class = cl$cluster; B = 100
library(stats)
class.tab <- table(class)
nclus <- length(class.tab)
if (min(class.tab)==1) stop("Singleton clusters not allowed")
if(!(length(class)==nrow(data))) stop("Length of class vector differs from nrow of data")
data <- as.matrix(data)
data <- scale(data, center = TRUE, scale = FALSE)
temp1 <- log(sum(by(data,factor(class),intern <- function(x)sum(dist(x)/ncol(x))/2)))
veigen <- svd(data)$v
x1 <- crossprod(t(data),veigen) # project data to pc-dimensions #
z1 <- matrix(data = NA, nrow = nrow(x1), ncol = ncol(x1))
tots <- vector(length = B)
for (k in 1:B){
   for (j in 1:ncol(x1)){
         min.x <- min(x1[,j]);max.x <- max(x1[,j])
         z1[,j] <- runif(nrow(x1), min = min.x, max = max.x) # x1[sample(1:nrow(x1),nrow(x1),replace=TRUE),j]
      }
z <- crossprod(t(z1),t(veigen))
new.clus <- cluster.func(data = z, k = nclus)
new.class <- new.clus$cluster
tots[k] <- log(sum(by(z, factor(new.class),intern <- function(x) sum(dist(x)/ncol(x))/2)))
}
out <- c(mean(tots)-temp1, sqrt(1+1/B)*sd(tots));names(out) <- c("Gap statistic", "one SE of simulation")
return(out)
}
