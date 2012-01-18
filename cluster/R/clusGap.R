#### Originally from orphaned package SLmisc
#### (Version: 1.4.1, 2007-04-12, Maintainer: Matthias Kohl <kohl@sirs-lab.com>)
#### License: GPL (version 2 or later)
####
#### which said
####  "function corresponds to function gap in package SAGx"

## MM: SAGx is now in Bioconductor --- 1.10.1{devel} or 1.11.1{release}
##     had gap() *corrected* to re-cluster using FUNcluster --> see ./gap-SAGx.R.~orig~
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
