\name{silhouette}
\alias{silhouette}
\alias{silhouette.default}
\alias{silhouette.partition}
\title{Compute or Extract Silhouette Information from Clustering}
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
silhouette(x, \dots)
silhouette.default  (x, dist, \dots)
silhouette.partition(x, \dots)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{a (clustering) \R objects.}
  \item{dist}{ ~~Describe \code{dist} here~~ }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
}
\details{
  ~~ If necessary, more details than the __description__  above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of `comp1'}
  \item{comp2 }{Description of `comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\seealso{\code{\link{partition.object}}, \code{\link{plot.partition}}.
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	     or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(x, dist, ...) {
    if(missing(dist)) stop("Need a dissimilarity (as 2nd argument)")
    dist <- as.dist(dist) # hopefully
    if(!is.null(cl <- x$clustering)) x <- cl
    n <- length(x)
    if(n != attr(dist, "Size"))
        stop("clustering `x' and dissimilarity `dist' are incompatible")
    if(!all(x == round(x))) stop("`x' must only have integer codes")
    k <- length(clid <- sort(unique(x)))
    if(k <= 1 || k >= n)
        return(NA)
    dmat <- as.matrix(dist)
    wds <- matrix(NA,n,3)
    nms <- character(n)
    c.avg <- numeric(k)
    iw <- 0
    for(j in 1:k) {
        Nj <- sum(iC <- x == clid[j])
        a.i <- colSums(dmat[iC, iC])/(Nj - 1)
        diC <- apply(dmat[!iC, iC], 2, function(r) tapply(r, x[!iC], mean))
        minC <- max.col(-t(diC))
        b.i <- diC[cbind(minC, seq(minC))]
        s.i <- (b.i - a.i) / pmax(b.i, a.i)
        ii <- sort.list(s.i, decreasing = TRUE)
        wds[iw + 1:Nj, ] <- cbind(j, clid[-j][minC], s.i)[ ii, ]
        nms[iw + 1:Nj] <- names(x[iC])[ii]
        iw <- iw + Nj
        c.avg[j] <- mean(s.i)
    }
    dimnames(wds) <- list(nms, c("cluster","neighbor","sil_width"))
    list(widths = wds, clus.avg.widths = c.avg,
         avg.width = mean(wds[,"sil_width"]))
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
