\name{clusGap}
\title{Gap Statistic for Estimating the Number of Clusters}
\alias{clusGap}
\alias{print.clusGap}
\alias{plot.clusGap}
\description{
  Calculates a goodness of clustering measure, the \dQuote{gap}
  statistic.  For each number of clusters \eqn{k}, it compares
  \eqn{\log(W(k))} with \eqn{E^*[\log(W(k))]} where the latter is
  defined via bootstrapping, i.e. simulating from a reference
  distribution.
}
\usage{
clusGap(x, FUNcluster, K.max, B = 500, verbose=1, \dots)

\S3method{print}{clusGap}(x, SE.factor = 1, \dots)
}
\arguments{
  \item{x}{numeric matrix or \code{\link{data.frame}}.}
  \item{FUNcluster}{a \code{\link{function}} which accepts as first
    argument a (data) matrix like \code{x}, second argument, say
    \eqn{k, k\geq 2}{k, k >= 2}, the number of clusters desired,
    and returns a \code{\link{list}} with a component named (or shortened to)
    \code{cluster} which is a vector of length \code{n = nrow(x)} of
    integers in \code{1:k} determining the clustering or grouping of the
    \code{n} observations.}
  \item{K.max}{the maximum number of clusters to consider, must be at
    least two.}
  \item{B}{integer, number of Monte Carlo (\dQuote{bootstrap}) samples.}
  \item{verbose}{integer or logical, determining if \dQuote{progress}
    output should be printed.  The default prints one bit per bootstrap
    sample.}
  \item{\dots}{optionally further arguments for \code{FUNcluster()}, see
    \code{kmeans} example below.}
  \item{SE.factor}{for the \code{print()} method only: When determining
    the optimal number of clusters, Tibshirani et al. proposed the
    \dQuote{1 S.E.}-rule.  Using an \code{SE.factor} \eqn{f}, the
    \dQuote{f S.E.}-rule is used, more generally.}
}
\details{
  The main result \code{<res>$Tab[,"gap"]} of course is from
  bootstrapping aka Monte Carlo simulation and hence random, or
  equivalently, depending on the initial random seed (see
  \code{\link{set.seed}()}).
  On the other hand, in our experience, using \code{B = 500} gives
  quite precise results such that the gap plot is basically unchanged
  after an another run.
}
\value{
  an object of S3 class \code{"clusGap"}, basically a list with
  components
  \item{Tab}{a matrix with \code{K.max} rows and 4 columns, named
    "logW", "E.logW", "gap", and "SE.sim",
    where \code{gap = E.logW - logW}
    ...... %% FIXME
  }
  \item{n}{number of observations, i.e., \code{nrow(x)}.}
  \item{B}{input \code{B}}
  \item{FUNcluster}{input function \code{FUNcluster}}
}
\references{
  T. Hastie, R. Tibshirani and G. Walther (2001)
  Estimating the number of data clusters via the Gap statistic.
  \emph{Journal of the Royal Statistical Society B}, \bold{63}, 411--423.

  Tibshirani, R., Walther, G. and Hastie, T. (2000).
  Estimating the number of clusters in a dataset via the Gap statistic.
  Technical Report. Stanford.

  Per Broberg (2006). SAGx: Statistical Analysis of the GeneChip.
  R package version 1.9.7.% moved to Bioconductor sometime after 2006
  \url{http://home.swipnet.se/pibroberg/expression_hemsida1.html}
}
\author{
  This function is originally based on the functions \code{gap} of
  (Bioconductor) package \pkg{SAGx} by Per Broberg,
  \code{gapStat()} from former package \pkg{SLmisc} by Matthias Kohl
  and ideas from \code{gap()} and its methods of package \pkg{lga} by
  Justin Harrington.

  The current implementation is by Martin Maechler.
}
\seealso{
  \code{\link{silhouette}} for a much simpler less sophisticated
  goodness of clustering measure.

  \code{\link[fpc]{cluster.stats}()} in package \pkg{fpc} for
  alternative measures.

  %\code{\link[SGAx]{gap}} in Bioconductor package \pkg{SGAx}.
}
\examples{
## ridiculously nicely separated clusters in 3 D :
x <- rbind(matrix(rnorm(150,           sd = 0.1), ncol = 3),
           matrix(rnorm(150, mean = 1, sd = 0.1), ncol = 3),
           matrix(rnorm(150, mean = 2, sd = 0.1), ncol = 3),
           matrix(rnorm(150, mean = 3, sd = 0.1), ncol = 3))

## Note we use  B = 60 in the following examples to keep them "speedy".
## ---- rather keep the default B = 500 for your analysis!

## note we can  pass 'nstart = 20' to kmeans() :
gskmn <- clusGap(x, FUN = kmeans, nstart = 20, K.max = 8, B = 60)
gskmn #-> its print() method
plot(gskmn, main = "clusGap(., FUN = kmeans, n.start=20, B= 60)")
set.seed(12)
system.time(
  gsPam0 <- clusGap(x, FUN = pam, K.max = 8, B = 60)
)
## Slightly faster way to use pam:
pam1 <- function(x,k) list(cluster = pam(x,k, cluster.only=TRUE))
set.seed(12)
system.time(
  gsPam1 <- clusGap(x, FUN = pam1, K.max = 8, B = 60)
)
## and show that it gives the same:
stopifnot(identical(gsPam1[-4], gsPam0[-4]))
gsPam1

gs.pam.RU <- clusGap(ruspini, FUN = pam1, K.max = 8, B = 60)
gs.pam.RU
plot(gs.pam.RU, main = "Gap statistic for the 'ruspini' data")
mtext("k = 4 is best .. and  k = 5  pretty close")

\donttest{## This takes a minute..
## No clustering ==> k = 1 ("one cluster") should be optimal:
Z <- matrix(rnorm(256*3), 256,3)
gsP.Z <- clusGap(Z, FUN = pam1, K.max = 8, B = 200)
plot(gsP.Z)
gsP.Z
}% dont
}
\keyword{cluster}
