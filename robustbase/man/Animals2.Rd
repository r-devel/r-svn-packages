\name{Animals2}
\alias{Animals2}
\title{Brain and Body Weights for 65 Species of Land Animals}
\description{
  A data frame with average brain and body weights for 62 species
  of land mammals and three others.

  Note that this is simply the union of \code{\link[MASS]{Animals}}
  and \code{\link[MASS]{mammals}}.
}
\usage{
Animals2
}
\format{
  \describe{
    \item{\code{body}}{body weight in kg}
    \item{\code{brain}}{brain weight in g}
  }
}
\source{
  Weisberg, S. (1985)
  \emph{Applied Linear Regression.}
  2nd edition.
  Wiley, pp. 144--5.

  P. J. Rousseeuw  and A. M. Leroy (1987)
  \emph{Robust Regression and Outlier Detection.}
  Wiley, p. 57.
}
\references{
  Venables, W. N. and Ripley, B. D. (2002)
  \emph{Modern Applied Statistics with S.} Forth Edition. Springer.
}
\note{
  After loading the \code{MASS} package, the data set is simply constructed by
  \code{Animals2 <- local({D <- rbind(Animals, mammals); unique(D[order(D$body,D$brain),])})}
}
\examples{
data(Animals2)
## Sensible Plot needs doubly logarithmic scale
plot(Animals2, log = "xy")
}
\keyword{datasets}
