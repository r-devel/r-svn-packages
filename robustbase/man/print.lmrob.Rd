\name{print.lmrob}
\alias{print.lmrob}
\alias{print.lmrob.null}
\title{Print Method for Objects of Class "lmrob"}
\description{
  Print method for elements of class lmrob
}
\usage{
\method{print}{lmrob}(x, digits = max(3, getOption("digits") - 3), \dots)
\method{print}{lmrob.null}(x, digits = max(3, getOption("digits") - 3), \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an \R object of class \code{lmrob}, typically created by
    \code{\link{lmrob}}.}
  \item{digits}{number of digits for printing, see \code{digits} in
    \code{\link{options}}.}
  \item{\dots}{potentially more arguments passed to methods.}
}
\seealso{\code{\link{lmrob}}, \code{\link{summary.lmrob}},
    \code{\link{print}} and \code{\link{summary}}.
}
\examples{
data(starsCYG, package = "robustbase")
plot(starsCYG)
 lmST <-    lm(log.light ~ log.Te, data = starsCYG)
RlmST <- lmrob(log.light ~ log.Te, data = starsCYG)
RlmST
abline(lmST, col = "red")
abline(RlmST, col = "blue")
}
\keyword{robust}
\keyword{regression}
