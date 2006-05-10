\name{lmrob.S}
\alias{lmrob.S}
\title{ S-regression estimators }
\description{
  Computes an S-estimator for linear regression,
  using the \dQuote{fast S} algorithm.% of Matias Salibian & Victor Yohai ...
}
\usage{
lmrob.S(x, y, control, trace.lev = 0)
}
\arguments{
  \item{x}{ design matrix }
  \item{y}{ response vector }
  \item{control}{ list as returned by \code{\link{lmrob.control}} }
  \item{trace.lev}{integer indicating if the progress of the algorithm
    should be traced (increasingly); default \code{trace.lev = 0} does
    no tracing.}
}
\details{
  This function is used by \code{\link{lmrob.fit.MM}} and not
  intended to be used on its own.
}
\value{
	A list containing the following elements:
  \item{coef    }{Vector of S-regression coefficient estimates }
  \item{cov     }{Covariance matrix of the coefficient estimators }
  \item{control }{The control list as given to it }
  \item{scale   }{The S-scale residual estimate }
  \item{seed    }{The value of \code{seed} }
}
\seealso{\code{\link{lmrob}}, also for references.
}
\author{ Matias Salibian-Barrera }
\keyword{robust}
\keyword{regression}
