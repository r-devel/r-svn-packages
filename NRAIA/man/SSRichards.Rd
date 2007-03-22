\name{SSRichards}
\title{Richards Growth Model}
\usage{
SSRichards(input, Asym, xmid, scal, lpow)
}
\alias{SSRichards}
\arguments{
 \item{input}{a numeric vector of values at which to evaluate the model.}
 \item{Asym}{a numeric parameter representing the asymptote.}
 \item{xmid}{a numeric parameter representing the \code{x} value at the
   inflection point of the curve.  The value of \code{SSlogis} will be
   \code{Asym/2} at \code{xmid}.}
 \item{scal}{a numeric scale parameter on the \code{input} axis.}
 \item{lpow}{the natural logarithm of the inverse of the power to which
   the denominator is raised.}
}
\description{
  This \code{selfStart} model evaluates the Richards growth model
  function and its gradient.  It has an \code{initial} attribute that
  creates initial estimates of the parameters \code{Asym},
  \code{xmid}, \code{scal} and \code{lpow}.
}
\value{
  a numeric vector of the same length as \code{input}.  It is the value of
  the expression \code{Asym*(1+exp((xmid-input)/scal))^(-exp(-lpow))}.
  If all of the arguments \code{Asym}, \code{xmid}, \code{scal} and
  \code{lpow} are names of objects, the gradient matrix with respect to
  these names is attached as an attribute named \code{gradient}.
}
%\author{Jose Pinheiro and Douglas Bates}
\seealso{\code{\link{nls}}, \code{\link{selfStart}}
}
\examples{
summary(fm1 <- nls(Length ~ SSRichards(Time, Asym, xmid, scal, lpow), Leaves))
}
\keyword{models}
