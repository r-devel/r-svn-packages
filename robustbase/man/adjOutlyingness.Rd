\name{adjOutlyingness}
\alias{adjOutlyingness}
\title{Compute Skewness-adjusted Multivariate Outlyingness}
\description{
  For an \eqn{n \times p}{n * p} data matrix (or data frame) \code{x},
  compute the \dQuote{\emph{outlyingness}} of all \eqn{n} observations.
  Outlyingness here is a generalization of the Donoho-Stahel
  outlyingness measure, where skewness is taken into account via the
  medcouple, \code{\link{mc}()}.
}
\usage{
adjOutlyingness(x, ndir = 250, clower = 3, cupper = 4)
}
\arguments{
  \item{x}{a numeric \code{\link{matrix}} or \code{\link{data.frame}}.}
  \item{ndir}{positive integer specifying the number of directions that
    should be searched.}
  \item{clower, cupper}{the constant to be used for the lower and upper
    tails, in order to transform the data towards symmetry.}
%%% MM: How can I modify it to the  *NON*-adjusted case ?
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{

}
\author{Guy Brys; help page and improvements by Martin Maechler}
\seealso{the adjusted boxplot, \code{\link{adjbox}} and the medcouple,
  \code{\link{mc}}.
}
\examples{
dim(longley)
adjOutlyingness(longley)

dim(hbk)
adjOutlyingness(hbk)
}
\keyword{robust}
\keyword{multivariate}
