\name{plotfit}
\alias{plotfit}
\alias{plotfit.nls}
\title{Plot x-y data and a fitted model}
\description{
  Extract the original data from a fitted model and plot these data
  along with a smooth curve of the fitted model function.  The fitted
  model must be for one covariate only.
}
\usage{
plotfit(fm, \dots)
\method{plotfit}{nls}(fm, \dots)
}
\arguments{
  \item{fm}{A fitted model object, typically a \code{nls} fitted model.}
  \item{\dots}{Arguments to be passed to \code{\link[lattice]{xyplot}}
    or \code{\link[lattice]{panel.curve}}.}
}
%\details{}
\value{
  A lattice object.
}
%\references{}
%\author{}
%\note{}
\seealso{\code{\link[lattice]{xyplot}}, \code{\link[lattice]{panel.curve}}}
\examples{
fm1 <- nls(rate ~ conc/(K + conc), Puromycin, c(K = 0.05),
           subset = state == "treated", alg = "plinear")
plotfit(fm1, from = 0)
}
\keyword{hplot}
\keyword{models}
\keyword{regression}
\keyword{nonlinear}
