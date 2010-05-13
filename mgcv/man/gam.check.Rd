\name{gam.check}
\alias{gam.check}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Some diagnostics for a fitted gam model}
\description{ Takes a fitted \code{gam} object produced by \code{gam()} and produces some diagnostic information
about the fitting procedure and results. The default is to produce 4 residual
plots, and some information about the convergence of the smoothness selection optimization.
}
\usage{
gam.check(b,old.style=FALSE,...)
}
%- maybe also `usage' for other objects documented here.
\arguments{ 
\item{b}{ a fitted \code{gam} object as produced by \code{gam()}.}
\item{old.style}{If you want old fashioned plots, exactly as in Wood, 2006, set to \code{TRUE}.}
\item{...}{extra graphics parameters to pass to plotting functions.}
}

\value{A vector of reference quantiles for the residual distribution, if these can be computed.}

\details{ This function plots 4 standard diagnostic plots, and some other
convergence diagnostics. Usually the 4 plots are various residual plots. The
printed information relates to the optimization used to select smoothing
parameters. For the default optimization methods the information is summarized in a
readable way, but for other optimization methods, whatever is returned by way of
convergence diagnostics is simply printed.

The QQ plot produced is usually created by a call to \code{\link{qq.gam}}, and plots deviance residuals 
against approximate theoretical quantilies of the deviance residual distribution, according to the fitted model. 
If this looks odd then investigate further using \code{\link{qq.gam}}. Note that raw binary residuals contain no 
information useul for model checking. 

}

\references{

Wood S.N. (2006) Generalized Additive Models: An Introduction with R. Chapman
and Hall/CRC Press.


\url{http://www.maths.bath.ac.uk/~sw283/}

}
\author{ Simon N. Wood \email{simon.wood@r-project.org}}


\seealso{ \code{\link{choose.k}},  \code{\link{gam}}, \code{\link{mgcv}}, \code{\link{magic}}}

\examples{
library(mgcv)
set.seed(0)
dat <- gamSim(1,n=200)
b<-gam(y~s(x0)+s(x1)+s(x2)+s(x3),data=dat)
plot(b,pages=1)
gam.check(b,pch=19,cex=.3)
}
\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ...






