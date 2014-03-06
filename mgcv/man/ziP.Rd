\name{ziP}
\alias{ziP}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{GAM zero-inflated Poisson regression family}
\description{Family for use with \code{\link{gam}}, implementing regression for zero inflated Poisson data
when the logit of the zero inflation rate is linearly dependent on the Poisson mean.
}

\usage{
ziP(theta = NULL, link = "log")
}
\arguments{
\item{theta}{the 2 parameters controlling the slope and intercept of the 
linear transform of the mean controlling the zero inflation rate. }

\item{link}{The link function: only the \code{"log"} is currently supported.}

}
\value{
An object of class \code{extended.family}.
}

\details{The probability of a zero count is given by \eqn{p + (1-p)\exp(-\mu)}{p + (1-p)exp(-mu)}, whereas the probability of
count \eqn{y>0}{y>0} is \eqn{(1-p)\exp(-\mu)\mu^y/y!}{(1-p)exp(-mu)mu^y/y!}. The linear predictor of the model controls the 
Poisson mean \eqn{\mu}{\mu}, while \eqn{logit(p) = \theta_1 + \theta_2 \mu }{logit(p) = theta_1 + theta_2 mu}. The \code{theta} 
parameters are estimated alongside the smoothing parameters.
}

%- maybe also `usage' for other objects documented here.

\author{ Natalya Pya (nyp20@bath.ac.uk)
}

\references{
}

\seealso{\code{\link{ziplss}}}

\examples{
rzip <- function(lambda,theta= c(2,2)) {
## generate zero inflated Poisson random variables, according 
## to model above
  if (length(theta)!=2) stop("theta must be of length 2")
  if (sum(lambda<0)) stop("lambda must be non negative")
  if (theta[2]< 0) stop("the second parameter must be non-negative")
  n.sim <- length(lambda)
  p <- exp(theta[1]- theta[2]*lambda)/(1 + exp(theta[1]- theta[2]*lambda)) 

  y <- rep(0,n.sim)
  good <- runif(n.sim) >= p
  y[good] <- rpois(sum(good),lambda[good])
  y
}

library(mgcv)
## Simulate some beta data...
set.seed(1);n<-400
dat <- gamSim(1,n=n)
mu <- exp(dat$f/4-2)
dat$y <- rzip(mu)

b <- gam(y~s(x0)+s(x1)+s(x2)+s(x3),family=ziP(),data=dat)

b
plot(b,pages=1)
}
\keyword{models} \keyword{regression}%-- one or more ..


