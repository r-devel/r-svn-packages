\name{ziP}
\alias{ziP}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{GAM zero-inflated Poisson regression family}
\description{Family for use with \code{\link{gam}}, implementing regression for zero inflated Poisson data
when the logit of the zero inflation rate is linearly dependent on the log of the Poisson mean. Use with care, 
it is relatively easy to produce models with an indefinite log likelihood using this family, by the very fact of 
providing two ways of generating a zero.
}

\usage{
ziP(theta = NULL, link = "identity")
}
\arguments{
\item{theta}{the 2 parameters controlling the slope and intercept of the 
linear transform of the mean controlling the zero inflation rate. If supplied (and second element is positive) 
then treated as fixed parameters (\eqn{\theta_1}{theta_1} and \eqn{\exp(\theta_2)}{exp(theta_2)}), otherwise estimated. If supplied and second element is negative then treated as starting values (with sign of second element reversed).}

\item{link}{The link function: only the \code{"identity"} is currently supported.}

}
\value{
An object of class \code{extended.family}.
}

\details{The probability of a zero count is given by \eqn{1-p + p \exp(-\mu)}{1- p + pexp(-mu)}, whereas the probability of
count \eqn{y>0}{y>0} is \eqn{p\exp(-\mu)\mu^y/y!}{(pexp(-mu)mu^y/y!}. The linear predictor of the model controls the 
log Poisson mean \eqn{\log \mu}{log(mu)}, while \eqn{logit(p) = \theta_1 - \exp(\theta_2) \log \mu }{logit(p) = theta_1 - exp(theta_2) log(mu)}. The \code{theta} parameters are estimated alongside the smoothing parameters.

For data with large areas of covariate space over which the response is zero it is advisable to use low order penalties to 
avoid problems. For 1D smooths uses e.g. \code{s(x,m=1)} and for isotropic smooths use \code{\link{Duchon.spline}}s in place of thin plaste terms with order 1 penalties, e.g \code{s(x,z,m=c(1,.5))} --- such smooths penalize towards constants, thereby avoiding extreme estimates when the data are uninformative. Uninformativeness arises because zero observations contain little or no information about the Poisson mean when the zero inflation is high, and contain only rather limited information when there is any zero inflation at all. 

Note that the fitted values for this model are the log of the Poisson mean. Use the \code{predict} function with \code{type=="response"} to get the predicted expected response.

}

%- maybe also `usage' for other objects documented here.

\author{ Simon N. Wood \email{simon.wood@r-project.org}
}

\section{WARNINGS }{
Data for which the response is zero (or binary) over substantial areas of the covariate space lead to identifiability problems that are very
difficult to handle automatically. It is possible to get parameter divergence in such cases which make it impossible to 
compute the likelihood and its derivatives accurately. Use of first order penalties on the smooths can help to regularize matters. 

Zero inflated models are often over-used. Having lots of zeroes in the data does not in itself imply zero inflation. Having too many zeroes *given the model mean* may imply zero inflation. 
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
b$outer.info ## check convergence!!
b
plot(b,pages=1)
## plot deviance residuals against expected response
plot(predict(b,type="response"),residuals(b))

}
\keyword{models} \keyword{regression}%-- one or more ..


