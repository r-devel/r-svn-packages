\name{multinom}
\alias{multinom}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{GAM multinomial logistic regression}
\description{Family for use with \code{\link{gam}}, implementing regression for categorical response data. 
Categories must be coded 0 to K, where K is a positive integer. \code{\link{gam}} should be called with a list of K formulae, one for each category except category zero (extra formulae for shared terms may also be supplied). The first formula also specifies the response variable.      
}

\usage{
multinom(K=1)
}
\arguments{
\item{K}{There are K+1 categories and K linear predictors. }

}
\value{
An object of class \code{general.family}.
}

\details{

Use \code{predict} with \code{type="response"} to get the predicted probabilities in each category.
}

%- maybe also `usage' for other objects documented here.

\author{ Simon N. Wood \email{simon.wood@r-project.org}
}

\seealso{\code{\link{ocat}}}

\examples{
library(mgcv)
## Simulate some ordered categorical data...
set.seed(3);n<-400
dat <- gamSim(1,n=n)
dat$f <- dat$f - mean(dat$f)

alpha <- c(-Inf,-1,0,5,Inf)
R <- length(alpha)-1
y <- dat$f
u <- runif(n)
u <- dat$f + log(u/(1-u)) 
for (i in 1:R) {
  y[u > alpha[i]&u <= alpha[i+1]] <- i
}
dat$y <- y

## plot the data...
par(mfrow=c(2,2))
with(dat,plot(x0,y));with(dat,plot(x1,y))
with(dat,plot(x2,y));with(dat,plot(x3,y))

## fit ocat model to data...
b <- gam(y~s(x0)+s(x1)+s(x2)+s(x3),family=ocat(R=R),data=dat)
b
plot(b,pages=1)
gam.check(b)
summary(b)
b$family$getTheta(TRUE) ## the estimated cut points

## predict probabilities of being in each category
predict(b,dat[1:2,],type="response",se=TRUE)
}
\keyword{models} \keyword{regression}%-- one or more ..


