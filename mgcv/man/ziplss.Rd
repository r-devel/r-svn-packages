\name{ziplss}
\alias{ziplss}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Zero inflated Poisson location-scale model family}
\description{The \code{ziplss} family implements a zero inflated Poisson model in which one linear predictor 
controls the probability of `potential presence' and the other controls the mean given `potential presence'.
For example, in animal abundance studies `potential presence' might mean `habitat suitability', with the number of animals observed being given by a Poisson random variable, with some mean, if the habitat is suitable. (The slightly involved description is because the model has two ways of generating zeroes/absences: Poisson zeroes, and `unsuitability zeroes'.) Useable 
only with \code{\link{gam}}, the linear predictors are specified via a list of formulae.
}

\usage{
ziplss(link=list("log","logit"))
}
\arguments{
\item{link}{two item list specifying the link for the Poisson mean and the probability of `potential presence'.}
}
\value{
 An object inheriting from class \code{general.family}.
}

\details{Used with \code{\link{gam}} to fit 2 stage zero inflated Poisson models. \code{gam} is called with 
a list containing 2 formulae, the first specifies the response on the left hand side and the structure of the linear predictor for the Poisson mean on the right hand side. The second is one sided, specifying the linear predictor for the probability of potential presence on the right hand side. 

Link functions \code{"log"}, \code{"identity"}  and \code{"sqrt"} are available for the Poisson mean. 
For the probability of potential presence the \code{"logit"} and \code{"probit"} links can be used.

The fitted values for this family will be a two column matrix. The first column is the Poisson mean, 
and the second column is the probability of potential presence. 
Predictions using \code{\link{predict.gam}} will also produce 2 column matrices for \code{type} 
\code{"link"} and \code{"response"}.

}

\references{

}


\examples{
library(mgcv)
## simulate some data...
f0 <- function(x) 2 * sin(pi * x); f1 <- function(x) exp(2 * x)
f2 <- function(x) 0.2 * x^11 * (10 * (1 - x))^6 + 10 * 
            (10 * x)^3 * (1 - x)^10
n <- 500;set.seed(5)
x0 <- runif(n); x1 <- runif(n)
x2 <- runif(n); x3 <- runif(n)

## Simulate probability of potential presence...
eta1 <- f0(x0) + f1(x1) - 3
p <- binomial()$linkinv(eta1) 
y <- as.numeric(runif(n)<p) ## 1 for presence, 0 for absence

## Simulate y given potentially present...
ind <- y>0
eta2 <- f2(x2[ind])/3
y[ind] <- rpois(exp(eta2),exp(eta2))

## Fit ZIP model... 
b <- gam(list(y~s(x2)+s(x3),~s(x0)+s(x1)),family=ziplss())

summary(b) 
plot(b,pages=1)
}
\keyword{models} \keyword{regression}%-- one or more ..


