\name{bam}
\alias{bam}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{Generalized additive models for very large datasets}

\description{ Fits a generalized additive model (GAM) to a very large
data set, the term `GAM' being taken to include any quadratically penalized GLM.  
The degree of smoothness of model terms is estimated as part of
fitting. In use the function is much like \code{\link{gam}}, except that the numerical methods
are designed for datasets containing upwards of several 10's of thousands of data.

}
\usage{
bam(formula,family=gaussian(),data=list(),weights=NULL,subset=NULL,na.action,
                offset=NULL,method="REML",control=gam.control(),scale=0,gamma=1,knots=NULL,
                sp=NULL,min.sp=NULL,paraPen=NULL,freq.cov=FALSE,chunk.size=10000,...)
}
%- maybe also `usage' for other objects documented here.

\arguments{ 
\item{formula}{ A GAM formula (see \code{\link{formula.gam}} and also \code{\link{gam.models}}). 
This is exactly like the formula for a GLM except that smooth terms, \code{s} and \code{te} can be added 
to the right hand side to specify that the linear predictor depends on smooth functions of predictors 
(or linear functionals of these).
} 

\item{family}{
This is a family object specifying the distribution and link to use in
fitting etc. See \code{\link{glm}} and \code{\link{family}} for more
details. A negative binomial family is provided: see \code{\link{negbin}}.
} 

\item{data}{ A data frame or list containing the model response variable and 
covariates required by the formula. By default the variables are taken 
from \code{environment(formula)}: typically the environment from 
which \code{gam} is called.} 

\item{weights}{ prior weights on the data.}

\item{subset}{ an optional vector specifying a subset of observations to be
          used in the fitting process.}

\item{na.action}{ a function which indicates what should happen when the data
          contain `NA's.  The default is set by the `na.action' setting
          of `options', and is `na.fail' if that is unset.  The
          ``factory-fresh'' default is `na.omit'.}

\item{offset}{Can be used to supply a model offset for use in fitting. Note
that this offset will always be completely ignored when predicting, unlike an offset 
included in \code{formula}: this conforms to the behaviour of
\code{lm} and \code{glm}.}



\item{method}{The smoothing parameter estimation method. \code{"GCV.Cp"} to use GCV for unknown scale parameter and
Mallows' Cp/UBRE/AIC for known scale. \code{"GACV.Cp"} is equivalent, but using GACV in place of GCV. \code{"REML"} 
for REML estimation, including of unknown scale, \code{"P-REML"} for REML estimation, but using a Pearson estimate 
of the scale. \code{"ML"} and \code{"P-ML"} are similar, but using maximum likelihood in place of REML. }

\item{control}{A list of fit control parameters returned by 
\code{\link{gam.control}}.}


\item{scale}{ If this is positive then it is taken as the known scale parameter. Negative signals that the 
scale paraemter is unknown. 0 signals that the scale parameter is 1  for Poisson and binomial and unknown otherwise. 
Note that (RE)ML methods can only work with scale parameter 1 for the Poisson and binomial cases.    
} 

\item{gamma}{It is sometimes useful to inflate the model degrees of 
freedom in the GCV or UBRE/AIC score by a constant multiplier. This allows 
such a multiplier to be supplied. } 


\item{knots}{this is an optional list containing user specified knot values to be used for basis construction. 
For most bases the user simply supplies the knots to be used, which must match up with the \code{k} value
supplied (note that the number of knots is not always just \code{k}). 
See \code{\link{tprs}} for what happens in the \code{"tp"/"ts"} case. 
Different terms can use different numbers of knots, unless they share a covariate.
}

\item{sp}{A vector of smoothing parameters can be provided here.
 Smoothing parameters must be supplied in the order that the smooth terms appear in the model 
formula. Negative elements indicate that the parameter should be estimated, and hence a mixture 
of fixed and estimated parameters is possible. If smooths share smoothing parameters then \code{length(sp)} 
must correspond to the number of underlying smoothing parameters.}

\item{min.sp}{Lower bounds can be supplied for the smoothing parameters. Note
that if this option is used then the smoothing parameters \code{full.sp}, in the 
returned object, will need to be added to what is supplied here to get the 
 smoothing parameters actually multiplying the penalties. \code{length(min.sp)} should 
always be the same as the total number of penalties (so it may be longer than \code{sp},
if smooths share smoothing parameters).}


\item{paraPen}{optional list specifying any penalties to be applied to parametric model terms. 
\code{\link{gam.models}} explains more.}

\item{chunk.size}{The model matrix is created in chunks of this size, rather than ever being formed whole.}


\item{...}{further arguments for 
passing on e.g. to \code{gam.fit} (such as \code{mustart}). }

}

\value{ 
An object of class \code{"gam"} as described in \code{\link{gamObject}}.
}


\details{ 

}


\references{


\url{http://www.maths.bath.ac.uk/~sw283/}
}
\author{ Simon N. Wood \email{simon.wood@r-project.org}


}

\section{WARNINGS }{

You must have more unique combinations of covariates than the model has total
parameters. (Total parameters is sum of basis dimensions plus sum of non-spline 
terms less the number of spline terms). 

This routine is less stable than `gam' for a the same dataset.

} 

\seealso{\code{\link{mgcv-package}}, \code{\link{gamObject}}, \code{\link{gam.models}}, \code{\link{smooth.terms}},
\code{\link{linear.functional.terms}}, \code{\link{s}},
\code{\link{te}} \code{\link{predict.gam}},
\code{\link{plot.gam}}, \code{\link{summary.gam}}, \code{\link{gam.side}},
\code{\link{gam.selection}},\code{\link{mgcv}}, \code{\link{gam.control}}
\code{\link{gam.check}}, \code{\link{linear.functional.terms}} \code{\link{negbin}}, \code{\link{magic}},\code{\link{vis.gam}}
 }

\examples{
## to come
}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..
\concept{Varying coefficient model}
\concept{Functional linear model}
\concept{Penalized GLM}
\concept{Generalized Additive Model}
\concept{Penalized regression}
\concept{Spline smoothing}
\concept{Penalized regression spline}
\concept{Generalized Cross Validation}
\concept{Smoothing parameter selection}
\concept{tensor product smoothing}
\concept{thin plate spline}
\concept{P-spline}
\concept{Generalized ridge regression}






