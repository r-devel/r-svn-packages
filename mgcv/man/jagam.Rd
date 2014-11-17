\name{jagam}
\alias{jagam}
\alias{sim2jam}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Just Another Gibb's Additive Modeller: JAGS support for mgcv.}

\description{Facilities to auto-generate model specification code and associated data to simulate with GAMs in JAGS (or BUGS). This is useful for inference about models with complex random effects structure best coded in JAGS. It is a very innefficient approach to making inferences about standard GAMs. The idea is that template JAGS code is produced for the smooth part of the model, which is the directly edited to include other stochastic components. After simulation with the resulting model, facilities are provided for plotting and prediction with the model smooth components. 
}
\usage{
jagam(formula,family=gaussian,data=list(),file,weights=NULL,na.action,
offset=NULL,knots=NULL,sp=NULL,drop.unused.levels=TRUE,control=gam.control(),
centred=TRUE,sp.prior = "gamma",diagonalize=FALSE)

sim2jam(sam,pregam,edf.type=2)
}
%- maybe also `usage' for other objects documented here.

\arguments{ 
\item{formula}{ A GAM formula (see \code{\link{formula.gam}} and also \code{\link{gam.models}}). 
This is exactly like the formula for a GLM except that smooth terms,  \code{\link{s}}, \code{\link{te}}, \code{\link{ti}} 
and \code{\link{t2}} can be added to the right hand side to specify that the linear predictor depends on smooth functions of predictors (or linear functionals of these).
} 

\item{family}{
This is a family object specifying the distribution and link to use in
fitting etc. See \code{\link{glm}} and \code{\link{family}} for more
details. A negative binomial family is provided: see \code{\link{negbin}}, but
only the known theta case is supported by \code{bam}.
} 

\item{data}{ A data frame or list containing the model response variable and 
covariates required by the formula. By default the variables are taken 
from \code{environment(formula)}: typically the environment from 
which \code{gam} is called.} 

\item{file}{Name of the file to which JAGS model specification code should be written. See \code{\link{setwd}} for setting and querying the current working directory.}

\item{weights}{ prior weights on the data.}

\item{na.action}{ a function which indicates what should happen when the data
          contain `NA's.  The default is set by the `na.action' setting
          of `options', and is `na.fail' if that is unset.  The
          ``factory-fresh'' default is `na.omit'.}

\item{offset}{Can be used to supply a model offset for use in fitting. Note
that this offset will always be completely ignored when predicting, unlike an offset 
included in \code{formula}: this conforms to the behaviour of
\code{lm} and \code{glm}.}


\item{control}{A list of fit control parameters to replace defaults returned by 
\code{\link{gam.control}}. Any control parameters not supplied stay at their default values.}


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

\item{drop.unused.levels}{by default unused levels are dropped from factors before fitting. For some smooths 
involving factor variables you might want to turn this off. Only do so if you know what you are doing.}

\item{centred}{Should centring constraints be applied to the smooths, as is usual with GAMS? Only set 
this to \code{FALSE} if you know exactly what you are doing. If \code{FALSE} there is a (usually global) 
intercept for each smooth.}

\item{sp.prior}{\code{"gamma"} or \code{"log.uniform"} prior for the smoothing parameters. Do check that the 
default parameters are appropriate for your model in the JAGS code.}

\item{diagonalize}{Should smooths be re-parameterized to have i.i.d. Gaussian priors? For Gaussian data this 
allows efficient conjugate samplers to be used, but otherwise it is better to update smoothers blockwise, and not do this.}

\item{sam}{jags sample object, containing at least fields \code{b} (coefficients) and \code{rho} (log 
smoothing parameters). May also contain field \code{mu} containing monitored expected response.}

\item{edf.type}{Since EDF is not uniquely defined and may be affected by the stochastic structure added to the JAGS model file, 3 options are offered. See details.}

}

\value{ 
For \code{jagam} a three item list containing 
\item{pregam}{standard \code{mgcv} GAM setup data.}
\item{jags.data}{list of arguments to be supplied to JAGS containing information referenced in model specification.}
\item{jags.ini}{initialization data for smooth coefficients and smoothing parameters.}

For \code{sim2jam} an object of class \code{"jam"}: a partial version of an \code{mgcv} \code{\link{gamObject}}, 
}


\details{

}


\references{
}

\author{ Simon N. Wood \email{simon.wood@r-project.org}


}

\section{WARNINGS }{
Gibb's sampling is a very slow inferential method for standard GAMs. It is only likely to be worthwhile when complex random effects structures are required above what is possible with direct GAMM methods.

Check the parameters of the priors. 
} 

\seealso{\code{\link{gam}}
 }

\examples{
require(mgcv)
## See help("mgcv-parallel") for using bam in parallel


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






