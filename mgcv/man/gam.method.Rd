\name{gam.method}
\alias{gam.method}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{GAM fitting methods}
\description{Smoothness slection criteria and numerical methods}
\section{Smoothness selection criteria}{

\code{\link{gam}} estimates the smoothing parameters for the terms of a GAM, or other penalized GLM
by optimizing one of several smoothness selection criteria. 

\describe{
\item{GCV}{Generalized Cross Validation. This is the number of data times the deviance over the square of the residual effective 
degrees of freedom. }

\item{UBRE}{ This criteria is used only when the scale parameter is known, as the alternative to GCV or GACV. 
In the additive model case it is Mallows' Cp and it is also equivalent to a scaled AIC 
(given that the scale parameter is known). This is basically the deviance plus twice the effective degrees of freedom. }

\item{GACV}{Generalized Approximate Cross Validation. Essentially a rescaled AIC, but with a Pearson estimator used for 
the scale parameter. This is easier to justify in the Generalized case than AIC or UBRE, but both are asymptotically 
equivalent to GACV.}

\item{REML}{Restricted Maximum Likelihood. If the smooths are treated as random effects (with partially improper 
distributions), then they can be integrated out of the likelihood to obtain a REML score to be maximized to 
find smoothing parameters. In the generalized case a Laplace approximation to the integral is used. }

\item{ML}{Maximum Likelihood. As REML, but only the wiggly components of the smooths are treated as random effects 
to be integrated out, resulting in smoother estimates.}

}

The smoothness selection criteria is controlled by the \code{method} argument of \code{\link{gam}}.

}

\section{Optimization methods}{

The numerical optimization method used to optimize the smoothness selection criterion is chosen using the 
\code{optimizer} argument to \code{\link{gam}}. 

 The default methods used by \code{gam} are based on Newton type optimization of
GCV/UBRE/AIC/REML scores with respect to smoothing parameters, as described in Wood (2004) 
and Wood (2008) for the additive and generalized additive model cases respectively. 
The smoothing criteria (GCV etc) are evaluated for the fitted model itself (rather than some working 
approximate model). Since this involves optimizing the criteria `outside' the PIRLS or LS method used 
for fitting, it is referred to as `outer' iteration. 

In the generalized case several alternative optimisation methods can be used for outer
optimization. Usually the fastest and most
reliable approach is to use a modified Newton optimizer with exact first and
second derivatives, and this is the default. However if there are large numbers of smoothing parameters
then it can be faster to use a hybrid Newton-BFGS quasi-Newton approach, which avoids the expense of 
frequent full second derivative evaluation.

 \code{nlm} can be used with finite differenced first derivatives. This is not ideal theoretically, since 
it is possible for the finite difference estimates of derivatives to be very
badly in error on rare
occasions when the P-IRLS convergence tolerance is close to being matched
exactly, so that two components of a finite differenced derivative require
different numbers of iterations of P-IRLS in their evaluation. An alternative
is provided in which \code{nlm} uses numerically exact first derivatives, this
is faster and less problematic than the other scheme. A further alternative is to use a quasi-Newton
scheme with exact derivtives, based on \code{optim}. In practice this usually
seems to be slower than the \code{nlm} method. 


The alternative approach of `performance oriented iteration' was suggested by Gu 
(and is rather similar to the PQL method in generalized linear mixed modelling). At each step of the
P-IRLS (penalized iteratively reweighted least squares) iteration, 
by which a gam is fitted, the smoothing parameters are
estimated by GCV or UBRE applied to the working penalized linear modelling
problem. In most cases, this process converges and gives smoothness estimates
that perform well (see e.g. Wood 2004). 

The performance iteration has two disadvantages. (i) in the presence of
co-linearity or concurvity (a frequent problem when spatial smoothers are
included in a model with other covariates) then the process can fail to
converge. Suppose we start with some coefficient and smoothing parameter
estimates, implying a working penalized linear model: the optimal smoothing
parameters and coefficients for this working model may in turn imply a working 
model for which the original estimates are better than the most recent
estimates. This sort of effect can prevent convergence.

Secondly it is often possible to find a set of smoothing parameters that
result in a lower GCV or UBRE score, for the final working model, than the
final score that results from the performance iterations. This is because the
performance iteration is only approximately optimizing this score (since
optimization is only performed on the working model). The disadvantage here is
not that the model with lower score would perform better (it usually doesn't),
but rather that it makes model comparison on the basis of GCV/UBRE score
rather difficult. 

In summary: performance iteration can fail to converge. It may occasionally be 
faster than outer iteration, but since the Wood (2008) method `outer' iteration 
is faster in most cases.
}


\references{

Gu and Wahba (1991) Minimizing GCV/GML scores with multiple smoothing parameters via
the Newton method. SIAM J. Sci. Statist. Comput. 12:383-398

Wood, S.N. (2004) Stable and efficient multiple smoothing parameter estimation for
generalized additive models. J. Amer. Statist. Ass.

Wood, S.N. (2008) Fast stable direct fitting and smoothness selection for generalized
additive models. J.R.Statist.Soc.B 70(3):495-518

\url{http://www.maths.bath.ac.uk/~sw283/}


}
\author{ Simon N. Wood \email{simon.wood@r-project.org}}


\seealso{\code{\link{gam.control}}   \code{\link{gam}}, \code{\link{gam.fit}}, \code{\link{glm.control}} }

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..




