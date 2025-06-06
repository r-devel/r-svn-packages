\name{coxph.object}
\alias{coxph.object}
\alias{extractAIC.coxph.penal}
\alias{print.coxph}
\alias{summary.coxph}
\title{
Proportional Hazards Regression Object
}
\description{
This class of objects is returned by the \code{coxph} class of functions
to represent a fitted proportional hazards model.


Objects of this class have methods for the functions \code{print},
\code{summary}, \code{residuals}, \code{predict} and \code{survfit}.
}
\section{COMPONENTS}{
The following components must be included in a legitimate \code{coxph} object.
\describe{
\item{coefficients}{
the coefficients of the linear predictor, which multiply the columns of the
model matrix.  If the model is over-determined there will be missing
values in the vector corresponding to the redundant columns in the model
matrix.
}
\item{var}{
the variance matrix of the coefficients.  Rows and columns corresponding to
any missing coefficients are set to zero.
}
\item{naive.var}{
this component will be present only if the \code{robust} option was true.  If so,
the \code{var} component will contain the robust estimate of variance, and this
component will contain the ordinary estimate.
}
\item{loglik}{
a vector of length 2 containing the log-likelihood with the initial values and
with the final values of the coefficients.
}
\item{score}{
value of the efficient score test, at the initial value of the coefficients.
}
\item{rscore}{
the robust log-rank statistic, if a robust variance was requested.
}
\item{wald.test}{
the Wald test of whether the final coefficients differ from the initial values.
}
\item{iter}{
number of iterations used.
}
\item{linear.predictors}{
the vector of linear predictors, one per subject.
}
\item{residuals}{
the martingale residuals.
}
\item{means}{
vector of column means of the X matrix.  Subsequent survival curves are
adjusted to this value.
}
\item{n}{
the number of observations used in the fit.
}
\item{weights}{
the vector of case weights, if one was used.
}
\item{method}{
the computation method used.
}
\item{na.action}{
the \code{na.action} attribute, if any, that was returned by the \code{na.action} routine.
}

The object will also contain the following, for documentation see the \code{lm}
object: \code{terms}, \code{assign}, \code{formula}, \code{call}, and, optionally, \code{x}, \code{y},
and/or \code{frame}.
}}
\seealso{
\code{\link{coxph}}, \code{\link{coxph.detail}}, \code{\link{cox.zph}}, \code{\link{survfit}}, \code{\link{residuals.coxph}}, \code{\link{survreg}}
}
\keyword{survival}
% Converted by Sd2Rd version 0.3-2.
