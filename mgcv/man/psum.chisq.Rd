\name{psum.chisq}
\alias{psum.chisq}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Evaluate the c.d.f. of a weighted sum of chi-squared deviates}
\description{ Evaluates the c.d.f. of a weighted sum of chi-squared radom variables
by the method of Davies (1980). That is it computes \deqn{P(q< \sum_{i=1}^r \lambda_i X_i + \sigma_z Z)}{P(q < sum_j lb[j] X_j + sigz Z)} where \eqn{X_i}{X_i} is a chi-squared random variable with \code{df[i]} (integer) degrees of freedom and non-centrality parameter \code{nc[i]}, while \eqn{Z}{Z} is a standard normal deviate. 
}
\usage{psum.chisq(q,lb,df=rep(1,length(lb)),nc=rep(0,length(lb)),sigz=0,
                       lower.tail=FALSE,tol=2e-5,nlim=100000,trace=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{q}{ is the vector of quantile values at which to evaluate.}
\item{lb}{ contains \eqn{\lambda_i}{lb[i]}, the weight for deviate code{i}.}
\item{df}{ is the integer vector of chi-squared degrees of freedom.}
\item{nc}{ is the vector of non-centrality parameters for the chi-squared deviates.}
\item{sigz}{ is the multiplier for the standard normal deviate. None positive to exclude this term. }
\item{lower.tail}{ indicates whether lower of upper tail probabilities are required.}
\item{tol}{ is the numerical tolerance to work to.}
\item{nlim}{ is the maximum number of integration steps to allow}
\item{trace}{ can be set to \code{TRUE} to return some trace information and a fault code as attributes.}
}

\details{This calls a C translation of the original Algol60 code from Davies (1980). Some modifications have been made to remove \code{goto} statements and global variables, to use a slightly more efficient ordering of \code{lb} and to use R functions for \code{log(1+x)}. In addition the integral and associate error are accumulated in single terms, rather than each being split into 2, since only their sums are ever used. If \code{q} is a vector then the algorithm has to be called separately for each \code{q[i]}.

If the Davies algorithm returns with an error then an attempt will be made to use the approximation of Liu et al (2009) and a warning will be issued. If that is not possible then an \code{NA} is returned. A warning will also be issued if the algorithm detects that round off errors may be significant.

If \code{trace} is set to \code{TRUE} then the result with have two attributes. \code{"ifault"} is 0 for no problem, 1 if the desired accuracy can not be obtained, 2 if round-off error may be significant, 3 is invalid parameters have been supplied or 4 if integration parameters can not be located. \code{"trace"} is a 7 element vector: 1. absolute value sum; 2. total number of integration terms; 3. number of integrations; 4. integration interval in main integration; 5. truncation point in initial integration; 6. sd of convergence factor term; 7. number of cycles to locate integration parameters. See Davies (1980) for more details. Note that for vector \code{q} these attributes relate to the final element of \code{q}.
}
 
\references{
Davies, R. B. (1980) Algorithm AS 155: The Distribution of a Linear Combination of Chi-squared Random Variables, J. R. Statist. Soc. C
29, 323-333

Liu, H.; Tang, Y. & Zhang, H. H (2009) A new chi-square approximation to the distribution of non-negative definite quadratic forms in non-central normal variables, Computational Statistics & Data Analysis 53,853-856
}

\author{ Simon N. Wood \email{simon.wood@r-project.org}}

\examples{
  require(mgcv)
  lb <- c(4.1,1.2,1e-3,-1)
  df <- c(2,1,1,1)
  psum.chisq(c(1,6,20),lb,df)
}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..









