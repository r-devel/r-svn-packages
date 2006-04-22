\name{lmrob.control}
\alias{lmrob.control}
\title{ Tuning parameters for lmrob }
\description{
  Tuning parameters for, \code{\link{lmrob}}, the MM-regression
  estimator and the associated S-estimator.
}
\usage{
lmrob.control(seed = 37, nResample = 500,
	tuning.chi = 1.54764, bb = 0.5, tuning.psi = 4.685061,
	max.it = 50, groups = 5, n.group = 400,
        best.r.s = 2, k.fast.s = 1, compute.rd = FALSE)
}
\arguments{
  \item{seed}{random seed for the re-samples used in obtaining candiates
    for the initial S-estimator.  The default, \code{37} used to be
    frozen in the underlying C code.}
  \item{nResample}{number of re-sampling candidates to be
    used to find the initial S-estimator.  Currently defaults to 500
    which works well in most situations (see references).}
  \item{tuning.chi}{tuning constant for the S-estimator.
    The default, \code{1.54764}, yields a 50\% breakdown estimator.}
  \item{bb}{expected value under the normal model of the
    "chi" function with tuning constant equal to
    \code{tuning.chi}.  This is used to compute the S-estimator.}
  \item{tuning.psi}{tuning constant for the re-descending M-estimator.
    The choice \code{4.685061} yields an estimator with asymptotic
    efficiency of 95\% for normal errors.}
  \item{max.it}{integer specifying the maximum number of IRWLS iterations.}
  \item{groups}{(for the fast-S algorithm): Number of
    random subsets to use when the data set is large.}
  \item{n.group}{(for the fast-S algorithm): Size of each of the
    \code{groups} above.  Note that this must be at least \eqn{p}.}
  \item{best.r.s}{(for the fast-S algorithm): Number of
    of best candidates to be iterated further (i.e.,
    \dQuote{\emph{\bold{r}efined}}); is denoted \eqn{t} in
    Salibian-Barrera \& Yohai(2006).}
  \item{k.fast.s}{(for the fast-S algorithm): Number of
    local improvement steps (\dQuote{\emph{I-steps}}) for each
    re-sampling candidate.}
  \item{compute.rd}{logical indicating if robust distances (based on
    the MCD robust covariance estimator \code{\link{covMcd}}) are to be
    computed for the robust diagnostic plots.  This may take some
    time to finish, particularly for large data sets, and can lead to
    singularity problems when there are \code{\link{factor}} explanatory
    variables (with many levels, or levels with \dQuote{few}
    observations). Hence, is \code{FALSE} by default.}
}
\author{ Matias Salibian-Barrera }
\seealso{ \code{\link{lmrob}}, also for references and examples.
}
\keyword{robust}
\keyword{regression}

