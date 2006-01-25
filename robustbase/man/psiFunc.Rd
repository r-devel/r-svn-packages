\name{psiFunc}
\alias{psiFunc}
\title{Constructor for psi Function Objects}
\description{
  Construct a psi function object, i.e., an object of class
  \code{psi_func}.  \code{psiFunc()} is just a slightly more convenient
  interface to the (more lower level) \code{new("psi_func", .....)}.
}
\usage{
psiFunc(rho, psi, wgt, Dpsi, Erho = NULL, Epsi2 = NULL, EDpsi = NULL, ...)
}
\arguments{
  \item{rho}{a \code{\link{function}} of \code{x} (say) and (optionally,
    but typically) of tuning parameters.}
  \item{psi}{the derivative of \code{rho(x)}.}
  \item{wgt}{the \code{\link{function}} computing \code{psi(x) / x}.}
  \item{Dpsi}{the derivative of \code{psi(x)} (i.e., 2nd derivative of
    \code{rho(x)}).}
%% FIXME: These should be part of
  \item{Erho, Epsi2, EDpsi}{Functionals, i.e. functions of only the
    tuning parameters, computing expected values (under \eqn{X \sim
      \mathcal{N}(0,1)}{X \~ N(0,1)}) of \eqn{\rho(X)}, \eqn{\psi^2(X)}
    and \eqn{\psi'(X)}.}
  \item{\dots}{specifications of tuning parameters and their defaults.}
}
\value{
  an object of class \code{\link{psi_func-class}}.
}
\author{Martin Maechler}
\seealso{\code{\link{psi_func-class}} for the class description.}
\examples{
## Trivial psi Function: Classical (Gaussian / Least Squares):
F1 <- function(x) rep.int(1, length(x))
cPsi <- psiFunc(rho = function(x) x^2 / 2, psi = function(x) x,
                wgt = F1, Dpsi = F1,
                Erho = function(x) rep.int(1/2, length(x)),
                Epsi2 = F1, EDpsi = F1)
                ##Epsi2 = function(...)1, EDpsi=function(...)1)
}
\keyword{classes}
\keyword{robust}
