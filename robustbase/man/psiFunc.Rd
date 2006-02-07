\name{psiFunc}
\alias{psiFunc}
\alias{huberPsi}
\title{Constructor for Objects "Psi Function" Class}
%% Rd Problem:  'psi_func' cannot be used in title!
\description{
  \code{psiFunc(..)} is a convenience interface to
  \code{new("psi_func",..)}, i.e. for constructing objects of class
  \code{"psi_func"}.
}
\usage{
psiFunc(rho, psi, wgt, Dpsi, Erho = NULL, Epsi2 = NULL, EDpsi = NULL, ...)
}
\arguments{
  \item{rho, psi, wgt, Dpsi}{each a \code{\link{function}} of \code{x} and
    tuning parameters typically.}
  %% FIXME: see ./psi_func-class.Rd
  \item{Erho}{ ~~Describe \code{Erho} here~~ }
  \item{Epsi2}{ ~~Describe \code{Epsi2} here~~ }
  \item{EDpsi}{ ~~Describe \code{EDpsi} here~~ }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
}
% \details{
%   ~~ If necessary, more details than the description above ~~
% }
% \value{
%   ~Describe the value returned
%   If it is a LIST, use
%   \item{comp1 }{Description of 'comp1'}
%   \item{comp2 }{Description of 'comp2'}
%   ...
% }
%%\references{ ~put references to the literature/web site here ~ }
\author{Martin Maechler}
\seealso{\code{\link{psi_func-class}} for the class description.}
\examples{
## classical {trivial, not interesting}:
F1 <- function(x) rep.int(1, length(x))
cPsi <- psiFunc(rho = function(x) x^2 / 2, psi = function(x) x,
                wgt = F1, Dpsi = F1,
                Erho = function(x) rep.int(1/2, length(x)),
                Epsi2 = F1, EDpsi = F1)
}
\keyword{classes}
\keyword{robust}
