\name{summary.glmrob}
\alias{summary.glmrob}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarizing Robust Fits of Discrete Generalized Linear Models}
\description{
  \code{summary} method for class \code{"glmrob"}.
}
\usage{
\method{summary}{glmrob}(object, correlation = FALSE, symbolic.cor = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{"glmrob"}, usually, a result of
    a call to \code{\link{glmrob}}.}
  \item{correlation}{logical; if \code{TRUE}, the correlation matrix of
    the estimated parameters is returned and printed.}
  \item{symbolic.cor}{logical. If \code{TRUE}, print the correlations in
    a symbolic form (see \code{\link{symnum}}) rather than as numbers.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\details{
  \code{\link{summary.glmrob}} tries to be smart about formatting the
  coefficients, standard errors, etc. and additionally gives
  "significance stars" if \code{signif.stars} is \code{TRUE}. ???
}

\value{
  The function \code{\link{summary.glmrob}} computes and returns a list
  of summary statistics of the robustly fitted linear model given in
  \code{object}. The following elements are in the list:
  \item{FixIt }{???}
}

%%\references{ ~put references to the literature/web site here ~ }
\author{Andreas Ruckstuhl}

%%\note{ ~~further notes~~ }
%%
%% ~Make other sections like Warning with \section{Warning }{....} ~
%%
%%\seealso{ ~~objects to See Also as \code{\link{~~fun~~}}, ~~~ }

\examples{
## To do
}
\keyword{models}
\keyword{regression}
\keyword{nonlinear}
\keyword{robust}
