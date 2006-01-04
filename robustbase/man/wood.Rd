\name{wood}
\alias{wood}
% \alias{wood.x}
% \alias{wood.y}
\docType{data}
\title{Modified Data on Wood Specific Gravity}
\description{
  The original data are from Draper and Smith (1966) and
  were used to determine the influence of anatomical factors on wood
  specific gravity, with five explanatory variables and an
  intercept. These data were contaminated by replacing a few
  observations with outliers.
}
\usage{data(wood)}
\format{
  A data frame with 20 observations on the following 6 variables.
  \describe{
    \item{x1, x2, x3, x4, x5}{explanatory \dQuote{anatomical} wood
      variables.}
    \item{y}{wood specific gravity, the target variable.}
  }
%  For convenience, the data sets \code{wood.x}, a matrix with the first 5 (independent)
%  variables of the data frame, and \code{wood.y}, the numeric vector giving the sixth
%  (dependent) variable, are provided as well.
}
\source{
  Draper and Smith (1966, p.227)

  Peter J. Rousseeuw and Annick M. Leroy (1987)
  \emph{Robust Regression and Outlier Detection} Wiley, p.243, table 8.
}
\examples{
data(wood)
plot(wood)
%%plot(wood.x)
%%covMcd(wood.x)
summary( lm.wood <-        lm(y ~ ., data = wood))
summary(rlm.wood <- MASS::rlm(y ~ ., data = wood))
}
\keyword{datasets}
