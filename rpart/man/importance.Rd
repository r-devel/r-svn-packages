\name{importance}
\alias{importance}
\title{Variable Importance for an rpart Fit}
\description{
Variable importance for an \code{rpart} fit.
}
\usage{importance(fit, rescale = TRUE)}
\arguments{
  \item{fit}{the result of an \code{rpart} call}
  \item{rescale}{if TRUE rescale the result to add to 100}
}
\details{
  A variable may appear many times in an rpart tree.  This function sums
  up the \dQuote{impact} of the variable over all splits for which it is the
  determining variable.

  To account for redundant variables the sum also includes surrogate
  splits.
  Say for instance that weight was the most important predictor, but the
  data set had two weight variables kg_monday and kg_friday, which are
  almost but not exactly identical.
  The fit might include one these for half the weight
  splits and the other for the rest, diluting the importance values for
  both.  If a primary split has an importance of z, then the primary split
  will be credited with z and the surrogates with z time the adjusted
  agreement for that surrogate.

  Agreement is the number of observations sent the same direction, left or
  right, as the primary split.
  Adjusted agreement is the impromvement on the naive \dQuote{go with the
  majority} surrogate.  So for instance if the primary split sent 80 left
  and 15 right and the surrogate agreed with the primary on 90 of the 95,
  then the agreement would be 90/95 and adjusted agreement would be
  (90-80)/(95-80).
}
\value{
  A vector of variable importance values.
}

\author{Terry Therneau}
\seealso{\code{\link{rpart}}}
\examples{
fit <- rpart(HP ~ ., data = car.test.frame)
importance(fit)
}
\keyword{tree}

