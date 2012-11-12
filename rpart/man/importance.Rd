\name{importance}
\alias{importance}
\title{Variable Importance for an rpart Fit}
\description{
Variable importance for an \code{rpart} fit.
}
\usage{
importance(fit, rescale = TRUE)
}
\arguments{
  \item{fit}{the result of an \code{rpart} call.}
  \item{rescale}{if \code{TRUE} (the default) rescale the results
    to sum to 100.}
}
\details{
  A variable may appear many times in an \code{rpart} tree.  This
  function sums up the \dQuote{impact} of the variable over all splits
  for which it is the determining variable.

  To account for redundant variables the sum also includes surrogate
  splits.  Say for instance that weight was the most important
  predictor, but the data set had two weight variables \code{kg_monday}
  and \code{kg_friday}, which are almost but not exactly identical.  The
  fit might include one of these for half the weight splits and the
  other for the rest, diluting the importance values for both.  If a
  primary split has an importance of \eqn{z}, then the primary split
  will be credited with \eqn{z} and the surrogates with \eqn{z} times
  the adjusted agreement for that surrogate.

  Agreement is the number of case sent in the same direction (left or
  right) as the primary split.  Adjusted agreement is the impromvement
  on the naive \dQuote{go with the majority} surrogate.  So for instance
  if the primary split sent 80 cases left and 15 cases right and the
  surrogate agreed with the primary on 90 of the 95, then the agreement
  would be 90/95 and adjusted agreement would be (90-80)/(95-80).
}
\value{
  A vector of variable importance values in decreasing order, labelled
  by the names of the variables which are used in \code{fit}.
}

\note{
  There are functions called \code{importance} in several packages.  Use
  namespace imports to make sure you get the one you intend.
}

\author{Terry Therneau}
\seealso{\code{\link{rpart}}}
\examples{
fit <- rpart(HP ~ ., data = car90)
importance(fit)
}
\keyword{tree}

