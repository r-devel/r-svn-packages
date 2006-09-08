\name{xmp13.06}
\alias{xmp13.06}
\non_function{}
\title{data from Example 13.6}
\description{
  The \code{xmp13.06} data frame has 16 rows and 2 columns on the yield
  of paddy (a grain farmed in India) versus the time of harvest.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{days}{
      date of harvesting (number of days after flowering.
    }
    \item{yield}{
      yield (kg/ha) of paddy
    }
  }
}
\source{
  Devore, J. L. (2000) \emph{Probability and Statistics for Engineering
    and the Sciences (5th ed)}, Duxbury

  (1975), ``Determination of biological maturity and effect of harvesting
  and drying conditions on milling quality of paddy'',
  \emph{J. Agricultural Eng. Research}, 353-361
}
\examples{
data(xmp13.06)
plot(yield ~ days, data = xmp13.06,
   main = "Compare to Figure 13.10, page 564")
fm1 <- lm(yield ~ days + I(days^2), data = xmp13.06)
summary(fm1)
anova(fm1)
predict(fm1, list(days = 25), interval = "conf")
predict(fm1, list(days = 25), interval = "pred")
}
\keyword{datasets}
