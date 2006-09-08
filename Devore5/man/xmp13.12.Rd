\name{xmp13.12}
\alias{xmp13.12}
\non_function{}
\title{data from Example 13.12}
\description{
  The \code{xmp13.12} data frame has 9 rows and 5 columns of data on
  characteristics of concrete.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{x1}{
      the \% limestone powder.
    }
    \item{x2}{
      the water-cement ratio.
    }
    \item{x1x2}{
      the interaction of limestone powder and water-cement ratio.
    }
    \item{strength}{
      the 28-day compressive strength (MPa).
    }
    \item{absorbability}{
      the absobability (\%).
    }
  }
}
\source{
  Devore, J. L. (2000) \emph{Probability and Statistics for Engineering
    and the Sciences (5th ed)}, Duxbury
}
\examples{
data(xmp13.12)
fm1 <- lm(strength ~ x1 * x2, data = xmp13.12)
summary(fm1)
}
\keyword{datasets}
