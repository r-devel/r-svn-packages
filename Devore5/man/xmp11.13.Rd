\name{xmp11.13}
\alias{xmp11.13}
\non_function{}
\title{data from Example 11.13}
\description{
  The \code{xmp11.13} data frame has 16 rows and 4 columns.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{Strength}{
      a numeric vector
    }
    \item{age}{
      a numeric vector
    }
    \item{tempture}{
      a numeric vector
    }
    \item{soil}{
      a numeric vector
    }
  }
}
\details{
}
\source{
  Devore, J. L. (2000) \emph{Probability and Statistics for Engineering
    and the Sciences (5th ed)}, Duxbury
}
\examples{
data(xmp11.13)
fm1 <- lm(Strength ~ age * tempture * soil, xmp11.13)
anova(fm1)   # compare with Table 11.12, page 471
}
\keyword{datasets}
