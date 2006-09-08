\name{xmp09.08}
\alias{xmp09.08}
\non_function{}
\title{data from Example 9.8}
\description{
  The \code{xmp09.08} data frame has 6 rows and 2 columns.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{bottom}{
      a numeric vector
    }
    \item{surface}{
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
data(xmp09.08)
boxplot(xmp09.08, main = "Boxplot of data from Example 9.8")
attach(xmp09.08)
boxplot(bottom-surface, main = "Boxplot of differences from Example 9.8")
t.test(bottom, surface, alt = "greater", paired = TRUE)
detach()
}
\keyword{datasets}
