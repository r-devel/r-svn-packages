\name{xmp07.15}
\alias{xmp07.15}
\title{data from Example 7.15}
\description{
  The \code{xmp07.15} data frame has 17 rows and 1 columns.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{voltage}{
      breakdown voltage of electrically stressed circuits
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
data(xmp07.15)
boxplot(xmp07.15, main = "Breakdown voltage")
with(xmp07.15, qqnorm(voltage, main = "Breakdown voltage"))
with(xmp07.15, qqline(voltage))
attach(xmp07.15)
var(voltage) * (length(voltage) - 1)/
  qchisq(c(0.975, 0.025), df = length(voltage) - 1)
detach()
}
\keyword{datasets}
