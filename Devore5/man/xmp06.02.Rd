\name{xmp06.02}
\alias{xmp06.02}
\non_function{}
\title{data from Example 6.2}
\description{
  The \code{xmp06.02} data frame has 20 rows and 1 columns.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{Voltage}{
      the dielectric breakdown voltage for pieces of expoxy resin
    }
  }
}
\details{
  This is the same data set as \code{xmp04.29}.
}
\source{
  Devore, J. L. (2000) \emph{Probability and Statistics for Engineering
    and the Sciences (5th ed)}, Duxbury
}
\examples{
data(xmp06.02)
summary(xmp06.02)                # gives mean, median, etc.
attach(xmp06.02)
mean(range(Voltage))             # average of the extremes
mean(Voltage, trim = 0.1)        # trimmed mean
}
\keyword{datasets}
