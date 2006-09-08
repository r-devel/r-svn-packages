\name{xmp01.15}
\alias{xmp01.15}
\non_function{}
\title{data from Example 1.15}
\description{
  The \code{xmp01.15} data frame has 20 rows and 1 columns.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{lifetime}{
      a numeric vector of lifetimes (hr) of a certain type of
      incandescent light bulb.
    }
  }
}
\source{
  Devore, J. L. (2000) \emph{Probability and Statistics for Engineering
    and the Sciences (5th ed)}, Duxbury
}
\examples{
data(xmp01.15)
attach(xmp01.15)
summary(lifetime)              # produces mean, median, etc.
mean(lifetime, trim = 0.1)     # 10\% trimmed mean
mean(lifetime, trim = 0.2)     # 20\% trimmed mean
dotchart(lifetime)
hist(lifetime)                 # display a histogram
rug(lifetime)                  # add the data
abline(v = median(lifetime), col = 2, lty = 2)
abline(v = mean(lifetime), col = 3, lty = 2)
abline(v = mean(lifetime, trim = 0.1), col = 4, lty = 2)
abline(v = mean(lifetime, trim = 0.2), col = 5, lty = 2)
legend(600, 6,
  c("median", "mean", "10\% trimmed mean", "20\% trimmed mean"),
  col = 2:5, lty = 2)
}
\keyword{datasets}
