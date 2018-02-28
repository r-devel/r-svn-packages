\name{Saccharin}
\alias{Saccharin}
\non_function{}
\title{Elimination of saccharin}
\description{
The \code{Saccharin} data frame has 9 rows and 3 columns from an
experiment on the elimination of saccharin.
}
\format{
This data frame contains the following columns:
  \describe{
    \item{start}{
      start of collection period, measured from the start of the
      experiment (hr).
    }
    \item{length}{
      duration of the collection period (hr).
    }
    \item{sacch}{
      exhaled radioactively tagged CO2 (g).
    }
}
}
\source{
  Bates and Watts (1998), \emph{Nonlinear Regression Analysis and Its
    Applications}, Wiley (Appendix A1.5).
}
\examples{
xyplot(cumsum(sacch) ~ I(start + length), Saccharin,
     type = c("g", "S"), scales = list(x = list(log = 2)),
     xlab = "Time (hr)", ylab = "Cumulative saccharin")
}
\keyword{datasets}
