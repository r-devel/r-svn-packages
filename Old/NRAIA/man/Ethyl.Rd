\name{Ethyl}
\alias{Ethyl}
\non_function{}
\title{Ethyl Acrylate data}
\description{
  The \code{Ethyl} data frame has 12 rows and 3 columns of ethyl acrylate
  concentrations in exhalate of a rat.
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
    \item{CO2}{
      exhaled radioactively tagged CO\eqn{_2}{2} (g).
    }
  }
}
\source{
  Bates and Watts (1998), \emph{Nonlinear Regression Analysis and Its
    Applications}, Wiley (Appendix A1.10).

  Originally from Watts, deBethizy, and Stiratelli (1986), Technical
  Report, Rohm and Haas Co.
}
\examples{
str(Ethyl)
xyplot(cumsum(CO2) ~ I(start + length), # compare to Figure 3.10, page 97
     data = Ethyl, type = c("g","S"), scales=list(x=list(log =2)),
     xlab = "Time (hr)",
     ylab = expression(plain("Normalized cumulative CO")[2]*(g)))
}
\keyword{datasets}
