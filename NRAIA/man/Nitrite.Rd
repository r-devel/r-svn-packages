\name{Nitrite}
\alias{Nitrite}
\non_function{}
\title{Nitrite utilization by bean leaves}
\description{
  The \code{Nitrite} data frame has 48 rows and 3 columns giving nitrite
  utilization by bean leaves under different light conditions on two
  different days.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{light}{
      a numeric vector of light intensities.
    }
    \item{utilization}{
      a numeric vector of nitrite utilizations.
    }
    \item{day}{
      a factor with levels
      \code{1} 
      \code{2} 
    }
  }
}
\source{
  Bates and Watts (1998), \emph{Nonlinear Regression Analysis and Its
    Applications}, Wiley (Appendix A1.10).
}
\examples{
str(Nitrite)
xyplot(utilization ~ light, Nitrite, groups = day,
     xlab = "Light intensity", type = c("a", "g", "p"),
     ylab = "Nitrite utilization",
     auto.key = list(space = "top", columns = 2))
}
\keyword{datasets}
