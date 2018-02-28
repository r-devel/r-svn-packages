\name{Pinene2}
\alias{Pinene2}
\non_function{}
\title{Alpha-pinene and by-products}
\description{
  The \code{Pinene} data frame has 8 rows and 6 columns giving the
  proportion of alpha-pinene and four of its by-products over time from
  a second experiment.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{time}{
      a numeric vector giving the time at which the observation is made. (hr)
    }
    \item{a.pin}{
      a numeric vector giving the alpha-pinene proportion. (\%)
    }
    \item{dipen}{
      a numeric vector giving the dipene proportion. (\%)
    }
    \item{alloo}{
      a numeric vector giving the allo-ocimene proportion. (\%)
    }
    \item{pyron}{
      a numeric vector giving the pyronene proportion. (\%)
    }
    \item{dimer}{
      a numeric vector giving the dimer proportion. (\%)
    }
  }
}
\source{
  Bates and Watts (1998), \emph{Nonlinear Regression Analysis and Its
    Applications}, Wiley (Appendix A1.6).
}
\seealso{\code{\link{Pinene}}}
\examples{
str(Pinene2)
xyplot(a.pin + dipen + alloo + pyron + dimer ~ time, Pinene2,
  type = c("b", "g"), aspect = 'xy',
  xlab = "Time (hr)", ylab = "Proportion (\%)",
  auto.key = list(space = "right", lines = TRUE))
}
\keyword{datasets}
