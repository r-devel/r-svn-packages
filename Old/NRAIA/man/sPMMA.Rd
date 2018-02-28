\name{sPMMA}
\alias{sPMMA}
\non_function{}
\title{syndiotactic poly-methyl-methacrylate}
\description{
The \code{sPMMA} data frame has 23 rows and 3 columns from an experiment
on the dielectric response of a sample of syndiotactic
poly-methyl-methacrylate (sPMMA).
}
\format{This data frame contains the following columns:
  \describe{
    \item{freq}{
      a numeric vector giving the frequency of the driving signal (Hz).
    }
    \item{real}{
      a numeric vector giving the real component of the dielectric response.
    }
    \item{imag}{
      a numeric vector giving the imaginary component of the dielectric response.
    }
  }
}
\source{
  Bates and Watts (1998), \emph{Nonlinear Regression Analysis and Its
    Applications}, Wiley (Appendix A1.5).
}
\examples{
str(sPMMA)
xyplot(imag ~ real, data = sPMMA, xlab = "Real component",
     ylab = "Imaginary component", type = c("p", "g"))
xyplot(real + imag ~ freq, data = sPMMA, xlab = "Frequency (Hz)",
     ylab = "Component of dielectric response",
     main = "sPMMA data", type = c("g", "p"), aspect = 'xy',
     auto.key = list(space = "right"))
}
\keyword{datasets}
