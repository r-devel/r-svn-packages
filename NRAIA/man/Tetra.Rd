\name{Tetra}
\alias{Tetra}
\non_function{}
\title{Tetracycline concentrations}
\description{
  The \code{Tetra} data frame has 9 rows and 2 columns from an experiment
  on the pharmacokinetics of tetracycline.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{time}{
      a numeric vector of time since drug administration (hr).
    }
    \item{conc}{
      a numeric vector of tetracycline concentrations.
    }
  }
}
\source{
  Bates and Watts (1998), \emph{Nonlinear Regression Analysis and Its
    Applications}, Wiley (Appendix A1.5).
}
\examples{
plot(conc ~ time, data = Tetra, las = 1,
     xlab = "Time since drug administration (hr)",
     ylab = "concentration",
     main = "Tetracycline data with fitted curve")
fm1 <- nls(conc ~ SSfol(1.0, time, lKe, lKa, lCl), Tetra)
summary(fm1)
usr <- par("usr")
xx <- seq(usr[1], usr[2], len = 101)
lines(xx, predict(fm1, list(time = xx)))
}
\keyword{datasets}
