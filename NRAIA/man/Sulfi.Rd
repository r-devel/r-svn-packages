\name{Sulfi}
\alias{Sulfi}
\non_function{}
\title{Pharmacokinetics of sulfisoxazole}
\description{
  The \code{Sulfi} data frame has 12 rows and 2 columns from an experiment
  on the pharmacokinetics of sulfisoxazole.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{time}{
      a numeric vector giving the time since drug administration (min).
    }
    \item{conc}{
      a numeric vector giving the observed concentration of
      sulfisoxazole (\eqn{\mu}{u}g/ml).
    }
  }
}
\source{
  Bates and Watts (1998), \emph{Nonlinear Regression Analysis and Its
    Applications}, Wiley (Appendix A1.5).
}
\examples{
plot(conc ~ time, data = Sulfi, las = 1,
     xlab = "Time since drug administration (min)",
     ylab = expression(plain("Sulfisoxazole concentration (")
               *mu*plain("g/ml)")),
     main = "Sulfisoxazole data and fitted curve")
fm1 <- nls(conc ~ SSbiexp(time, A1, lrc1, A2, lrc2),
    data = Sulfi)
summary(fm1)
usr <- par("usr")
xx <- seq(usr[1], usr[2], len = 50)
lines(xx, predict(fm1, list(time = xx)))
plot(conc ~ time, data = Sulfi, las = 1, log = "y",
     xlab = "Time since drug administration (min)",
     ylab = expression(plain("Sulfisoxazole concentration (")
               *mu*plain("g/ml)")),
     main = "Sulfisoxazole data (log scale)")
lines(xx, predict(fm1, list(time = xx)))
}
\keyword{datasets}
