\name{Sacch2}
\alias{Sacch2}
\non_function{}
\title{Pharmacokinetics of saccharin}
\description{
  The \code{Sacch2} data frame has 10 rows and 2 columns from an
  experiment on the pharmacokinetics of saccharin.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{time}{
      a numeric vector giving the time since drug administration (min).
    }
    \item{conc}{
      a numeric vector giving the observed concentration of saccharin.
    }
  }
}
\source{
  Bates and Watts (1998), \emph{Nonlinear Regression Analysis and Its
    Applications}, Wiley (Appendix A1.5).
}
\examples{
str(Sacch2)
xyplot(conc ~ time, Sacch2, type = c("g", "b"),
     xlab = "Time since drug administration (min)",
     ylab = "Saccharin concentration", aspect = "xy")
xyplot(conc ~ time, data = Sacch2, type = c("g", "b"),
     scales = list(y = list(log = 2)), aspect = 'xy',
     xlab = "Time since drug administration (min)",
     ylab = "Saccharin concentration")
\dontrun{
fm1 <- nls(conc ~ SSfol(1.0, time, lKe, lKa, lCl), data = Sacch2)
summary(fm1)
xpred <- seq(0, 140, len = 51)
ypred <- predict(fm1, list(time = xpred, Dose = rep(1.0, length(xpred))))
lines(xpred, ypred)
}
}
\keyword{datasets}
