\name{Lipo}
\alias{Lipo}
\non_function{}
\title{Lipoprotein concentrations}
\description{
  The \code{Lipo} data frame has 12 rows and 2 columns of lipoprotein
  concentrations over time.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{time}{
      a numeric vector giving the time of the concentration measurement (hr)
    }
    \item{conc}{
      a numeric vector of concentrations.
    }
  }
}
\source{
  Bates and Watts (1998), \emph{Nonlinear Regression Analysis and Its
    Applications}, Wiley (Appendix A4.1).
}
\examples{
plot(conc ~ time, data = Lipo,
     xlab = "Time since drug administration (hr)",
     ylab = "Lipoprotein concentration",
     main = "Lipo data and fitted curve", las = 1)
fm1 <- nls(conc ~ SSbiexp(time, A1, lrc1, A2, lrc2), data = Lipo)
summary(fm1)
usr <- par("usr")
xx <- seq(usr[1], usr[2], len = 51)
lines(xx, predict(fm1, list(time = xx)))
title(sub = deparse(fm1$call$formula))
}
\keyword{datasets}
