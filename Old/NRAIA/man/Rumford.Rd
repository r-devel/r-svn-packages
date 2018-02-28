\name{Rumford}
\alias{Rumford}
\non_function{}
\title{Count Rumford's cooling data}
\description{
The \code{Rumford} data frame has 13 rows and 2 columns from an
experiment by Count Rumford on the rate of cooling.
}
\format{
This data frame contains the following columns:
\describe{
\item{time}{
a numeric vector giving the time since the beginning of the experiment (hr).
}
\item{temp}{
a numeric vector giving the temperature (degrees Fahrenheit) of the cannon.
}
}
}
\source{
  Bates and Watts (1998), \emph{Nonlinear Regression Analysis and Its
    Applications}, Wiley (Appendix A1.5).
}
\examples{
plot(temp ~ time, Rumford, las = 1,
     xlab = "Time since beginning of experiment (hr)",
     ylab = "Temperature of cannon (degrees F)",
     main = "Rumford data and fitted curve")
fm1 <- nls(temp ~ 60 + 70 * exp(tc * time), data = Rumford,
           start = c(tc = -0.01))
summary(fm1)
usr <- par("usr")
xx <- seq(usr[1], usr[2], len = 50)
lines(xx, predict(fm1, list(time = xx)))
}
\keyword{datasets}
