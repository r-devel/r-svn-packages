\name{Chloride}
\alias{Chloride}
\non_function{}
\title{Chloride ion concentrations}
\description{
  The \code{Chloride} data frame has 54 rows and 2 columns representing 
  measurements of the chloride ion concentration in blood cells suspended
  in a salt solution.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{conc}{
      A numeric vector giving the chloride ion concentration (\%).
    }
    \item{Time}{
      A numeric vector giving the time of the concentration
      measurement (min).
    }
  }
}
\details{
  There is a noticeable serial correlation in the residuals from a fit of
  a simple asymptotic regression model to these data, as described in
  section 3.8 of Bates and Watts (Wiley, 1988).  We can either use a
  modified model with \code{nls}, as described in the reference, or model
  it directly with the \code{gnls} function.
}
\source{
  Bates and Watts (1998), \emph{Nonlinear Regression Analysis and Its
    Applications}, Wiley (Appendix A4.1).

  Originally from Sredni (1970), Ph.D. Thesis, University of Wisconsin - Madison.
}
\examples{
xyplot(conc ~ Time, Chloride, aspect = 'xy',
     type = c("g", "p"), xlab = "Time (min)",
     ylab = "Chloride ion concentration (\%)")
## Fit nonlinear regression model
(fm1 <- nls(conc ~ Asym*(1 - prop*exp(-exp(lrc)*Time)), data = Chloride,
           start = c(Asym = 50, prop = 0.6, lrc = log(0.25)), trace = TRUE))
# usr <- par("usr")
# xx <- seq(usr[1], usr[2], len = 50)
# lines(xx, predict(fm1, list(Time = xx)))
plot(resid(fm1) ~ fitted(fm1),  # plot shows patterns in the residuals
     xlab = "Least squares fitted values",
     ylab = "Least squares residuals",
     main = "Residuals from least squares fit")
abline(h = 0, lty = 2, lwd = 0)
plot(resid(fm1)[-1], resid(fm1)[-length(resid(fm1))],  # lag plot
     xlab = "Least squares residual",
     ylab = "Lagged least squares residual",
     main = "Lag plot of least squares residuals")
abline(h = 0, lty = 2, lwd = 0)
abline(v = 0, lty = 2, lwd = 0)
## Use conditional linearity in another formulation of the model
(fm2 <- nls(conc ~ cbind(1, exp(-exp(lrc)*Time)), data = Chloride,
           start = c(lrc = log(0.25)), algorithm = "plinear", trace = TRUE))
## Using a self-Starting model in yet another formulation
(fm3 <- nls(conc ~ SSasympOff(Time, Asym, c0, lrc), data = Chloride))
}
\keyword{datasets}
