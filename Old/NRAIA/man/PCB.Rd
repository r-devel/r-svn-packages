\name{PCB}
\alias{PCB}
\non_function{}
\title{PCB concentrations in Lake Cayuga fish}
\description{
  The \code{PCB} data frame has 28 rows and 2 columns
  giving the concentration of polychlorinated biphenyls (PCB's) in fish
  caught in Lake Cayuga in northern New York state.  The fish had been
  tagged as yearlings so their ages could be accurately determined.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{age}{
      a numeric vector giving the age of the fish in years.
    }
    \item{conc}{
      a numeric vector giving the concentration of PCB's in the fish.
    }
  }
}
\source{
  Bates and Watts (1998), \emph{Nonlinear Regression Analysis and Its
    Applications}, Wiley (Appendix A1.1).
}
\examples{
str(PCB)
# compare to Figure 1.1 (p. xxx)
xyplot(conc ~ age, PCB, aspect = 'xy',
     type = c("g", "p"), xlab = "Age of fish (yr)",
     ylab = "PCB concentration")
# compare to Figure 1.2 (p. xxx)
xyplot(conc ~ age, PCB, scales = list(y = list(log = 2)),
     aspect = 'xy', type = c("g", "p", "smooth"),
     xlab = "Age of fish (yr)",
     ylab = "PCB concentration")
# linear model in cube root of age
summary(fm1 <- lm(log(conc) ~ I(age^(1/3)), data = PCB))
xyplot(log(conc) ~ I(age^(1/3)), data = PCB, aspect = 'xy',
     xlab = "Cube root of age (yr)", type = c("g", "p", "r"),
     ylab = "log(PCB concentration)",
     main = "Transformed PCB data and fitted line",
     sub = deparse(fm1$call$formula))
# diagnostic plots
opar <- par(mfrow = c(2, 2))
plot(fm1, which = 1:4, las = 1)
par(opar)
}
\keyword{datasets}
