\name{Isom}
\alias{Isom}
\non_function{}
\title{Isomerization data}
\description{
  The \code{Isom} data frame has 24 rows and 4 columns from an
  isomerization experiment.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{hyd}{
      partial pressure of hydrogen (psia).
    }
    \item{n.pent}{
      partial pressure of n-pentane (psia).
    }
    \item{iso.pen}{
      partial pressure of isopentane (psia).
    }
    \item{rate}{
      reaction rate for isomerization of n-pentane to isopentane (1/hr).
    }
  }
}
\details{
  These are data on the reaction rate of catalytic isomerization of
  n-pentane to isopentane.
}
\source{
  Bates and Watts (1998), \emph{Nonlinear Regression Analysis and Its
    Applications}, Wiley (Appendix A1.5).

  Originally from Carr (1960), \emph{Industrial and Engineering
    Chemistry}, 52, pp. 391-396.
}
\examples{
str(Isom)
splom(Isom, main = "Isom data")
xyplot(rate ~ hyd, data = Isom, type = c("g", "p"),
     xlab = "Hydrogen partial pressure (psia)",
     ylab = expression(plain("Reaction rate (hr")^{-1}*plain(")")))
xyplot(rate ~ n.pent, data = Isom, type = c("g", "p"),
     xlab = "n-pentane partial pressure (psia)",
     ylab = expression(plain("Reaction rate (hr")^{-1}*plain(")")))
xyplot(rate ~ iso.pen, data = Isom, type = c("g", "p"),
     xlab = "iso-pentane partial pressure (psia)",
     ylab = expression(plain("Reaction rate (hr")^{-1}*plain(")")))
## Note - the model is mis-stated on page 272, x2 and x3 are reversed
fm1 <-
 nls(rate ~ b3*(n.pent - iso.pen/1.632)/(1+b2*hyd+b3*n.pent+b4*iso.pen),
         data = Isom, start = c(b2 = 0.1, b3 = 0.1, b4 = 0.1),
         algorithm = "plinear", trace = TRUE)
summary(fm1)        # compare to Table 2.2, page 56
}
\keyword{datasets}
