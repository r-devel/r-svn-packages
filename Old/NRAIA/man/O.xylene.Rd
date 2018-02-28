\name{O.xylene}
\alias{O.xylene}
\non_function{}
\title{O-xylene reaction data}
\description{
  The \code{O.xylene} data frame has 57 rows and 4 columns of data on a
  chemical engineering experiment.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{oxygen}{
      a numeric vector giving the oxygen partial pressure.
    }
    \item{o.xyl}{
      a numeric vector giving the O.xylene partial pressure.
    }
    \item{temp}{
      a numeric vector giving the temperature (Kelvins).
    }
    \item{rate}{
      a numeric vector giving the rate of reaction.
    }
  }
}
\source{
  Bates and Watts (1998), \emph{Nonlinear Regression Analysis and Its
    Applications}, Wiley (Appendix A1.6).
}
\examples{
str(O.xylene)
tempf <- with(O.xylene, ordered(temp))
levels(tempf) <- paste(levels(tempf), "K")
Oxygen <- with(O.xylene, equal.count(oxygen))
O.xyl <- with(O.xylene, equal.count(o.xyl))
xyplot(rate ~ o.xyl | Oxygen * tempf, O.xylene,
       type = c("g", "p"), aspect = 'xy',
       xlab = "O.xylene partial pressure",
       ylab = "Reaction rate")
xyplot(rate ~ o.xyl | Oxygen, O.xylene, groups = tempf,
       type = c("g", "p"), aspect = 2.5,
       xlab = "O.xylene partial pressure",
       ylab = "Reaction rate",
       auto.key = list(space = "right"))
xyplot(rate ~ oxygen | O.xyl, O.xylene, groups = tempf,
       type = c("g", "p"), aspect = 2.5,
       xlab = "Oxygen partial pressure",
       ylab = "Reaction rate",
       auto.key = list(space = "right"))
}
\keyword{datasets}
