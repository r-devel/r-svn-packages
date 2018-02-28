\name{Oilshale}
\alias{Oilshale}
\non_function{}
\title{Pyrolysis of Oil Shale}
\description{
    The \code{Oilshale} data frame has 64 rows and 4 columns from an
    experiment on the pyrolysis of oil shale.
}
\format{
    This data frame contains the following columns:
    \describe{
	\item{time}{
	    a numeric vector giving the time since the beginning of the
	    experiment.
	}
	\item{bitumen}{
	    a numeric vector giving the proportion of bitumen.
	}
	\item{oil}{
	    a numeric vector giving the proportion of oil.
	}
	\item{temp}{
	    a numeric vector giving the temperature of the run.
	}
    }
}
\source{
    Bates and Watts (1998), \emph{Nonlinear Regression Analysis and Its
	Applications}, Wiley (Appendix A1.5).
}
\examples{
str(Oilshale)
}
\keyword{datasets}
