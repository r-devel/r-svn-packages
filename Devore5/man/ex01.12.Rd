\name{ex01.12}
\alias{ex01.12}
\non_function{}
\title{data from exercise 1.12}
\description{
    The \code{ex01.12} data frame has 129 rows and 1 columns.
}
\format{
    This data frame contains the following columns:
    \describe{
	\item{rate}{
	    a numeric vector of shower-flow rates (L/min)
	}
    }
}
\details{
    The shower-flow rates for a sample of houses in Perth, Australia.
}
\source{
    Devore, J. L. (2000) \emph{Probability and Statistics for Engineering
	and the Sciences (5th ed)}, Duxbury 
    
    James, I. R. and Knuiman, M. W. (1987) An application of Bayes
    methodology to the analysis of diary records in a water use
    study. \emph{Journal of the American  Statistical Association},
    \bold{82}, 705--711. 
}
\examples{
data(ex01.12)
attach(ex01.12)
boxplot(rate, main = "Boxplot of rate", ylab = "rate")
summary(rate)
stem(rate)   # stem-and-leaf plot
hist(rate, xlab = "Shower-flow rate (L/min)")
hist(rate, prob = TRUE,  # histogram on probability scale
     xlab = "Shower-flow rate (L/min)")
lines(density(rate, bw = 1), col = "blue")
rug(rate)    # add the observations as a rug
detach()
}
\keyword{datasets}
