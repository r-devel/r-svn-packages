\name{ex01.11}
\alias{ex01.11}
\non_function{}
\title{data from exercise 1.11}
\description{
    The \code{ex01.11} data frame has 79 rows and 1 columns.
}
\format{
    This data frame contains the following columns:
    \describe{
	\item{octane}{
	    a numeric vector
	}
    }
}
\details{
    Motor octane rating for various gasoline blends.
}
\source{
    Devore, J. L. (2000) \emph{Probability and Statistics for
	Engineering and the Sciences (5th ed)}, Duxbury

    Snee, R. D. (1977) Validation of regression models: methods and
    examples, \emph{Technometrics}, 415--428.
}
\examples{
data(ex01.11)
attach(ex01.11)
stem(octane)              # compact
stem(octane, scale = 2)   # expanded
summary(octane)
hist(octane)              # standard histogram
hist(octane, prob = TRUE)
lines(density(octane), col = "blue")
rug(octane)
detach()
}
\keyword{datasets}
