\name{ex04.80}
\alias{ex04.80}
\non_function{}
\title{data from exercise 4.80}
\description{
    The \code{ex04.80} data frame has 10 rows and 1 columns of bearing
    lifetimes.
}
\format{
    This data frame contains the following columns:
    \describe{
	\item{lifetime}{
	    a numeric vector
	}
    }
}
\source{
    Devore, J. L. (2000) \emph{Probability and Statistics for
	Engineering and the Sciences (5th ed)}, Duxbury

    (1985), "Modified Moment Estimation for the three-parameter
    Log-normal distribution", \emph{Journal of Quality Technology},
    92--99.
}
\examples{
data(ex04.80)
attach(ex04.80)
boxplot(lifetime, ylab = "Lifetime (hr)",
        main = "Bearing lifetimes from exercise 4.80",
        col = "lightgray")
## Normal probability plot on the original time scale
qqnorm(lifetime, ylab = "Lifetime (hr)", las = 1)
qqline(lifetime)
## Try normal probability plot of the log(lifetime)
qqnorm(log(lifetime), ylab = "log(lifetime) (log(hr))",
       las = 1)
qqline(log(lifetime))
detach()
}
\keyword{datasets}
