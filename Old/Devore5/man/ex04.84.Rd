\name{ex04.84}
\alias{ex04.84}
\title{data from exercise 4.84}
\description{
    The \code{ex04.84} data frame has 20 rows and 1 columns of bearing
    load lifetimes.
}
\format{
    This data frame contains the following columns:
    \describe{
	\item{loadlife}{
	    a numeric vector of bearing load life (million revs) for
	    bearings tested at a 6.45 kN load.
	}
    }
}
\source{
    Devore, J. L. (2000) \emph{Probability and Statistics for
	Engineering and the Sciences (5th ed)}, Duxbury

    (1984) "The load-life relationship for M50 bearings with silicon
    nitrate ceramic balls", \emph{Lubrication Engineering}, 153--159.
}
\examples{
data(ex04.84)
attach(ex04.84)
boxplot(loadlife, ylab = "Load-Life (million revs)",
        main = "Bearing load-lifes from exercise 4.84",
        col = "lightgray")
## Normal probability plot
qqnorm(loadlife, ylab = "Load-life (million revs)")
qqline(loadlife)
## Weibull probability plot
plot(
  log(-log(1 - (seq(along = loadlife) - 0.5)/length(loadlife))),
  log(sort(loadlife)), xlab = "Theoretical Quantiles",
  ylab = "log(load-life) (log(million revs))",
  main = "Weibull Q-Q Plot", las = 1)
detach()
}
\keyword{datasets}
