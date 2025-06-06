\name{Gasoline}
\alias{Gasoline}
\title{Refinery yield of gasoline}
\description{
  The \code{Gasoline} data frame has 32 rows and 6 columns.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{yield}{
      a numeric vector giving the percentage of crude oil converted to
      gasoline after distillation and fractionation
    }
    \item{endpoint}{
        a numeric vector giving the temperature (degrees F) at which all
	the gasoline is vaporized
    }
    \item{Sample}{
      the inferred crude oil sample number - a factor with levels
      \code{A} to \code{J}
    }
    \item{API}{
      a numeric vector giving the crude oil gravity (degrees API)
    }
    \item{vapor}{
	a numeric vector giving the vapor pressure of the crude oil
	\eqn{(\mathrm{lbf}/\mathrm{in}^2)}{(lbf/in^2)}
    }
    \item{ASTM}{
      a numeric vector giving the crude oil 10\% point ASTM---the
      temperature at which 10\% of the crude oil has become vapor. 
    }
  }
}
\details{
  Prater (1955) provides data on crude oil properties and
  gasoline yields.  Atkinson (1985)
  uses these data to illustrate the use of diagnostics in multiple
  regression analysis.  Three of the covariates---\textsf{API},
  \textsf{vapor}, and \textsf{ASTM}---measure characteristics of the
  crude oil used to produce the gasoline.  The other covariate ---
  \textsf{endpoint}---is a characteristic of the refining process.
  Daniel and Wood (1980) notice that the covariates characterizing
  the crude oil occur in only ten distinct groups and conclude that the
  data represent responses measured on ten different crude oil samples.
}
\source{
  Prater, N. H. (1955), Estimate gasoline yields from crudes,
  \emph{Petroleum Refiner}, \bold{35} (5).

  Atkinson, A. C. (1985), \emph{Plots, Transformations, and
    Regression}, Oxford Press, New York.

  Daniel, C. and Wood, F. S. (1980), \emph{Fitting Equations to Data},
  Wiley, New York

  Venables, W. N. and Ripley, B. D. (1999) \emph{Modern Applied
    Statistics with S-PLUS (3rd ed)}, Springer, New York.
}
\examples{
str(Gasoline)
xyplot(yield ~ endpoint | Sample, Gasoline, aspect = 'xy',
       main = "Gasoline data", xlab = "Endpoint (degrees F)",
       ylab = "Percentage yield",
       type = c("g", "p", "r"),
       index.cond = function(x,y) coef(lm(y~x))[2],
       layout = c(5,2))
print(m1 <- lmer(yield ~ endpoint + (1|Sample), Gasoline), corr = FALSE)
m2 <- lmer(yield ~ endpoint + (endpoint|Sample), Gasoline, verbose = 1)
print(m2)
Gasoline$endptC <- with(Gasoline, endpoint - mean(endpoint))
m3 <- lmer(yield ~ endpoint + (endptC|Sample), Gasoline, verbose = 1)
print(m3)
xyplot(endptC ~ `(Intercept)`, ranef(m3)[[1]], type = c("g", "p", "r"),
       aspect = 1)
}
\keyword{datasets}
