\name{Leaves}
\alias{Leaves}
\non_function{}
\title{Growth of Leaves}
\description{
  The \code{Leaves} data frame has 15 rows and 2 columns of leaf length
  over time.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{Time}{
      time from initial emergence (days).
    }
    \item{Length}{
      leaf length (cm).
    }
  }
}
\source{
  Bates and Watts (1998), \emph{Nonlinear Regression Analysis and Its
    Applications}, Wiley (Appendix A4.5).

  Originally from Heyes and Brown (1956) in F.L. Milthorpe (ed),
  \emph{The Growth of Leaves}, London: Butterworth. 
}
\examples{
plot(Length ~ Time, data = Leaves, las = 1,
     xlab = "Time from initial emergence (days)",
     ylab = "Leaf length (cm)",
     main = "Leaves data and fitted curve")
## first fit a logistic model
fm1 <- nls(Length ~ SSlogis(Time, Asym, xmid, scal), data = Leaves)
summary(fm1)
usr <- par("usr")
xx <- seq(usr[1], usr[2], len = 51)
lines(xx, predict(fm1, list(Time = xx)))
title(sub = deparse(fm1$call$formula))
## compare with Richards growth model
\dontrun{
fm2 <- nls(Length ~ SSRichards(Time, Asym, xmid, scal, pow),
           data = Leaves)
anova( fm2, fm1 )
}
}
\keyword{datasets}
