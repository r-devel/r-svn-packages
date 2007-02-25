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
options(show.signif.stars = FALSE)
## first fit a logistic model
fm1 <- nls(Length ~ SSlogis(Time, Asym, xmid, scal), data = Leaves)
summary(fm1)
plotfit(fm1, xlab = "Time from initial emergence (days)",
        ylab = "Leaf length (cm)", main = "Logistic growth model")
## compare with Richards growth model
fm2 <- nls(Length ~ Asym/(1+exp(-(Time - xmid)/scal))^exp(-lpow), Leaves,
           c(coef(fm1),c(lpow = 0)))
summary(fm2)
anova(fm1, fm2)
plotfit(fm2, xlab = "Time from initial emergence (days)",
        ylab = "Leaf length (cm)", main = "Richards growth model")
pm1 <- profile(fm1)
plot(pm1, aspect = 'xy', layout = c(3,1))
plot(pm1, absVal = FALSE, aspect = 'xy', layout = c(3,1))
pm2 <- profile(fm2, alpha = 0.05)
plot(pm2, aspect = 'xy', layout = c(4,1))
plot(pm2, absVal = FALSE, aspect = 'xy', layout = c(4,1))
}
\keyword{datasets}
