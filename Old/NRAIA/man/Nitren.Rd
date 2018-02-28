\name{Nitren}
\alias{Nitren}
\non_function{}
\title{Nitrendipene data}
\description{
  The \code{Nitren} data frame has 26 rows and 5 columns from an experiment
  in cardiology.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{log.NIF}{
      a numeric vector giving the log of the NIF concentration
    }
    \item{tiss1}{
      a numeric vector giving the reaction in tissue 1.
    }
    \item{tiss2}{
      a numeric vector giving the reaction in tissue 2.
    }
    \item{tiss3}{a
      a numeric vector giving the reaction in tissue 3.
    }
    \item{tiss4}{
      a numeric vector giving the reaction in tissue 4.
    }
  }
}
\source{
  Bates and Watts (1998), \emph{Nonlinear Regression Analysis and Its
    Applications}, Wiley (Appendix A1.5).
}
\examples{
matplot(Nitren[, 1], Nitren[, -1], las = 1,
    xlab = "log(NIF concentration)",
    ylab = "reaction level",
    main = paste("Nitren data and fitted curves,",
         "omitting zero concentration data"))
## without the data at NIF concentration of zero
options(na.action = na.omit)
fm1 <- nls(tiss1 ~ SSfpl(log.NIF, A, B, xmid, scal), Nitren)
fm1
fm2 <- nls(tiss2 ~ SSfpl(log.NIF, A, B, xmid, scal), Nitren)
fm2
fm3 <- nls(tiss3 ~ SSfpl(log.NIF, A, B, xmid, scal), Nitren)
fm3
fm4 <- nls(tiss4 ~ SSfpl(log.NIF, A, B, xmid, scal), Nitren)
fm4
usr <- par("usr")
xx <- seq(usr[1], usr[2], len = 50)
lines(xx, predict(fm1, list(log.NIF = xx)), col = 1, lty = 2)
lines(xx, predict(fm2, list(log.NIF = xx)), col = 2, lty = 2)
lines(xx, predict(fm3, list(log.NIF = xx)), col = 3, lty = 2)
lines(xx, predict(fm4, list(log.NIF = xx)), col = 4, lty = 2)
title(sub = deparse(fm1$call$formula))
## replacing the data at NIF concentration of zero by a very small value
log.NIF <- Nitren[, 1]
log.NIF[ is.na(log.NIF) ] <- -18
Nitren[, 1] <- log.NIF
matplot(Nitren[, 1], Nitren[, -1], las = 1,
    xlab = "log(NIF concentration)",
    ylab = "reaction level",
    main = paste("Nitren data and fitted curves",
       "- zero concentration recoded as -18"))
fm1 <- nls(tiss1 ~ SSfpl(log.NIF, A, B, xmid, scal), Nitren)
fm1
fm2 <- nls(tiss2 ~ SSfpl(log.NIF, A, B, xmid, scal), Nitren)
fm2
fm3 <- nls(tiss3 ~ SSfpl(log.NIF, A, B, xmid, scal), Nitren)
fm3
fm4 <- nls(tiss4 ~ SSfpl(log.NIF, A, B, xmid, scal), Nitren)
fm4
usr <- par("usr")
xx <- seq(usr[1], usr[2], len = 50)
lines(xx, predict(fm1, list(log.NIF = xx)), col = 1, lty = 2)
lines(xx, predict(fm2, list(log.NIF = xx)), col = 2, lty = 2)
lines(xx, predict(fm3, list(log.NIF = xx)), col = 3, lty = 2)
lines(xx, predict(fm4, list(log.NIF = xx)), col = 4, lty = 2)
title(sub = deparse(fm1$call$formula))
}
\keyword{datasets}

