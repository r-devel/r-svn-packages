\name{gtk}
\alias{gtk}
\alias{GTK}
\title{GTK+ Graphics Device}
\description{
  A variant of the X11 graphics device.
}
\usage{
gtk(display="", width=7, height=7, pointsize=12)
GTK(display="", width=7, height=7, pointsize=12)
}
\arguments{
  \item{display}{the display on which the graphics window will appear.
    The default is to use the value in the user's environment variable
    \code{DISPLAY}.}
  \item{width}{the width of the plotting window in inches.}
  \item{height}{the height of the plotting window in inches.}
  \item{pointsize}{the default pointsize to be used.}
}
\author{Lyndon Drake} 
\references{
  \url{http://www.gtk.org} for the \code{GTK+} (GIMP Tool Kit) libraries.
}
\seealso{
  \code{\link{x11}},
  \code{\link{Devices}}.
}
\examples{
  gtk()
  plot(1:10)
  gtk()
  par(mfrow=c(2,2))
  x <- rnorm(100)
  plot(x)
  hist(x)
  boxplot(x)  
}
\keyword{device}
