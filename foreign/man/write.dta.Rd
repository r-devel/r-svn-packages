\name{write.dta}
\alias{write.dta}
\title{Write files in Stata binary format}
\usage{
write.dta(dataframe, file, version = 6,convert.dates=TRUE,tz="GMT")
}
\arguments{
  \item{dataframe}{a data frame.}
  \item{file}{character string giving filename.}
  \item{version}{Stata version: 6 and 7 are supported.}
  \item{convert.dates}{Convert \code{POSIXt} objects to Stata dates?}
  \item{tz}{timezone for date conversion}
  } 
\description{
  Writes the data frame to file in the Stata version 6.0 or 7.0 binary
  format.  Does not write matrix variables.
}
\details{
  The columns in the data frame become variables in the Stata data set.
  Missing values are correctly handled.  Nothing is done with factor
  levels, which should end up as variable labels.  The major difference
  between versions is that 7.0 allows 32-character variable names.

  Optionally, R date/time objects (\code{POSIXt} classes) are converted
  into the Stata format. This loses information -- Stata dates are in
  days since 1960-1-1.  \code{POSIXct} objects can be written without
  conversion but will not be understood as dates by Stata;
  \code{POSIXlt} objects cannot be written.
  }
\value{
  \code{NULL}
}
\references{
  Stata 6.0 Users Manual, Stata 7.0 Programming manual describe the file
  formats.
} 
\author{Thomas Lumley}
\seealso{
  \code{\link{read.dta}},
  \code{\link{attributes}}
  \code{\link{DateTimeClasses}}
}
\examples{
data(swiss)
write.dta(swiss,swissfile<-tempfile())
read.dta(swissfile)
}
\keyword{file}
