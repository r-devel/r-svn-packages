\name{write.dbf}
\alias{write.dbf}
\title{Write a DBF File}
\description{
  The function tries to write a data frame to a DBF file.
}
\usage{
write.dbf(dataframe, file, factor2char = TRUE)
}
\arguments{
  \item{dataframe}{a data frame object.}
  \item{file}{a file name to be written to.}
  \item{factor2char}{logical, default \code{TRUE}, convert factor columns to
    character: otherwise they are written as the internal integer codes.}
}
\details{
  The function calls code from \code{shapelib} to write a DBF format file.  At
  present, \code{NA}s will cause an error and no file will be written.  If
  factors should be represented as character fields, the data frame must
  either be pre-processed to convert those colums using
  \code{I(as.character())}, or the \code{factor2char} argument set to
  \code{TRUE}.  In addition logical columns are converted to integer; maximum
  precision (number of digits including minus sign and decimal sign) for
  numeric is 19 - scale (digits after the decimal sign) is calculated
  internally based on the number of digits before the decimal sign.
}
\value{
  no return value.
}
\references{\url{http://shapelib.maptools.org/}}

\author{
  Nicholas J. Lewin-Koh, modified by Roger Bivand; shapelib by Frank
  Warmerdam.
}

\seealso{
  \code{\link{read.dbf}}
}
\examples{
str(warpbreaks)
try1 <- paste(tempfile(), ".dbf", sep = "")
write.dbf(warpbreaks, try1, factor2char = FALSE)
in1 <- read.dbf(try1)
str(in1)
try2 <- paste(tempfile(), ".dbf", sep = "")
write.dbf(warpbreaks, try2, factor2char = TRUE)
in2 <- read.dbf(try2)
str(in2)
}
\keyword{file}
