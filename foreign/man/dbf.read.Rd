\name{dbf.read}
\alias{dbf.read}
\title{Read a DBF File}
\description{
  The function reads a DBF file into a data frame, converting character
  fields to factors, and trying to respect NULL fields.
}
\usage{
dbf.read(filen)
}
\arguments{
  \item{filen}{name of input file}
}

\value{
  A data frame of data from the DBF file; note that the field names are
  adjusted to use in R using \code{make.names()}, and so will permit the
  underscore character.
}

\references{\url{http://shapelib.maptools.org/}}

\author{
  Nicholas Lewin-Koh and  Roger Bivand; shapelib by Frank Warmerdam
}

\seealso{
  \code{\link{dbf.write}}
}
\examples{
x <- dbf.read(system.file("files/sids.dbf", package="foreign")[1])
str(x)
summary(x)
}
\keyword{file}
