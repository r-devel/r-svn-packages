% File nlme/man/collapse.Rd
% Part of the nlme package for R
% Distributed under GPL 2 or later: see nlme/LICENCE.note

\name{collapse}
\title{Collapse According to Groups}
\usage{
collapse(object, \dots)
}
\alias{collapse}
\arguments{
 \item{object}{an object to be collapsed, usually a data frame.}
 \item{\dots}{some methods for the generic may require additional
   arguments.}
}
\description{
  This function is generic; method functions can be written to handle
  specific classes of objects. Currently, only a \code{groupedData}
  method is available. 
}
\value{
  will depend on the method function used; see the appropriate documentation.
}
\author{José Pinheiro and Douglas Bates \email{bates@stat.wisc.edu}}

\seealso{\code{\link{collapse.groupedData}}}

\keyword{models}
