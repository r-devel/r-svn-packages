\name{sqlite.data.frame}
\alias{sqlite.data.frame}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Creates an sqlite data frame from ordinary data frame.}
\description{
This function creates an sqlite data frame from the as.data.frame()-d argument.
}
\usage{
sqlite.data.frame(x)
}
%- maybe also 'usage' for other objects documented here.
\details{

}
\value{
  \item{obj}{The object to be converted into a data frame and then stored in a SQLite database.}
}
\references{ ~put references to the literature/web site here ~ }
\author{Miguel A. R. Manese}
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
    iris.sdf <- sqlite.data.frame(iris)
    #names(iris.sdf)
    #apply(iris.sdf[1:4], 2, mean)

}
