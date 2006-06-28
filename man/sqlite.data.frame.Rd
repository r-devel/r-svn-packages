\name{sqlite.data.frame}
\alias{sqlite.data.frame}
\alias{sdf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{SQLite Data Frames}
\description{ Creates an sqlite data frame (SDF) from ordinary data frame.  }
\usage{
sqlite.data.frame(x, name=NULL)
}
\details{

}
\value{
  \item{x}{The object to be converted into a data frame and then stored in a 
SQLite database. }
  \item{name}{The name of the SDF. If none is provided, a generic
name \emph{data}<n> is used. Names will be appended with numbers to prevent
duplicates. E.g. if name arg is \emph{iris}, and it already exists, then the
new SDF will have a name \emph{iris1}. If it still exists, then the name
will be \emph{iris2}, and so on.
}
\references{ ~put references to the literature/web site here ~ }
\author{Miguel A. R. Manese}
\note{
A SDF is created by creating a sqlite database, with the specified name
as the filename with a db extension (e.g. data1.db). Each SDF is contained in
a separate database. The name of the stored in the SDF.
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
    iris.sdf <- sqlite.data.frame(iris)
    #names(iris.sdf)
    #apply(iris.sdf[1:4], 2, mean)

}
