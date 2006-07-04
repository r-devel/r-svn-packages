\name{getSdf}
\alias{getSdf}
\title{ Get a handle to an SDF in the workspace. }
\description{
Gets a handle to an SDF registered in the workspace. The returned
handle is a sqlite.data.frame object which can be used with the most 
common operators for data frames.
}
\usage{
sdf <- getSdf(name)
}
\arguments{
  \item{name}{ The name of the SDF. }
}
\details{
When SQLiteDF is loaded, it creates a SQLite data base in the file
workspace.db in the working directory (note that databases in SQLite
are all contained in one file). This database contains a list of SQLite
databases containing SDF data that has been attached.
}
\value{
Returns an sqlite.data.frame object.
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ 
    \code{\link[SQLiteDF]{sqlite.data.frame}}
    \code{\link[SQLiteDF]{lsSdf}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
#function (name) 
#.Call("sdf_get_sdf", name)
}
