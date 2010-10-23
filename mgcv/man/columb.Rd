\name{columb}
\alias{columb}
\alias{columb.polys}
\docType{data}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Reduced version of Columbus OH crime data}

\description{By district crime data from Columbus OH, together with polygons describing district shape. Useful for 
illustrating use of simple Markov Random Field smoothers.
}
\usage{
data(columb)
data(columb.polys)
}
%- maybe also `usage' for other objects documented here.

\format{ \code{columb} is a 49 row data frame with the following columns 
\describe{
\item{area}{land area of district}
\item{home.value}{housing value in 1000USD.}
\item{income}{household income in 1000USD.}
\item{crime}{residential burglaries and auto thefts per 1000 households.}
\item{open.space}{measure of open space in district.}
\item{district}{code identifying district, and matching \code{names(columb.polys)}. }
}
\code{columb.polys} contains the polygons defining the areas in the format described below.
} 


\details{The data frame \code{columb} relates to the districts whose boundaries are coded in \code{columb.polys}.
\code{columb.polys[[i]]} contains a list, each element of which describes a polygon delineating part of district i, 
which has identification code \code{names(columb.polys)[i]}. The jth polygon for the ith district is described by the list \code{columb.polys[[i]][[j]]}, which has elements: \code{coords} a two column matrix of co-ordinates of 
polygon nodes; \code{hole} which is \code{FALSE} if this polygon represents a hole in the district; \code{plotpt} which is a point suitable for treating as a typical location in the district for plotting purposes.
}

\source{
The data are adapted from the \code{columbus} example in the \code{spdep} package, where the original source is given as: 

Anselin, Luc. 1988. Spatial econometrics: methods and models. Dordrecht: Kluwer Academic, Table 12.1 p. 189.

}



\examples{
## see ?mrf help files
}








