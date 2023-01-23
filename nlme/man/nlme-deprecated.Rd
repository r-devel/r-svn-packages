% File nlme/man/nlme-deprecated.Rd
% Part of the nlme package for R
% Distributed under GPL 2 or later: see nlme/LICENCE.note

\name{nlme-deprecated}
\alias{nlme-deprecated}
\title{Deprecated Functions in Package \pkg{nlme}}

%% said to be deprecated in 2018-08, with warning since 2023-01
\alias{nfGroupedData}
\alias{nmGroupedData}

\usage{
## internal functions wrapped by groupedData()
nfGroupedData(formula, data, order.groups, FUN, outer, inner,
              labels, units)
nmGroupedData(formula, data, order.groups, FUN, outer, inner,
              labels, units)
}

\description{
  These functions are provided for compatibility with older versions of
  \pkg{nlme} only, and may be defunct as soon as the next release.
}

\keyword{internal}
