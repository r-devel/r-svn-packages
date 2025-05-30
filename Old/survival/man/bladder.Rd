\name{bladder}
\docType{data}
\alias{bladder}
\alias{bladder2}
\title{Bladder Cancer Recurrences}
\usage{
bladder
bladder2
}
\description{Data on recurrences of bladder cancer, used by many people
  to demonstrate methodology for recurrent event modelling.}
\format{
  bladder
  \tabular{ll}{
    id:\tab Patient id\cr
    rx:\tab Treatment 1=placebo  2=thiotepa\cr
    number:\tab Initial number of tumours (8=8 or more)\cr
    size:\tab size (cm) of largest initial tumour\cr
    stop:\tab recurrence or censoring time\cr
    enum:\tab which recurrence (up to 4)\cr
  }
  bladder2 
  \tabular{ll}{
    id:\tab Patient id\cr
    rx:\tab Treatment 1=placebo  2=thiotepa\cr
    number:\tab Initial number of tumours (8=8 or more)\cr
    size:\tab size (cm) of largest initial tumour\cr
    start:\tab start of interval (0 or previous recurrence time)\cr
    stop:\tab recurrence or censoring time\cr
    enum:\tab which recurrence (up to 4)\cr
  }
  
}
\source{
  LJ Wei, DY Lin, L Weissfeld (1989),
  Regression analysis of multivariate incomplete failure time data by
  modeling marginal distributions.
  \emph{Journal of the American Statistical Association},
  \bold{84}.
}
\keyword{datasets}
\keyword{survival}
