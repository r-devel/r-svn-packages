\name{survexp.fit}
\alias{survexp.fit}
\title{
Work Function to Compute Expected Survival
}
\description{
Compute expected survival. This function is not to be called by the user.
}
\usage{
survexp.fit(x, y, times, death, ratetable)
}
\arguments{
\item{x}{
a matrix. The first column contains the group,
an integer value that divides the subjects into subsets.
Remaining columns must match the dimensions of the \code{ratetable},
in the correct order.
}
\item{y}{
the follow up time for each subject.
}
\item{times}{
the vector of times at which a result will be computed.
}
\item{death}{death indicator}
\item{ratetable}{
a rate table, such as \code{survexp.uswhite}.
}}
\value{
A list containing the number of subjects and the expected survival(s)
at each time point.
If there are multiple groups, these will be
matrices with one column per group.
}
\details{
  For cohort survival it must be the potential censoring time for
each subject, ignoring death.


For an exact estimate \code{times} should be a superset of \code{y}, so that each
subject at risk is at risk for the entire sub-interval of time.
For a large data set, however, this can use an inordinate amount of
storage and/or compute time.  If the \code{times} spacing is more coarse than
this, an actuarial approximation is used which should, however, be extremely
accurate as long as all of the returned values are > .99.


For a subgroup of size 1 and \code{times} > \code{y},
the conditional method reduces to exp(-h) where
h is the expected cumulative hazard for the subject over his/her
observation time.  This is used to compute individual expected survival.
}
\section{WARNING}{
Most users will call the higher level routine \code{survexp}.
Consequently, this function has very few error checks on its input arguments.
}
\seealso{
\code{\link{survexp}}, \code{\link{survexp.us}}
}
\keyword{survival}
% Converted by Sd2Rd version 0.3-2.
