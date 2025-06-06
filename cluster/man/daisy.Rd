\name{daisy}
\alias{daisy}
\title{Dissimilarity Matrix Calculation}
\concept{Gower's formula}
\concept{Gower's distance}
\concept{Gower's coefficient}% FIXME: see ../TODO-MM
\description{
  Compute all the pairwise dissimilarities (distances) between observations
  in the data set.  The original variables may be of mixed types.  In
  that case, or whenever \code{metric = "gower"} is set, a
  generalization of Gower's formula is used, see \sQuote{Details}
  below.
}
\usage{
daisy(x, metric = c("euclidean", "manhattan", "gower"),
      stand = FALSE, type = list(), weights = rep.int(1, p),
      warnBin = warnType, warnAsym = warnType, warnConst = warnType,
      warnType = TRUE)
}
\arguments{
  \item{x}{
    numeric matrix or data frame, of dimension \eqn{n\times p}{n x p},
    say.  Dissimilarities will be computed
    between the rows of \code{x}.  Columns of mode \code{numeric}
    (i.e. all columns when \code{x} is a matrix) will be recognized as
    interval scaled variables, columns of class \code{factor} will be
    recognized as nominal variables, and columns of class \code{ordered}
    will be recognized as ordinal variables.  Other variable types
    should be specified with the \code{type} argument.  Missing values
    (\code{\link{NA}}s) are allowed.
  }
  \item{metric}{
    character string specifying the metric to be used.
    The currently available options are \code{"euclidean"} (the default),
    \code{"manhattan"} and \code{"gower"}.\cr
    Euclidean distances are root sum-of-squares of differences, and
    manhattan distances are the sum of absolute differences.

    \dQuote{Gower's distance} is chosen by metric \code{"gower"}
    or automatically if some columns of \code{x} are not numeric. Also
    known as Gower's coefficient (1971),
    expressed as a dissimilarity, this implies that a particular
    standardisation will be applied to each variable, and the
    \dQuote{distance} between two units is the sum of all the
    variable-specific distances, see the details section.
  }
  \item{stand}{logical flag: if TRUE, then the measurements in \code{x}
    are standardized before calculating the
    dissimilarities.  Measurements are standardized for each variable
    (column), by subtracting the variable's mean value and dividing by
    the variable's mean absolute deviation.

    If not all columns of \code{x} are numeric, \code{stand} will
    be ignored and Gower's standardization (based on the
    \code{\link{range}}) will be applied in any case, see argument
    \code{metric}, above, and the details section.
  }
  \item{type}{list for specifying some (or all) of the types of the
    variables (columns) in \code{x}.  The list may contain the following
    components:
    \describe{
      \item{\code{"asymm"}  }{\bold{A}symmetric binary variable, aka
	\code{"A"} in result \code{Types}, see \code{\link{dissimilarity.object}}.}
      \item{\code{"symm"}   }{\bold{S}ymmetric  binary variable, aka \code{"S"}.}
      \item{\code{"factor"} }{\bold{N}ominal -- the default for \code{\link{factor}}
	variables, aka \code{"N"}.  When the factor has 2 levels, this is
	equivalent to \code{type = "S"} for a (symmetric) binary variable.}
      \item{\code{"ordered"}}{\bold{O}rdinal -- the default for \code{\link{ordered}}
	(factor) variables, aka \code{"O"}, see \code{\link{dissimilarity.object}}.}
      \item{\code{"logratio"}}{ratio scaled numeric variables that are to
	be logarithmically transformed (\code{\link{log10}}) and then
	treated as numeric (\code{"I"}): must be \emph{positive} numeric variable.}
      \item{\code{"ordratio"}}{\dQuote{ra\bold{T}io}-like
	variable to be treated as \code{\link{ordered}} (using the factor
	codes \code{unclass(\link{as.ordered}(x[,j]))}), aka \code{"T"}.}
      \item{\code{"numeric"}/\code{"integer"}}{\bold{I}nterval
	scaled -- the \bold{default} for all numeric (incl \code{integer})
	columns of \code{x}, aka \code{"I"} in result \code{Types}, see
	\code{\link{dissimilarity.object}}.}
    }
    Each component is a (character or numeric) vector, containing either
    the names or the numbers of the corresponding columns of \code{x}.

    Variables not mentioned in \code{type} are interpreted as usual, see
    argument \code{x}, and also \sQuote{default} above.  Consequently,
    the default \code{type = list()} may often be sufficient.
  }
  \item{weights}{an optional numeric vector of length \eqn{p}(=\code{ncol(x)}); to
    be used in \dQuote{case 2} (mixed variables, or \code{metric = "gower"}),
    specifying a weight for each variable (\code{x[,k]}) instead of
    \eqn{1} in Gower's original formula.}
  \item{warnBin, warnAsym, warnConst}{logicals indicating if the
    corresponding type checking warnings should be signalled (when found).}
  \item{warnType}{logical indicating if \emph{all} the type checking
    warnings should be active or not.}
}% end{arg..}

\value{
  an object of class \code{"dissimilarity"} containing the
  dissimilarities among the rows of \code{x}.  This is typically the
  input for the functions \code{pam}, \code{fanny}, \code{agnes} or
  \code{diana}.  For more details, see \code{\link{dissimilarity.object}}.
}
\details{
  The original version of \code{daisy} is fully described in chapter 1
  of Kaufman and Rousseeuw (1990).
  Compared to \code{\link{dist}} whose input must be numeric
  variables, the main feature of \code{daisy} is its ability to handle
  other variable types as well (e.g. nominal, ordinal, (a)symmetric
  binary) even when different types occur in the same data set.

  The handling of nominal, ordinal, and (a)symmetric binary data is
  achieved by using the general dissimilarity coefficient of Gower
  (1971).  If \code{x} contains any columns of these
  data-types, both arguments \code{metric} and \code{stand} will be
  ignored and Gower's coefficient will be used as the metric.  This can
  also be activated for purely numeric data by \code{metric = "gower"}.
  With that, each variable (column) is first standardized by dividing
  each entry by the range of the corresponding variable, after
  subtracting the minimum value; consequently the rescaled variable has
  range \eqn{[0,1]}, exactly.
  %% FIXME: Use something like  "gowerRob" which uses *robust* rescaling

  Note that setting the type to \code{symm} (symmetric binary) gives the
  same dissimilarities as using \emph{nominal} (which is chosen for
  non-ordered factors) only when no missing values are present, and more
  efficiently.

  Note that \code{daisy} signals a warning when 2-valued numerical
  variables do not have an explicit \code{type} specified, because the
  reference authors recommend to consider using \code{"asymm"}; the
  warning may be silenced by \code{warnBin = FALSE}.

  In the \code{daisy} algorithm, missing values in a row of x are not
  included in the dissimilarities involving that row.  There are two
  main cases,
  \enumerate{
    \item If all variables are interval scaled (and \code{metric} is
      \emph{not} \code{"gower"}), the metric is "euclidean", and
      \eqn{n_g} is the number of columns in which
      neither row i and j have NAs, then the dissimilarity d(i,j) returned is
      \eqn{\sqrt{p/n_g}}{sqrt(p/n_g)} (\eqn{p=}ncol(x)) times the
      Euclidean distance between the two vectors of length \eqn{n_g}
      shortened to exclude NAs.  The rule is similar for the "manhattan"
      metric, except that the coefficient is \eqn{p/n_g}.  If \eqn{n_g = 0},
      the dissimilarity is NA.

    \item When some variables have a type other than interval scaled, or
      if \code{metric = "gower"} is specified, the
      dissimilarity between two rows is the weighted mean of the contributions of
      each variable.  Specifically,
      \deqn{d_{ij} = d(i,j) = \frac{\sum_{k=1}^p w_k \delta_{ij}^{(k)} d_{ij}^{(k)}}{
	  \sum_{k=1}^p w_k \delta_{ij}^{(k)}}.
      }{d_ij = d(i,j) =	sum(k=1:p; w_k delta(ij;k) d(ij,k)) / sum(k=1:p; w_k delta(ij;k)).}
      In other words, \eqn{d_{ij}}{d_ij} is a weighted mean of
      \eqn{d_{ij}^{(k)}}{d(ij,k)} with weights \eqn{w_k \delta_{ij}^{(k)}}{w_k delta(ij;k)},
      where \eqn{w_k}\code{= weigths[k]},
      \eqn{\delta_{ij}^{(k)}}{delta(ij;k)} is 0 or 1, and
      \eqn{d_{ij}^{(k)}}{d(ij,k)}, the k-th variable contribution to the
      total distance, is a distance between \code{x[i,k]} and \code{x[j,k]},
      see below.

      The 0-1 weight \eqn{\delta_{ij}^{(k)}}{delta(ij;k)} becomes zero
      when the variable \code{x[,k]} is missing in either or both rows
      (i and j), or when the variable is asymmetric binary and both
      values are zero.  In all other situations it is 1.

      The contribution \eqn{d_{ij}^{(k)}}{d(ij,k)} of a nominal or binary variable to the total
      dissimilarity is 0 if both values are equal, 1 otherwise.
      The contribution of other variables is the absolute difference of
      both values, divided by the total range of that variable.  Note
      that \dQuote{standard scoring} is applied to ordinal variables,
      i.e., they are replaced by their integer codes \code{1:K}.  Note
      that this is not the same as using their ranks (since there
      typically are ties).
      % contrary to what Kaufman and Rousseeuw write in their book, and
      % the original help page.

      As the individual contributions \eqn{d_{ij}^{(k)}}{d(ij,k)} are in
      \eqn{[0,1]}, the dissimilarity \eqn{d_{ij}}{d_ij} will remain in
      this range.
      If all weights \eqn{w_k \delta_{ij}^{(k)}}{w_k delta(ij;k)} are zero,
      the dissimilarity is set to \code{\link{NA}}.
  }
}
\section{Background}{
  Dissimilarities are used as inputs to cluster analysis and
  multidimensional scaling.  The choice of metric may have a
  large impact.
}
\references{
  Gower, J. C. (1971)
  A general coefficient of similarity and some of its properties,
  \emph{Biometrics} \bold{27}, 857--874.

  Kaufman, L. and Rousseeuw, P.J. (1990)
  \emph{Finding Groups in Data: An Introduction to Cluster Analysis}.
  Wiley, New York.

  Struyf, A., Hubert, M. and Rousseeuw, P.J. (1997)
  Integrating Robust Clustering Techniques in S-PLUS,
  \emph{Computational Statistics and Data Analysis} \bold{26}, 17--37.
}
\author{
  Anja Struyf, Mia Hubert, and Peter and Rousseeuw, for the original
  version.
  \cr
  Martin Maechler improved the \code{\link{NA}} handling and
  \code{type} specification checking, and extended functionality to
  \code{metric = "gower"} and the optional \code{weights} argument.
}
\seealso{
  \code{\link{dissimilarity.object}}, \code{\link{dist}},
  \code{\link{pam}}, \code{\link{fanny}}, \code{\link{clara}},
  \code{\link{agnes}}, \code{\link{diana}}.
}
\examples{% NB: >>> ../tests/daisy-ex.R  <<<<
data(agriculture)
## Example 1 in ref:
##  Dissimilarities using Euclidean metric and without standardization
d.agr <- daisy(agriculture, metric = "euclidean", stand = FALSE)
d.agr
as.matrix(d.agr)[,"DK"] # via as.matrix.dist(.)
## compare with
as.matrix(daisy(agriculture, metric = "gower"))

## Example 2 in reference, extended  ---  different ways of "mixed" / "gower":

example(flower) # -> data(flower) *and* provide 'flowerN'

summary(d0    <- daisy(flower))  # -> the first 3 {0,1} treated as *N*ominal
summary(dS123 <- daisy(flower,  type = list(symm = 1:3))) # first 3 treated as *S*ymmetric
stopifnot(dS123 == d0) # i.e.,  *S*ymmetric <==> *N*ominal {for 2-level factor}
summary(dNS123<- daisy(flowerN, type = list(symm = 1:3)))
stopifnot(dS123 == d0)
## by default, however ...
summary(dA123 <- daisy(flowerN)) # .. all 3 logicals treated *A*symmetric binary (w/ warning)
summary(dA3  <- daisy(flower, type = list(asymm = 3)))
summary(dA13 <- daisy(flower, type = list(asymm = c(1, 3), ordratio = 7)))
## Mixing variable *names* and column numbers (failed in the past):
summary(dfl3 <- daisy(flower, type = list(asymm = c("V1", "V3"), symm= 2,
                                          ordratio= 7, logratio= 8)))

## If we'd treat the first 3 as simple {0,1}
Nflow <- flower
Nflow[,1:3] <- lapply(flower[,1:3], function(f) as.integer(as.character(f)))
summary(dN <- daisy(Nflow)) # w/ warning: treated binary .. 1:3 as interval
## Still, using Euclidean/Manhattan distance for {0-1} *is* identical to treating them as "N" :
stopifnot(dN == d0)
stopifnot(dN == daisy(Nflow, type = list(symm = 1:3))) # or as "S"
}
\keyword{cluster}
