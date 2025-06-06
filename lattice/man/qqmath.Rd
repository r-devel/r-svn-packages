\name{B_04_qqmath}
\alias{qqmath}
\alias{qqmath.formula}
\alias{qqmath.numeric}
\title{ Q-Q Plot with Theoretical Distribution }
\description{
  Quantile-Quantile plot of a sample and a theoretical distribution
}
\usage{
qqmath(x, data, \dots)

\method{qqmath}{formula}(x,
       data,
       allow.multiple = is.null(groups) || outer,
       outer = !is.null(groups),
       distribution = qnorm,
       f.value = NULL,
       auto.key = FALSE,
       aspect = "fill",
       panel = lattice.getOption("panel.qqmath"),
       prepanel = NULL,
       scales, strip, groups,
       xlab, xlim, ylab, ylim,
       drop.unused.levels = lattice.getOption("drop.unused.levels"),
       \dots,
       lattice.options = NULL,
       default.scales = list(),
       subscripts,
       subset)
\method{qqmath}{numeric}(x, data = NULL, ylab, \dots)
}
\arguments{
  \item{x}{
    The object on which method dispatch is carried out.

    For the \code{"formula"} method, a formula of the form \code{~ x |
      g1 * g2 * \ldots}, where \code{x} must be a numeric.  For the
    \code{"numeric"} method, a numeric vector.
  }
  \item{data}{
    For the \code{formula} method, an optional data frame in which
    variables in the formula (as well as \code{groups} and
    \code{subset}, if any) are to be evaluated.  Usualll ignored with a
    warning in other methods.
  }
  \item{distribution}{ a quantile function that takes a vector of
    probabilities as argument and produces the corresponding
    quantiles. Possible values are \code{qnorm}, \code{qunif} etc.
    Distributions with other required arguments need to be passed in as
    user defined functions.
  }
  \item{f.value}{ optional numeric vector of probabilities, quantiles
    corresponding to which should be plotted.  Can also be a function of
    a single integer (representing sample size) that returns such a
    numeric vector.  The typical value for this argument is the function
    \code{ppoints}, which is also the S-PLUS default.  If specified, the
    probabilities generated by this function is used for the plotted
    quantiles, using the \code{quantile} function for the sample, and
    the function specified as the \code{distribution} argument for the
    theoretical distribution.  

    \code{f.value} defaults to \code{NULL}, which has the effect of
    using \code{ppoints} for the quantiles of the theoretical
    distribution, but the exact data values for the sample.  This is
    similar to what happens for \code{qqnorm}, but different from the
    S-PLUS default of \code{f.value=ppoints}.

    For large \code{x}, this argument can be useful in plotting a
    smaller set of quantiles, which is usually enough to capture the
    pattern.
  }
  \item{panel}{
    The panel function to be used.  Unlike in older versions, the
    default panel function does most of the actual computations and has
    support for grouping.  See \code{\link{panel.qqmath}} for details.
  }
  \item{allow.multiple, outer, auto.key, aspect, prepanel, scales,
    strip, groups, xlab, xlim, ylab, ylim, drop.unused.levels,
    lattice.options, default.scales, subscripts, subset}{ See
    \code{\link{xyplot}} }
  \item{\dots}{
    Further arguments.  See corresponding entry in \code{\link{xyplot}}
    for non-trivial details.
  }
}
\value{

  An object of class \code{"trellis"}. The
  \code{\link[lattice:update.trellis]{update}} method can be used to
  update components of the object and the
  \code{\link[lattice:print.trellis]{print}} method (usually called by
  default) will plot it on an appropriate plotting device.

}
\details{
  \code{qqmath} produces a Q-Q plot of the given sample and a
  theoretical distribution.  The default behaviour of \code{qqmath} is
  different from the corresponding S-PLUS function, but is similar to
  \code{qqnorm}. See the entry for \code{f.value} for specifics.

  The implementation details are also different from S-PLUS.  In
  particular, all the important calculations are done by the panel (and
  prepanel function) and not \code{qqmath} itself.  In fact, both the
  arguments \code{distribution} and \code{f.value} are passed unchanged
  to the panel and prepanel function.  This allows, among other things,
  display of grouped Q-Q plots, which are often useful.  See the help
  page for \code{\link{panel.qqmath}} for further details.

  This and all other high level Trellis functions have several
  arguments in common. These are extensively documented only in the
  help page for \code{xyplot}, which should be consulted to learn more
  detailed usage.
}
\author{Deepayan Sarkar \email{Deepayan.Sarkar@R-project.org}}
\seealso{\code{\link{xyplot}}, \code{\link{panel.qqmath}},
  \code{\link{panel.qqmathline}}, \code{\link{prepanel.qqmathline}},
  \code{\link{Lattice}}, \code{\link{quantile}}
}
\examples{
qqmath(~ rnorm(100), distribution = function(p) qt(p, df = 10))
qqmath(~ height | voice.part, aspect = "xy", data = singer,
       prepanel = prepanel.qqmathline,
       panel = function(x, ...) {
          panel.qqmathline(x, ...)
          panel.qqmath(x, ...)
       })
vp.comb <-
    factor(sapply(strsplit(as.character(singer$voice.part), split = " "),
                  "[", 1),
           levels = c("Bass", "Tenor", "Alto", "Soprano"))
vp.group <-
    factor(sapply(strsplit(as.character(singer$voice.part), split = " "),
                  "[", 2))
qqmath(~ height | vp.comb, data = singer,
       groups = vp.group, auto.key = list(space = "right"),
       aspect = "xy",
       prepanel = prepanel.qqmathline,
       panel = function(x, ...) {
          panel.qqmathline(x, ...)
          panel.qqmath(x, ...)
       })
}
\keyword{dplot}

