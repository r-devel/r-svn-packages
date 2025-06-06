\name{plot.survfit}
\alias{plot.survfit}
\title{
Plot Method for `survfit'
}
\usage{
\method{plot}{survfit}(x, conf.int=, mark.time=TRUE,
 mark=3, col=1, lty=1, lwd=1, cex=1, log=FALSE, xscale=1, yscale=1, 
 firstx=0, firsty=1, xmax, ymin=0, fun,
 xlab="", ylab="", xaxs="S", bty=NULL,
legend.text=NULL,legend.pos=0,legend.bty="n",main=NULL,\dots)
}
\arguments{
\item{x}{
an object of class \code{survfit}, usually returned by the \code{survfit} function.
}
\item{conf.int}{
determines whether confidence intervals will be plotted.  The default is to
do so if there is only 1 curve, i.e., no strata.
}
\item{mark.time}{
controls the labeling of the curves.  If set to \code{FALSE}, no labeling is done.
If \code{TRUE}, then curves are marked at each censoring time which is not also 
a death time.  If \code{mark.time} is a
numeric vector, then curves are marked at the specified time points.
}
\item{mark}{
vector of mark parameters, which will be used to label the curves.
The \code{lines} help file contains examples of the possible marks.
The vector is reused cyclically if it is shorter than the number of curves.
}
\item{col}{
a vector of integers specifying colors for each curve.
The default value is 1.
}
\item{lty}{
a vector of integers specifying line types for each curve.
The default value is 1.
}
\item{lwd}{
a vector of numeric values for line widths. The default value is 1.
}
\item{cex}{
a numeric value specifying the size of the marks.
Not a vector; all marks have the same size.
}
\item{log}{
a logical value, if TRUE the y axis wll be on a log scale.
Alternately, one of the standard character strings "x", "y", or "xy"
can be given to specific logarithmic horizontal and/or vertical axes.
}
\item{xscale}{
scale the x-axis values before plotting.
If time were in days, then a value of
365.25 will give labels in years instead of the original days.
}
\item{yscale}{
will be used to multiply the labels on the y axis.
A value of 100, for instance, would be used to give a percent scale.
Only the labels are
changed, not the actual plot coordinates, so that adding a curve with
"\code{lines(surv.exp(\dots{}))}", say, 
will perform as it did without the \code{yscale} argument.
}
\item{firstx,firsty}{
the starting point for the survival curves.  If either of these is set to
\code{NA} the plot will start at the first time point of the curve.
}
\item{xmax}{
the maximum horizontal plot coordinate.  This can be used to shrink the range
of a plot.  It shortens the curve before plotting it, so that unlike using the
\code{xlim} graphical parameter, warning messages about out of bounds points are
not generated.
}
\item{ymin}{
lower boundary for y values.  Survival curves are most often drawn in the
range of 0-1, even if none of the curves approach zero.
The parameter is ignored if the \code{fun} argument is present, or if it has been
set to \code{NA}.
}
\item{fun}{
an arbitrary function defining a transformation of the survival curve.
For example \code{fun=log} is an alternative way to draw a log-survival curve
(but with the axis labeled with log(S) values),
and \code{fun=sqrt} would generate a curve on square root scale.
Four often used transformations can be specified with a character
argument instead: \code{"log"} is the same as using the \code{log=T} option,
\code{"event"} plots cumulative events (f(y) = 1-y), 
\code{"cumhaz"} plots the cumulative hazard function (f(y) = -log(y)), and
\code{"cloglog"} creates a complimentary log-log survival plot (f(y) =
log(-log(y)) along with log scale for the x-axis). 
}
\item{xlab}{
label given to the x-axis.
}
\item{ylab}{
label given to the y-axis.
}
\item{xaxs}{
either \code{"S"} for a survival curve or a standard x axis style as listed in \code{par}.
Survival curves are usually displayed with the curve touching the y-axis,
but not touching the bounding box of the plot on the other 3 sides.
Type \code{"S"} accomplishes this by manipulating the plot range and then using
the \code{"i"} style internally.
}\item{bty}{see \code{\link{par}}}\item{legend.pos}{position for the legend: a vector of length 2, or
  0 to put the legend in the lower left, 1 to put it in the upper
  right. These may not work well with transformed curves.}
\item{legend.bty}{Box type, see \code{\link{legend}}}
\item{legend.text}{Text for legend, see \code{\link{legend}}}
\item{main}{Plot title}
\item{...}{other graphical parameters}
}
\value{
a list with components \code{x} and \code{y}, containing the coordinates of the last point
on each of the curves (but not the confidence limits).  
This may be useful for labeling.
}
\description{
A plot of survival curves is produced, one curve for each strata.
The \code{log=T} option does extra work to avoid log(0), and to try to create a
pleasing result.  If there are zeros, they are plotted by default at
0.8 times the smallest non-zero value on the curve(s).

If \code{legend.text} is supplied a legend is created.
}
\section{BUGS}{
Survival curve objects created from a \code{coxph} model does not include
the censoring times. Therefore, specifying \code{mark.time=T} does not work.
If you want to mark censoring times on the curve(s) resulting 
from a \code{coxph} fit, provide a vector of times as the \code{mark.time} argument 
in the call to \code{plot.survfit} or \code{lines.survfit}. 
}
\seealso{
  \code{\link{par}}, \code{\link{survfit}}, \code{\link{lines.survfit}}
\code{\link{legend}}.
}
\examples{
leukemia.surv <- survfit(Surv(time, status) ~ x, data = aml)
plot(leukemia.surv, lty = 2:3)
legend(100, .9, c("Maintenance", "No Maintenance"), lty = 2:3)
title("Kaplan-Meier Curves\nfor AML Maintenance Study")


lsurv2 <- survfit(Surv(time, status) ~ x, aml, type='fleming')
plot(lsurv2, lty=2:3, fun="cumhaz",
	xlab="Months", ylab="Cumulative Hazard")

plot(leukemia.surv,lty=1:2,legend.pos=0,col=c("Red","Blue"),legend.text=c("Maintenance",
"No Maintenance"))
if (interactive()){
cat("Click on the graphics device to place a legend\n")
plot(leukemia.surv,lty=2:3,legend.pos=locator(1),legend.text=c("Maintenance",
"No Maintenance"))
}
}
\keyword{survival}
\keyword{hplot}
% Converted by Sd2Rd version 0.3-2.
