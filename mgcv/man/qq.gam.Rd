\name{qq.gam}
\alias{qq.gam}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{QQ plots for gam model residuals}
\description{ Takes a fitted \code{gam} object produced by \code{gam()} and produces QQ plots of its residuals (conditional 
on the fitted model coefficients and scale parameter). If the model distributional assumptions are met then usually these plots should be close to a 
straight line (although discrete data can yield marked random departures from this line). 
}
\usage{
qq.gam(object,rep=0,level=.9,type="deviance",pch=".",rl.col=2,rep.col="gray80",...)
}
%- maybe also `usage' for other objects documented here.
\arguments{ 
\item{object}{ a fitted \code{gam} object as produced by \code{gam()} (or a \code{glm} object).}
\item{rep}{How many replicate datasets to generate to simulate quantiles of the residual distribution. 0 
           results in an efficient simulation free method for direct calculation, if this is possible for the 
object family. }
\item{level}{If simulation is used for the quantiles, then reference intervals can be provided for the QQ-plot, this specifies the level. 
             0 or less for no intervals, 1 or more to simply plot the QQ plot for each replicate generated.}
\item{type}{What sort of residuals should be plotted?}
\item{pch}{plot character to use. 19 is good.}
\item{rl.col}{color for the reference line on the plot.}
\item{rep.col}{color for reference bands or replicate reference plots.}
\item{...}{extra graphics parameters to pass to plotting functions.}
}
\details{QQ-plots of the the model residuals can be produced in one of two ways. The cheapest method generates reference quantiles by 
ordering the deviance residuals, associating a quantile of the uniform distribution with the corresponding datum according to the ordering, 
and then feeding these uniform quantiles into the quantile function associated with each datum. The resulting quantiles are then used in place of 
each datum to generate residuals. When sorted, these give approximate quantiles for the residual distribution expected for the fitted model.

The second method is to simulate. For each replicate, data are simulated from the fitted model, and the corresponding residuals computed. 
Quantiles are readily obtained from the empirical distribution of residuals so obtained. From this method reference bands are also computable.  

Even if \code{rep} is set to zero, the routine will attempt to simulate quantiles if no quantile function is available for the family. If no random 
deviate generating function family is available (e.g. for the quasi families), then a normal QQ-plot is produced. The routine conditions on the 
fitted model coefficents and the scale parameter estimate. 

The plots are similar in spirit to those proposed in Ben and Yohai (2004), but are substantially cheaper to produce.  

It is important to be aware that plots for raw residuals from fits to binary data contain no useful information 
about model fit (for example the deviance residuals for binary data depend only on the 
fitted values and not on the data in any other way). To check distributional 
assumptions about binary residuals you have to be able to group the data somehow. Binomial models other than binary are ok.   
}

\references{

M.G. Ben and V.J. Yohai (2004) JCGS 13(1), 36-47. (but note that binary residuals actually contain no information for model checking.)

\url{http://www.maths.bath.ac.uk/~sw283/}

}
\author{ Simon N. Wood \email{simon.wood@r-project.org}}


\seealso{ \code{\link{choose.k}},  \code{\link{gam}}, \code{\link{mgcv}}, \code{\link{magic}}}

\examples{

library(mgcv)
## simulate binomial data...
set.seed(0)
n.samp <- 400
dat <- gamSim(1,n=n.samp,dist="binary",scale=.33)
p <- binomial()$linkinv(dat$f) ## binomial p
n <- sample(c(1,3),n.samp,replace=TRUE) ## binomial n
dat$y <- rbinom(n,n,p)
dat$n <- n

lr.fit <- gam(y/n~s(x0)+s(x1)+s(x2)+s(x3)
             ,family=binomial,data=dat,weights=n,method="REML")

par(mfrow=c(2,2))
## normal QQ-plot of deviance residuals
qqnorm(residuals(lr.fit),pch=19,cex=.3)
## Quick QQ-plot of deviance residuals
qq.gam(lr.fit,pch=19,cex=.3)
## Simulation based QQ-plot with reference bands 
qq.gam(lr.fit,rep=100,level=.9)
## Simulation based QQ-plot, Pearson resids, all
## simulated reference plots shown...  
qq.gam(lr.fit,rep=100,level=1,type="pearson",pch=19,cex=.2)

## Now fit the wrong model and check....

pif <- gam(y~s(x0)+s(x1)+s(x2)+s(x3)
             ,family=poisson,data=dat,method="REML")
par(mfrow=c(2,2))
qqnorm(residuals(pif),pch=19,cex=.3)
qq.gam(pif,pch=19,cex=.3)
qq.gam(pif,rep=100,level=.9)
qq.gam(pif,rep=100,level=1,type="pearson",pch=19,cex=.2)


}
\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ...






