\name{gam.models}
\alias{gam.models}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Specifying generalized additive models}
\description{ This page is intended to provide some more information on
 how to specify GAMs. A GAM is a GLM in which the linear predictor depends, 
in part, on a sum of smooth functions of predictors and (possibly) linear 
functionals of smooth functions of (possibly dummy) predictors.

Specifically let \eqn{y_i}{y_i} denote an independent random variable 
with mean \eqn{\mu_i}{mu_i} and an exponential family distribution, or failing 
that a known mean variance relationship suitable for use of quasi-likelihood methods. 
Then the the linear predictor of a GAM has a structure something like
 
\deqn{g(\mu_i) = {\bf X}_i{\beta} + f_1(x_{1i},x_{2i}) + f_2(x_{3i}) + L_i f_3(x_4) + 
\ldots}{g(mu_i)=X_ib + f_1(x_1i,x_2i) + f_2(x_3i) + L_i f_3(x_4) +}

where \eqn{g}{g} is a known smooth monotonic `link' function, \eqn{{\bf X}_i\beta}{X_ib} 
is the parametric part of the linear predictorm the \eqn{x_j}{x_j} are predictor variables,
the \eqn{f_j}{f_j} are smooth functions and \eqn{L_i}{L_i} is some linear functional of 
\eqn{f_3}{f_3}. There may of course be multiple linear functional terms, or none.

The key idea here is that the
dependence of the response on the predictors can be represented as a
parametric sub-model plus the sum of some (functionals of) smooth functions of one or
more of the predictor variables. Thus the model is quite flexible
relative to strictly parametric linear or generalized linear models,
but still has much more structure than the completely general model
that says that the response is just some smooth function of all the
covariates.

Note one important point. In order for the model to be identifiable
the smooth functions have to be constrained to have zero mean (usually
taken over the set of covariate values). Such constraints are always
applied by \code{gam}.

The following sections discuss specifying model structures for \code{gam}. 
Specification of the distribution and link function is done using the \code{\link{family}} 
 argument to \code{\link{gam}} and works in the same way as for \code{\link{glm}}. 
This page therefore concentrates on the model formula for \code{gam}.

}

\section{Models with simple smooth terms}{
Consider the example model.
\deqn{\mu_i = \beta_0 + \beta_1 x_{1i} + \beta_2 x_{2i} + f_1(x_{3i}) + 
f_2(x_{4i},x_{5i})}{ mu_i = b_0 + b_1 x_1i + b_2 x_2i + f1(x_3i) + f2(x_4i,x_5i)}
where the (independent) response variables have distribution 
\eqn{y_i \sim N(\mu_i,\sigma^2)}{y_i ~ N(mu_i,s^2)}.

The \code{gam} formula for this would be \cr
  \code{y~x1+x2+s(x3)+s(x4,x5)}. \cr
  This would use the default basis for the smooths (a thin plate
  regression spline basis for each), with automatic selection of the
  effective degrees of freedom for both smooths. The dimension of the
  smoothing basis is given a default value as well (the dimension of the
  basis sets an upper limit on the maximum possible degrees of
  freedom for the basis - the limit is typically one less than basis
  dimension). Full details of how to control smooths are given in
  \code{\link{s}} and \code{\link{te}}, and further discussion of basis
  dimension choice can be found in \code{\link{choose.k}}. 
  For the moment suppose that we would like to change
  the basis of the first smooth to a cubic regression spline basis with
  a dimension of 20, while fixing the second term at 25 degrees of
  freedom. The appropriate formula would be:\cr
  \code{y~x1+x2+s(x3,bs="cr",k=20)+s(x4,x5,k=26,fx=TRUE)}.

The above assumes that \eqn{x_{4}}{x_4} and \eqn{x_5}{x_5} are naturally on 
similar scales (e.g. they might be co-ordinates), so that isotropic smoothing 
is appropriate. If this assumption is false then tensor product smoothing might be 
more appropriate (see \code{\link{te}}). Then\cr
\code{y~x1+x2+s(x3)+te(x4,x5)}\cr
would generate a tensor product smooth for \eqn{x_{4}}{x_4} and \eqn{x_5}{x_5}. 
By default this would have basis dimension 25 and use cubic regression spline marginals. 
Varying the default is easy. For example\cr
\code{y~x1+x2+s(x3)+te(x4,x5,bs=c("cr","ps"),k=c(6,7))}\cr
specifyies that the tensor product should use a rank 6 cubic regression spline marginal
and a rank 7 P-spline marginal to create a smooth with basis dimension 42.
}

\section{nested terms}{
 Another situation that occurs quite often is the one in which we
  would like to find out if the model:
  \deqn{E(y_i) = f(x_i,z_i)}{E(y)=f(x,z)}
  is really necessary or whether:
  \deqn{E(y_i)=f_1(x_i)+f_2(z_i)}{E(y)=f1(x)+f2(z)}
  wouldn't do just as well. One way to do this is to look at the results
  of fitting:\cr
  \code{y~s(x)+s(z)+s(x,z)}.\cr
  \code{gam} automatically generates side conditions to make this model
  identifiable. You can also estimate `overlapping' models like:\cr
  \code{y~s(x,z)+s(z,v)}.  
}

\section{`by' variables}{
 \code{by} variables are the means for constructing `varying-coefficient models' (geographic regression models) and 
for letting smooths `interact' with factors or parametric terms.

The \code{\link{s}} and \code{\link{te}} terms used to specify smooths in a model accept an argument \code{by}, 
which is a numeric or factor variable of the same dimension as the covariates of the smooth. 
If a \code{by} variable is numeric, then its \eqn{i^{th}}{ith} element multiples the  \eqn{i^{th}}{ith}
row of the model matrix corresponding to the smooth term concerned. If a \code{by} variable is a factor then it generates 
an indicator vector for each level of the factor. The model matrix for the smooth term is then replicated for each factor level,
and each copy has its rows multiplied by the corresponding rows of its indicator variable. In short a different smooth is generated
for each factor level (the \code{id} argument to \code{\link{s}} and \code{\link{te}} can be used to force all 
such smooths to have the same smoothing parameter).  

 
As an example, consider the model
\deqn{E(y_i) = \beta_0+ f(x_i)z_i}{E(y_i)=b_0+f(x_i)z_i}
where \eqn{f}{f} is a smooth function, and \eqn{z_i}{z_i} is a numeric variable.
The appropriate formula is:\cr
\code{y~z+s(x,by=z)}\cr
- the \code{by} argument ensures that the smooth function gets multiplied by
covariate \code{z}, but GAM smooths are centred (average value zero),
so the \code{z+} term is needed as well (\eqn{f}{f} is being
represented by a constant plus a centred smooth). If we'd wanted:
\deqn{E(y_i) = f(x_i)z_i}{E(y_i)=f(x_i)z_i}
then the appropriate formula would be:\cr
\code{y~z+s(x,by=z)-1}.\cr
The example code below also includes factor \code{by} variables.

\code{by} variables may be supplied as numeric matrices as part of specifying general linear functional terms.
}

\section{Linking smooths with `id'}{
It is sometimes desirable to insist that different smooth terms have the same degree of smoothness. 
This can be done by using the \code{id} argument to \code{\link{s}} or \code{\link{te}} terms. Smooths 
which share an \code{id} will have the same smoothing parameter. Really this only makes sense if the 
smooths use the same basis functions, and the default behaviour is to force this to happen: all smooths 
sharing an \code{id} have the same basis functions as the first smooth occurring with that \code{id}. Note 
that if you want exactly the same function for each smooth, then this is best achieved by making use of the 
summation convention covered under `linear functional terms'. 
}

\section{Linear functional terms}{

General linear functional terms have a long history in the spline literature including in the penalized 
GLM context (see e.g. Wahba 1990). Such terms encompass variable coefficient models/ geographic 
regression, functional GLMs (i.e. GLMs with functional predictors), GLASS models, etc, and allow 
smoothing with respect to aggregated covariate values, for example. 

}

%- maybe also `usage' for other objects documented here.

\references{
Wahba (1990) Spline Models of Observational Data SIAM.
}

\author{ Simon N. Wood \email{simon.wood@r-project.org}}

\examples{
set.seed(10)
## simulate date from y = f(x2)*x1 + error
dat <- gamSim(3,n=400)

b<-gam(y~x1-1+s(x2,by=x1),data=dat)
plot(b,pages=1)
summary(b)

## Factor `by' variable example (with a spurious covariate x0)
## simulate data...

dat <- gamSim(4)
attach(dat)

## fit model...
b <- gam(y ~ fac+s(x2,by=fac)+s(x0))
plot(b,pages=1)
summary(b)

## note that the preceding fit is the same as....
b1<-gam(y ~ fac+s(x2,by=as.numeric(fac==1))+s(x2,by=as.numeric(fac==2))+
            s(x2,by=as.numeric(fac==3))+s(x0))
plot(b1,pages=1)
summary(b1)

## repeat, but force all the s(x2) terms to have the same smoothing param
## (not a very good idea for these data!)
b2 <- gam(y ~ fac+s(x2,by=fac,id=1)+s(x0))
plot(b2,pages=1)
summary(b2)

detach(dat);rm(dat)
}
\keyword{models} \keyword{regression}%-- one or more ..


