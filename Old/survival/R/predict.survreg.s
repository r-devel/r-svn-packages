# SCCS @(#)predict.survreg.s	4.11 02/06/99
predict.survreg <-
    function(object, newdata, type=c('response', "link", 'lp', 'linear',
				     'terms', 'quantile','uquantile'),
				se.fit=FALSE,  terms=NULL,
	                        p=c(.1, .9),...)
    {
#
# What do I need to do predictions ?
#   
#  linear predictor: exists
#           +se    : X matrix
#          newdata : new X matrix
#
#  response -- same as lp, +transform, from distribution
#  
#  p --  density function from distribution
#          scale(s) -- if multiple I need the strata
#          +se : variance matrix
#	   newdata: new X
#
    ripley<-FALSE ##obsolete undocumented option
    type <-match.arg(type)
    if (type=='link') type<- 'lp'  #true until their are link functions
    if (type=='linear') type<- 'lp'
    n <- length(object$linear.predictors)
    Terms <- object$terms
    if(!inherits(Terms, "terms"))
	    stop("invalid terms component of  object")

    strata <- attr(Terms, 'specials')$strata
    Terms <- delete.response(Terms)
    coef <- object$coefficients
    intercept <- attr(Terms, "intercept")
    nvar <- length(object$coefficients)
    vv <- object$var[1:nvar, 1:nvar]
    fixedscale <- (nvar == ncol(object$var)) || ripley

    if (missing(newdata) && (type=='terms' || se.fit)) need.x <- TRUE
    else  need.x <- FALSE

    if (length(strata) && (type=='quantile' || type=='uquantile') &&
	      !fixedscale) {
	#
	# We need to reconstruct the "strata" variable
	#
	if (is.null(object$model)) m <- model.frame(object)
	else m <- object$model
	temp <- untangle.specials(Terms, 'strata', 1)
	dropx <- temp$terms
	if (length(temp$vars)==1) strata.keep <- m[[temp$vars]]
	else strata.keep <- strata(m[,temp$vars], shortlabel=TRUE)
	strata <- as.numeric(strata.keep)
	nstrata <- max(strata)
	    
	if (missing(newdata) && need.x){
	    x <- object$x
	    if (is.null(x)) x <- model.matrix(Terms[-dropx], m)
	    }

	else if (!missing(newdata)) {
	    newframe <- model.frame(Terms, newdata, na.action=function(x)x)
	    if (length(temp$vars)==1) newstrat <- newframe[[temp$vars]]
	    else newstrat <- strata(newframe[,temp$vars], shortlabel=TRUE)
	    strata <- match(newstrat, levels(strata.keep))
	    x <- model.matrix(Terms[-dropx], newframe)
	    offset <- model.extract(newframe, 'offset')
	    }
	}

    else {  # per subject strata not needed
	temp <- untangle.specials(Terms, 'strata', 1)
	if (length(temp$terms)) Terms <- Terms[-temp$terms]
	strata <- rep(1,n); nstrata<- 1
	if (missing(newdata) && need.x) {
	    x <- object$x
	    if (is.null(x)) {
		if (is.null(object$model)) 
			x <- model.matrix(Terms, model.frame(object))
		else    x <- model.matrix(Terms, object$model)
		}
	    }

	else if (!missing(newdata)) {
	    x <- model.matrix(Terms, newdata)
	    offset <- 0
	    strata <- rep(1, nrow(x))
	    }
	}
    scale <- object$scale[strata]
    #center x if terms are to be computed
    if(type=='p' || (type == "terms" && intercept)) 
	    x <- sweep(x, 2, object$means)

    #
    # Grab the distribution
    #
    if (is.character(object$dist)) dd <- survreg.distributions[[object$dist]]
    else dd <- object$dist
    if (is.null(dd$itrans)) {
	itrans <- function(x) x
        dtrans <- function (x) 1 ## bug fix from TMT, 2002-17-6
	}
    else {
	itrans <- dd$itrans
	dtrans <- dd$dtrans
	}
    if (!is.null(dd$dist)) dd <- survreg.distributions[[dd$dist]]

    #
    # Now, lay out the code one case at a time.
    #  There is some repetition this way, but otherwise the code just gets
    #    too complicated.
    #
    if (type=='lp' || type=='response') {
	if (missing(newdata)) {
 	    pred <- object$linear.predictors
#	    names(pred) <- names(object$residuals)
	    }
	else  pred <- x %*% coef  + offset
	if (se.fit) se <- sqrt(diag(x %*% vv %*% t(x)))

	if (type=='response') {
	    pred <- itrans(pred)
	    if (se.fit) se <- se/ dtrans(pred)
	    }
	}
    else if (type=='quantile' || type=='uquantile') {
	if (missing(newdata)) pred <- object$linear.predictors
	else  pred <- x %*% coef 
	# "pred" is the mean of the distribution,
	#   now add quantiles and then invert
	qq <- dd$quantile(p, object$parm)
	if (length(qq)==1 || length(pred)==1) {
	    pred <- pred + qq*scale
	    if (se.fit && fixedscale) {
		var <- ((x %*% vv) * x) %*% rep(1., ncol(x))
		se <- rep(sqrt(drop(var)), length(qq))
		}
	    else if (se.fit) {
		x.strata <- outer(strata, 1:nstrata, 
				  function(x,y) 1*(x==y))
		se <- matrix(0, ncol=length(qq), nrow=nrow(x))
		for (i in 1:(length(qq))) {
		    temp <- cbind(x, (qq[i]*scale)* x.strata)
		    var <- ((temp %*% object$var) *temp) %*% rep(1, ncol(temp))
		    se[,i] <- sqrt(drop(var))
		    }
		se <- drop(se)
		}
	    }
	else {
	    pred <- c(pred) + outer(scale, qq)
	    if (se.fit && fixedscale) {
		var <- ((x %*% vv) * x) %*% rep(1., ncol(x))
		if (length(qq) >1) {
		    se <- rep(sqrt(drop(var)), length(qq))
		    se <- matrix(se, ncol=length(qq))
		    }
		else se <- sqrt(drop(var))
		}
	    else if (se.fit) {
		x.strata <- outer(strata, 1:nstrata, 
				  function(x,y) 1*(x==y))
		se <- pred
		nc <- rep(1., ncol(object$var))
		for (i in 1:length(qq)) {
		    temp <- cbind(x, (qq[i]*scale)*x.strata)
		    var <- ((temp %*% object$var)* temp) %*% nc
		    se[,i] <- sqrt(drop(var))
		    }
		se <- drop(se)
		}
	    }
	pred <- drop(pred)
	if (type == 'quantile') {
	    pred <- itrans(pred)
	    if (se.fit) se <- se/dtrans(pred)
	    }
	}

    else {  #terms
	asgn <- attrassign(x,Terms)
        hasintercept<-attr(Terms,"intercept")>0
        if (hasintercept)
          asgn$"(Intercept)"<-NULL
        nterms<-length(asgn)
        pred<-matrix(ncol=nterms,nrow=NROW(x))
        dimnames(pred)<-list(rownames(x),names(asgn))
        if (se.fit){
          se<-matrix(ncol=nterms,nrow=NROW(x))
          dimnames(se)<-list(rownames(x),names(asgn))
          R<-object$var
          ip <- real(NROW(x))
        }
        for (i in 1:nterms){
          ii<-asgn[[i]]
          pred[,i]<-x[,ii,drop=FALSE]%*%(coef[ii])
          if (se.fit){
            for(j in (1:NROW(x))){
              xi<-x[j,ii,drop=FALSE]*(coef[ii])
              vci<-R[ii,ii]
              se[j,i]<-sqrt(sum(xi%*% vci %*%t( xi)))
            }
          }
        }
        if (!is.null(terms)){
          pred<-pred[,terms,drop=FALSE]
          if (se.fit)
            se<-se[,terms,drop=FALSE]
        }
      }
    if (type=='terms') attr(pred, 'constant') <- if (hasintercept) coef(object)["(Intercept)"] else 0
    #Expand out the missing values in the result
    # But only if operating on the original dataset
    if (missing(newdata) && !is.null(object$na.action)) {
	pred <- naresid(object$na.action, pred)
	if(se.fit) se <- naresid(object$na.action, se)

	
	}
    if (se.fit) list(fit=pred, se.fit=se)
    else pred
    }
