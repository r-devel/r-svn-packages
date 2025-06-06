# SCCS @(#)coxph.getdata.s	4.3 09/08/94
#
# Reconstruct the Cox model data.  This is done in so many routines
#  that I extracted it out.
#
# The "stratax" name is to avoid conflicts with the strata() function, but
#   still allow users to type "strata" as an arg.
#
coxph.getdata <- function(fit, y=TRUE, x=TRUE, stratax=TRUE, offset=FALSE) {
    ty <- fit$y
    tx <- fit$x
    strat <- fit$strata
    Terms <- fit$terms
    if (is.null(attr(Terms, 'offset'))) offset <- FALSE
    if (offset) x<- TRUE
    if (!inherits(Terms, 'terms'))
	    stop("invalid terms component of fit")
    strats <- attr(Terms, "specials")$strata
    if (length(strats)==0) stratax <- FALSE

    if ( (y && is.null(ty)) || (x && is.null(tx)) ||
	     (stratax && is.null(strat)) || offset) {
	# get the model frame
	m <- fit$model
	if (is.null(m)) m <- model.frame(fit)

	# Pull things out
	if (y && is.null(ty)) ty <- model.extract(m, 'response')

	if (offset) toff <- model.extract(m, 'offset')

	# strata was saved in the fit if and only if x was
	if (x && is.null(tx)) {
	    dropx <- untangle.specials(Terms, 'cluster')$terms
	    if (stratax) {
		temp <- untangle.specials(Terms, 'strata', 1)
		tx <- model.matrix(Terms[-c(dropx,temp$terms)], m)[,-1,drop=FALSE]
		strat <- strata(m[temp$vars], shortlabel=TRUE)
		}
	    else {
		if (length(dropx)) tx <- model.matrix(Terms[-dropx], m)[,-1,drop=FALSE]
		else               tx <- model.matrix(Terms, m)[,-1,drop=FALSE]
		}
	    }
	}
    else if (offset)
       toff <- fit$linear.predictors -(c(tx %*% fit$coef) - sum(fit$means*fit$coef))

    temp <- NULL
    if (y) temp <- c(temp, "y=ty")
    if (x) temp <- c(temp, "x=tx")
    if (stratax)  temp <- c(temp, "strata=strat")
    if (offset)  temp <- c(temp, "offset=toff")

    eval(parse(text=paste("list(", paste(temp, collapse=','), ")")))
    }
