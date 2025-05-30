#SCCS  @(#)pyears.s	5.8 12/30/01
pyears <- function(formula=formula(data), data=sys.parent(),
	weights, subset, na.action,
	ratetable=survexp.us, scale=365.25,  expect=c('event', 'pyears'),
	model=FALSE, x=FALSE, y=FALSE) {

    expect <- match.arg(expect)
    call <- match.call()
    m <- match.call(expand.dots=FALSE)
    m$ratetable <- m$model <- m$x <- m$y <- m$scale<- m$expect <- NULL

    Terms <- if(missing(data)) terms(formula, 'ratetable')
	     else              terms(formula, 'ratetable',data=data)

    rate <- attr(Terms, "specials")$ratetable
    if (length(rate) >1 )
	stop ("Can have only 1 ratetable() call in a formula")
    else if (length(rate) == 0 && !missing(ratetable)) {
	# add a 'ratetable' call to the internal formula
        # The dummy function stops an annoying warning message "Looking for
        #  'formula' of mode function, ignored one of mode ..."
	xx <- function(x) formula(x)
    
	if(is.ratetable(ratetable))   varlist <- attr(ratetable, "dimid")
	else stop("Invalid rate table")

	ftemp <- deparse(substitute(formula))
	formula <- xx( paste( ftemp, "+ ratetable(",
			  paste( varlist, "=", varlist, collapse = ","), ")"))
	Terms <- if (missing(data)) terms(formula, "ratetable")
	         else               terms(formula, "ratetable", data = data)
	rate <- attr(Terms, "specials")$ratetable
	}

    if (any(attr(Terms, 'order') >1))
	    stop("Pyears cannot have interaction terms")
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())

    Y <- model.extract(m, 'response')
    if (is.null(Y)) stop ("Follow-up time must appear in the formula")
    if (!is.Surv(Y)){
	if (any(Y <0)) stop ("Negative follow up time")
	Y <- as.matrix(Y)
	if (ncol(Y) >2) stop("Y has too many columns")
	if (ncol(Y)==2 && any(Y[,2] <= Y[,1]))
	    stop("stop time must be > start time")
	}
    n <- nrow(Y)
    if (is.null(n) || n==0) stop("Data set has 0 observations")

    weights <- model.extract(m, 'weights')
    if (is.null(weights)) weights <- rep(1.0, n)

    if (length(rate)==1) {
	ovars <- (dimnames(attr(Terms, 'factors'))[[1]])[-c(1, rate)]
	rtemp <- match.ratetable(m[,rate], ratetable)
	R <- rtemp$R
	if (!is.null(rtemp$call)) {#need to drop some dimensions from ratetable
	    ratetable <- eval(parse(text=rtemp$call))
	    }
	}
    else ovars <- (dimnames(attr(Terms, 'factors'))[[1]])[-1]

    # Now process the other (non-ratetable) variables
    if (length(ovars)==0)  {
	# no categories!
	X <- rep(1,n)
	ofac <- odim <- odims <- ocut <- 1
	}
    else {
	odim <- length(ovars)
	ocut <- NULL
	odims <- ofac <- double(odim)
	X <- matrix(0, n, odim)
	outdname <- vector("list", odim)
	for (i in 1:odim) {
	    temp <- m[[ovars[i]]]
	    ctemp <- oldClass(temp)
	    if (!is.null(ctemp) && ctemp=='tcut') {
		X[,i] <- temp
		temp2 <- attr(temp, 'cutpoints')
		odims[i] <- length(temp2) -1
		ocut <- c(ocut, temp2)
		ofac[i] <- 0
		outdname[[i]] <- attr(temp, 'labels')
		}
	    else {
		temp2 <- factor(temp)
		X[,i] <- temp2
		temp3 <- levels(temp2)
		odims[i] <- length(temp3)
		ofac[i] <- 1
		outdname[[i]] <- temp3
		}
	    }
	}

    # Now do the computations
    ocut <-c(ocut,0)   #just in case it were of length 0
    osize <- prod(odims)
    if (length(rate)) {  #include expected
	atts <- attributes(ratetable)
	cuts <- atts$cutpoints
	rfac <- atts$factor
	us.special <- (rfac >1)
	if (any(us.special)) {  #special handling for US pop tables
	    if (sum(us.special) >1)
		stop("Two columns marked for special handling as a US rate table")
	    #slide entry date so that it appears that they were born on Jan 1
	    cols <- match(c("age", "year"), atts$dimid)
	    if (any(is.na(cols))) stop("Ratetable does not have expected shape")
            temp <- date.mdy(R[, cols[2]] - R[, cols[1]])
            R[, cols[2]] <- R[, cols[2]] - mdy.date(temp$month,
                temp$day, 1960)
	    # Doctor up "cutpoints"
	    temp <- (1:length(rfac))[us.special]
	    nyear <- length(cuts[[temp]])
	    nint <- rfac[temp]       #intervals to interpolate over
	    cuts[[temp]] <- round(approx(nint*(1:nyear), cuts[[temp]],
					    nint:(nint*nyear))$y - .0001)
	    }

	temp <- .C("pyears1",
			as.integer(n),
			as.integer(ncol(Y)),
			as.integer(is.Surv(Y)),
			as.double(Y),
		        as.double(weights),
			as.integer(length(atts$dim)),
			as.integer(rfac),
			as.integer(atts$dim),
			as.double(unlist(cuts)),
			as.double(ratetable),
			as.double(R),
			as.integer(odim),
			as.integer(ofac),
			as.integer(odims),
			as.double(ocut),
			as.integer(expect=='event'),
			as.double(X),
			pyears=double(osize),
			pn    =double(osize),
			pcount=double(if(is.Surv(Y)) osize else 1),
			pexpect=double(osize),
			offtable=double(1), PACKAGE="survival")[18:22]
	}
    else {
	temp <- .C('pyears2',
			as.integer(n),
			as.integer(ncol(Y)),
			as.integer(is.Surv(Y)),
			as.double(Y),
		        as.double(weights),
			as.integer(odim),
			as.integer(ofac),
			as.integer(odims),
			as.double(ocut),
			as.double(X),
			pyears=double(osize),
			pn    =double(osize),
			pcount=double(if(is.Surv(Y)) osize else 1),
			offtable=double(1), PACKAGE="survival") [11:14]
	}

    if (prod(odims) ==1) {  #don't make it an array
	out <- list(call=call, pyears=temp$pyears/scale, n=temp$pn,
		    offtable=temp$offtable/scale)
	if (length(rate)) {
	    out$expected <- temp$pexpect
	    if (expect=='pyears') out$expected <- out$expected/scale
	    if (!is.null(rtemp$summ)) out$summary <- rtemp$summ
	    }
	if (is.Surv(Y)) out$event <- temp$pcount
	}
    else {
	out <- list(call = call,
		pyears= array(temp$pyears/scale, dim=odims, dimnames=outdname),
		n     = array(temp$pn,     dim=odims, dimnames=outdname),
		offtable = temp$offtable/scale)
	if (length(rate)) {
	    out$expected <- array(temp$pexpect, dim=odims, dimnames=outdname)
	    if (expect=='pyears') out$expected <- out$expected/scale
	    if (!is.null(rtemp$summ)) out$summary <- rtemp$summ
	    }
	if (is.Surv(Y))
		out$event <- array(temp$pcount, dim=odims, dimnames=outdname)
	}
    na.action <- attr(m, "na.action")
    if (length(na.action))  out$na.action <- na.action
    if (model) out$model <- m
    else {
	if (x) out$x <- cbind(X, R)
	if (y) out$y <- Y
	}
    oldClass(out) <- 'pyears'
    out
    }

