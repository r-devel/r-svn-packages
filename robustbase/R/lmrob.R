### FIXME:
### ----- MM wants to change
### 1) compute.rd = FALSE  {no robust Mahalanobis distances;
###			   are expensive, and only needed for some plotting
###    will want 'x = !compute.rd' (need X matrix for Maha.Dist)
### 2) 'seed': By default always same seed --> same result
###	       even though algorithm is random.
###	INSTEAD: I want to use R's .Random.seed!
### 3) allow the 'control' entries to enter via "..." as well

### 4) lmrob() should really behave like lm() {in R; not S !}
###		--> 'subset' etc

### 5) There are still quite a few things hard-coded in C.
###    In particular,  'Nres = 500' cannot be changed.
lmrob <-
    function(formula, data = list(), weights, na.action,
	     model = TRUE, x = FALSE, y = FALSE,
	     singular.ok = TRUE, contrasts = NULL, offset = NULL,
	     control = lmrob.control())
{
    ret.x <- x
    ret.y <- y
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    mf$singular.ok <- mf$model <- NULL
    mf$x <- mf$y <- mf$contrasts <- mf$control <- NULL
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    na.act <- attr(mf, "na.action")
    xvars <- as.character(attr(mt, "variables"))[-1]
    x <- as.matrix(mf[-1])
    if ((yvar <- attr(mt, "response")) > 0)
	xvars <- xvars[-yvar]
    xlev <- if (length(xvars) > 0) {
	xlev <- lapply(mf[xvars], levels)
	xlev[!sapply(xlev, is.null)]
    }
    if (!singular.ok)
	warning("only `singular.ok = TRUE' is currently implemented.")
    y <- model.response(mf, "numeric")
    w <- model.weights(mf)
    offset <- model.offset(mf)
    if (!is.null(offset) && length(offset) != NROW(y))
	stop(paste("Number of offsets is", length(offset), ", should equal",
		   NROW(y), "(number of observations)"))
    if (is.empty.model(mt)) {
	xx <- NULL
	z <- list(coefficients = numeric(0), residuals = y, fitted.values = 0 *
		  y + offset, weights = w, rank = 0, df.residual = length(y))
	class(z) <- c("lmrob.null", "lmrob")
    }
    else {
	xx <- model.matrix(mt, mf, contrasts)
	z <- if (is.null(w))
	    lmrob.fit.MM(xx, y, control = control)
	else stop("Weights are not implemented for this estimator")
	class(z) <- c("lmrob")
    }
    if (!is.null(na.act))
	z$na.action <- na.act
    z$offset <- offset
    z$contrasts <- attr(xx, "contrasts")
    z$xlevels <- xlev
    z$call <- cl
    z$terms <- mt
    if( control$compute.rd ) {
	rob <- MASS::cov.rob(x, method = 'mcd')
	z$MD <- sqrt( mahalanobis(x, rob$center, rob$cov) )
    }
    if (model)
	z$model <- mf
    if (ret.x) ## Careful!  returns 'xx' whereas Maha.Dist use 'x'
	z$x <- xx
    if (ret.y)
	z$y <- y
    z$control <- control
    z
}
