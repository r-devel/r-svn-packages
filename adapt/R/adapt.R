.First.lib <- function(lib, pkg) {
    library.dynam("adapt", pkg, lib)
}

adapt <- function (ndim, lower, upper, minpts, maxpts, functn, eps, ...)
{
    keep.trying <- is.null(maxpts)
    if (is.null(maxpts)) maxpts <- max(minpts,500)

    if (ndim == 1) { ## fudge for 1-d functions
	warning("Using integrate() from base package for 1-d integration")
	return(integrate(functn,lower,upper,subdivisions=maxpts,rel.tol=eps,...))
    }
    else ff <-
	if(length(list(...)) && length(formals(functn)) > 1)
	    function(x) functn(x, ...)
	else functn # .Alias

    ## Check to make sure that upper and lower are reasonable lengths
    ## Both the upper and lower limits should be at least of length ndim
    ##
    if (length(lower) < ndim || length(upper) < ndim)#MM: dropped 'at least':
	stop(paste("The lower and upper vectors need to have ndim elements\n",
		   "Your parameters are:  ndim", ndim, ", length(lower)",
		   length(lower), ", length(upper)", length(upper), "\n"))

    if (minpts >= maxpts) {
	warning(paste("maxpts must be > minpts.\n",
		      "Maxpts has be increased to  minpts + 1"))
	maxpts <- minpts + 1
    }
    ## rulcls and lenwrk are mandated in the adapt source
    ##
    repeat{
	rulcls <- 2^ndim + 2*ndim^2 + 6*ndim + 1
	lenwrk <- (2*ndim + 3) * (1 + maxpts/rulcls)/2
	## maxpts should be large enough.  Prefer 10*rulclc, but use 2*rulclc.
	##
	if (maxpts < 2 * rulcls) {
	    warning(paste("You have maxpts (= ", maxpts, ") too small\n",
			  "It needs to be at least 2 times 2^ndim + 2*ndim^2 + 6*ndim+1\n",
			  "It has been reset to ", 2 * rulcls, "\n", sep=""))
	    maxpts <- 2 * rulcls
	}

	x <- .C("cadapt",
		as.integer(ndim),
		as.double(lower),
		as.double(upper),
		minpts = as.integer(minpts),
		maxpts = as.integer(maxpts),
		## now pass ff and current environment
		ff, rho = environment(),
		as.double(eps),
		relerr = double(1),
		lenwrk = as.integer(lenwrk),
		value = double(1),	# will contain the value of the integral
		ifail = integer(1))[
                c("value","relerr","minpts", "lenwrk", "ifail")]

	if (x$ifail == 1 && keep.trying)
	    maxpts <- maxpts*2
	else
	    break
    }
    if(x$ifail)
	warning(x$warn <-
		c("Ifail=1, maxpts was too small. Check the returned relerr!",
		  paste("Ifail=2, lenwrk was too small. -- fix adapt() !\n",
			"Check the returned relerr!"),
		  "Ifail=3: ndim > 20 -- rewrite the fortran code ;-) !",
		  "Ifail=4, minpts > maxpts; should not happen!",
		  "Ifail=5, internal non-convergence; should not happen!",
		  )[x$ifail])

    class(x) <- "integration"
    x
}

print.integration <- function(x, ...) {
    print(noquote(sapply(x, format, ...)),...)
    invisible(x)
}

