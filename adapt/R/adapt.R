.First.lib <- function(lib, pkg) {
    library.dynam("adapt", pkg, lib)
}


adapt <- function (ndim, lower, upper, minpts = 100, maxpts = NULL,
                   functn, eps = 0.01, xtra.lenwrk = 0:0, ...)
{
    keep.trying <- is.null(maxpts)
    minpts <- as.integer(minpts)

    if (ndim == 1) { ## fudge for 1-d functions
	warning("Using integrate() from base package for 1-d integration")
        if (keep.trying) maxpts <- minpts
	return(integrate(functn,lower,upper,subdivisions=maxpts,rel.tol=eps,...))
    }
    ## else ndim >= 2 :

    ## Check to make sure that upper and lower are reasonable lengths
    ## Both the upper and lower limits should be at least of length ndim
    if (length(lower) < ndim || length(upper) < ndim)#MM: dropped 'at least':
	stop("The lower and upper vectors need to have ndim elements\n",
	     "Your parameters are:  ndim ", ndim, ", length(lower) ",
	     length(lower), ", length(upper) ", length(upper))
    ff <-
	if(length(list(...)) && length(formals(functn)) > 1)
	    function(x) functn(x, ...)
	else functn # .Alias

    rulcls <- 2^ndim + 2*ndim^2 + 6*ndim + 1 #-> ../src/adapt.f

    ## maxpts should be large enough.  Prefer 10*rulclc, but use 2*rulclc.
    if (keep.trying)
	maxpts <- max(minpts, 500, 2 * rulcls)
    else {
	if (minpts >= maxpts) {
	    warning("maxpts must be > minpts.\n",
		    "Maxpts has be increased to	 minpts + 1")
	    maxpts <- minpts + 1
	}
	##
	if (maxpts < 2 * rulcls) {
	    warning("You have maxpts (= ", maxpts, ") too small.\n It needs to",
		    "be at least 2 times 2^ndim + 2*ndim^2 + 6*ndim+1\n",
		    " and has been reset to ", 2 * rulcls, "\n")
            maxpts <- 2 * rulcls
        }
    }

    rgnstr <- 2*ndim + 3 # as in Fortran
    it <- 0
    repeat {
        it <- it+1
        ## lenwrk= length work array{allocated in C} mandated in adapt source:
	lenwrk <- as.integer(rgnstr * (1 + maxpts/rulcls)/2)
        if(verbose) cat(it,": lenwrk.orig=", lenwrk,", ")
        lenwrk <- as.integer(lenwrk + xtra.lenwrk)

	x <- .C("cadapt",
		as.integer(ndim),
		as.double(lower),
		as.double(upper),
		minpts = minpts,
		maxpts = as.integer(maxpts),
		## now pass ff and current environment
		ff, rho = environment(),
		eps = as.double(eps),
		relerr = double(1),
		lenwrk = lenwrk,
		value = double(1),	# will contain the value of the integral
		ifail = integer(1),
                PACKAGE = "adapt")[
                c("value","relerr","minpts", "lenwrk", "ifail", "eps")]

        if(verbose)
            cat("  needed (lenwrk equivalent of)", x$eps+ 2*(ndim + 3)+ 2,"\n")
	if (x$ifail == 1 && keep.trying)
	    maxpts <- maxpts*2
	else
	    break
    }
    if(x$ifail)
	warning(x$warn <-
		c("Ifail=1, maxpts was too small. Check the returned relerr!",
		  paste("Ifail=2, lenwrk=",lenwrk," was too small, needed",
                        ## sbrgns + rgnstr + 2 (<- ../src/adapt.f )  =
                        x$eps + 2*(ndim + 3) + 2,
                        ".\n ***> fix adapt(); check the returned relerr!"),
		  "Ifail=3: ndim > 20 -- rewrite the fortran code ;-) !",
		  "Ifail=4, minpts > maxpts; should not happen!",
		  "Ifail=5, internal non-convergence; should not happen!",
		  )[x$ifail])

    x$eps <- NULL
    x$iter <- it
    class(x) <- "integration"
    x
}

print.integration <- function(x, ...) {
    print(noquote(sapply(x, format, ...)),...)
    invisible(x)
}

