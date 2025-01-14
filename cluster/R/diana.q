### $Id$

diana <- function(x, diss = inherits(x, "dist"),
		  metric = "euclidean", stand = FALSE,
		  stop.at.k = FALSE,
                  keep.diss = n < 100, keep.data = !diss, trace.lev = 0)
{
    if((diss <- as.logical(diss))) {
	## check type of input vector
	if(anyNA(x)) stop("NA values in the dissimilarity matrix not allowed.")
	if(data.class(x) != "dissimilarity") { # try to convert to
	    if(!is.null(dim(x))) {
		x <- as.dist(x) # or give an error
	    } else {
		## possibly convert input *vector*
		if(!is.numeric(x) || is.na(n <- sizeDiss(x)))
		    stop("'x' is not and cannot be converted to class \"dissimilarity\"")
		attr(x, "Size") <- n
	    }
	    class(x) <- dissiCl
	    if(is.null(attr(x,"Metric"))) attr(x, "Metric") <- "unspecified"
	}
	n <- as.integer(attr(x, "Size"))
	dv <- x[lower.to.upper.tri.inds(n)] # <==> n >= 2 or error;  *slow* [c * n^2 ; but large c]  in large cases
	## prepare arguments for the Fortran call
	jp <- 1L
	mdata <- FALSE
	ndyst <- 0L
	x2 <- double(1)
    }
    else {
	## check input matrix and standardize, if necessary
	x <- data.matrix(x)
	if(!is.numeric(x)) stop("x is not a numeric dataframe or matrix.")
	x2 <- if(stand) scale(x, scale = apply(x, 2, meanabsdev)) else x
	ndyst <- if(metric == "manhattan") 2L else 1L
	n <- nrow(x2)
	jp <- ncol(x2)
        if(!jp) stop("x has zero columns") # easier to read than later error
        stopifnot(is.integer(n), is.integer(jp)) # always currently
	if((mdata <- any(inax <- is.na(x2)))) { # TRUE if x[] has any NAs
	    jtmd <- integer(jp)
	    jtmd[apply(inax, 2L, any)] <- -1L
	    ## VALue for MISsing DATa
            ## __ FIXME __ now have C and R only, could use true NA (double | int.) or 'Inf'
            ##    =====   the following fails e.g. when max(x2) == double.xmax
	    valmisdat <- 1.1* max(abs(range(x2, na.rm=TRUE)))
	    x2[inax] <- valmisdat
	}
	dv <- double((n * (n - 1))/2) # FIXME: an allocation waste for large n !
    }
    stopifnot(is.logical(stop.at.k) ||
	      (is.numeric(stop.at.k) && 1 <= stop.at.k && stop.at.k <= n))
    stopifnot(length(trace.lev <- as.integer(trace.lev)) == 1)
    C.keep.diss <- keep.diss && !diss
    res <- .C(twins,
		    n,
		    jp,
		    as.double(x2),
		    dv, # input, w/ length >= 1 \\  output:
		    dis = double(if(C.keep.diss) length(dv) else 1L),
		    jdyss = if(C.keep.diss) diss + 10L else as.integer(diss),
		    if(mdata && jp) rep(valmisdat, jp) else double(1L),
		    if(mdata) jtmd else integer(jp),
		    ndyst,
		    2L,# jalg = 2 <==> DIANA
                    as.integer(stop.at.k),# 'method'; default = 0L  :  do *not* stop early
		    integer(n),
		    ner = integer(n),
		    ban = double(n),
		    dc = double(1L), # coef
		    double(1L), # { unused for diana() }
		    merge = matrix(0L, n - 1L, 2L), # integer or error if(n == 0) !
                    trace = trace.lev)[c("dis", "jdyss", "ner", "ban", "dc", "merge")]
    if(!diss) {
	## give warning if some dissimilarities are missing.
	if(res$jdyss == -1)
	    stop("No clustering performed, NA values in the dissimilarity matrix.")
        if(keep.diss) {
            ## adapt Fortran output to S:
            ## convert lower matrix, read by rows, to upper matrix, read by rows.
            disv <- res$dis
            disv[disv == -1] <- NA
            disv <- disv[upper.to.lower.tri.inds(n)] # <==> n >= 2 or error
            class(disv) <- dissiCl
            attr(disv, "Size") <- nrow(x)
            attr(disv, "Metric") <- metric
            attr(disv, "Labels") <- dimnames(x)[[1]]
        }
	## add labels to Fortran output
	if(length(dimnames(x)[[1]]) != 0)
	    order.lab <- dimnames(x)[[1]][res$ner]
    }
    else {
        if(keep.diss) disv <- x
	## add labels to Fortran output
	if(length(attr(x, "Labels")) != 0)
	    order.lab <- attr(x, "Labels")[res$ner]
    }
    clustering <- list(order = res$ner, height = res$ban[-1], dc = res$dc,
		       merge = res$merge, diss = if(keep.diss)disv,
                       call = match.call())
    if(exists("order.lab"))
	clustering$order.lab <- order.lab
    if(keep.data && !diss) {
	if(mdata) x2[x2 == valmisdat] <- NA
	clustering$data <- x2 # NB: in non-NA cases, x2 may still be "integer"
    }
    class(clustering) <- c("diana", "twins")
    clustering
}

print.diana <- function(x, ...)
{
    cat("Merge:\n")
    print(x$merge, ...)
    cat("Order of objects:\n")
    print(if (length(x$order.lab) != 0) x$order.lab else x$order,
	  quote = FALSE, ...)
    cat("Height:\n")
    print(x$height, ...)
    cat("Divisive coefficient:\n")
    print(x$dc, ...)
    cat("\nAvailable components:\n")
    print(names(x), ...)
    invisible(x)
}

summary.diana <- function(object, ...)
{
    class(object) <- "summary.diana"
    object
}

print.summary.diana <- function(x, ...)
{
    cat("Merge:\n");			print(x$merge, ...)
    cat("Order of objects:\n")
    print(if(length(x$order.lab)) x$order.lab else x$order, quote = FALSE, ...)
    cat("Height:\n");			print(x$height, ...)
    cat("Divisive coefficient:\n");	print(x$dc, ...)
    if(!is.null(x$diss)) { ## Dissimilarities:
	cat("\n");			print(summary(x$diss, ...))
    }
    cat("\nAvailable components:\n");	print(names(x), ...)
    invisible(x)
}
