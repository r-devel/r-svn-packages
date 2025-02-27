#### Toplevel ``virtual'' class "Matrix"


### Virtual coercions -- via smart "helpers" (-> ./Auxiliaries.R)

setAs("Matrix", "sparseMatrix", function(from) as(from, "CsparseMatrix"))
setAs("Matrix", "CsparseMatrix", function(from) as_Csparse(from))
setAs("Matrix", "denseMatrix",  function(from) as_dense(from))

## Maybe TODO:
## setAs("Matrix", "nMatrix", function(from) ....)

## Most of these work; this is a last resort:
setAs(from = "Matrix", to = "matrix", # do *not* call base::as.matrix() here:
      function(from) .bail.out.2("coerce", class(from), class(to)))
setAs(from = "matrix", to = "Matrix", function(from) Matrix(from))

## ## probably not needed eventually:
## setAs(from = "ddenseMatrix", to = "matrix",
##       function(from) {
## 	  if(length(d <- dim(from)) != 2) stop("dim(.) has not length 2")
## 	  array(from@x, dim = d, dimnames = dimnames(from))
##       })

## should propagate to all subclasses:
setMethod("as.matrix", signature(x = "Matrix"), function(x) as(x, "matrix"))
## for 'Matrix' objects, as.array() should be equivalent:
setMethod("as.array",  signature(x = "Matrix"), function(x) as(x, "matrix"))

## head and tail apply to all Matrix objects for which subscripting is allowed:
setMethod("head", signature(x = "Matrix"), utils::head.matrix)
setMethod("tail", signature(x = "Matrix"), utils::tail.matrix)

setMethod("drop", signature(x = "Matrix"),
	  function(x) if(all(dim(x) != 1)) x else drop(as(x, "matrix")))

## slow "fall back" method {subclasses should have faster ones}:
setMethod("as.vector", signature(x = "Matrix", mode = "missing"),
	  function(x) as.vector(as(x, "matrix")))

## mainly need these for "dMatrix" or "lMatrix" respectively, but why not general:
setMethod("as.numeric", signature(x = "Matrix"),
	  function(x, ...) as.numeric(as.vector(x)))
setMethod("as.logical", signature(x = "Matrix"),
	  function(x, ...) as.logical(as.vector(x)))

setMethod("cov2cor", signature(V = "Matrix"),
	  function(V) { ## was as(cov2cor(as(V, "matrix")), "dpoMatrix"))
	      r <- V
	      p <- (d <- dim(V))[1]
	      if(p != d[2]) stop("'V' is not a square matrix")
	      Is <- sqrt(1/diag(V)) # diag( 1/sigma_i )
	      if(any(!is.finite(Is)))
		  warning("diag(.) had 0 or NA entries; non-finite result is doubtful")
              Is <- Diagonal(x = Is)
              r <- Is %*% V %*% Is
	      r[cbind(1L:p,1L:p)] <- 1 # exact in diagonal
	      as(forceSymmetric(r), "dpoMatrix")
          })

## "base" has an isSymmetric() S3-generic since R 2.3.0
setMethod("isSymmetric", signature(object = "symmetricMatrix"),
	  function(object, ...) TRUE)
setMethod("isSymmetric", signature(object = "triangularMatrix"),
	  ## TRUE iff diagonal:
	  function(object, ...) isDiagonal(object))

setMethod("isTriangular", signature(object = "matrix"), isTriMat)

setMethod("isDiagonal", signature(object = "matrix"), .is.diagonal)

## The "catch all" methods -- far from optimal:
setMethod("symmpart", signature(x = "Matrix"),
	  function(x) as((x + t(x))/2, "symmetricMatrix"))
setMethod("skewpart", signature(x = "Matrix"),
	  function(x) (x - t(x))/2)

## FIXME: do this (similarly as for "ddense.." in C
setMethod("symmpart", signature(x = "matrix"), function(x) (x + t(x))/2)
setMethod("skewpart", signature(x = "matrix"), function(x) (x - t(x))/2)




setMethod("dim", signature(x = "Matrix"),
	  function(x) x@Dim, valueClass = "integer")

setMethod("length", "Matrix", function(x) prod(dim(x)))

setMethod("dimnames", signature(x = "Matrix"), function(x) x@Dimnames)


## not exported but used more than once for "dimnames<-" method :
## -- or do only once for all "Matrix" classes ??
dimnamesGets <- function (x, value) {
    d <- dim(x)
    if (!is.list(value) || length(value) != 2 ||
	!(is.null(v1 <- value[[1]]) || length(v1) == d[1]) ||
	!(is.null(v2 <- value[[2]]) || length(v2) == d[2]))
	stop(sprintf("invalid dimnames given for '%s' object", class(x)))
    x@Dimnames <- list(if(!is.null(v1)) as.character(v1),
		       if(!is.null(v2)) as.character(v2))
    x
}
setMethod("dimnames<-", signature(x = "Matrix", value = "list"),
	  dimnamesGets)

setMethod("unname", signature("Matrix", force="missing"),
	  function(obj) { obj@Dimnames <- list(NULL,NULL); obj})

setMethod("all", signature(x = "Matrix"),
	  function(x, ..., na.rm)
	  callGeneric(as(x, "lMatrix"), ..., na.rm=na.rm))

setMethod("any", signature(x = "Matrix"),
	  function(x, ..., na.rm)
	  callGeneric(as(x, "lMatrix"), ..., na.rm=na.rm))

## NOTE:  "&" and "|"  are now in group "Logic" c "Ops" --> ./Ops.R
##        "!" is in ./not.R


Matrix <- function (data = NA, nrow = 1, ncol = 1, byrow = FALSE,
                    dimnames = NULL, sparse = NULL, forceCheck = FALSE)
{
    sparseDefault <- function(m) prod(dim(m)) > 2*sum(isN0(as(m, "matrix")))

    i.M <- is(data, "Matrix")
    if(!i.M && inherits(data, "table")) # special treatment
	class(data) <- "matrix" # "matrix" first for S4 dispatch
    if(is.null(sparse1 <- sparse) && (i.M || is(data, "matrix")))
	sparse <- sparseDefault(data)
    sM <- FALSE
    doDN <- TRUE
    if (i.M) {
        if(!missing(nrow) || !missing(ncol)|| !missing(byrow))
            warning("'nrow', 'ncol', etc, are disregarded when 'data' is \"Matrix\" already")
	sM <- is(data,"sparseMatrix")
	if(!forceCheck && ((sparse && sM) || (!sparse && !sM)))
	    return(data)
	## else : convert  dense <-> sparse -> at end
    }
    else if (!is.matrix(data)) { ## cut & paste from "base::matrix" :
	if (missing(nrow))
	    nrow <- ceiling(length(data)/ncol)
	else if (missing(ncol))
	    ncol <- ceiling(length(data)/nrow)
	if(length(data) == 1 && is0(data) && !identical(sparse, FALSE)) {
	    ## Matrix(0, ...) : always sparse unless "sparse = FALSE":
	    if(is.null(sparse)) sparse1 <- sparse <- TRUE
	    i.M <- sM <- TRUE
            isSym <- nrow == ncol
	    ## will be sparse: do NOT construct full matrix!
	    data <- new(paste(if(is.numeric(data)) "d" else
			      if(is.logical(data)) "l" else
			      stop("invalid 'data'"),
			      if(isSym) "s" else "g", "CMatrix", sep=''),
			p = rep.int(0L, ncol+1L),
			Dim = as.integer(c(nrow,ncol)),
			Dimnames = if(is.null(dimnames)) list(NULL,NULL)
			else dimnames)
	} else { ## normal case - using .Internal() to avoid more copying
	    if(getRversion() >= "2.7.0")
		data <- .Internal(matrix(data, nrow, ncol, byrow, dimnames))
	    else {
		data <- .Internal(matrix(data, nrow, ncol, byrow))
		dimnames(data) <- dimnames
	    }
	    if(is.null(sparse))
		sparse <- sparseDefault(data)
	}
        doDN <- FALSE
    } else if(!missing(nrow) || !missing(ncol)|| !missing(byrow))
	warning("'nrow', 'ncol', etc, are disregarded for matrix 'data'")

    ## 'data' is now a "matrix" or "Matrix"
    if (doDN && !is.null(dimnames))
	dimnames(data) <- dimnames

    ## check for symmetric / triangular / diagonal :
    isSym <- isSymmetric(data)
    if((isTri <- !isSym))
	isTri <- isTriangular(data)
    isDiag <- isSym # cannot be diagonal if it isn't symmetric
    if(isDiag)
	isDiag <- !isTRUE(sparse1) && isDiagonal(data)

    ## try to coerce ``via'' virtual classes
    if(isDiag) { ## diagonal is preferred to sparse !
	data <- as(data, "diagonalMatrix")
	isSym <- FALSE
    } else if(sparse && !sM)
	data <- as(data, "sparseMatrix")
    else if(!sparse) {
	if(i.M) { ## data is 'Matrix'
	    if(!is(data, "denseMatrix"))
		data <- as(data, "denseMatrix")
	} else { ## data is "matrix" (and result "dense" -> go via "general"
	    ctype <- typeof(data)
	    if (ctype == "complex")
		stop("complex matrices not yet implemented in Matrix package")
	    if (ctype == "integer") ## integer Matrices not yet implemented
		storage.mode(data) <- "double"
	    data <- new(paste(.M.kind(data), "geMatrix", sep=''),
			Dim = dim(data),
			Dimnames = .M.DN(data),
			x = c(data))
	}
    }

    if(isTri && !is(data, "triangularMatrix")) {
	data <- if(attr(isTri,"kind") == "L") tril(data) else triu(data)
					#was as(data, "triangularMatrix")
    } else if(isSym && !is(data, "symmetricMatrix"))
	data <- forceSymmetric(data) #was as(data, "symmetricMatrix")

    data
}

## Methods for operations where one argument is numeric

## Using as.matrix() and rbind()
## in order to get dimnames from names {at least potentially}:

setMethod("%*%", signature(x = "Matrix", y = "numeric"),
	  function(x, y) callGeneric(x, as.matrix(y)))
setMethod("%*%", signature(x = "numeric", y = "Matrix"),
	  function(x, y) callGeneric(matrix(x, nrow = 1, byrow=TRUE), y))

setMethod("%*%", signature(x = "Matrix", y = "matrix"),
	  function(x, y) callGeneric(x, Matrix(y)))
setMethod("%*%", signature(x = "matrix", y = "Matrix"),
	  function(x, y) callGeneric(Matrix(x), y))


setMethod("crossprod", signature(x = "Matrix", y = "numeric"),
	  function(x, y = NULL) callGeneric(x, as.matrix(y)))
setMethod("crossprod", signature(x = "numeric", y = "Matrix"),
	  function(x, y = NULL)	 callGeneric(as.matrix(x), y))

setMethod("crossprod", signature(x = "Matrix", y = "matrix"),
	  function(x, y = NULL) callGeneric(x, Matrix(y)))
setMethod("crossprod", signature(x = "matrix", y = "Matrix"),
	  function(x, y = NULL) callGeneric(Matrix(x), y))

## The as.matrix() promotion seems illogical to MM,
## but is according to help(tcrossprod, package = "base") :
setMethod("tcrossprod", signature(x = "Matrix", y = "numeric"),
	  function(x, y = NULL) callGeneric(x, as.matrix(y)))
setMethod("tcrossprod", signature(x = "numeric", y = "Matrix"),
	  function(x, y = NULL)	 callGeneric(as.matrix(x), y))
setMethod("tcrossprod", signature(x = "Matrix", y = "matrix"),
	  function(x, y = NULL) callGeneric(x, Matrix(y)))
setMethod("tcrossprod", signature(x = "matrix", y = "Matrix"),
	  function(x, y = NULL) callGeneric(Matrix(x), y))

## maybe not 100% optimal, but elegant:
setMethod("solve", signature(a = "Matrix", b = "missing"),
	  function(a, b, ...) solve(a, Diagonal(nrow(a))))

setMethod("solve", signature(a = "Matrix", b = "numeric"),
	  function(a, b, ...) callGeneric(a, Matrix(b)))
setMethod("solve", signature(a = "Matrix", b = "matrix"),
	  function(a, b, ...) callGeneric(a, Matrix(b)))
setMethod("solve", signature(a = "matrix", b = "Matrix"),
	  function(a, b, ...) callGeneric(Matrix(a), b))

## when no sub-class method is found, bail out
setMethod("solve", signature(a = "Matrix", b = "Matrix"),
	  function(a, b, ...) .bail.out.2("solve", class(a), class(b)))

## bail-out methods in order to get better error messages
setMethod("%*%", signature(x = "Matrix", y = "Matrix"),
	  function (x, y)
          stop(gettextf('not-yet-implemented method for <%s> %%*%% <%s>',
                        class(x), class(y))))

setMethod("crossprod", signature(x = "Matrix", y = "ANY"),
	  function (x, y = NULL) .bail.out.2(.Generic, class(x), class(y)))
setMethod("crossprod", signature(x = "ANY", y = "Matrix"),
	  function (x, y = NULL) .bail.out.2(.Generic, class(x), class(y)))
setMethod("tcrossprod", signature(x = "Matrix", y = "ANY"),
	  function (x, y = NULL) .bail.out.2(.Generic, class(x), class(y)))
setMethod("tcrossprod", signature(x = "ANY", y = "Matrix"),
	  function (x, y = NULL) .bail.out.2(.Generic, class(x), class(y)))

## cheap fallbacks
setMethod("crossprod", signature(x = "Matrix", y = "Matrix"),
	  function(x, y = NULL) t(x) %*% y)
setMethod("tcrossprod", signature(x = "Matrix", y = "Matrix"),
	  function(x, y = NULL) x %*% t(y))

## There are special sparse methods; this is a "fall back":
setMethod("kronecker", signature(X = "Matrix", Y = "ANY",
				 FUN = "ANY", make.dimnames = "ANY"),
	  function(X, Y, FUN, make.dimnames, ...) {
	      if(is(X, "sparseMatrix"))
		  warning("using slow kronecker() method")
	      X <- as(X, "matrix") ; Matrix(callGeneric()) })

setMethod("kronecker", signature(X = "ANY", Y = "Matrix",
				 FUN = "ANY", make.dimnames = "ANY"),
	  function(X, Y, FUN, make.dimnames, ...) {
	      if(is(Y, "sparseMatrix"))
		  warning("using slow kronecker() method")
	      Y <- as(Y, "matrix") ; Matrix(callGeneric()) })


## FIXME: All of these should never be called
setMethod("chol", signature(x = "Matrix"),
	  function(x, pivot = FALSE, ...) .bail.out.1(.Generic, class(x)))
setMethod("determinant", signature(x = "Matrix"),
	  function(x, logarithm = TRUE, ...) .bail.out.1(.Generic, class(x)))

setMethod("diag", signature(x = "Matrix"),
	  function(x, nrow, ncol) .bail.out.1(.Generic, class(x)))
setMethod("t", signature(x = "Matrix"),
	  function(x) .bail.out.1(.Generic, class(x)))

setMethod("norm", signature(x = "Matrix", type = "character"),
	  function(x, type, ...) .bail.out.1(.Generic, class(x)))
setMethod("rcond", signature(x = "Matrix", norm = "character"),
	  function(x, norm, ...) .bail.out.1(.Generic, class(x)))


## for all :
setMethod("norm", signature(x = "ANY", type = "missing"),
	  function(x, type, ...) norm(x, type = "O", ...))
setMethod("rcond", signature(x = "ANY", norm = "missing"),
	  function(x, norm, ...) rcond(x, norm = "O", ...))





## MM: More or less "Cut & paste" from
## --- diff.default() from  R/src/library/base/R/diff.R :
setMethod("diff", signature(x = "Matrix"),
	  function(x, lag = 1, differences = 1, ...) {
	      if (length(lag) > 1 || length(differences) > 1 ||
		  lag < 1 || differences < 1)
		  stop("'lag' and 'differences' must be integers >= 1")
	      xlen <- nrow(x)
	      if (lag * differences >= xlen)
		  return(x[,FALSE][0])	# empty of proper mode

	      i1 <- -1:-lag
	      for (i in 1:differences)
		  x <- x[i1, , drop = FALSE] -
		      x[-nrow(x):-(nrow(x)-lag+1), , drop = FALSE]
	      x
	  })

setMethod("image", "Matrix",
	  function(x, ...) { # coercing to sparse is not inefficient,
	      ##	       since we need 'i' and 'j' for levelplot()
	      x <- as(as(x, "sparseMatrix"), "dMatrix")
	      callGeneric()
	  })


## Group Methods

##-> see ./Ops.R
##         ~~~~~
## For all  non-dMatrix objects, and note that  "all" and "any" have their own
setMethod("Summary", signature(x = "Matrix", na.rm = "ANY"),
	  function(x, ..., na.rm)
	  callGeneric(as(x,"dMatrix"), ..., na.rm = na.rm))


### --------------------------------------------------------------------------
###
### Subsetting "["  and
### SubAssign  "[<-" : The "missing" cases can be dealt with here, "at the top":

## Using "index" for indices should allow
## integer (numeric), logical, or character (names!) indices :

## "x[]":
setMethod("[", signature(x = "Matrix",
			 i = "missing", j = "missing", drop = "ANY"),
	  function (x, i, j, ..., drop) x)

## missing 'drop' --> 'drop = TRUE'
##                     -----------
## select rows __ or __ vector indexing:
setMethod("[", signature(x = "Matrix", i = "index", j = "missing",
			 drop = "missing"),
	  function(x,i,j, ..., drop) {
	      if(nargs() == 2) { ## e.g. M[0] , M[TRUE],  M[1:2]
		  if(any(as.logical(i)) || prod(dim(x)) == 0)
                      ## FIXME: for *large sparse*, use sparseVector !
		      as.vector(x)[i]
		  else ## save memory (for large sparse M):
		      as.vector(x[1,1])[FALSE]
	      } else {
		  callGeneric(x, i=i, , drop=TRUE)
		  ##		      ^^
	      }
	  })

## select columns
setMethod("[", signature(x = "Matrix", i = "missing", j = "index",
			 drop = "missing"),
	  function(x,i,j, ..., drop) callGeneric(x, j=j, drop= TRUE))
setMethod("[", signature(x = "Matrix", i = "index", j = "index",
                         drop = "missing"),
	  function(x,i,j, ..., drop) callGeneric(x, i=i, j=j, drop= TRUE))

## bail out if any of (i,j,drop) is "non-sense"
setMethod("[", signature(x = "Matrix", i = "ANY", j = "ANY", drop = "ANY"),
	  function(x,i,j, ..., drop)
          stop("invalid or not-yet-implemented 'Matrix' subsetting"))

## logical indexing, such as M[ M >= 7 ] *BUT* also M[ M[,1] >= 3,],
## The following is *both* for    M [ <logical>   ]
##                 and also for   M [ <logical> , ]
.M.sub.i.logical <- function (x, i, j, ..., drop)
{
    nA <- nargs()
    if(nA == 2) { ##  M [ M >= 7 ]
	## FIXME: when both 'x' and 'i' are sparse, this can be very inefficient
	if(is(x, "sparseMatrix"))
	    message("<sparse>[ <logic> ] : .M.sub.i.logical() maybe inefficient")
	toC <- geClass(x)
	if(canCoerce(x, toC)) as(x, toC)@x[as.vector(i)]
	else as(as(as(x, "generalMatrix"), "denseMatrix"), toC)@x[as.vector(i)]
	## -> error when lengths don't match
    } else if(nA == 3) { ##  M [ M[,1, drop=FALSE] >= 7, ]
	stop("not-yet-implemented 'Matrix' subsetting") ## FIXME

    } else stop("nargs() = ", nA,
		".  Extraneous illegal arguments inside '[ .. ]' (i.logical)?")
}
setMethod("[", signature(x = "Matrix", i = "lMatrix", j = "missing",
			 drop = "ANY"),
	  .M.sub.i.logical)
setMethod("[", signature(x = "Matrix", i = "logical", j = "missing",
			 drop = "ANY"),
	  .M.sub.i.logical)


subset.ij <- function(x, ij) {
    m <- nrow(ij)
    if(m > 3) {
        cld <- getClassDef(class(x))
	sym.x <- extends(cld, "symmetricMatrix")
	if(sym.x) {
	    W <- if(x@uplo == "U") # stored only [i,j] with i <= j
		ij[,1] > ij[,2] else ij[,1] < ij[,2]
	    if(any(W))
		ij[W,] <- ij[W, 2:1]
        }
        if(extends(cld, "sparseMatrix")) {
	    ## do something smarter:
	    nr <- nrow(x)
	    if(!extends(cld, "CsparseMatrix")) {
		x <- as(x, "CsparseMatrix") # simpler; our standard
		cld <- getClassDef(class(x))
	    }
	    tri.x <- extends(cld, "triangularMatrix")
	    if(tri.x) {
		## need these for the 'x' slot in any case
		if (x@diag == "U") x <- .Call(Csparse_diagU2N, x)
		## slightly more efficient than non0.i() or non0ind():
		ij.x <- .Call(compressed_non_0_ij, x, isC=TRUE)
	    } else { ## symmetric / general : for symmetric, only "existing"b
		ij.x <- non0.i(x, cld)
	    }

	    mi <- match(encodeInd(ij.x,	  nr),
			encodeInd(ij -1L, nr), nomatch=0)
	    mmi <- mi != 0
	    ## Result:
	    ans <- vector(mode = .type.kind[.M.kindC(cld)], length = m)
	    ## those that are *not* zero:
	    ans[mi[mmi]] <-
		if(extends(cld, "nsparseMatrix")) TRUE else x@x[mmi]
	    ans

        } else { ## non-sparse : dense
            ##---- NEVER happens:  'denseMatrix' has its own setMethod(.) !
            message("m[ <ij-matrix> ]: inefficiently indexing single elements")
            i1 <- ij[,1]
            i2 <- ij[,2]
            ## very inefficient for large m
            unlist(lapply(seq_len(m), function(j) x[i1[j], i2[j]]))
        }
    } else { # 1 <= m <= 3
        i1 <- ij[,1]
        i2 <- ij[,2]
        unlist(lapply(seq_len(m), function(j) x[i1[j], i2[j]]))
    }
}

## A[ ij ]  where ij is (i,j) 2-column matrix -- but also when that is logical mat!
.M.sub.i.2col <- function (x, i, j, ..., drop)
{
    nA <- nargs()
    if(nA == 2) { ##  M [ cbind(ii,jj) ] or M [ <logical matrix> ]
	if(!is.integer(nc <- ncol(i)))
	    stop(".M.sub.i.2col(): 'i' has no integer column number;\n",
		 "should never happen; please report")
	if(is.logical(i))
	    return(.M.sub.i.logical(x, i=i)) # call with 2 args!
	else if(!is.numeric(i) || nc != 2)
	    stop("such indexing must be by logical or 2-column numeric matrix")
	m <- nrow(i)
        if(m == 0) return(vector(mode = .type.kind[.M.kind(x)]))
        ## else
        subset.ij(x, i)

    } else stop("nargs() = ", nA,
		".  Extraneous illegal arguments inside '[ .. ]' (i.2col)?")
}
setMethod("[", signature(x = "Matrix", i = "matrix", j = "missing"),# drop="ANY"
	  .M.sub.i.2col)


### "[<-" : -----------------

## x[] <- value :
setReplaceMethod("[", signature(x = "Matrix", i = "missing", j = "missing",
                                value = "ANY"),## double/logical/...
	  function (x, value) {
	      ## Fails for 'nMatrix' ... FIXME : make sure have method there
	      x@x <- rep(value, length = length(x@x))
	      validObject(x)# check if type and lengths above match
	      x
          })

## A[ ij ] <- value,  where ij is (i,j) 2-column matrix :
## ----------------
## The cheap general method --- FIXME: provide special ones; done for Tsparse..
## NOTE:  need '...' below such that setMethod() does
##	  not use .local() such that nargs() will work correctly:
.M.repl.i.2col <- function (x, i, j, ..., value)
{
    nA <- nargs()
    if(nA == 3) { ##  M [ cbind(ii,jj) ] <- value  or M [ Lmat ] <- value
	if(!is.integer(nc <- ncol(i)))
	    stop(".M.repl.i.2col(): 'i' has no integer column number;\n",
		 "should never happen; please report")
	else if(!is.numeric(i) || nc != 2)
	    stop("such indexing must be by logical or 2-column numeric matrix")
	if(is.logical(i)) {
	    message(".M.repl.i.2col(): drop 'matrix' case ...")
	    ## c(i) : drop "matrix" to logical vector
	    return( callGeneric(x, i=c(i), value=value) )
	}
	if(!is.integer(i)) storage.mode(i) <- "integer"
	if(any(i < 0))
	    stop("negative values are not allowed in a matrix subscript")
	if(any(is.na(i)))
	    stop("NAs are not allowed in subscripted assignments")
	if(any(i0 <- (i == 0))) # remove them
            i <- i[ - which(i0, arr.ind = TRUE)[,"row"], ]
        ## now have integer i >= 1
	m <- nrow(i)
	## mod.x <- .type.kind[.M.kind(x)]
	if(length(value) > 0 && m %% length(value) != 0)
	    warning("number of items to replace is not a multiple of replacement length")
	## recycle:
	value <- rep(value, length = m)
	i1 <- i[,1]
	i2 <- i[,2]
	if(m > 2)
	    message("m[ <ij-matrix> ] <- v: inefficiently treating single elements")
	## inefficient -- FIXME -- (also loses "symmetry" unnecessarily)
	for(k in seq_len(m))
	    x[i1[k], i2[k]] <- value[k]

	x
    } else stop("nargs() = ", nA,
		".  Extraneous illegal arguments inside '[ .. ]' ?")
}

setReplaceMethod("[", signature(x = "Matrix", i = "matrix", j = "missing",
				value = "replValue"),
	  .M.repl.i.2col)


setReplaceMethod("[", signature(x = "Matrix", i = "missing", j = "ANY",
				value = "Matrix"),
		 function (x, i, j, ..., value)
		 callGeneric(x=x, , j=j, value = as.vector(value)))

setReplaceMethod("[", signature(x = "Matrix", i = "ANY", j = "missing",
				value = "Matrix"),
		 function (x, i, j, ..., value)
		 callGeneric(x=x, i=i, , value = as.vector(value)))

setReplaceMethod("[", signature(x = "Matrix", i = "ANY", j = "ANY",
				value = "Matrix"),
		 function (x, i, j, ..., value)
		 callGeneric(x=x, i=i, j=j, value = as.vector(value)))

setReplaceMethod("[", signature(x = "Matrix", i = "missing", j = "ANY",
				value = "matrix"),
		 function (x, i, j, ..., value)
		 callGeneric(x=x, , j=j, value = c(value)))

setReplaceMethod("[", signature(x = "Matrix", i = "ANY", j = "missing",
				value = "matrix"),
		 function (x, i, j, ..., value)
		 callGeneric(x=x, i=i, , value = c(value)))

setReplaceMethod("[", signature(x = "Matrix", i = "ANY", j = "ANY",
				value = "matrix"),
		 function (x, i, j, value)
		 callGeneric(x=x, i=i, j=j, value = c(value)))

## (ANY,ANY,ANY) is used when no `real method' is implemented :
setReplaceMethod("[", signature(x = "Matrix", i = "ANY", j = "ANY",
                                value = "ANY"),
	  function (x, i, j, value) {
              if(!is.atomic(value))
		  stop(sprintf("RHS 'value' (class %s) matches 'ANY', but must match matrix class %s",
			       class(value),class(x)))
              else stop("not-yet-implemented 'Matrix[<-' method")
          })
