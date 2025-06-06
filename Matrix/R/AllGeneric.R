#### Define those generics that we need, if they don't exist;
#### not all will be exported

    setGeneric("expand", function(x, ...) standardGeneric("expand"))

## tcrossprod() is now in R's base

## these two are *not* exported (yet) :

    setGeneric("isDiagonal", function(object, ...)
               standardGeneric("isDiagonal"))

    setGeneric("isTriangular", function(object, ...) ## 'upper = NA'
               standardGeneric("isTriangular"))

## isSymmetric is "S3 generic" in R's  base/R/eigen.R

    setGeneric("facmul",
               function(x, factor, y, transpose, left, ...)
               standardGeneric("facmul"))

setGeneric("BunchKaufman", function(x, ...) standardGeneric("BunchKaufman"))
setGeneric("lu", function(x, ...) standardGeneric("lu"))

##NB ## do not redefine the "base signature"

##NB setGeneric("chol", def = function(x, pivot= FALSE,...) standardGeneric("chol"),
##NB            useAsDefault= function(x, pivot= FALSE,...) base::chol(x, pivot, ...))

##NB     setGeneric("qr", def =   function(x, tol=1e-7,...) standardGeneric("qr"),
##NB                useAsDefault= function(x, tol=1e-7,...) base::qr(x, tol, ...))

    setGeneric("norm", function(x, type, ...) standardGeneric("norm"))

## base::rcond() has additional argument 'triangular'
## which should *not* be part of the signature
setGeneric("rcond", function(x, norm, ...) standardGeneric("rcond"),
	   useAsDefault =
           ## for now:
           if(existsFunction("rcond", where=baseenv()))
           ## later: if(exists(getRversion() >= "2.7.0")
	   function(x, norm, ...) base::rcond(x, norm=norm, ...))


    setGeneric("Schur", function(x, vectors, ...) standardGeneric("Schur"))

    setGeneric("unpack", function(x, ...) standardGeneric("unpack"))

##-     setGeneric("%p%", function(a, b) standardGeneric("%p%"))

    setGeneric("expm", function(x) standardGeneric("expm"))

##     setGeneric("writeHB", function(obj, file, ...)
##                standardGeneric("writeHB"))

    setGeneric("writeMM", function(obj, file, ...)
               standardGeneric("writeMM"))

##     setGeneric("qqmath", function(x, data, ...)
##                standardGeneric("qqmath"))

    setGeneric("tril", function(x, k = 0, ...)
               standardGeneric("tril"))

    setGeneric("triu", function(x, k = 0, ...)
               standardGeneric("triu"))

    setGeneric("band", function(x, k1, k2, ...)
               standardGeneric("band"))

    setGeneric("Cholesky",
	       function(A, perm = TRUE, LDL = !super, super = FALSE,
			Imult = 0, ...)
               standardGeneric("Cholesky"))

setGeneric("symmpart", function(x) standardGeneric("symmpart"))
setGeneric("skewpart", function(x) standardGeneric("skewpart"))

## A version of coercion to  "symmetric" which does *NOT* check,
## but just takes the ## upper (or lower) values and
## ``declares'' the symmetric:
setGeneric("forceSymmetric",
	   function(x, uplo) standardGeneric("forceSymmetric"))



###---- Group Generics ----

if (getRversion() < "2.6.0") {
## The following are **WORKAROUND** s currently needed for all non-Primitives:

## [The following is more future-proof than direct  setGeneric(.) calls:
## FIX (in R!) : "trunc" should really be in Math, but we try both for the time

## "Math"
for(fname in intersect(getGroupMembers("Math"),
		       c("log", "log2", "log10", "logb", "log1p", "expm1",
			 "gamma", "lgamma", "digamma", "trigamma",
			 "cummax", "cummin", "trunc")))
    if(!is.primitive(get(fname))) setGeneric(fname, group="Math")

## "Math2"
for(fname in intersect(getGroupMembers("Math2"),
		       c("round", "signif", "trunc")))
    if (!is.primitive(get(fname))) setGeneric(fname, group="Math2")

## "Summary"

    ## --- some hoop jumping that is needed for R versions <= 2.5.x

.max_def <- function(x, ..., na.rm = FALSE) base::max(x, ..., na.rm = na.rm)
.min_def <- function(x, ..., na.rm = FALSE) base::min(x, ..., na.rm = na.rm)
.range_def <- function(x, ..., na.rm = FALSE) base::range(x, ..., na.rm = na.rm)
.prod_def <- function(x, ..., na.rm = FALSE) base::prod(x, ..., na.rm = na.rm)
.sum_def <- function(x, ..., na.rm = FALSE) base::sum(x, ..., na.rm = na.rm)
.any_def <- function(x, ..., na.rm = FALSE) base::any(x, ..., na.rm = na.rm)
.all_def <- function(x, ..., na.rm = FALSE) base::all(x, ..., na.rm = na.rm)

setGeneric("max", function(x, ..., na.rm = FALSE) standardGeneric("max"),
           useAsDefault = .max_def, group = "Summary")
setGeneric("min", function(x, ..., na.rm = FALSE) standardGeneric("min"),
           useAsDefault = .min_def, group="Summary")
setGeneric("range", function(x, ..., na.rm = FALSE) standardGeneric("range"),
           useAsDefault = .range_def, group="Summary")
setGeneric("prod", function(x, ..., na.rm = FALSE) standardGeneric("prod"),
           useAsDefault = .prod_def, group="Summary")
setGeneric("sum", function(x, ..., na.rm = FALSE) standardGeneric("sum"),
           useAsDefault = .sum_def, group="Summary")
setGeneric("any", function(x, ..., na.rm = FALSE) standardGeneric("any"),
           useAsDefault = .any_def, group="Summary")
setGeneric("all", function(x, ..., na.rm = FALSE) standardGeneric("all"),
           useAsDefault = .all_def, group="Summary")

} # end{hoop jumping etc for R versions <= 2.5.1 }


## Add '...' so our methods can add  'sparseResult':
setGeneric("colSums",
	   def = function(x, na.rm = FALSE, dims = 1, ...) standardGeneric("colSums"),
	   useAsDefault = function(x,  na.rm = FALSE, dims = 1, ...)
	   base::colSums(x, na.rm=na.rm, dims=dims))
setGeneric("colMeans",
	   def = function(x, na.rm = FALSE, dims = 1, ...) standardGeneric("colMeans"),
	   useAsDefault = function(x,  na.rm = FALSE, dims = 1, ...)
	   base::colMeans(x, na.rm=na.rm, dims=dims))
setGeneric("rowSums",
	   def = function(x, na.rm = FALSE, dims = 1, ...) standardGeneric("rowSums"),
	   useAsDefault = function(x,  na.rm = FALSE, dims = 1, ...)
	   base::rowSums(x, na.rm=na.rm, dims=dims))
setGeneric("rowMeans",
	   def = function(x, na.rm = FALSE, dims = 1, ...) standardGeneric("rowMeans"),
	   useAsDefault = function(x,  na.rm = FALSE, dims = 1, ...)
	   base::rowMeans(x, na.rm=na.rm, dims=dims))
