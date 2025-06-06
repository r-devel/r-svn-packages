### Coercion and Methods for Dense Numeric Symmetric Matrices

setAs("dgeMatrix", "dsyMatrix",
     function(from) {
	if(isSymmetric(from))# < with tolerance!
	    .Call(dense_to_symmetric, from, "U", FALSE)
	else
	    stop("not a symmetric matrix; consider forceSymmetric() or symmpart()")
     })
## NB: The alternative, 'zero tolerance' { <=> isSymmetric(*, tol=0) }
##     breaks too much previous code -- though it would be much faster --
## setAs("dgeMatrix", "dsyMatrix",
##       function(from) .Call(dense_to_symmetric, from, "U", TRUE))


setAs("matrix", "dsyMatrix",
      function(from) as(as(from, "dgeMatrix"), "dsyMatrix"))

setAs("dsyMatrix", "matrix",
      function(from) .Call(dsyMatrix_as_matrix, from, TRUE))

setAs("dsyMatrix", "dspMatrix",
      function(from) .Call(dsyMatrix_as_dspMatrix, from))

setAs("dsyMatrix", "dsTMatrix",
      function(from) { # 'dsT': only store upper *or* lower
	  ## working via "matrix" : not very efficient	(FIXME)
	  m <- .Call(dsyMatrix_as_matrix, from, FALSE) # no dimnames!
	  ij <- which(m != 0, arr.ind = TRUE)
	  uplo <- from@uplo
	  ij <- ij[if(uplo == "U") ij[,1] <= ij[,2] else ij[,1] >= ij[,2] ,
		   , drop = FALSE]
	  new("dsTMatrix", i = ij[,1] - 1L, j = ij[,2] - 1L,
	      x = as.vector(m[ij]), uplo = uplo,
	      Dim = from@Dim, Dimnames = from@Dimnames)
      })

setAs("dsyMatrix", "dsCMatrix",
      function(from) as(as(from, "dsTMatrix"), "dsCMatrix"))


## Note: Just *because* we have an explicit  dtr -> dge coercion,
##       show( <ddenseMatrix> ) is not okay, and we need our own:
setMethod("show", "dsyMatrix", function(object) prMatrix(object))


setMethod("rcond", signature(x = "dsyMatrix", norm = "character"),
          function(x, norm, ...)
          .Call(dsyMatrix_rcond, x, norm),
          valueClass = "numeric")

setMethod("rcond", signature(x = "dsyMatrix", norm = "missing"),
          function(x, norm, ...)
          .Call(dsyMatrix_rcond, x, "O"),
          valueClass = "numeric")

setMethod("%*%", signature(x = "dsyMatrix", y = "ddenseMatrix"),
          function(x, y) .Call(dsyMatrix_matrix_mm, x, y, FALSE))
setMethod("%*%", signature(x = "dsyMatrix", y = "matrix"),
          function(x, y) .Call(dsyMatrix_matrix_mm, x, y, FALSE))

setMethod("%*%", signature(x = "ddenseMatrix", y = "dsyMatrix"),
          function(x, y) .Call(dsyMatrix_matrix_mm, y, x, TRUE))
setMethod("%*%", signature(x = "matrix", y = "dsyMatrix"),
          function(x, y) .Call(dsyMatrix_matrix_mm, y, x, TRUE))

setMethod("solve", signature(a = "dsyMatrix", b = "missing"),
          function(a, b, ...) .Call(dsyMatrix_solve, a),
          valueClass = "dsyMatrix")

setMethod("solve", signature(a = "dsyMatrix", b = "matrix"),
          function(a, b, ...) .Call(dsyMatrix_matrix_solve, a, b),
          valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dsyMatrix", b = "ddenseMatrix"),
          function(a, b, ...) .Call(dsyMatrix_matrix_solve, a, b),
          valueClass = "dgeMatrix")

setMethod("norm", signature(x = "dsyMatrix", type = "character"),
          function(x, type, ...) .Call(dsyMatrix_norm, x, type),
          valueClass = "numeric")

setMethod("norm", signature(x = "dsyMatrix", type = "missing"),
          function(x, type, ...) .Call(dsyMatrix_norm, x, "O"),
          valueClass = "numeric")

## Should this create the opposite storage format - i.e. "U" -> "L"
## and vice-versa?
## MM: I think yes, since the other part can be filled arbitrarily (wrongly)
##WAS setMethod("t", signature(x = "dsyMatrix"), function(x) x)
setMethod("t", signature(x = "dsyMatrix"), t_trMatrix,
          valueClass = "dsyMatrix")

setMethod("BunchKaufman", signature(x = "dsyMatrix"),
	  function(x) .Call(dsyMatrix_trf, x))

## The following has the severe effect of making
## "dsyMatrix" a subclass of "dpoMatrix" and since the reverse is
## by definition of "dpoMatrix", the class-hierarchy gets a *cycle* !
##
setIs("dsyMatrix", "dpoMatrix",
      test = function(obj)
          "try-error" != class(try(.Call(dpoMatrix_chol, obj), silent=TRUE)),
      replace = function(obj, value) { ## copy all slots (is needed)
          for(n in slotNames(obj)) slot(obj, n) <- slot(value, n)
          obj
      })

## Now that we have "chol", we can define  "determinant" methods,
## exactly like in ./dsCMatrix.R
## DB - Probably figure out how to use the BunchKaufman decomposition instead
## {{FIXME: Shouldn't it be possible to have "determinant" work by
## default automatically for "Matrix"es  when there's a "chol" method available?
## ..> work with ss <- selectMethod("chol", signature("dgCMatrix"))
## -- not have to define showMethod("determinant", ...) for all classes

