#### Triangular Packed Matrices -- Coercion and Methods

setAs("dtpMatrix", "dtrMatrix",
      function(from) .Call(dtpMatrix_as_dtrMatrix, from))

## Is this needed?  already have coercion to "TsparseMatrix" {FIXME}
setAs("dtpMatrix", "dtTMatrix",
      function(from) {
	  x <- as(from, "TsparseMatrix")
          cld <- getClassDef(class(x))
	  if(extends(cld, "dtTMatrix"))
	      x
	  else { ## triangularity lost: should not have happened
	      warning("inefficient coercion (lost triangularity); please report")
	      gT2tT(as(x, "dgTMatrix"), uplo = from@uplo, diag = from@diag,
		    cl = "dgTMatrix", toClass = "dtTMatrix", cld = cld)
	  }
      })

setAs("dtpMatrix", "matrix",
      function(from) as(as(from, "dtrMatrix"), "matrix"))
setAs("matrix", "dtpMatrix",
      function(from) as(as(from, "dtrMatrix"), "dtpMatrix"))

setAs("pCholesky", "lMatrix",
      function(from) as(as(from, "dtpMatrix"), "lMatrix"))
setAs("pBunchKaufman", "lMatrix",
      function(from) as(as(from, "dtpMatrix"), "lMatrix"))


setMethod("%*%", signature(x = "dtpMatrix", y = "ddenseMatrix"),
	  function(x, y) .Call(dtpMatrix_matrix_mm, x, y))
setMethod("%*%", signature(x = "dgeMatrix", y = "dtpMatrix"),
	  function(x, y) .Call(dgeMatrix_dtpMatrix_mm, x, y))
## DB: I don't think this is needed any more
## %*% should always work for  <fooMatrix> %*% <fooMatrix>
## setMethod("%*%", signature(x = "dtpMatrix", y = "dtpMatrix"),
##           function(x, y)
##           ## FIXME: this is cheap; could we optimize chosing the better of
##           ## callGeneric(x, as(y, "dgeMatrix"))  and
##           ## callGeneric(as(x "dgeMatrix"), y))  depending on their 'uplo' ?
##           callGeneric(x, as(y, "dgeMatrix")))

## dtpMatrix <-> matrix : will be used by the "numeric" one
setMethod("%*%", signature(x = "dtpMatrix", y = "matrix"),
          function(x, y) .Call(dtpMatrix_matrix_mm, x, y))
setMethod("%*%", signature(x = "matrix", y = "dtpMatrix"),
          function(x, y) callGeneric(as(x, "dgeMatrix"), y))

## dtpMatrix <-> numeric : the auxiliary functions are R version specific!
##setMethod("%*%", signature(x = "dtpMatrix", y = "numeric"), .M.v)
##setMethod("%*%", signature(x = "numeric", y = "dtpMatrix"), .v.M)


setMethod("determinant", signature(x = "dtpMatrix", logarithm = "missing"),
	  function(x, logarithm, ...) determinant(x, TRUE))

setMethod("determinant", signature(x = "dtpMatrix", logarithm = "logical"),
	  function(x, logarithm, ...) {
	      dg <- diag(x)
	      if (logarithm) {
		  modulus <- sum(log(abs(dg)))
		  sgn <- prod(sign(dg))
	      } else {
		  modulus <- prod(dg)
		  sgn <- sign(modulus)
		  modulus <- abs(modulus)
	      }
	      attr(modulus, "logarithm") <- logarithm
	      val <- list(modulus = modulus, sign = sgn)
	      class(val) <- "det"
	      val
	  })

setMethod("diag", signature(x = "dtpMatrix"),
	  function(x, nrow, ncol) .Call(dtpMatrix_getDiag, x),
	  valueClass = "numeric")

setMethod("norm", signature(x = "dtpMatrix", type = "character"),
	  function(x, type, ...) .Call(dtpMatrix_norm, x, type),
	  valueClass = "numeric")

setMethod("norm", signature(x = "dtpMatrix", type = "missing"),
	  function(x, type, ...) .Call(dtpMatrix_norm, x, "O"),
	  valueClass = "numeric")

setMethod("rcond", signature(x = "dtpMatrix", norm = "character"),
	  function(x, norm, ...)
	  .Call(dtpMatrix_rcond, x, norm),
	  valueClass = "numeric")

setMethod("rcond", signature(x = "dtpMatrix", norm = "missing"),
	  function(x, norm, ...)
	  .Call(dtpMatrix_rcond, x, "O"),
	  valueClass = "numeric")

setMethod("solve", signature(a = "dtpMatrix", b="missing"),
	  function(a, b, ...) .Call(dtpMatrix_solve, a),
	  valueClass = "dtpMatrix")

setMethod("solve", signature(a = "dtpMatrix", b="ddenseMatrix"),
	  function(a, b, ...) .Call(dtpMatrix_matrix_solve, a, b),
	  valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dtpMatrix", b="matrix"),
	  function(a, b, ...) .Call(dtpMatrix_matrix_solve, a, b),
	  valueClass = "dgeMatrix")

setMethod("t", signature(x = "dtpMatrix"),
          function(x) as(t(as(x, "dtrMatrix")), "dtpMatrix"),
          valueClass = "dtpMatrix")

setMethod("unpack", signature(x = "dtpMatrix"),
          function(x, ...) as(x, "dtrMatrix"),
          valueClass = "dtrMatrix")
###
