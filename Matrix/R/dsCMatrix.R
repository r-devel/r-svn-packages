#### Symmetric Sparse Matrices in compressed column-oriented format

setAs("dgCMatrix", "dsCMatrix",
      function(from) {
	  if(!exists(".warn.dsC")) { ## now only warn *once* ..
	      warning("as(.,\"dsCMatrix\") is deprecated; do use as(., \"symmetricMatrix\")")
	      assign(".warn.dsC", "DONE", envir = .GlobalEnv)
	  }
	  as(from, "symmetricMatrix")
      })

## Specific conversions, should they be necessary.  Better to convert as
## as(x, "TsparseMatrix") or as(x, "denseMatrix")

## Moved to ./Csparse.R
## setAs("dsCMatrix", "dsTMatrix",
##       function(from) .Call(Csparse_to_Tsparse, from, FALSE))

setAs("dsCMatrix", "dgTMatrix", # needed for show(), image()
      function(from)
      ## pre-Cholmod -- FIXME: get rid of
      .Call(dsCMatrix_to_dgTMatrix, from))

setAs("dsCMatrix", "dgeMatrix",
      function(from) as(as(from, "dgTMatrix"), "dgeMatrix"))

setAs("dsCMatrix", "matrix",
      function(from) as(as(from, "generalMatrix"), "matrix"))
setAs("matrix", "dsCMatrix",
      function(from) as(as(from, "CsparseMatrix"), "symmetricMatrix"))

setAs("dsCMatrix", "lsCMatrix",
      function(from) new("lsCMatrix", i = from@i, p = from@p, uplo = from@uplo,
                         x = as.logical(from@x),
                         Dim = from@Dim, Dimnames = from@Dimnames))
setAs("dsCMatrix", "nsCMatrix",
      function(from) new("nsCMatrix", i = from@i, p = from@p, uplo = from@uplo,
                         Dim = from@Dim, Dimnames = from@Dimnames))

setAs("dsCMatrix", "dgCMatrix",
      function(from) .Call(Csparse_symmetric_to_general, from))

setAs("dsCMatrix", "dsyMatrix",
      function(from) as(from, "denseMatrix"))

## have rather tril() and triu() methods than
## setAs("dsCMatrix", "dtCMatrix", ....)
setMethod("tril", "dsCMatrix",
	  function(x, k = 0, ...) {
	      if(x@uplo == "L" && k == 0)
		  ## same internal structure (speedup potential !?)
		  new("dtCMatrix", uplo = x@uplo, i = x@i, p = x@p,
		      x = x@x, Dim = x@Dim, Dimnames = x@Dimnames)
	      else tril(as(x, "dgCMatrix"), k = k, ...)
	  })

setMethod("triu", "dsCMatrix",
	  function(x, k = 0, ...) {
	      if(x@uplo == "U" && k == 0)
		  ## same internal structure (speedup potential !?)
		  new("dtCMatrix", uplo = x@uplo, i = x@i, p = x@p,
		      x = x@x, Dim = x@Dim, Dimnames = x@Dimnames)
	      else triu(as(x, "dgCMatrix"), k = k, ...)
	  })

setMethod("solve", signature(a = "dsCMatrix", b = "ddenseMatrix"),
          function(a, b, ...) {
              if (class(b) != "dgeMatrix")
                  b <- .Call(dup_mMatrix_as_dgeMatrix, b)
              .Call(dsCMatrix_matrix_solve, a, b)
          },
          valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dsCMatrix", b = "matrix"),
          function(a, b, ...)
          .Call(dsCMatrix_matrix_solve, a,
                .Call(dup_mMatrix_as_dgeMatrix, b)),
          valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dsCMatrix", b = "numeric"),
          function(a, b, ...)
          .Call(dsCMatrix_matrix_solve, a,
                .Call(dup_mMatrix_as_dgeMatrix, b)),
          valueClass = "dgeMatrix")

## `` Fully-sparse'' solve() :
setMethod("solve", signature(a = "dsCMatrix", b = "dsparseMatrix"),
	  function(a, b, ...) {
	      if (!is(b, "CsparseMatrix"))
		  b <- as(b, "CsparseMatrix")
	      if (is(b, "symmetricMatrix")) ## not supported (yet) by cholmod_spsolve
		  b <- as(b, "dgCMatrix")
	      .Call(dsCMatrix_Csparse_solve, a, b)
	  })


setMethod("chol", signature(x = "dsCMatrix"),
	  function(x, pivot = FALSE, ...) .Call(dsCMatrix_chol, x, pivot),
	  valueClass = "dtCMatrix")

setMethod("Cholesky", signature(A = "dsCMatrix"),
          function(A, perm = TRUE, LDL = !super, super = FALSE, Imult = 0, ...)
          .Call(dsCMatrix_Cholesky, A, perm, LDL, super, Imult))


setMethod("t", signature(x = "dsCMatrix"),
          function(x) .Call(Csparse_transpose, x, FALSE),
          valueClass = "dsCMatrix")

setMethod("determinant", signature(x = "dsCMatrix", logarithm = "missing"),
          function(x, logarithm, ...) determinant(x, TRUE))

.diag.dsC <- function(x, Chx = Cholesky(x, LDL=TRUE), res.kind = "diag") {
    force(Chx)
    stopifnot(is.integer(Chx@p), is.double(Chx@x))
    .Call(diag_tC, Chx@p, Chx@x, Chx@perm, res.kind)
}

## FIXME:  kind = "diagBack" is not yet implemented
##	would be much more efficient, but there's no CHOLMOD UI (?)
##
## Note: for det(), permutation is unimportant;
##       for diag(), apply *inverse* permutation
##    	q <- p ; q[q] <- seq_along(q); q

ldet1.dsC <- function(x, ...) .Call(CHMfactor_ldetL2, Cholesky(x, ...))
ldet2.dsC <- function(x, ...) {
    Ch <- Cholesky(x, super = FALSE, ...)
    .Call(diag_tC, Ch@p, Ch@x, Ch@perm, "sumLog")
}
ldet3.dsC <- function(x, perm = TRUE)
    .Call(dsCMatrix_LDL_D, x, perm=perm, "sumLog")

setMethod("determinant", signature(x = "dsCMatrix", logarithm = "logical"),
	  function(x, logarithm, ...)
      {
          ## Chx <- Cholesky(x, LDL=TRUE)
          ## ldet <- .Call(diag_tC, Chx@p, Chx@x, Chx@perm, res.kind = "sumLog")
          ## or
          ## ldet <- .Call("CHMfactor_ldetL2", Chx) # which would also work
          ##                                 when Chx <- Cholesky(x, super=TRUE)

### FIXME: not okay when the matrix is *NOT* pos.def.
          ldet <- .Call(dsCMatrix_LDL_D, x, perm=TRUE, "sumLog")
	  modulus <- if (logarithm) ldet else exp(ldet)
	  attr(modulus, "logarithm") <- logarithm
	  structure(list(modulus = modulus, sign = as.integer(1)),
		    class = "det")
      })

## setMethod("writeHB", signature(obj = "dsCMatrix"),
##           function(obj, file, ...) {
##               .Deprecated("writeMM")
##               .Call(Matrix_writeHarwellBoeing,
##                     if (obj@uplo == "U") t(obj) else obj,
##                     as.character(file), "DSC")
##           })
