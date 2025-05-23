#### Logical Sparse Matrices in Compressed column-oriented format

### contains = "lsparseMatrix"

## Can use CsparseMatrix methods for all of these

## setMethod("%*%", signature(x = "lgCMatrix", y = "lgCMatrix"),
##           function(x, y)
##           .Call(lgCMatrix_lgCMatrix_mm, x, y),
##           valueClass = "lgCMatrix")

## setMethod("t", signature(x = "lgCMatrix"),
##           function(x) .Call(lgCMatrix_trans, x),
##           valueClass = "lgCMatrix")


## setMethod("diag", signature(x = "lgCMatrix"),
## 	  function(x, nrow, ncol) .Call(lgCMatrix_diag, x))

## setMethod("crossprod", signature(x = "lgCMatrix", y = "missing"),
## 	  function(x, y = NULL)
##           .Call(lgCMatrix_crossprod, x, TRUE, NULL),
## 	  valueClass = "lsCMatrix")

## setMethod("tcrossprod", signature(x = "lgCMatrix", y = "missing"),
## 	  function(x, y = NULL)
##           .Call(lgCMatrix_crossprod, x, FALSE, NULL),
## 	  valueClass = "lsCMatrix")

setAs("lgCMatrix", "dgCMatrix",
      function(from) new("dgCMatrix", i = from@i, p = from@p,
                         x = as.double(from@x),
                         Dim = from@Dim, Dimnames = from@Dimnames))

setAs("lgCMatrix", "lgTMatrix",
      function(from) new("lgTMatrix", i = from@i, x = from@x,
                         j = .Call(Matrix_expand_pointers, from@p),
                         Dim = from@Dim, Dimnames = from@Dimnames))

setAs("lgCMatrix", "lgeMatrix",
      function(from)
	  new("lgeMatrix", x = c(as(from, "matrix")), # is fast,
	      Dim = from@Dim, Dimnames = from@Dimnames))

setAs("lgCMatrix", "matrix",
      function(from) .Call(lcsc_to_matrix, from))
## not this: .Call(Csparse_to_matrix, from)), since it goes via dense -> double precision


## TODO (maybe): write  matrix_to_lcsc()  in ../src/lgCMatrix.c
setAs("matrix", "lgCMatrix",
      function(from) as(as(from, "lgTMatrix"), "lgCMatrix"))


setMethod("image", "lgCMatrix",
          function(x, ...) {
              x <- as(as(x, "dgCMatrix"), "dgTMatrix")
              callGeneric()
          })

