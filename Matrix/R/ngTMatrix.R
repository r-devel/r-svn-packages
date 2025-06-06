#### Nonzero Pattern Sparse Matrices in triplet format

### contains = "nsparseMatrix"
###             ============= ---> superclass methods in ./nsparseMatrix.R


setAs("ngTMatrix", "ngeMatrix",
      function(from) .Call(lgTMatrix_to_lgeMatrix, as(from,"lgTMatrix")))

setAs("ngTMatrix", "matrix",
      function(from) .Call(lgTMatrix_to_matrix, as(from, "lgTMatrix")))
## setAs("ngTMatrix", "matrix", # go via fast C code:
##       function(from) as(as(from, "ngCMatrix"), "matrix"))

setAs("matrix", "ngTMatrix",
      function(from) {
	  if(!is.logical(from))
	      storage.mode(from) <- "logical"
	  if(any(is.na(from)))
	      warning("'NA's coerced to 'FALSE' in coercion to logical sparse")
          dn <- dimnames(from)
          if(is.null(dn))
              dn <- list(NULL,NULL)
          else dimnames(from) <- NULL # such that which(.) does not see any:
	  ij <- which(from, arr.ind = TRUE) - 1L
	  if(length(ij) == 0) ij <- matrix(ij, 0, 2)
	  new("ngTMatrix",
	      i = ij[,1],
	      j = ij[,2],
	      Dim = as.integer(dim(from)),
	      Dimnames = dn)
	  })

setAs("matrix", "nMatrix", function(from) as(from, "ngTMatrix"))


setAs("ngTMatrix", "dgTMatrix",
      function(from)
      ## more efficient than
      ## as(as(as(sM, "ngCMatrix"), "dgCMatrix"), "dgTMatrix")
      new("dgTMatrix", i = from@i, j = from@j,
	  x = rep.int(1, length(from@i)),
	  ## cannot copy factors, but can we use them?
	  Dim = from@Dim, Dimnames= from@Dimnames))

setAs("ngTMatrix", "lgTMatrix",
      function(from)
      new("lgTMatrix", i = from@i, j = from@j,
	  x = rep.int(TRUE, length(from@i)),
	  ## cannot copy factors, but can we use them?
	  Dim = from@Dim, Dimnames= from@Dimnames))

setAs("ngTMatrix", "ntTMatrix",
      function(from) check.gT2tT(from, cl = "ngTMatrix", toClass = "ntTMatrix"))


setMethod("t", signature(x = "ngTMatrix"),
	  function(x) new("ngTMatrix", i = x@j, j = x@i,
			  Dim = x@Dim[2:1],
			  Dimnames= x@Dimnames[2:1]),
	  valueClass = "ngTMatrix")
