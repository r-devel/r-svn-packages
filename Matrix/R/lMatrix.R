setAs("matrix", "lMatrix",
      function(from) { storage.mode(from) <- "logical" ; Matrix(from) })

setAs("lMatrix", "nMatrix",
      function(from) {
	  if(any(is.na(from@x)))
	      stop("\"lMatrix\" object with NAs cannot be coerced to \"nMatrix\"")
	  ## i.e. from@x are only TRUE (or FALSE in dense case)
	  cl <- class(from)
	  if(extends(cl, "diagonalMatrix")) # have no "ndi*" etc class
	      cl <- class(from <- as(from, "sparseMatrix"))
	  nCl <- sub("^l", "n", cl)
	  sNams <- slotNames(if(extends(cl, "sparseMatrix")) .sp.class(cl) else cl)
	  r <- new(nCl)# default => no validity check; and copy slots:
	  for(nm in sNams)
	      slot(r, nm) <- slot(from, nm)
	  r
      })

## and the reverse as well :

setAs("nMatrix", "lMatrix",
      function(from) {
	  cl <- class(from)
	  nCl <- sub("^n", "l", cl)
	  r <- new(nCl)# default => no validity check; and copy slots:
	  ## result is "same", for sparse just with an 'x' slot
	  if(extends(cl, "sparseMatrix"))
	      slot(r, "x") <- rep.int(TRUE, nnzero(from))
	  for(nm in slotNames(cl))
	      slot(r, nm) <- slot(from, nm)
	  r
      })

setAs("dMatrix", "lMatrix",
      function(from) {
	  r <- new(class2(class(from), "l"))# default => no validity
	  r@x <- as.logical(from@x)
	  sNams <- slotNames(r)
	  for(nm in sNams[sNams != "x"])
	      slot(r, nm) <- slot(from, nm)
	  r
      })

setAs("lMatrix", "dMatrix",
      function(from) {
	  r <- new(sub("^l", "d", class(from)))
	  r@x <- as.double(from@x)
	  sNams <- slotNames(r)
	  for(nm in sNams[sNams != "x"])
	      slot(r, nm) <- slot(from, nm)
	  r
      })

## needed at least for lsparse* :
setAs("lMatrix", "dgCMatrix",
      function(from) as(as(from, "lgCMatrix"), "dgCMatrix"))

## all() methods ---> ldenseMatrix.R and lsparseMatrix.R

setMethod("any", signature(x = "lMatrix"),
	  function(x, ..., na.rm = FALSE)
	  ## logical unit-triangular has TRUE diagonal:
	  (prod(dim(x)) >= 1 && is(x, "triangularMatrix") && x@diag == "U") ||
	  any(x@x, ..., na.rm = na.rm))
