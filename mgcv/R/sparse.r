## (c) Simon N. Wood 2011
## functions for sparse smoothing.



nearest <- function(k,X,get.a=FALSE,balanced=FALSE,cut.off=5) {
## The rows of X contain coordinates of points.
## For each point, this routine finds its k nearest 
## neighbours, returning a list of 2, n by k matrices:
## ni - ith row indexes the rows of X containing 
##      the k nearest neighbours of X[i,]
## dist - ith row is the distances to the k nearest 
##        neighbours.
## a - area associated with each point, if get.a is TRUE 
## ties are broken arbitrarily.
## if balanced = TRUE, then some nearest neighbours are sacrificed 
## for neighbours chosen to be on either side of the box in each 
## direction in this case k>2*ncol(X). These neighbours are only used
## if closer than cut.off*max(k nearest distances).
  require(mgcv)
  Xu <- uniquecombs(X);ind <- attr(Xu,"index") ## Xu[ind,] == X
  n <- nrow(Xu)
  d <- ncol(Xu)
  dist <- matrix(0,n,k)
  if (get.a) a <- 1:n else a=1
  if (balanced) {
    kk <- k - 2*ncol(X)
    if (kk<0) stop("k too small for balanced neighbours")
    oo <- .C(C_kba_nn,Xu=as.double(Xu),dist=as.double(dist),a=as.double(a),ni=as.integer(dist),
                    n=as.integer(n),d=as.integer(d),k=as.integer(kk),get.a=as.integer(get.a),
                    as.double(cut.off))
  } else {
    oo <- .C(C_k_nn,Xu=as.double(Xu),dist=as.double(dist),a=as.double(a),ni=as.integer(dist),
                    n=as.integer(n),d=as.integer(d),k=as.integer(k),get.a=as.integer(get.a))
  }
  dist <- matrix(oo$dist,n,k)[ind,]
  rind <- 1:n
  rind[ind] <- 1:n
  ni <- matrix(rind[oo$ni+1],n,k)[ind,]
  list(ni=ni,dist=dist,a=oo$a[ind])
}


sparse.pen <- function(X,area.weight=TRUE) {
## Function to generate a sparse approximate TPS penalty matrix
## NOTE: duplicates not yet handled!!
  require(Matrix)
  n <- nrow(X)
  d <- ncol(X)
  if (d!=2) stop("only 2D case available so far")
  m <- 3
  D <- matrix(0,n,6*m)
  k <- 5
  ni <- matrix(0,n,k)
 
  ## Get the sqrt penalty matrix entries...

  oo <- .C(C_sparse_penalty,as.double(X),as.integer(n),as.integer(d),D=as.double(D),
           ni=as.integer(ni),as.integer(k),as.integer(m),as.integer(area.weight))

  ## Now make put the entries into sparse matrices...

  ## First sort out indexing...
  ##ni <- cbind(1:n,matrix(oo$ni,n,k)) ## neighbour indices, including self
  ii <- rep(1:n,k+1) ## row index
  jj <- c(1:n,oo$ni+1) ## col index
  
  ni <- length(ii)
  Kx <- sparseMatrix(i=ii,j=jj,x=oo$D[1:ni],dims=c(n,n))
  Kz <- sparseMatrix(i=ii,j=jj,x=oo$D[1:ni+ni],dims=c(n,n))
  Kxz <- sparseMatrix(i=ii,j=jj,x=oo$D[1:ni+2*ni],dims=c(n,n))
  list(Kx=Kx,Kz=Kz,Kxz=Kxz)
}