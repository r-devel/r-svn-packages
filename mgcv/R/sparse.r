## (c) Simon N. Wood 2011
## functions for sparse smoothing.



## Efficient stable full rank cubic spline routines, based on    
## deHoog and Hutchinson, 1987 and Hutchinson and deHoog,
## 1985....

setup.spline <- function(x,w=rep(1,length(x)),lambda=1,tol=1e-9) { 
## setup a cubic smoothing spline given data locations in x. 
## ties will be treated by removing duplicate x's, and averaging corresponding
## y's. Averaging is \sum_i y_i w_i^2 / \sum_i w_i^2, and the weight
## assigned to this average is then w_a^2 = \sum_i w_i^2... 
## spline object has to record duplication information, as well as 
## rotations defining spline.  
   n <- length(x)
   ind <- order(x)
   x <- x[ind] ## sort x
   w <- w[ind]
   U <- V <- rep(0,4*n)
   diagA <- rep(0,n)
   lb <- rep(0,2*n)
   oo <- .C(C_sspl_construct,as.double(lambda),x=as.double(x),w=as.double(w),U=as.double(U),V=as.double(V),
                          diagA=as.double(diagA),lb=as.double(lb),n=as.integer(n),tol=as.double(tol))
   list(trA=sum(oo$diagA), ## trace of influence matrix
        U=oo$U,V=oo$V,     ## spline defining Givens rotations
        lb=oo$lb,          ## final lower band
        x=x,               ## original x sequence, ordered
        ind=ind,           ## x0 <- x; x0[ind] <- x, puts original ordering in x0 
        w=w,               ## original weights
        ns=oo$n,           ## number of unique x values (maximum spline rank)
        tol=tol)           ## tolerance used to judge tied x values
}

apply.spline <- function(spl,y) {
## Use cubic spline object spl, from setup.spline, to smooth data in y.
  if (is.matrix(y)) { 
    m <- ncol(y) 
    y <- y[spl$ind,] ## order as x
  } else {
    m <- 1
    y <- y[spl$ind]
  }
  n <- length(spl$x)
  oo <- .C(C_sspl_mapply,f = as.double(y),x=as.double(spl$x),as.double(spl$w),
            U=as.double(spl$U),as.double(spl$V),n=as.integer(spl$ns),
            nf=as.integer(n),tol=as.double(spl$tol),m=as.integer(m))
  if (is.matrix(y)) {
    y <- matrix(oo$f,n,m)
    y[spl$ind,] <- y ## original order
  } else {
    y[spl$ind] <- oo$f
  }
  y
}



## Routines for sparse thin plate splines...

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

tieMatrix <- function(x) {
## takes matrix x, and produces sparse matrix P that maps list of unique 
## rows to full set. Matrix of unique rows is returned in xu.
## If a smoothing penalty matrix, S, is set up based on rows of xu,
## then P%*%solve(t(P)%*%P + S,t(P)) is hat matrix. 
  x <- as.matrix(x)
  n <- nrow(x)
  xu <- uniquecombs(x)
  if (nrow(xu)==nrow(x)) return(NULL)
  ind <- attr(xu,"index")
  x <- as.matrix(x)
  n <- nrow(x)
  P <- sparseMatrix(i=1:n,j=ind,x=rep(1,n),dims=c(n,nrow(xu)))
  return(list(P=P,xu=xu))
}


## sparse smooths must be initialized with...
## 1. a set of variable names, a blocking factor and a type.  


spasm.construct.spatps <- function(object,data) {
## entry object inherits from "spatps" & contains:
## * terms, the name of 2 arguments of the smooth
## * block, the name of a blocking factor. Can be NULL.
## return object also has...
## * nblock - number of blocks.
## * ind, list where ind[[i]] gives rows to which block i applies.
## * S, a list where S[[i]] is the ith penalty coefficient matrix.
## * Ri, an empty list to contain cholski factors later.
  if (is.null(object$area.weight)) object$area.weight <- FALSE
  dat <- list()
  d <- length(object$terms)
  if (d != 2) stop("sparse tps only deals with 2D data so far") 
  for (i in 1:length(object$terms)) 
    dat[[object$term[i]]] <- get.var(object$term[i],data)
  ind <- list()
  n <- length(dat[[1]])
  ## if there is a blocking factor then set up indexes 
  ## indicating which data go with which block...
  if (!is.null(object$block)) {
    block <- as.factor(get.var(object$block,data))
    nb <- length(levels(block))
    for (i in 1:nb) { 
      ind[[i]] <- (1:n)[block==levels(block)[i]]
    }
  } else { ## all one block
    nb <- 1
    ind[[1]] <- 1:n
  }
  object$nblock <- nb
  object$ind <- ind
  ## so ind[[i]] indexes the elements operated on by the ith smoother.
  object$P <- object$S <- list()
  for (i in 1:nb) {
    X <- cbind(dat[[1]][ind[[i]]],dat[[2]][ind[[i]]])
    dup <- tieMatrix(X) ## strip out any dumplicates...
    if (is.null(dup)) {     
      um <- sparse.pen(X,object$area.weight)
    } else {
      object$P[[i]] <- dup$P ## maps this block's data to unique set
      um <- sparse.pen(dup$xu,object$area.weight)
    }
    object$S[[i]] <- t(um$Kx)%*%um$Kx + t(um$Kz)%*%um$Kz + 2* t(um$Kxz)%*%um$Kxz
  }
  object$Ri <- list()
  class(object) <- "spatps"
  object
}

spasm.sp.spatps <- function(object,sp,get.trH=FALSE) { 
## Set up smooth, given new smoothing parameter.
## In particular, construct sparse choleski factor of inverse 
## smoother matrix, for each block. Optionally returns tr(H) and
## always log|H|. 
  if (is.null(object$S)) stop("object not fully initialized")
  trH <- ldetH <- 0
  for (i in 1:object$nblock) {
    n <- length(object$ind[[i]])
    if (length(object$P) < i || is.null(object$P[[i]])) {
      Hi <- Diagonal(n) + object$S[[i]] * sp 
    } else {
      Hi <- t(object$P[[i]])%*%object$P[[i]] + object$S[[i]] * sp 
    }
    object$Ri[[i]] <- chol(Hi,pivot=TRUE)
    ## next line is wrong and doesn't deal with P
    ldetH <- ldetH - 2*sum(log(diag(object$Ri[[i]]))) ## WRONG det! need |I-H|
    if (get.trH) {
      if (length(object$P) < i || is.null(object$P[[i]])) {
        R <- solve(object$Ri[[i]]) ## most costly part of routine
      } else {
        R <- solve(object$Ri[[i]],t(object$P[[i]]))
      }
      trH <- trH + sum(R^2)
    }
  }
  if (get.trH) object$trH <- trH
  object$ldetH <- ldetH 
  object
}

spasm.smooth.spatps <- function(object,X,residual=FALSE,block=0) {
## apply smooth, or its residual operation to X.
## if block == 0 then apply whole thing, otherwise X must have the correct 
## number of rows for the smooth block.
  if (block>0) { 
    piv <- attr(object$Ri[[block]],"pivot")
    n <- length(object$ind[[block]])
    ipiv <- piv
    ipiv[piv] <- 1:n
    if (length(object$P) < block || is.null(object$P[[block]])) {
      P <- Diagonal(n)
    } else { P <- object$P[[block]] }
    if (is.matrix(X)) {
      X1 <- solve(t(object$Ri[[block]]),t(P)%*%X[piv,])
      X1 <- solve(object$Ri[[block]],X1)[ipiv,] 
      if (residual) X <- X - P%*%X1 else X <- P%*%X1;
    } else {
      X1 <- solve(t(object$Ri[[block]]),t(P)%*%X[piv])
      X1 <- solve(objectRi[[block]],X1)[ipiv]
      if (residual) X <- X - as.numeric(P%*%X1) else X <- as.numeric(P%*%X1);
    }
  } else for (i in 1:object$nblock) { ## work through all blocks
    piv <- attr(object$Ri[[i]],"pivot")
    n <- length(object$ind[[i]])
    ipiv <- piv
    ipiv[piv] <- 1:n
    if (length(object$P) < i || is.null(object$P[[i]])) 
    if (is.matrix(X)) {
      Xi <- t(P)%*%X[object$ind[[i]],]
      X1 <- solve(t(object$Ri[[i]]),Xi[piv,])
      X1 <- solve(object$Ri[[i]],X1)[ipiv,]
      if (residual) X1 <- X[object$ind[[i]],] - P%*%X1
      X[object$ind[[i]],] <- P%*%X1 
    } else {
      Xi <- t(P)%*%X[object$ind[[i]]]
      X1 <- solve(t(object$Ri[[i]]),Xi[piv])
      X1 <- solve(object$Ri[[i]],X1)[ipiv]
      if (residual) X1 <- X[object$ind[[i]]] - as.numeric(P%*%X1)
      X[object$ind[[i]]] <- as.numeric(P%*%X1) 
    }    
  }
  X
} 



## generics for sparse smooth classes...

spasm.construct <- function(object,data,knots) UseMethod("spasm.construct")
spasm.sp <- function(object,sp) UseMethod("spasm.sp")




