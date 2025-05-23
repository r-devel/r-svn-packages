"update.bigqr" <-
function(bigQR, X, y, w=NULL,
                       singcheck=FALSE, add.intercept=FALSE){
  if (NCOL(X)+add.intercept!=length(bigQR$D))
    stop("Wrong number of columns")
  if (length(y)!=NROW(X))
    stop("Wrong number of rows")
  if (is.null(w)) w<-rep(1.0, length(y))
  if (length(y)!=length(w))
    stop("`weights' has wrong length")
  bigQR<-.Call("updateQR",X, y, w, bigQR, add.intercept)
  
  if (singcheck)
    bigQR<-.Call("singcheckQR",bigQR);

  bigQR
}

