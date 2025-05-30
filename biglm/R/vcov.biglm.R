"vcov.biglm" <-
function(object,...){
   if(!object$qr$checked)
      object$qr<-singcheck.bigqr(object$qr)
   p<-length(object$qr$D)
   R<-diag(p)
   R[row(R)>col(R)]<-object$qr$rbar
   R<-t(R)
   R<-sqrt(object$qr$D)*R
   ok<-object$qr$D!=0
   R[ok,ok]<-chol2inv(R[ok,ok])
   R[!ok,]<-R[,!ok]<-NA
   dimnames(R)<-list(object$names, object$names)
   
   if(!is.null(object$sandwich)){
     xyqr<-singcheck.bigqr(object$sandwich$xy)
     rxy<-diag(p*(p+1))
     rxy[row(rxy)>col(rxy)]<-xyqr$rbar
     rxy<-t(rxy)
     rxy<-sqrt(xyqr$D)*rxy
     M<-t(rxy)%*%rxy
     
     beta<-coef(object)
     beta[!ok]<-0
     bbeta<-kronecker(diag(p),beta)
     ##FIXME: singularities in beta
     Vcenter<-M[1:p,1:p] + t(bbeta)%*%M[-(1:p),-(1:p)]%*%bbeta -
       t(bbeta)%*%M[-(1:p),1:p] - M[1:p,-(1:p)]%*%bbeta
     
     V<-matrix(NA,p,p)
     V[ok,ok]<-R[ok,ok]%*%Vcenter[ok,ok]%*%R[ok,ok]
     dimnames(V)<-list(object$names, object$names)
     attr(V,"model-based")<-R*object$qr$ss/(object$n-p+sum(!ok))
   } else {
     V<-R*object$qr$ss/(object$n-p+sum(!ok))
   }

   V
}

