## (c) Simon N. Wood (2013, 2014) mvn model extended family. 
## Released under GPL2 ...

mvn <- function() { 
## Extended family object for multivariate normal additive model.
 
  env <- new.env(parent = .GlobalEnv)
  validmu <- function(mu) all(is.finite(mu))

   
  aic <- function(y, mu, theta=NULL, wt, dev) {
    
  }
    

  ## initialization has to add in the extra parameters of 
  ## the cov matrix...
  
    preinitialize <- expression({
    ## code to evaluate in estimate.gam...
      ydim <- ncol(G$y) ## dimension of response
      nbeta <- ncole(G$X)
      ntheta <- ydim*(ydim+1)/2 ## numer of cov matrix factor params
      lpi <- attr(G$X,"lpi")
      G$X <- cbind(G$X,matrix(0,nrow(G$X),ntheta)) ## add dummy columns to G$X
      attr(G$X,"lpi") <- lpi
      G$family.data <- list(ydim = ydim,nbeta=nbeta)
      ## now get initial parameters
      for (k in 1:ydim) {
        
        magic(G$y[,k],G$X[,lpi[[k]]],sp,S,off)
      }
    })
    
    postproc <- expression({
    ## code to evaluate in estimate.gam, to do with estimated factor of
    ## precision matrix, etc...
      ydim <- object$family.data$ydim
      object$family.data$R <- matrix(0,ydim,ydim)
      ind <- object$family.data$nbeta + 1:(ydim*(ydim+1)/2);
      theta <- object$coefficients[ind]
      k <- 1;for (i in 1:ydim) for (j in i:ydim) {
        if (i==j) R[i,j] <- exp(theta[k]) else R[i,j] <- theta[k]
        k <- k + 1
      } 
      ##object$fitted.values <-
      ## object$null.deviance <-
    })
    
    initialize <- expression({
      ## Ideally fit separate models to each component and
      ## extract initial coefs, s.p.s and variances this way
    })


    residuals <- function(object,type=c("deviance","martingale")) {
      type <- match.arg(type)
      w <- object$prior.weights;log.s <- log(object$fitted.values)
      res <- w + log.s ## martingale residuals
      if (type=="deviance") res <- sign(res)*sqrt(-2*(res + w * log(-log.s)))
      res 
    }


    predict <- function(family,se=FALSE,eta=NULL,y=NULL,
               X=NULL,beta=NULL,off=NULL,Vb=NULL,family.data=NULL) {
      ## prediction function.
      ii <- order(y,decreasing=TRUE) ## C code expects non-increasing
      n <- nrow(X)
      oo <- .C("coxpred",as.double(X[ii,]),t=as.double(y[ii]),as.double(beta),as.double(Vb),
                a=as.double(family.data$a),h=as.double(family.data$h),q=as.double(family.data$q),
                tr = as.double(family.data$tr),
                n=as.integer(n),p=as.integer(ncol(X)),nt = as.integer(family.data$nt),
                s=as.double(rep(0,n)),se=as.double(rep(0,n)),PACKAGE="mgcv")
      s <- sef <- oo$s
      s[ii] <- oo$s
      sef[ii] <- oo$se    
      if (se) return(list(fit=s,se.fit=sef)) else return(list(fit=s))
    }

    rd <- qf <- NULL ## these functions currently undefined for Cox PH

    ll <- function(y,X,coef,wt,family,deriv=0,d1b=0,d2b=0,Hp=NULL,rank=0,fh=NULL,D=NULL) {
    ## function defining the cox model log lik.
    ## Calls C code "coxlpl"
    ## deriv codes: 0 - eval; 1 - grad and Hessian
    ##              2 - d1H (diagonal only)
    ##              3 - d1H; 4 d2H (diag)
    ## Hp is the preconditioned penalized hessian of the log lik
    ##    which is of rank 'rank'.
    ## fh is a factorization of Hp - either its eigen decomposition 
    ##    or its Choleski factor
    ## D is the diagonal pre-conditioning matrix used to obtain Hp
    ##   if Hr is the raw Hp then Hp = D*t(D*Hr)
      ##tr <- sort(unique(y),decreasing=TRUE)
      tr <- unique(y)
      r <- match(y,tr)
      p <- ncol(X)
      deriv <- deriv - 1
      mu <- X%*%coef
      g <- rep(0,p);H <- rep(0,p*p)
      if (deriv > 0) {
        M <- ncol(d1b)
        d1H <- if (deriv==1) rep(0,p*M) else rep(0,p*p*M)
      } else M <- d1H <- 0
      if (deriv > 2) {
        d2H <- rep(0,p*M*(M+1)/2)
        #X <- t(forwardsolve(t(L),t(X)))
        #d1b <- L %*% d1b; d2b <- L %*% d2b
        if (is.list(fh)) {
          ev <- fh
        } else  { ## need to compute eigen here
          ev <- eigen(Hp,symmetric=TRUE)
          if (rank < p) ev$values[(rank+1):p] <- 0
        } 
        X <- X%*%(ev$vectors*D)
        d1b <- t(ev$vectors)%*%(d1b/D); d2b <- t(ev$vectors)%*%(d2b/D)
      } else trHid2H <- d2H <- 0
      ## note that the following call can not use .C(C_coxlpl,...) since the ll
      ## function is not in the mgcv namespace.
      oo <- .C("coxlpl",as.double(mu),as.double(X),as.integer(r),as.integer(wt),
            as.double(tr),n=as.integer(length(y)),p=as.integer(p),nt=as.integer(length(tr)),
            lp=as.double(0),g=as.double(g),H=as.double(H),
            d1b=as.double(d1b),d1H=as.double(d1H),d2b=as.double(d2b),d2H=as.double(d2H),
            n.sp=as.integer(M),deriv=as.integer(deriv),PACKAGE="mgcv");
      if (deriv==1) d1H <- matrix(oo$d1H,p,M) else
      if (deriv>1) {
        ind <- 1:(p^2)
        d1H <- list()
        for (i in 1:M) { 
          d1H[[i]] <- matrix(oo$d1H[ind],p,p)
          ind <- ind + p^2
        }
      } 
      if (deriv>2) { 
        d2H <- matrix(oo$d2H,p,M*(M+1)/2)
        #trHid2H <- colSums(d2H)
        d <- ev$values
        d[d>0] <- 1/d[d>0];d[d<=0] <- 0
        trHid2H <- colSums(d2H*d)
      }
      assign(".log.partial.likelihood", oo$lp, envir=environment(sys.function()))
      list(l=oo$lp,lb=oo$g,lbb=matrix(oo$H,p,p),d1H=d1H,d2H=d2H,trHid2H=trHid2H)
    }

    # environment(dev.resids) <- environment(aic) <- environment(getTheta) <- 
    # environment(rd)<- environment(qf)<- environment(variance) <- environment(putTheta) 
    environment(aic) <- environment(ll) <- env
    structure(list(family = "Cox PH", link = linktemp, linkfun = stats$linkfun,
        linkinv = stats$linkinv, ll=ll,
        aic = aic, mu.eta = stats$mu.eta, 
        initialize = initialize,preinitialize=preinitialize,postproc=postproc,
        hazard=hazard,predict=predict,residuals=residuals,
        validmu = validmu, valideta = stats$valideta, 
        rd=rd,qf=qf,drop.intercept = TRUE,
        ls=1 ## signal ls not needed
        ),
        class = c("general.family","extended.family","family"))
} ## cox.ph

