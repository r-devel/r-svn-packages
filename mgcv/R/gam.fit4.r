
## (c) Simon N. Wood (2013). Provided under GPL 2.
## Routines for gam estimation beyond exponential family.

dDeta <- function(y,mu,wt,theta,fam,deriv=0) {
## What is available directly from the family are derivatives of the 
## deviance and link w.r.t. mu. This routine converts these to the
## required derivatives of the deviance w.r.t. eta.
## deriv is the order of derivative of the smoothing parameter score 
## required.
   r <- fam$Dd(y, mu, theta, wt, level=deriv)  
   d <- list(Deta=0,Dth=0,Dth2=0,Deta2=0,EDeta2=0,Detath=0,
             Deta3=0,Deta2th=0,Detath2=0,
             Deta4=0,Deta3th=0,Deta2th2=0)
   if (fam$link=="identity") { ## don't waste time on transformation
      d$Deta <- r$Dmu;d$Deta2 <- r$Dmu2
      d$EDeta2 <- r$EDmu2
      if (deriv>0) {
        d$Dth <- r$Dth; d$Detath <- r$Dmuth
        d$Deta3 <- r$Dmu3; d$Deta2th <- r$Dmu2th
      }
      if (deriv>1) {
        d$Deta4 <- r$Dmu4; d$Dth2 <- r$Dth2; d$Detath2 <- r$Dmuth2
        d$Deta2th2 <- r$Dmu2th2; d$Deta3th <- r$Dmu3th
      }
      return(d)
   }

   ig1 <- fam$mu.eta(fam$linkfun(mu)) 
   g2 <- fam$d2link(mu)
   g22 <- g2^2
   ig12 <- ig1^2;ig13 <- ig12 * ig1

   d$Deta <- r$Dmu * ig1
   d$Deta2 <- (r$Dmu2 - r$Dmu*g2*ig1)*ig12
   d$EDeta2 <- r$EDmu2*ig12
   if (deriv>0) {
      d$Dth <- r$Dth 
      d$Detath <- r$Dmuth * ig1
      g3 <- fam$d3link(mu)
      d$Deta3 <- (r$Dmu3 - 3*r$Dmu2 * g2 * ig1 + r$Dmu * (3*g22*ig12 - g3 * ig1))*ig13
      d$Deta2th <- (r$Dmu2th - r$Dmuth*g2*ig1)*ig12
   }
   if (deriv>1) {
     g4 <- fam$d4link(mu)
     d$Deta4 <- ig12^2*(r$Dmu4 - 6*r$Dmu3*ig1*g2 + r$Dmu2*(15*g22*ig12-4*g3*ig1) - 
                       r$Dmu*(15*g2^3*ig13-10*g2*g3*ig12  +g4*ig1))
     d$Dth2 <- r$Dth2
     d$Detath2 <- r$Dmuth2 * ig1 
     d$Deta2th2 <- ig12*(r$Dmu2th2 - r$Dmuth2*g2*ig1)
     d$Deta3th <-  ig13*(r$Dmu3th - 3 *r$Dmu2th*g2*ig1 + r$Dmuth*(3*g22*ig12-g3*ig1))
   }
   d
} ## dDmu

fetad.test <- function(y,mu,wt,theta,fam,eps = 1e-7) {
## test family derivatives w.r.t. eta
  
  dd <- dDeta(y,mu,wt,theta,fam,deriv=2)
  dev <- fam$dev.resids(y, mu, wt,theta)
  mu1 <- fam$linkinv(fam$linkfun(mu)+eps)
  dev1 <- fam$dev.resids(y,mu1, wt,theta)
  Deta.fd <- (dev1-dev)/eps
  cat("Deta: rdiff = ",range(dd$Deta-Deta.fd)," cor = ",cor(dd$Deta,Deta.fd),"\n")
  for (i in 1:length(theta)) {
    th1 <- theta;th1[i] <- th1[i] + eps
    dev1 <- fam$dev.resids(y, mu, wt,th1)
    Dth.fd <- (dev1-dev)/eps
    cat("Dth[",i,"]: rdiff = ",range(dd$Dth-Dth.fd)," cor = ",cor(dd$Dth,Dth.fd),"\n")
  }
  ## second order up...
  dd1 <- dDeta(y,mu1,wt,theta,fam,deriv=2)
  Deta2.fd <- (dd1$Deta - dd$Deta)/eps
  cat("Deta2: rdiff = ",range(dd$Deta2-Deta2.fd)," cor = ",cor(dd$Deta2,Deta2.fd),"\n")
  Deta3.fd <- (dd1$Deta2 - dd$Deta2)/eps
  cat("Deta3: rdiff = ",range(dd$Deta3-Deta3.fd)," cor = ",cor(dd$Deta3,Deta3.fd),"\n")
  Deta4.fd <- (dd1$Deta3 - dd$Deta3)/eps
  cat("Deta4: rdiff = ",range(dd$Deta4-Deta4.fd)," cor = ",cor(dd$Deta4,Deta4.fd),"\n")
  ## and now the higher derivs wrt theta...
  for (i in 1:length(theta)) {
    th1 <- theta;th1[i] <- th1[i] + eps
    dd1 <- dDeta(y,mu,wt,th1,fam,deriv=2)
    Detath.fd <- (dd1$Deta - dd$Deta)/eps
    cat("Detath[",i,"]: rdiff = ",range(dd$Detath-Detath.fd)," cor = ",cor(dd$Detath,Detath.fd),"\n")
    Deta2th.fd <- (dd1$Deta2 - dd$Deta2)/eps
    cat("Deta2th[",i,"]: rdiff = ",range(dd$Deta2th-Deta2th.fd)," cor = ",cor(dd$Deta2th,Deta2th.fd),"\n")
    Deta3th.fd <- (dd1$Deta3 - dd$Deta3)/eps
    cat("Deta3th[",i,"]: rdiff = ",range(dd$Deta3th-Deta3th.fd)," cor = ",cor(dd$Deta3th,Deta3th.fd),"\n")
    ## now the 3 second derivative w.r.t. theta terms
    Dth2.fd <- (dd1$Dth - dd$Dth)/eps
    cat("Dth2[",i,",]: rdiff = ",range(dd$Dth2-Dth2.fd)," cor = ",cor(dd$Dth2,Dth2.fd),"\n")
    Detath2.fd <- (dd1$Detath - dd$Detath)/eps
    cat("Detath2[",i,",]: rdiff = ",range(dd$Detath2-Detath2.fd)," cor = ",cor(dd$Detath2,Detath2.fd),"\n")
    Deta2th2.fd <- (dd1$Deta2th - dd$Deta2th)/eps
    cat("Deta2th2[",i,",]: rdiff = ",range(dd$Deta2th2-Deta2th2.fd)," cor = ",cor(dd$Deta2th2,Deta2th2.fd),"\n")
  }
} ## fetad.test

fmud.test <- function(y,mu,wt,theta,fam,eps = 1e-7) {
## test family deviance derivatives w.r.t. mu
  dd <- fam$Dd(y, mu, theta, wt, level=2) 
  dev <- fam$dev.resids(y, mu, wt,theta)
  dev1 <- fam$dev.resids(y, mu+eps, wt,theta)
  Dmu.fd <- (dev1-dev)/eps
  cat("Dmu: rdiff = ",range(dd$Dmu-Dmu.fd)," cor = ",cor(dd$Dmu,Dmu.fd),"\n")
  for (i in 1:length(theta)) {
    th1 <- theta;th1[i] <- th1[i] + eps
    dev1 <- fam$dev.resids(y, mu, wt,th1)
    Dth.fd <- (dev1-dev)/eps
    cat("Dth[",i,"]: rdiff = ",range(dd$Dth-Dth.fd)," cor = ",cor(dd$Dth,Dth.fd),"\n")
  }
  ## second order up...
  dd1 <- fam$Dd(y, mu+eps, theta, wt, level=2)
  Dmu2.fd <- (dd1$Dmu - dd$Dmu)/eps
  cat("Dmu2: rdiff = ",range(dd$Dmu2-Dmu2.fd)," cor = ",cor(dd$Dmu2,Dmu2.fd),"\n")
  Dmu3.fd <- (dd1$Dmu2 - dd$Dmu2)/eps
  cat("Dmu3: rdiff = ",range(dd$Dmu3-Dmu3.fd)," cor = ",cor(dd$Dmu3,Dmu3.fd),"\n")
  Dmu4.fd <- (dd1$Dmu3 - dd$Dmu3)/eps
  cat("Dmu4: rdiff = ",range(dd$Dmu4-Dmu4.fd)," cor = ",cor(dd$Dmu4,Dmu4.fd),"\n")
  ## and now the higher derivs wrt theta 
  for (i in 1:length(theta)) {
    th1 <- theta;th1[i] <- th1[i] + eps
    dd1 <- fam$Dd(y, mu, th1, wt, level=2)
    Dmuth.fd <- (dd1$Dmu - dd$Dmu)/eps
    cat("Dmuth[",i,"]: rdiff = ",range(dd$Dmuth-Dmuth.fd)," cor = ",cor(dd$Dmuth,Dmuth.fd),"\n")
    Dmu2th.fd <- (dd1$Dmu2 - dd$Dmu2)/eps
    cat("Dmu2th[",i,"]: rdiff = ",range(dd$Dmu2th-Dmu2th.fd)," cor = ",cor(dd$Dmu2th,Dmu2th.fd),"\n")
    Dmu3th.fd <- (dd1$Dmu3 - dd$Dmu3)/eps
    cat("Dmu3th[",i,"]: rdiff = ",range(dd$Dmu3th-Dmu3th.fd)," cor = ",cor(dd$Dmu3th,Dmu3th.fd),"\n")
    ## now the 3 second derivative w.r.t. theta terms
    Dth2.fd <- (dd1$Dth - dd$Dth)/eps
    cat("Dth2[",i,",]: rdiff = ",range(dd$Dth2-Dth2.fd)," cor = ",cor(dd$Dth2,Dth2.fd),"\n")
    Dmuth2.fd <- (dd1$Dmuth - dd$Dmuth)/eps
    cat("Dmuth2[",i,",]: rdiff = ",range(dd$Dmuth2-Dmuth2.fd)," cor = ",cor(dd$Dmuth2,Dmuth2.fd),"\n")
    Dmu2th2.fd <- (dd1$Dmu2th - dd$Dmu2th)/eps
    cat("Dmu2th2[",i,",]: rdiff = ",range(dd$Dmu2th2-Dmu2th2.fd)," cor = ",cor(dd$Dmu2th2,Dmu2th2.fd),"\n")
  }
}




gam.fit4 <- function(x, y, sp, Eb,UrS=list(),
            weights = rep(1, nobs), start = NULL, etastart = NULL, 
            mustart = NULL, offset = rep(0, nobs),U1=diag(ncol(x)), Mp=-1, family = gaussian(), 
            control = gam.control(), deriv=2,
            scale=1,scoreType="REML",null.coef=rep(0,ncol(x)),...) {
## Routine for fitting GAMs beyond exponential family.
## Inputs as gam.fit3 except that family is of class "extended.family", while
## sp contains the vector of extended family parameters, followed by the log smoothing parameters,
## followed by the log scale parameter if scale < 0

  if (family$n.theta>0) { 
    ind <- 1:family$n.theta
    theta <- sp[ind] ## parameters of the family
    family$putTheta(theta)
    sp <- sp[-ind]   ## log smoothing parameters
  }

  if (scale>0) scale.known <- TRUE else {
    ## unknown scale parameter, trial value supplied as 
    ## final element of sp. 
    scale.known <- FALSE
    nsp <- length(sp)
    scale <- exp(sp[nsp])
    sp <- sp[-nsp]
  }
  
  x <- as.matrix(x)  
  nSp <- length(sp) 
  rank.tol <- .Machine$double.eps*100 ## tolerance to use for rank deficiency
  q <- ncol(x)
  nobs <- nrow(x)  
  
  xnames <- dimnames(x)[[2]]
  ynames <- if (is.matrix(y)) rownames(y) else names(y)
  ## Now a stable re-parameterization is needed....

  if (length(UrS)) {
      rp <- mgcv:::gam.reparam(UrS,sp,deriv)
      T <- diag(q)
      T[1:ncol(rp$Qs),1:ncol(rp$Qs)] <- rp$Qs
      T <- U1%*%T ## new params b'=T'b old params
    
      null.coef <- t(T)%*%null.coef  
     
      ## form x%*%T in parallel 
      x <- .Call("mgcv_pmmult2",x,T,0,0,control$nthreads,PACKAGE="mgcv")
      ## x <- .Call(C_mgcv_pmmult2,x,T,0,0,control$nthreads) ## within package version
      ## x <- x%*%T   ## model matrix 0(nq^2)
      rS <- list()
      for (i in 1:length(UrS)) {
        rS[[i]] <- rbind(rp$rS[[i]],matrix(0,Mp,ncol(rp$rS[[i]])))
      } ## square roots of penalty matrices in current parameterization
      Eb <- Eb%*%T ## balanced penalty matrix
      rows.E <- q-Mp
      Sr <- cbind(rp$E,matrix(0,nrow(rp$E),Mp))
      St <- rbind(cbind(rp$S,matrix(0,nrow(rp$S),Mp)),matrix(0,Mp,q))
  } else { 
      T <- diag(q); 
      St <- matrix(0,q,q) 
      rSncol <- sp <- rows.E <- Eb <- Sr <- 0   
      rS <- list(0)
      rp <- list(det=0,det1 = rep(0,0),det2 = rep(0,0),fixed.penalty=FALSE)
  }

  ## re-parameterization complete. Initialization....

  nvars <- ncol(x)
  if (nvars==0) stop("emtpy models not available")
  if (is.null(weights)) weights <- rep.int(1, nobs)
  if (is.null(offset)) offset <- rep.int(0, nobs)

  ## call the families initialization code...

  if (is.null(mustart)) {
    eval(family$initialize)
  } else {
    mukeep <- mustart
    eval(family$initialize)
    mustart <- mukeep
  }
  
  ## and now finalize initialization of mu and eta...

  coefold <- NULL
  eta <- if (!is.null(etastart)) etastart
         else if (!is.null(start)) 
              if (length(start) != nvars) 
                  stop("Length of start should equal ", nvars, 
                  " and correspond to initial coefs for ", deparse(xnames))
              else {
                  coefold <- start
                  offset + as.vector(if (NCOL(x) == 1) 
                  x * start
                  else x %*% start)
              }
              else family$linkfun(mustart)
 

   mu.eta <- family$mu.eta
   Dd <- family$Dd
   d2link <- family$d2link
   linkinv <- family$linkinv
   valideta <- family$valideta
   validmu <- family$validmu
   dev.resids <- family$dev.resids
 
   mu <- linkinv(eta)
     
   ## need an initial `null deviance' to test for initial divergence...
   null.eta <- as.numeric(x%*%null.coef + as.numeric(offset))
   old.pdev <- sum(dev.resids(y, linkinv(null.eta), weights,theta)) + t(null.coef)%*%St%*%null.coef 
   conv <-  boundary <- FALSE
 
   for (iter in 1:control$maxit) { ## start of main fitting iteration 
      dd <- dDeta(y,mu,weights,theta,family,0) ## derivatives of deviance w.r.t. eta
      good <- dd$Deta2 != 0
      w <- dd$Deta2[good] * .5
      z <- eta[good] - .5 * dd$Deta[good] / w
      oo <- .C("pls_fit1",   ##C_pls_fit1, reinstate for use in mgcv
               y=as.double(z),X=as.double(x[good,]),w=as.double(w),
                     E=as.double(Sr),Es=as.double(Eb),n=as.integer(sum(good)),
                     q=as.integer(ncol(x)),rE=as.integer(rows.E),eta=as.double(z),
                     penalty=as.double(1),rank.tol=as.double(rank.tol),
                     nt=as.integer(control$nthreads),PACKAGE="mgcv") ## delete PACKAGE in mgcv
      if (oo$n<0) { ## then problem is indefinite - switch to +ve weights for this step
        good <- dd$Deta2 > 0
        w <- dd$Deta2[good] * .5
        z <- eta[good] - .5 * dd$Deta[good] / w
        oo <- .C("pls_fit1", ##C_pls_fit1,
                  y=as.double(z),X=as.double(x[good,]),w=as.double(w),
                     E=as.double(Sr),Es=as.double(Eb),n=as.integer(sum(good)),
                     q=as.integer(ncol(x)),rE=as.integer(rows.E),eta=as.double(z),
                     penalty=as.double(1),rank.tol=as.double(rank.tol),
                     nt=as.integer(control$nthreads),PACKAGE="mgcv")
      }
      start <- oo$y[1:ncol(x)] ## current coefficient estimates
      penalty <- oo$penalty ## size of penalty
      eta <- drop(x%*%start) ## the linear predictor

      if (any(!is.finite(start))) { ## test for breakdown
          conv <- FALSE
          warning("Non-finite coefficients at iteration ", 
                  iter)
          break
      }        
     
      mu <- linkinv(eta <- eta + offset)
      dev <- sum(dev.resids(y, mu, weights,theta))

      ## now step halve under non-finite deviance...
      if (!is.finite(dev)) {
         if (is.null(coefold)) {
            if (is.null(null.coef)) 
              stop("no valid set of coefficients has been found:please supply starting values", 
                    call. = FALSE)
            ## Try to find feasible coefficients from the null.coef and null.eta
            coefold <- null.coef
            etaold <- null.eta
         }
         warning("Step size truncated due to divergence", 
                     call. = FALSE)
         ii <- 1
         while (!is.finite(dev)) {
               if (ii > control$maxit) 
                    stop("inner loop 1; can't correct step size")
               ii <- ii + 1
               start <- (start + coefold)/2
               eta <- (eta + etaold)/2               
               mu <- linkinv(eta)
               dev <- sum(dev.resids(y, mu, weights,theta))
         }
         boundary <- TRUE
         if (control$trace) 
                  cat("Step halved: new deviance =", dev, "\n")
      } ## end of infinite deviance correction

      ## now step halve if mu or eta are out of bounds... 
      if (!(valideta(eta) && validmu(mu))) {
         warning("Step size truncated: out of bounds", 
                  call. = FALSE)
         ii <- 1
         while (!(valideta(eta) && validmu(mu))) {
                  if (ii > control$maxit) 
                    stop("inner loop 2; can't correct step size")
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  eta <- (eta + etaold)/2 
                  mu <- linkinv(eta)
         }
         boundary <- TRUE
         dev <- sum(dev.resids(y, mu, weights))
         if (control$trace) 
                  cat("Step halved: new deviance =", dev, "\n")
      } ## end of invalid mu/eta handling

      ## now check for divergence of penalized deviance....
  
      pdev <- dev + penalty  ## the penalized deviance 
      if (control$trace) cat("penalized deviance =", pdev, "\n")
     
      div.thresh <- 10*(.1+abs(old.pdev))*.Machine$double.eps^.5

      if (pdev-old.pdev>div.thresh) { ## solution diverging
         ii <- 1 ## step halving counter
         if (iter==1) { ## immediate divergence, need to shrink towards zero 
               etaold <- null.eta; coefold <- null.coef
         }
         while (pdev -old.pdev > div.thresh)  { ## step halve until pdev <= old.pdev
           if (ii > 100) 
              stop("inner loop 3; can't correct step size")
           ii <- ii + 1
           start <- (start + coefold)/2 
           eta <- (eta + etaold)/2               
           mu <- linkinv(eta)
           dev <- sum(dev.resids(y, mu, weights,theta))
           pdev <- dev + t(start)%*%St%*%start ## the penalized deviance
           if (control$trace) 
                  cat("Step halved: new penalized deviance =", pdev, "\n")
        }
     } ## end of pdev divergence

     ## convergence testing...

     if (abs(pdev - old.pdev)/(0.1 + abs(pdev)) < control$epsilon) {
       if (max(abs(start-coefold))>control$epsilon*max(abs(start+coefold))/2) {
         old.pdev <- pdev  ## not converged quite enough
         coef <- coefold <- start
         etaold <- eta 
         muold <- mu
       } else { ## converged
         conv <- TRUE
         coef <- start
         break 
       }
     } else { ## not converged
       old.pdev <- pdev
       coef <- coefold <- start
       etaold <- eta 
     }
   } ## end of main loop
   
   ## so at this stage the model has been fully estimated
   coef <- as.numeric(T %*% coef)
   #coef ## final parameters

   ## now obtain derivatives, if these are needed...
   dd <- dDeta(y,mu,weights,theta,family,deriv)
   good <- dd$Deta2 != 0
   w <- dd$Deta2[good] * .5
   z <- eta[good] - .5 * dd$Deta[good] / w
   wf <- dd$EDeta2[good] * .5 ## Fisher type weights 

   residuals <- rep.int(NA, nobs)
   residuals[good] <- z - (eta - offset)[good]

   ntot <- length(theta) + length(sp)
   if (deriv>1) n2d <- ntot*(1+ntot)/2 else n2d <- 0 
   rSncol <- unlist(lapply(UrS,ncol))
   oo <- .C("gdi2", ## C_gdi2,
            X=as.double(x[good,]),E=as.double(Sr),Es=as.double(Eb),rS=as.double(unlist(rS)),
            U1 = as.double(U1),sp=as.double(exp(sp)),theta=as.double(theta),
            z=as.double(z),w=as.double(w),wf=as.double(wf),Dth=as.double(dd$Dth),Det=as.double(dd$Deta),
            Det2=as.double(dd$Deta2),Dth2=as.double(dd$Dth2),Det.th=as.double(dd$Detath),
            Det2.th=as.double(dd$Deta2th),Det3=as.double(dd$Deta3),Det.th2 = as.double(dd$Detath2),
            Det4 = as.double(dd$Deta4),Det3.th=as.double(dd$Deta3th), Deta2.th2=as.double(dd$Deta2th2),
            beta=as.double(coef),D1=as.double(rep(0,ntot)),D2=as.double(rep(0,ntot^2)),
            P=as.double(0),P1=as.double(rep(0,ntot)),P2 = as.double(rep(0,ntot^2)),
            ldet=as.double(0),ldet1 = as.double(rep(0,ntot)), ldet2 = as.double(rep(0,ntot^2)),
            rV=as.double(rep(0,ncol(x)^2)),
            rank.tol=as.double(.Machine$double.eps^.75),rank.est=as.integer(0),
	    n=as.integer(sum(good)),q=as.integer(ncol(x)),M=as.integer(length(UrS)),
            n.theta=as.integer(length(theta)), Mp=as.integer(Mp),Enrow=as.integer(rows.E),
            rSncol=as.integer(rSncol),deriv=as.integer(deriv),
	    fixed.penalty = as.integer(rp$fixed.penalty),nt=as.integer(control$nthreads),PACKAGE="mgcv")

   rV <- matrix(oo$rV,ncol(x),ncol(x)) ## rV%*%t(rV)*scale gives covariance matrix 
   rV <- T %*% rV   
   Kmat <- matrix(0,nrow(x),ncol(x)) 
   Kmat[good,] <- oo$X                    ## rV%*%t(K)%*%(sqrt(wf)*X) = F; diag(F) is edf array 

   D2 <- matrix(oo$D2,ntot,ntot); ldet2 <- matrix(oo$ldet2,ntot,ntot)
   bSb2 <- matrix(oo$P2,ntot,ntot)
   ## compute the REML score...
   ls <- family$ls(y,weights,n,theta,scale)
   nt <- length(theta)
   lsth1 <- ls$lsth1[1:nt];
   lsth2 <- as.matrix(ls$lsth2)[1:nt,1:nt] ## exclude any derivs w.r.t log scale here
   REML <- (dev+oo$P)/(2*scale) - ls$ls + (oo$ldet - rp$det)/2 - Mp * log(2*pi*scale)/2
   REML1 <- REML2 <- NULL
   if (deriv) {
     ind <- 1:nSp + length(theta)
     det1 <- oo$ldet1;det1[ind] <- det1[ind] - rp$det1
     REML1 <- (oo$D1+oo$P1)/(2*scale) - c(lsth1,rep(0,length(sp))) + (det1)/2
     if (deriv>1) {
       ls2 <- D2*0;ls2[1:nt,1:nt] <- lsth2 
       ldet2[ind,ind] <- ldet2[ind,ind] - rp$det2
       REML2 <- (D2+bSb2)/(2*scale) - ls2 + ldet2/2
     }
   } 

   if (!scale.known&&deriv) { ## need derivatives wrt log scale, too 
      Dp <- dev + oo$P
      dlr.dlphi <- -Dp/(2 *scale) - ls$lsth1[nt+1] - Mp/2
      d2lr.d2lphi <- Dp/(2*scale) - ls$lsth2[nt+1,nt+1] 
      d2lr.dspphi <- -(oo$D1+oo$P1)/(2*scale) 
      d2lr.dspphi[1:nt] <- d2lr.dspphi[1:nt] - ls$lsth2[nt+1,1:nt]
      REML1 <- c(REML1,dlr.dlphi)
      if (deriv==2) {
              REML2 <- rbind(REML2,as.numeric(d2lr.dspphi))
              REML2 <- cbind(REML2,c(as.numeric(d2lr.dspphi),d2lr.d2lphi))
      }
   }

   ## for the moment return the deviance and its derivatives, only...
   names(coef) <- xnames
   names(residuals) <- ynames
   wtdmu <- sum(weights * y)/sum(weights)
   nulldev <- sum(dev.resids(y, wtdmu, weights))
   n.ok <- nobs - sum(weights == 0)
   nulldf <- n.ok
   wt <- rep.int(0, nobs)
   wt[good] <- wf 

   aic.model <- family$aic(y, mu, theta, weights, dev) # note: incomplete 2*edf needs to be added
 
  ## fitted values might not just be mu....
  ## if (!is.null(family$fv)) mu <- family$fv(mu,theta) 
  ## .... actually, probably do not want to return fv values here

   list(coefficients = coef,residuals=residuals,fitted.values = mu,
        family=family, linear.predictors = eta,deviance=dev,
        null.deviance=nulldev,iterr=iter,
        weights=wf, ## note that these are Fisher type weights 
        prior.weights=weights,
        df.null = nulldf, y = y, converged = conv,
        boundary = boundary,
        REML=REML,REML1=REML1,REML2=REML2,
        rV=rV,
        scale.est=scale,reml.scale=scale,
        aic=aic.model,
        rank=oo$rank.est,
        K=Kmat,control=control
        ,D1=oo$D1,D2=D2,
        ldet=oo$ldet,ldet1=oo$ldet1,ldet2=ldet2,
        bSb=oo$P,bSb1=oo$P1,bSb2=bSb2,
        ls=ls$ls,ls1=ls$lsth1,ls2=ls$lsth2
       )
 
} ## gam.fit4






