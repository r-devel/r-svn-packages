## R routines for gam fitting with calculation of derivatives w.r.t. sp.s
## (c) Simon Wood 2004-2008

## This routine is for type 3 gam fitting. The basic idea is that a P-IRLS
## is run to convergence, and only then is a scheme for evaluating the 
## derivatives iterated to convergence. The advantage is that many key
## quantities are fixed at this stage, including the key decompositions
## In addition the R side work is simplified considerably.The routine
## evaluates first and second derivatives of the deviance and tr(A).


gam.fit2 <- function (x, y, sp, S=list(),rS=list(),off, H=NULL, 
            weights = rep(1, nobs), start = NULL, etastart = NULL, 
            mustart = NULL, offset = rep(0, nobs), family = gaussian(), 
            control = gam.control(), intercept = TRUE,deriv=2,use.svd=TRUE,
            gamma=1,scale=1,printWarn=TRUE,scoreType="REML",...) 
## deriv, sp, S, rS, H added to arg list. 
## need to modify family before call.
## This is essentially a backup of the version of gam.fit3 implementing a Wood (2008) based
## on derivative iteration, but with Newton extensions. This version will fail if 
## newton weights become negative, but could easily be modified to revert to fisher only.
{   if (family$link==family$canonical) fisher <- TRUE else fisher=FALSE ## Newton = Fisher, but Fisher cheaper!
    if (scale>0) scale.known <- TRUE else scale.known <- FALSE
    scale <- abs(scale)
    if (!deriv%in%c(0,1,2)) stop("unsupported order of differentiation requested of gam.fit3")
    x <- as.matrix(x)
    iter <- 0;coef <- rep(0,ncol(x))
    xnames <- dimnames(x)[[2]]
    ynames <- if (is.matrix(y)) 
        rownames(y)
    else names(y)
    conv <- FALSE
    n <- nobs <- NROW(y) ## n is just to keep codetools happy
    nvars <- ncol(x)
    EMPTY <- nvars == 0
    if (is.null(weights)) 
        weights <- rep.int(1, nobs)
    if (is.null(offset)) 
        offset <- rep.int(0, nobs)
    variance <- family$variance
    dev.resids <- family$dev.resids
    aic <- family$aic
    linkinv <- family$linkinv
    mu.eta <- family$mu.eta
    if (!is.function(variance) || !is.function(linkinv)) 
        stop("illegal `family' argument")
    valideta <- family$valideta
    if (is.null(valideta)) 
        valideta <- function(eta) TRUE
    validmu <- family$validmu
    if (is.null(validmu)) 
        validmu <- function(mu) TRUE
    if (is.null(mustart)) {
        eval(family$initialize)
    }
    else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
    }

    ## Added code
    if (family$family=="gaussian"&&family$link=="identity") strictly.additive <- TRUE else
      strictly.additive <- FALSE
    nSp <- length(S)
    if (nSp==0) deriv <- FALSE 
    St <- totalPenalty(S,H,off,sp,ncol(x))
    Sr <- mroot(St)

    ## end of added code

    D1 <- D2 <- P <- P1 <- P2 <- trA <- trA1 <- trA2 <- 
        GCV<- GCV1<- GCV2<- GACV<- GACV1<- GACV2<- UBRE <-
        UBRE1<- UBRE2<- REML<- REML1<- REML2 <-NULL

    if (EMPTY) {
        eta <- rep.int(0, nobs) + offset
        if (!valideta(eta)) 
            stop("Invalid linear predictor values in empty model")
        mu <- linkinv(eta)
        if (!validmu(mu)) 
            stop("Invalid fitted means in empty model")
        dev <- sum(dev.resids(y, mu, weights))
        w <- ((weights * mu.eta(eta)^2)/variance(mu))^0.5
        residuals <- (y - mu)/mu.eta(eta)
        good <- rep(TRUE, length(residuals))
        boundary <- conv <- TRUE
        coef <- numeric(0)
        iter <- 0
        V <- variance(mu)
        alpha <- dev
        trA2 <- trA1 <- trA <- 0
        if (deriv) GCV2 <- GCV1<- UBRE2 <- UBRE1<-trA1 <- rep(0,nSp)
        GCV <- nobs*alpha/(nobs-gamma*trA)^2
        UBRE <- alpha/nobs - scale + 2*gamma/n*trA
        scale.est <- alpha / (nobs - trA)
    } ### end if (EMPTY)
    else {
        coefold <- NULL
        eta <- if (!is.null(etastart)) 
            etastart
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
        etaold <- eta
        muold <- mu <- linkinv(eta)
        if (!(validmu(mu) && valideta(eta))) 
            stop("Can't find valid starting values: please specify some")
    
        boundary <- conv <- FALSE
        rV=matrix(0,ncol(x),ncol(x))   
        old.pdev <- 0     
        for (iter in 1:control$maxit) {
            good <- weights > 0
            varmu <- variance(mu)[good]
            if (any(is.na(varmu))) 
                stop("NAs in V(mu)")
            if (any(varmu == 0)) 
                stop("0s in V(mu)")
            mu.eta.val <- mu.eta(eta)
            if (any(is.na(mu.eta.val[good]))) 
                stop("NAs in d(mu)/d(eta)")
            good <- (weights > 0) & (mu.eta.val != 0)
            if (all(!good)) {
                conv <- FALSE
                warning("No observations informative at iteration ", 
                  iter)
                break
            }
            mevg<-mu.eta.val[good];mug<-mu[good];yg<-y[good]
            weg<-weights[good];var.mug<-variance(mug)
            if (fisher) { ## Conventional Fisher scoring
              z <- (eta - offset)[good] + (yg - mug)/mevg
              w <- sqrt((weg * mevg^2)/var.mug)
            } else { ## full Newton
              c <- yg - mug
              e <- mevg*(1 + c*(family$dvar(mug)/mevg+var.mug*family$d2link(mug))*mevg/var.mug)
              z <- (eta - offset)[good] + c/e ## offset subtracted as eta = X%*%beta + offset
              w <- sqrt(weg*e*mevg/var.mug)
              ## correct `good', `w' and `z' to remove any zero weights
              if (sum(w==0)) {
                wf <- weights*0
                wf[good] <- w
                good <- (wf!=0)&good
                ind <- w!=0
                w <- w[ind];z <- z[ind]
              }
            }

            ## Here a Fortran call has been replaced by update.beta call
           
            if (sum(good)<ncol(x)) stop("Not enough informative observations.")
     
            oo<-.C(C_pls_fit,y=as.double(z),as.double(x[good,]),as.double(w),as.double(Sr),as.integer(sum(good)),
            as.integer(ncol(x)),as.integer(ncol(Sr)),eta=as.double(z),penalty=as.double(1),
            as.double(.Machine$double.eps*100))
       
            start <- oo$y[1:ncol(x)];
            penalty <- oo$penalty
            eta <- drop(x%*%start)

            if (any(!is.finite(start))) {
                conv <- FALSE
                warning("Non-finite coefficients at iteration ", 
                  iter)
                break
            }        
     
           mu <- linkinv(eta <- eta + offset)
           dev <- sum(dev.resids(y, mu, weights))
          
           if (control$trace) 
                cat("Deviance =", dev, "Iterations -", iter, 
                  "\n")
            boundary <- FALSE
            
            if (!is.finite(dev)) {
                if (is.null(coefold)) 
                  stop("no valid set of coefficients has been found:please supply starting values", 
                    call. = FALSE)
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
                  dev <- sum(dev.resids(y, mu, weights))
                }
                boundary <- TRUE
                if (control$trace) 
                  cat("Step halved: new deviance =", dev, "\n")
            }
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
            }

            pdev <- dev + penalty  ## the penalized deviance 

            if (control$trace) 
                  cat("penalized deviance =", pdev, "\n")

            div.thresh <- 10*(.1+abs(old.pdev))*.Machine$double.eps^.5 
            ## ... threshold for judging divergence --- too tight and near
            ## perfect convergence can cause a failure here

            if (iter>1&&(pdev-old.pdev>div.thresh)) { ## solution diverging
             ii <- 1 ## step halving counter
             while (pdev -old.pdev > div.thresh)  
             { ## step halve until pdev <= old.pdev
                if (ii > 200) 
                   stop("inner loop 3; can't correct step size")
                ii <- ii + 1
                start <- (start + coefold)/2 
                eta <- (eta + etaold)/2               
                mu <- linkinv(eta)
                  dev <- sum(dev.resids(y, mu, weights))
                  pdev <- dev + t(start)%*%St%*%start ## the penalized deviance
                if (control$trace) 
                  cat("Step halved: new penalized deviance =", pdev, "\n")
              }
            } 
            
            if (strictly.additive) { conv <- TRUE;coef <- start;break;}

            if (abs(pdev - old.pdev)/(0.1 + abs(pdev)) < control$epsilon) {
               ## if (max(abs(start-coefold))>control$epsilon*max(abs(start+coefold))/2) {
                if (max(abs(mu-muold))>control$epsilon*max(abs(mu+muold))/2) {
                  old.pdev <- pdev
                  coef <- coefold <- start
                  etaold <- eta 
                  muold <- mu
                } else {
                  conv <- TRUE
                  coef <- start
                  break 
                }
            }
            else {  old.pdev <- pdev
                coef <- coefold <- start
                etaold <- eta 
            }
        } ### end main loop 
       
        dev <- sum(dev.resids(y, mu, weights)) 
       
        ## Now call the derivative calculation scheme. This requires the
        ## following inputs:
        ## z and w - the pseudodata and weights
        ## X the model matrix and E where EE'=S
        ## rS the single penalty square roots
        ## sp the log smoothing parameters
        ## y and mu the data and model expected values
        ## g1,g2,g3 - the first 3 derivatives of g(mu) wrt mu
        ## V,V1,V2 - V(mu) and its first two derivatives wrt mu
        ## on output it returns the gradient and hessian for
        ## the deviance and trA 

         good <- weights > 0
         varmu <- variance(mu)[good]
         if (any(is.na(varmu))) stop("NAs in V(mu)")
         if (any(varmu == 0)) stop("0s in V(mu)")
         mu.eta.val <- mu.eta(eta)
         if (any(is.na(mu.eta.val[good]))) 
                stop("NAs in d(mu)/d(eta)")
         good <- (weights > 0) & (mu.eta.val != 0)
   
         mevg <- mu.eta.val[good];mug <- mu[good];yg <- y[good]
         weg <- weights[good];etag <- eta[good]
         var.mug<-variance(mug)

         if (fisher) { ## Conventional Fisher scoring
              
              z <- (eta - offset)[good] + (yg - mug)/mevg
            
              w <- sqrt((weg * mevg^2)/var.mug)
         } else { ## full Newton
          
              c <- yg - mug
              e <- mevg*(1 + c*(family$dvar(mug)/mevg+var.mug*family$d2link(mug))*mevg/var.mug)
              z <- (eta - offset)[good] + c/e ## offset subtracted as eta = X%*%beta + offset
              w <- sqrt(weg*e*mevg/var.mug)

              ## correct `good', `w' and `z' to remove any zero weights
              if (sum(w==0)) {
                wf <- weights*0
                wf[good] <- w
                good <- (wf!=0)&good
                ind <- w!=0
                w <- w[ind];z <- z[ind]
                mevg <- mu.eta.val[good];mug <- mu[good];yg <- y[good]
                weg <- weights[good];etag <- eta[good]
                var.mug<-variance(mug)
              }
         
         }
        
         g1 <- 1/mevg
         g2 <- family$d2link(mug)
         g3 <- family$d3link(mug)

         V <- family$variance(mug)
         V1 <- family$dvar(mug)
         V2 <- family$d2var(mug)      
        
         if (fisher) {
           g4 <- V3 <- 0
         } else {
           g4 <- family$d4link(mug)
           V3 <- family$d3var(mug)
         }

         if (TRUE) { ### TEST CODE for derivative ratio based versions of code... 
           g2 <- g2/g1;g3 <- g3/g1;g4 <- g4/g1
           V1 <- V1/V;V2 <- V2/V;V3 <- V3/V
         }

         P1 <- D1 <- array(0,nSp);P2 <- D2 <- matrix(0,nSp,nSp) # for derivs of deviance/ Pearson
         trA1 <- array(0,nSp);trA2 <- matrix(0,nSp,nSp) # for derivs of tr(A)
         rV=matrix(0,ncol(x),ncol(x));
         dum <- 1
         if (control$trace) cat("calling gdi...")

       oo <- .C(C_gdi,X=as.double(x[good,]),E=as.double(Sr),rS = as.double(unlist(rS)),
           sp=as.double(exp(sp)),z=as.double(z),w=as.double(w),mu=as.double(mug),eta=as.double(etag),y=as.double(yg),
           p.weights=as.double(weg),g1=as.double(g1),g2=as.double(g2),g3=as.double(g3),g4=as.double(g4),V0=as.double(V),
           V1=as.double(V1),V2=as.double(V2),V3=as.double(V3),beta=as.double(coef),D1=as.double(D1),D2=as.double(D2),
           P=as.double(dum),P1=as.double(P1),P2=as.double(P2),trA=as.double(dum),
           trA1=as.double(trA1),trA2=as.double(trA2),rV=as.double(rV),rank.tol=as.double(.Machine$double.eps*100),
           conv.tol=as.double(control$epsilon),rank.est=as.integer(1),n=as.integer(length(z)),
           p=as.integer(ncol(x)),M=as.integer(nSp),Encol = as.integer(ncol(Sr)),
           rSncol=as.integer(unlist(lapply(rS,ncol))),deriv=as.integer(deriv),use.svd=as.integer(use.svd),
           REML = as.integer(scoreType=="REML"),fisher=as.integer(fisher))      
       
         if (control$trace) cat("done! (iteration took ",oo$deriv," steps)\n")
 
         rV <- matrix(oo$rV,ncol(x),ncol(x))
         coef <- oo$beta;
         trA <- oo$trA;
         scale.est <- dev/(nobs-trA)

        if (scoreType=="REML") {
          if (scale.known) { ## use Fisher-Laplace REML
             ls <- family$ls(y,weights,n,scale) ## saturated likelihood and derivatives
             REML <- (dev + oo$conv.tol)/(2*scale) - ls[1] + oo$rank.tol/2
             if (deriv) {
               REML1 <- oo$D1/(2*scale) + oo$trA1/2
               if (deriv==2) REML2 <- (matrix(oo$D2,nSp,nSp)/scale + matrix(oo$trA2,nSp,nSp))/2
               if (sum(!is.finite(REML2))) {
                 stop("Smoothing parameter derivate iteration diverging. Decrease fit tolerance! See `epsilon' in `gam.contol'")
               }
             }
           } else { ## scale unknown use Pearson-Fisher-Laplace REML
             phi <- oo$P ## REMLish scale estimate
             ls <- family$ls(y,weights,n,phi) ## saturated likelihood and derivatives
             phi1 <- oo$P1;phi2 <- matrix(oo$P2,nSp,nSp)
             Dp <- dev + oo$conv.tol
             Dp1 <- oo$D1
             Dp2 <- matrix(oo$D2,nSp,nSp)
             K <- oo$rank.tol/2
             K1 <- oo$trA1/2;K2 <- matrix(oo$trA2,nSp,nSp)/2             

             REML <- Dp/(2*phi) - ls[1] + K
             if (deriv) {
               REML1 <- Dp1/(2*phi) - phi1*(Dp/(2*phi^2) + ls[2]) + K1
               if (deriv==2) REML2 <- 
                      Dp2/(2*phi) - (outer(Dp1,phi1)+outer(phi1,Dp1))/(2*phi^2) +
                      (Dp/phi^3 - ls[3])*outer(phi1,phi1) -
                      (Dp/(2*phi^2)+ls[2])*phi2 + K2
             }
           } 
         } else { ## Not REML ....

           P <- oo$P
           
           delta <- nobs - gamma * trA
           delta.2 <- delta*delta           
  
           GCV <- nobs*dev/delta.2
           GACV <- dev/nobs + P * 2*gamma*trA/(delta * nobs) 

           UBRE <- dev/nobs - 2*delta*scale/nobs + scale
        
           if (deriv) {
             trA1 <- oo$trA1
           
             D1 <- oo$D1
             P1 <- oo$P1
          
             if (sum(!is.finite(D1))||sum(!is.finite(P1))||sum(!is.finite(trA1))) { 
                 stop("Smoothing parameter derivate iteration diverging. Decrease fit tolerance! See `epsilon' in `gam.contol'")}
         
             delta.3 <- delta*delta.2
  
             GCV1 <- nobs*D1/delta.2 + 2*nobs*dev*trA1*gamma/delta.3
             GACV1 <- D1/nobs + 2*P/delta.2 * trA1 + 2*gamma*trA*P1/(delta*nobs)

             UBRE1 <- D1/nobs + gamma * trA1 *2*scale/nobs
             if (deriv==2) {
               trA2 <- matrix(oo$trA2,nSp,nSp) 
               D2 <- matrix(oo$D2,nSp,nSp)
               P2 <- matrix(oo$P2,nSp,nSp)
              
               if (sum(!is.finite(D2))||sum(!is.finite(P2))||sum(!is.finite(trA2))) { 
                 stop("Smoothing parameter derivate iteration diverging. Decrease fit tolerance! See `epsilon' in `gam.contol'")}
             
               GCV2 <- outer(trA1,D1)
               GCV2 <- (GCV2 + t(GCV2))*gamma*2*nobs/delta.3 +
                      6*nobs*dev*outer(trA1,trA1)*gamma*gamma/(delta.2*delta.2) + 
                      nobs*D2/delta.2 + 2*nobs*dev*gamma*trA2/delta.3  
               GACV2 <- D2/nobs + outer(trA1,trA1)*4*P/(delta.3) +
                      2 * P * trA2 / delta.2 + 2 * outer(trA1,P1)/delta.2 +
                      2 * outer(P1,trA1) *(1/(delta * nobs) + trA/(nobs*delta.2)) +
                      2 * trA * P2 /(delta * nobs) 
               GACV2 <- (GACV2 + t(GACV2))*.5
               UBRE2 <- D2/nobs +2*gamma * trA2 * scale / nobs
             } ## end if (deriv==2)
           } ## end if (deriv)
        } ## end !REML
        # end of inserted code
        if (!conv&&printWarn) 
            warning("Algorithm did not converge")
        if (printWarn&&boundary) 
            warning("Algorithm stopped at boundary value")
        eps <- 10 * .Machine$double.eps
        if (printWarn&&family$family[1] == "binomial") {
            if (any(mu > 1 - eps) || any(mu < eps)) 
                warning("fitted probabilities numerically 0 or 1 occurred")
        }
        if (printWarn&&family$family[1] == "poisson") {
            if (any(mu < eps)) 
                warning("fitted rates numerically 0 occurred")
        }
 
        residuals <- rep.int(NA, nobs)
        residuals[good] <- z - (eta - offset)[good]
          
        names(coef) <- xnames 
    } ### end if (!EMPTY)
    names(residuals) <- ynames
    names(mu) <- ynames
    names(eta) <- ynames
    wt <- rep.int(0, nobs)
    wt[good] <- w^2
    names(wt) <- ynames
    names(weights) <- ynames
    names(y) <- ynames
   
    wtdmu <- if (intercept) 
        sum(weights * y)/sum(weights)
    else linkinv(offset)
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    n.ok <- nobs - sum(weights == 0)
    nulldf <- n.ok - as.integer(intercept)
   
    aic.model <- aic(y, n, mu, weights, dev) # note: incomplete 2*edf needs to be added

    list(coefficients = coef, residuals = residuals, fitted.values = mu, 
         family = family, linear.predictors = eta, deviance = dev, 
        null.deviance = nulldev, iter = iter, weights = wt, prior.weights = weights, 
        df.null = nulldf, y = y, converged = conv,
        boundary = boundary,D1=D1,D2=D2,P=P,P1=P1,P2=P2,trA=trA,trA1=trA1,trA2=trA2,
        GCV=GCV,GCV1=GCV1,GCV2=GCV2,GACV=GACV,GACV1=GACV1,GACV2=GACV2,UBRE=UBRE,
        UBRE1=UBRE1,UBRE2=UBRE2,REML=REML,REML1=REML1,REML2=REML2,rV=rV,
        scale.est=scale.est,aic=aic.model,rank=oo$rank.est)
} ## end of gam.fit2


gam.fit3 <- function (x, y, sp, S=list(),rS=list(),off, H=NULL, 
            weights = rep(1, nobs), start = NULL, etastart = NULL, 
            mustart = NULL, offset = rep(0, nobs), family = gaussian(), 
            control = gam.control(), intercept = TRUE,deriv=2,use.svd=TRUE,
            gamma=1,scale=1,printWarn=TRUE,scoreType="REML",...) 
## This version is designed to allow iterative weights to be negative. This means that 
## it deals with weights, rather than sqrt weights.
## deriv, sp, S, rS, H added to arg list. 
## need to modify family before call.
{   if (family$link==family$canonical) fisher <- TRUE else fisher=FALSE ## Newton = Fisher, but Fisher cheaper!
    if (scale>0) scale.known <- TRUE else scale.known <- FALSE
    scale <- abs(scale)
    if (!deriv%in%c(0,1,2)) stop("unsupported order of differentiation requested of gam.fit3")
    x <- as.matrix(x)
    iter <- 0;coef <- rep(0,ncol(x))
    xnames <- dimnames(x)[[2]]
    ynames <- if (is.matrix(y)) 
        rownames(y)
    else names(y)
    conv <- FALSE
    n <- nobs <- NROW(y) ## n is just to keep codetools happy
    nvars <- ncol(x)
    EMPTY <- nvars == 0
    if (is.null(weights)) 
        weights <- rep.int(1, nobs)
    if (is.null(offset)) 
        offset <- rep.int(0, nobs)
    variance <- family$variance
    dev.resids <- family$dev.resids
    aic <- family$aic
    linkinv <- family$linkinv
    mu.eta <- family$mu.eta
    if (!is.function(variance) || !is.function(linkinv)) 
        stop("illegal `family' argument")
    valideta <- family$valideta
    if (is.null(valideta)) 
        valideta <- function(eta) TRUE
    validmu <- family$validmu
    if (is.null(validmu)) 
        validmu <- function(mu) TRUE
    if (is.null(mustart)) {
        eval(family$initialize)
    }
    else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
    }

    ## Added code
    if (family$family=="gaussian"&&family$link=="identity") strictly.additive <- TRUE else
      strictly.additive <- FALSE
    nSp <- length(S)
    if (nSp==0) deriv <- FALSE 
    St <- totalPenalty(S,H,off,sp,ncol(x))
    Sr <- mroot(St)

    ## end of added code

    D1 <- D2 <- P <- P1 <- P2 <- trA <- trA1 <- trA2 <- 
        GCV<- GCV1<- GCV2<- GACV<- GACV1<- GACV2<- UBRE <-
        UBRE1<- UBRE2<- REML<- REML1<- REML2 <-NULL

    if (EMPTY) {
        eta <- rep.int(0, nobs) + offset
        if (!valideta(eta)) 
            stop("Invalid linear predictor values in empty model")
        mu <- linkinv(eta)
        if (!validmu(mu)) 
            stop("Invalid fitted means in empty model")
        dev <- sum(dev.resids(y, mu, weights))
        w <- (weights * mu.eta(eta)^2)/variance(mu)   ### BUG: incorrect for Newton
        residuals <- (y - mu)/mu.eta(eta)
        good <- rep(TRUE, length(residuals))
        boundary <- conv <- TRUE
        coef <- numeric(0)
        iter <- 0
        V <- variance(mu)
        alpha <- dev
        trA2 <- trA1 <- trA <- 0
        if (deriv) GCV2 <- GCV1<- UBRE2 <- UBRE1<-trA1 <- rep(0,nSp)
        GCV <- nobs*alpha/(nobs-gamma*trA)^2
        UBRE <- alpha/nobs - scale + 2*gamma/n*trA
        scale.est <- alpha / (nobs - trA)
    } ### end if (EMPTY)
    else {
        coefold <- NULL
        eta <- if (!is.null(etastart)) 
            etastart
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
        etaold <- eta
        muold <- mu <- linkinv(eta)
        if (!(validmu(mu) && valideta(eta))) 
            stop("Can't find valid starting values: please specify some")
    
        boundary <- conv <- FALSE
        rV=matrix(0,ncol(x),ncol(x))   
        old.pdev <- 0     
        for (iter in 1:control$maxit) { ## start of main fitting iteration
            good <- weights > 0
            varmu <- variance(mu)[good]
            if (any(is.na(varmu))) 
                stop("NAs in V(mu)")
            if (any(varmu == 0)) 
                stop("0s in V(mu)")
            mu.eta.val <- mu.eta(eta)
            if (any(is.na(mu.eta.val[good]))) 
                stop("NAs in d(mu)/d(eta)")
            good <- (weights > 0) & (mu.eta.val != 0)
            if (all(!good)) {
                conv <- FALSE
                warning("No observations informative at iteration ", 
                  iter)
                break
            }
            mevg<-mu.eta.val[good];mug<-mu[good];yg<-y[good]
            weg<-weights[good];var.mug<-variance(mug)
            if (fisher) { ## Conventional Fisher scoring
              z <- (eta - offset)[good] + (yg - mug)/mevg
              w <- (weg * mevg^2)/var.mug
            } else { ## full Newton
              c <- yg - mug
              alpha <- mevg*(1 + c*(family$dvar(mug)/mevg+var.mug*family$d2link(mug))*mevg/var.mug)
              z <- (eta - offset)[good] + c/alpha ## offset subtracted as eta = X%*%beta + offset
              w <- weg*alpha*mevg/var.mug
            }

            ## Here a Fortran call has been replaced by update.beta call
           
            if (sum(good)<ncol(x)) stop("Not enough informative observations.")
     
            oo<-.C(C_pls_fit,y=as.double(z),as.double(x[good,]),as.double(w),as.double(Sr),n=as.integer(sum(good)),
            as.integer(ncol(x)),as.integer(ncol(Sr)),eta=as.double(z),penalty=as.double(1),
            as.double(.Machine$double.eps*100))
       
            if (!fisher&&oo$n<0) { ## likelihood indefinite - switch to Fisher for this step
              z <- (eta - offset)[good] + (yg - mug)/mevg
              w <- (weg * mevg^2)/var.mug
              oo<-.C(C_pls_fit,y=as.double(z),as.double(x[good,]),as.double(w),as.double(Sr),n=as.integer(sum(good)),
                     as.integer(ncol(x)),as.integer(ncol(Sr)),eta=as.double(z),penalty=as.double(1),
                     as.double(.Machine$double.eps*100))
            }

            start <- oo$y[1:ncol(x)];
            penalty <- oo$penalty
            eta <- drop(x%*%start)

            if (any(!is.finite(start))) {
                conv <- FALSE
                warning("Non-finite coefficients at iteration ", 
                  iter)
                break
            }        
     
           mu <- linkinv(eta <- eta + offset)
           dev <- sum(dev.resids(y, mu, weights))
          
           if (control$trace) 
                cat("Deviance =", dev, "Iterations -", iter, 
                  "\n")
            boundary <- FALSE
            
            if (!is.finite(dev)) {
                if (is.null(coefold)) 
                  stop("no valid set of coefficients has been found:please supply starting values", 
                    call. = FALSE)
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
                  dev <- sum(dev.resids(y, mu, weights))
                }
                boundary <- TRUE
                if (control$trace) 
                  cat("Step halved: new deviance =", dev, "\n")
            }
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
            }

            pdev <- dev + penalty  ## the penalized deviance 

            if (control$trace) 
                  cat("penalized deviance =", pdev, "\n")

            div.thresh <- 10*(.1+abs(old.pdev))*.Machine$double.eps^.5 
            ## ... threshold for judging divergence --- too tight and near
            ## perfect convergence can cause a failure here

            if (iter>1&&(pdev-old.pdev>div.thresh)) { ## solution diverging
             ii <- 1 ## step halving counter
             while (pdev -old.pdev > div.thresh)  
             { ## step halve until pdev <= old.pdev
                if (ii > 200) 
                   stop("inner loop 3; can't correct step size")
                ii <- ii + 1
                start <- (start + coefold)/2 
                eta <- (eta + etaold)/2               
                mu <- linkinv(eta)
                  dev <- sum(dev.resids(y, mu, weights))
                  pdev <- dev + t(start)%*%St%*%start ## the penalized deviance
                if (control$trace) 
                  cat("Step halved: new penalized deviance =", pdev, "\n")
              }
            } 
            
            if (strictly.additive) { conv <- TRUE;coef <- start;break;}

            if (abs(pdev - old.pdev)/(0.1 + abs(pdev)) < control$epsilon) {
               ## if (max(abs(start-coefold))>control$epsilon*max(abs(start+coefold))/2) {
                if (max(abs(mu-muold))>control$epsilon*max(abs(mu+muold))/2) {
                  old.pdev <- pdev
                  coef <- coefold <- start
                  etaold <- eta 
                  muold <- mu
                } else {
                  conv <- TRUE
                  coef <- start
                  break 
                }
            }
            else {  old.pdev <- pdev
                coef <- coefold <- start
                etaold <- eta 
            }
        } ### end main loop 
       
        dev <- sum(dev.resids(y, mu, weights)) 
       
        ## Now call the derivative calculation scheme. This requires the
        ## following inputs:
        ## z and w - the pseudodata and weights
        ## X the model matrix and E where EE'=S
        ## rS the single penalty square roots
        ## sp the log smoothing parameters
        ## y and mu the data and model expected values
        ## g1,g2,g3 - the first 3 derivatives of g(mu) wrt mu
        ## V,V1,V2 - V(mu) and its first two derivatives wrt mu
        ## on output it returns the gradient and hessian for
        ## the deviance and trA 

         good <- weights > 0
         varmu <- variance(mu)[good]
         if (any(is.na(varmu))) stop("NAs in V(mu)")
         if (any(varmu == 0)) stop("0s in V(mu)")
         mu.eta.val <- mu.eta(eta)
         if (any(is.na(mu.eta.val[good]))) 
                stop("NAs in d(mu)/d(eta)")
         good <- (weights > 0) & (mu.eta.val != 0)
   
         mevg <- mu.eta.val[good];mug <- mu[good];yg <- y[good]
         weg <- weights[good];etag <- eta[good]
         var.mug<-variance(mug)

         if (fisher) { ## Conventional Fisher scoring
              z <- (eta - offset)[good] + (yg - mug)/mevg
              w <- (weg * mevg^2)/var.mug
         } else { ## full Newton
              c <- yg - mug
              alpha <- mevg*(1 + c*(family$dvar(mug)/mevg+var.mug*family$d2link(mug))*mevg/var.mug)
              z <- (eta - offset)[good] + c/alpha ## offset subtracted as eta = X%*%beta + offset
              w <- weg*alpha*mevg/var.mug
         }
        
         g1 <- 1/mevg
         g2 <- family$d2link(mug)
         g3 <- family$d3link(mug)

         V <- family$variance(mug)
         V1 <- family$dvar(mug)
         V2 <- family$d2var(mug)      
        
         if (fisher) {
           g4 <- V3 <- 0
         } else {
           g4 <- family$d4link(mug)
           V3 <- family$d3var(mug)
         }

         if (TRUE) { ### TEST CODE for derivative ratio based versions of code... 
           g2 <- g2/g1;g3 <- g3/g1;g4 <- g4/g1
           V1 <- V1/V;V2 <- V2/V;V3 <- V3/V
         }

         P1 <- D1 <- array(0,nSp);P2 <- D2 <- matrix(0,nSp,nSp) # for derivs of deviance/ Pearson
         trA1 <- array(0,nSp);trA2 <- matrix(0,nSp,nSp) # for derivs of tr(A)
         rV=matrix(0,ncol(x),ncol(x));
         dum <- 1
         if (control$trace) cat("calling gdi...")

       oo <- .C(C_gdi,X=as.double(x[good,]),E=as.double(Sr),rS = as.double(unlist(rS)),
           sp=as.double(exp(sp)),z=as.double(z),w=as.double(w),mu=as.double(mug),eta=as.double(etag),y=as.double(yg),
           p.weights=as.double(weg),g1=as.double(g1),g2=as.double(g2),g3=as.double(g3),g4=as.double(g4),V0=as.double(V),
           V1=as.double(V1),V2=as.double(V2),V3=as.double(V3),beta=as.double(coef),D1=as.double(D1),D2=as.double(D2),
           P=as.double(dum),P1=as.double(P1),P2=as.double(P2),trA=as.double(dum),
           trA1=as.double(trA1),trA2=as.double(trA2),rV=as.double(rV),rank.tol=as.double(.Machine$double.eps*100),
           conv.tol=as.double(control$epsilon),rank.est=as.integer(1),n=as.integer(length(z)),
           p=as.integer(ncol(x)),M=as.integer(nSp),Encol = as.integer(ncol(Sr)),
           rSncol=as.integer(unlist(lapply(rS,ncol))),deriv=as.integer(deriv),use.svd=as.integer(use.svd),
           REML = as.integer(scoreType=="REML"),fisher=as.integer(fisher))      
       
         if (control$trace) cat("done! (iteration took ",oo$deriv," steps)\n")
 
         rV <- matrix(oo$rV,ncol(x),ncol(x))
         coef <- oo$beta;
         trA <- oo$trA;
         scale.est <- dev/(nobs-trA)

        if (scoreType=="REML") {
          if (scale.known) { ## use Fisher-Laplace REML
             ls <- family$ls(y,weights,n,scale) ## saturated likelihood and derivatives
             REML <- (dev + oo$conv.tol)/(2*scale) - ls[1] + oo$rank.tol/2
             if (deriv) {
               REML1 <- oo$D1/(2*scale) + oo$trA1/2
               if (deriv==2) REML2 <- (matrix(oo$D2,nSp,nSp)/scale + matrix(oo$trA2,nSp,nSp))/2
               if (sum(!is.finite(REML2))) {
                 stop("Non finite derivatives. Try decreasing fit tolerance! See `epsilon' in `gam.contol'")
               }
             }
           } else { ## scale unknown use Pearson-Fisher-Laplace REML
             phi <- oo$P ## REMLish scale estimate
             ls <- family$ls(y,weights,n,phi) ## saturated likelihood and derivatives
             phi1 <- oo$P1;phi2 <- matrix(oo$P2,nSp,nSp)
             Dp <- dev + oo$conv.tol
             Dp1 <- oo$D1
             Dp2 <- matrix(oo$D2,nSp,nSp)
             K <- oo$rank.tol/2
             K1 <- oo$trA1/2;K2 <- matrix(oo$trA2,nSp,nSp)/2             

             REML <- Dp/(2*phi) - ls[1] + K
             if (deriv) {
               REML1 <- Dp1/(2*phi) - phi1*(Dp/(2*phi^2) + ls[2]) + K1
               if (deriv==2) REML2 <- 
                      Dp2/(2*phi) - (outer(Dp1,phi1)+outer(phi1,Dp1))/(2*phi^2) +
                      (Dp/phi^3 - ls[3])*outer(phi1,phi1) -
                      (Dp/(2*phi^2)+ls[2])*phi2 + K2
             }
           } 
         } else { ## Not REML ....

           P <- oo$P
           
           delta <- nobs - gamma * trA
           delta.2 <- delta*delta           
  
           GCV <- nobs*dev/delta.2
           GACV <- dev/nobs + P * 2*gamma*trA/(delta * nobs) 

           UBRE <- dev/nobs - 2*delta*scale/nobs + scale
        
           if (deriv) {
             trA1 <- oo$trA1
           
             D1 <- oo$D1
             P1 <- oo$P1
          
             if (sum(!is.finite(D1))||sum(!is.finite(P1))||sum(!is.finite(trA1))) { 
                 stop("Non-finite derivatives. Try decreasing fit tolerance! See `epsilon' in `gam.contol'")}
         
             delta.3 <- delta*delta.2
  
             GCV1 <- nobs*D1/delta.2 + 2*nobs*dev*trA1*gamma/delta.3
             GACV1 <- D1/nobs + 2*P/delta.2 * trA1 + 2*gamma*trA*P1/(delta*nobs)

             UBRE1 <- D1/nobs + gamma * trA1 *2*scale/nobs
             if (deriv==2) {
               trA2 <- matrix(oo$trA2,nSp,nSp) 
               D2 <- matrix(oo$D2,nSp,nSp)
               P2 <- matrix(oo$P2,nSp,nSp)
              
               if (sum(!is.finite(D2))||sum(!is.finite(P2))||sum(!is.finite(trA2))) { 
                 stop("Non-finite derivatives. Try decreasing fit tolerance! See `epsilon' in `gam.contol'")}
             
               GCV2 <- outer(trA1,D1)
               GCV2 <- (GCV2 + t(GCV2))*gamma*2*nobs/delta.3 +
                      6*nobs*dev*outer(trA1,trA1)*gamma*gamma/(delta.2*delta.2) + 
                      nobs*D2/delta.2 + 2*nobs*dev*gamma*trA2/delta.3  
               GACV2 <- D2/nobs + outer(trA1,trA1)*4*P/(delta.3) +
                      2 * P * trA2 / delta.2 + 2 * outer(trA1,P1)/delta.2 +
                      2 * outer(P1,trA1) *(1/(delta * nobs) + trA/(nobs*delta.2)) +
                      2 * trA * P2 /(delta * nobs) 
               GACV2 <- (GACV2 + t(GACV2))*.5
               UBRE2 <- D2/nobs +2*gamma * trA2 * scale / nobs
             } ## end if (deriv==2)
           } ## end if (deriv)
        } ## end !REML
        # end of inserted code
        if (!conv&&printWarn) 
            warning("Algorithm did not converge")
        if (printWarn&&boundary) 
            warning("Algorithm stopped at boundary value")
        eps <- 10 * .Machine$double.eps
        if (printWarn&&family$family[1] == "binomial") {
            if (any(mu > 1 - eps) || any(mu < eps)) 
                warning("fitted probabilities numerically 0 or 1 occurred")
        }
        if (printWarn&&family$family[1] == "poisson") {
            if (any(mu < eps)) 
                warning("fitted rates numerically 0 occurred")
        }
 
        residuals <- rep.int(NA, nobs)
        residuals[good] <- z - (eta - offset)[good]
          
        names(coef) <- xnames 
    } ### end if (!EMPTY)
    names(residuals) <- ynames
    names(mu) <- ynames
    names(eta) <- ynames
    wt <- rep.int(0, nobs)
    wt[good] <- w
    names(wt) <- ynames
    names(weights) <- ynames
    names(y) <- ynames
   
    wtdmu <- if (intercept) 
        sum(weights * y)/sum(weights)
    else linkinv(offset)
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    n.ok <- nobs - sum(weights == 0)
    nulldf <- n.ok - as.integer(intercept)
   
    aic.model <- aic(y, n, mu, weights, dev) # note: incomplete 2*edf needs to be added

    list(coefficients = coef, residuals = residuals, fitted.values = mu, 
         family = family, linear.predictors = eta, deviance = dev, 
        null.deviance = nulldev, iter = iter, weights = wt, prior.weights = weights, 
        df.null = nulldf, y = y, converged = conv,
        boundary = boundary,D1=D1,D2=D2,P=P,P1=P1,P2=P2,trA=trA,trA1=trA1,trA2=trA2,
        GCV=GCV,GCV1=GCV1,GCV2=GCV2,GACV=GACV,GACV1=GACV1,GACV2=GACV2,UBRE=UBRE,
        UBRE1=UBRE1,UBRE2=UBRE2,REML=REML,REML1=REML1,REML2=REML2,rV=rV,
        scale.est=scale.est,aic=aic.model,rank=oo$rank.est)
}


deriv.check <- function(x, y, sp, S=list(),rS=list(),off, H=NULL, 
            weights = rep(1, length(y)), start = NULL, etastart = NULL, 
            mustart = NULL, offset = rep(0, length(y)), family = gaussian(), 
            control = gam.control(), intercept = TRUE,deriv=2,use.svd=TRUE,
            gamma=1,scale=1,printWarn=TRUE,scoreType="REML",eps=1e-7,...)
## FD checking of derivatives: basically a debugging routine
{  if (!deriv%in%c(1,2)) stop("deriv should be 1 or 2")
   if (control$epsilon>1e-9) control$epsilon <- 1e-9 
   b<-gam.fit3(x=x, y=y, sp=sp, S=S,rS=rS,off=off, H=H,
      offset = offset,family = family,weights=weights,deriv=deriv,
      control=control,gamma=gamma,scale=scale,
      printWarn=FALSE,use.svd=use.svd,mustart=mustart,scoreType=scoreType,...)

   P0 <- b$P;fd.P1 <- P10 <- b$P1;  if (deriv==2) fd.P2 <- P2 <- b$P2 
   trA0 <- b$trA;fd.gtrA <- gtrA0 <- b$trA1 ; if (deriv==2) fd.htrA <- htrA <- b$trA2 
   dev0 <- b$deviance;fd.D1 <- D10 <- b$D1 ; if (deriv==2) fd.D2 <- D2 <- b$D2 

   if (scoreType=="REML") {
     score0 <- b$REML;grad0 <- b$REML1; if (deriv==2) hess <- b$REML2 
   } else if (scoreType=="GACV") {
     score0 <- b$GACV;grad0 <- b$GACV1;if (deriv==2) hess <- b$GACV2 
   } else if (scoreType=="UBRE"){
     score0 <- b$UBRE;grad0 <- b$UBRE1;if (deriv==2) hess <- b$UBRE2 
   } else { ## default to deviance based GCV
     score0 <- b$GCV;grad0 <- b$GCV1;if (deriv==2) hess <- b$GCV2
   }
  
   fd.grad <- grad0
   if (deriv==2) fd.hess <- hess
   for (i in 1:length(sp)) {
     sp1 <- sp;sp1[i] <- sp[i]+eps
     b<-gam.fit3(x=x, y=y, sp=sp1, S=S,rS=rS,off=off, H=H,
      offset = offset,family = family,weights=weights,deriv=deriv,
      control=control,gamma=gamma,scale=scale,
      printWarn=FALSE,use.svd=use.svd,mustart=mustart,scoreType=scoreType,...)
      
      if (scoreType!="REML") {
        P <- b$P; P1 <- b$P1
        trA <- b$trA;gtrA <- b$trA1
      }
      dev <- b$deviance;D1 <- b$D1

      if (scoreType=="REML") {
        score <- b$REML;if (deriv==2) grad <- b$REML1;
      } else if (scoreType=="GACV") {
        score <- b$GACV;if (deriv==2) grad <- b$GACV1;
      } else if (scoreType=="UBRE"){
        score <- b$UBRE;if (deriv==2) grad <- b$UBRE1; 
      } else { ## default to deviance based GCV
        score <- b$GCV;if (deriv==2) grad <- b$GCV1;
      }

      if (scoreType!="REML") {
        fd.P1[i] <- (P-P0)/eps
        fd.gtrA[i] <- (trA-trA0)/eps
      }
      fd.D1[i] <- (dev - dev0)/eps
     
      fd.grad[i] <- (score-score0)/eps
      if (deriv==2) { 
        fd.hess[,i] <- (grad-grad0)/eps
        if (scoreType!="REML") {
          fd.htrA[,i] <- (gtrA-gtrA0)/eps
          fd.P2[,i] <- (P1-P10)/eps
        }
        fd.D2[,i] <- (D1-D10)/eps
       
      }
   }
   
   if (scoreType!="REML") {
     cat("\n Pearson Statistic... \n")
     cat("grad    ");print(P10)
     cat("fd.grad ");print(fd.P1)
     if (deriv==2) {
       fd.P2 <- .5*(fd.P2 + t(fd.P2))
       cat("hess\n");print(P2)
       cat("fd.hess\n");print(fd.P2)
     }

     cat("\n\n tr(A)... \n")
     cat("grad    ");print(gtrA0)
     cat("fd.grad ");print(fd.gtrA)
     if (deriv==2) {
       fd.htrA <- .5*(fd.htrA + t(fd.htrA))
       cat("hess\n");print(htrA)
       cat("fd.hess\n");print(fd.htrA)
     }
   }

   cat("\n Deviance... \n")
   cat("grad    ");print(D10)
   cat("fd.grad ");print(fd.D1)
   if (deriv==2) {
     fd.D2 <- .5*(fd.D2 + t(fd.D2))
     cat("hess\n");print(D2)
     cat("fd.hess\n");print(fd.D2)
   }



   cat("\n\n The objective...\n")

   cat("grad    ");print(grad0)
   cat("fd.grad ");print(fd.grad)
   if (deriv==2) {
     fd.hess <- .5*(fd.hess + t(fd.hess))
     cat("hess\n");print(hess)
     cat("fd.hess\n");print(fd.hess)
   }
   NULL
}


newton <- function(lsp,X,y,S,rS,off,L,lsp0,H,offset,family,weights,
                   control,gamma,scale,conv.tol=1e-6,maxNstep=5,maxSstep=2,
                   maxHalf=30,printWarn=FALSE,scoreType="deviance",
                   use.svd=TRUE,mustart = NULL,...)
## Newton optimizer for GAM gcv/aic optimization that can cope with an 
## indefinite Hessian! Main enhancements are: i) always peturbs the Hessian
## to +ve definite ii) step halves on step 
## failure, without obtaining derivatives until success; (iii) carries start
## values forward from one evaluation to next to speed convergence.    
## L is the matrix such that L%*%lsp + lsp0 gives the logs of the smoothing 
## parameters actually multiplying the S[[i]]'s
{  
  

  ## sanity check L
  if (is.null(L)) L <- diag(length(S)) else {
    if (!inherits(L,"matrix")) stop("L must be a matrix.")
    if (nrow(L)<ncol(L)) stop("L must have at least as many rows as columns.")
    if (nrow(L)!=length(S)||ncol(L)!=length(lsp)) stop("L has inconsistent dimensions.")
  }
  if (is.null(lsp0)) lsp0 <- rep(0,ncol(L))

  ## code designed to be turned on during debugging...
  check.derivs <- FALSE 
  if (check.derivs) {
     deriv <- 2
     eps <- 1e-4
     deriv.check(x=X, y=y, sp=L%*%lsp+lsp0, S=S,rS=rS,off=off, H=H,
         offset = offset,family = family,weights=weights,deriv=deriv,
         control=control,gamma=gamma,scale=scale,
         printWarn=FALSE,use.svd=use.svd,mustart=mustart,
         scoreType=scoreType,eps=eps,...)
  }
  ## ... end of debugging code 


  ## initial fit
  b<-gam.fit3(x=X, y=y, sp=L%*%lsp+lsp0, S=S,rS=rS,off=off, H=H,
     offset = offset,family = family,weights=weights,deriv=2,
     control=control,gamma=gamma,scale=scale,
     printWarn=FALSE,use.svd=use.svd,mustart=mustart,scoreType=scoreType,...)

  mustart<-b$fitted.values

  if (scoreType=="REML") {
     old.score <- score <- b$REML;grad <- b$REML1;hess <- b$REML2 
  } else if (scoreType=="GACV") {
    old.score <- score <- b$GACV;grad <- b$GACV1;hess <- b$GACV2 
  } else if (scoreType=="UBRE"){
    old.score <- score <- b$UBRE;grad <- b$UBRE1;hess <- b$UBRE2 
  } else { ## default to deviance based GCV
    old.score <- score <- b$GCV;grad <- b$GCV1;hess <- b$GCV2
  }
  
  grad <- t(L)%*%grad
  hess <- t(L)%*%hess%*%L

  Slength <- maxSstep 
  score.scale <- b$scale.est + score;    
  uconv.ind <- abs(grad) > score.scale*conv.tol
  ## check for all converged too soon, and undo !
  if (!sum(uconv.ind)) uconv.ind <- uconv.ind | TRUE
  for (i in 1:200) {
    ## exclude apparently converged gradients from computation
    hess1 <- hess[uconv.ind,uconv.ind] 
    grad1 <- grad[uconv.ind]
    ## get the trial step ...
    eh <- eigen(hess1,symmetric=TRUE)
    d <- eh$values;U <- eh$vectors
    ind <- d < 0
    d[ind] <- -d[ind] ## see Gill Murray and Wright p107/8
    d <- 1/d
    
    Nstep <- 0 * grad
    Nstep[uconv.ind] <- -drop(U%*%(d*(t(U)%*%grad1))) # (modified) Newton direction
   
    Sstep <- -Slength * grad/max(abs(grad)) # steepest descent direction 
    
    ms <- max(abs(Nstep))
    if (ms>maxNstep) Nstep <- maxNstep * Nstep/ms

    ## try the step ...
    lsp1 <- lsp + Nstep
    b<-gam.fit3(x=X, y=y, sp=L%*%lsp1+lsp0, S=S,rS=rS,off=off, H=H,
       offset = offset,family = family,weights=weights,deriv=2,
       control=control,gamma=gamma,scale=scale,
       printWarn=FALSE,mustart=mustart,use.svd=use.svd,scoreType=scoreType,...)
    
    if (scoreType=="REML") {
      score1 <- b$REML
    } else if (scoreType=="GACV") {
      score1 <- b$GACV
    } else if (scoreType=="UBRE") {
      score1 <- b$UBRE
    } else score1 <- b$GCV
    ## accept if improvement, else step halve
    ii <- 0 ## step halving counter
    if (score1<score) { ## accept
      old.score <- score 
      mustart <- b$fitted.values
      lsp <- lsp1
      if (scoreType=="REML") {
          score <- b$REML;grad <- b$REML1;hess <- b$REML2 
      } else if (scoreType=="GACV") {
          score <- b$GACV;grad <- b$GACV1;hess <- b$GACV2
      } else if (scoreType=="UBRE") {
          score <- b$UBRE;grad <- b$UBRE1;hess <- b$UBRE2 
      } else { score <- b$GCV;grad <- b$GCV1;hess <- b$GCV2} 
      grad <- t(L)%*%grad
      hess <- t(L)%*%hess%*%L
    } else { ## step halving ...
      step <- Nstep ## start with the (pseudo) Newton direction
      while (score1>score && ii < maxHalf) {
        if (ii==3) { ## Newton really not working
          step <- Sstep ## use steepest descent direction
        } else step <- step/2
        if (ii>3) Slength <- Slength/2 ## keep track of SD step length
        lsp1 <- lsp + step
        b1<-gam.fit3(x=X, y=y, sp=L%*%lsp1+lsp0, S=S,rS=rS,off=off, H=H,
           offset = offset,family = family,weights=weights,deriv=0,
           control=control,gamma=gamma,scale=scale,
           printWarn=FALSE,mustart=mustart,use.svd=use.svd,
           scoreType=scoreType,...)
         
        if (scoreType=="REML") {       
          score1 <- b1$REML
        } else if (scoreType=="GACV") {
          score1 <- b1$GACV
        } else if (scoreType=="UBRE") {
          score1 <- b1$UBRE
        } else score1 <- b1$GCV

        if (score1 <= score) { ## accept
          b<-gam.fit3(x=X, y=y, sp=L%*%lsp1+lsp0, S=S,rS=rS,off=off, H=H,
             offset = offset,family = family,weights=weights,deriv=2,
             control=control,gamma=gamma,scale=scale,
             printWarn=FALSE,mustart=mustart,use.svd=use.svd,scoreType=scoreType,...)
          mustart <- b$fitted.values
          old.score <- score;lsp <- lsp1
         
          if (scoreType=="REML") {
            score <- b$REML;grad <- b$REML1;hess <- b$REML2 
          } else if (scoreType=="GACV") {
            score <- b$GACV;grad <- b$GACV1;hess <- b$GACV2
          } else if (scoreType=="UBRE") {
            score <- b$UBRE;grad <- b$UBRE1;hess <- b$UBRE2 
          } else { score <- b$GCV;grad <- b$GCV1;hess <- b$GCV2}
          grad <- t(L)%*%grad
          hess <- t(L)%*%hess%*%L
          if (ii>3) Slength <- min(Slength*2,maxSstep) ## try increasing SD step length
        }  # end of if (score1<= score )
        ii <- ii + 1
      } # end of step halving
    }
    ## test for convergence
    converged <- TRUE
    score.scale <- b$scale.est + score;    
    uconv.ind <- abs(grad) > score.scale*conv.tol
    if (sum(uconv.ind)) converged <- FALSE
    if (abs(old.score-score)>score.scale*conv.tol) { 
      if (converged) uconv.ind <- uconv.ind | TRUE ## otherwise can't progress
      converged <- FALSE      
    }
    if (ii==maxHalf) converged <- TRUE ## step failure
    if (converged) break
  } ## end of iteration loop
  if (ii==maxHalf) ct <- "step failed"
  else if (i==200) ct <- "iteration limit reached" 
  else ct <- "full convergence"
  list(score=score,lsp=lsp,lsp.full=L%*%lsp,grad=grad,hess=hess,iter=i,conv =ct,object=b)
}

bfgs <- function(lsp,X,y,S,rS,off,L,lsp0,H,offset,family,weights,
                   control,gamma,scale,conv.tol=1e-6,maxNstep=5,maxSstep=2,
                   maxHalf=30,printWarn=FALSE,scoreType="deviance",use.svd=TRUE,
                   mustart = NULL,...)
## This optimizer is experimental... The main feature is to alternate infrequent 
## Newton steps with BFGS Quasi-Newton steps. In theory this should be faster 
## than Newton, because of the cost of full Hessian calculation, but
## in practice the extra steps required by QN tends to mean that the advantage
## is not realized...
## Newton optimizer for GAM gcv/aic optimization that can cope with an 
## indefinite Hessian, and alternates BFGS and Newton steps for speed reasons
## Main enhancements are: i) always peturbs the Hessian
## to +ve definite ii) step halves on step 
## failure, without obtaining derivatives until success; (iii) carries start
## values forward from one evaluation to next to speed convergence.    
## L is the matrix such that L%*%lsp + lsp0 gives the logs of the smoothing 
## parameters actually multiplying the S[[i]]'s
{ ## sanity check L
  if (is.null(L)) L <- diag(length(S)) else {
    if (!inherits(L,"matrix")) stop("L must be a matrix.")
    if (nrow(L)<ncol(L)) stop("L must have at least as many rows as columns.")
    if (nrow(L)!=length(S)||ncol(L)!=length(lsp)) stop("L has inconsistent dimensions.")
  }
  if (is.null(lsp0)) lsp0 <- rep(0,ncol(L))
  ## initial fit
#  ptm <- proc.time()
  b<-gam.fit3(x=X, y=y, sp=L%*%lsp+lsp0, S=S,rS=rS,off=off, H=H,
     offset = offset,family = family,weights=weights,deriv=2,
     control=control,gamma=gamma,scale=scale,
     printWarn=FALSE,use.svd=use.svd,mustart=mustart,scoreType=scoreType,...)
#  ptm <- proc.time()-ptm
#  cat("deriv=2 ",ptm,"\n")

  mustart<-b$fitted.values

  QNsteps <- floor(length(S)/2) ## how often to Newton should depend on cost...

  if (scoreType=="REML") {
     score <- b$REML;grad <- b$REML1;hess <- b$REML2 
  } else if (scoreType=="GACV") {
    old.score <- score <- b$GACV;grad <- b$GACV1;hess <- b$GACV2 
  } else if (scoreType=="UBRE"){
    old.score <- score <- b$UBRE;grad <- b$UBRE1;hess <- b$UBRE2 
  } else { ## default to deviance based GCV
    old.score <- score <- b$GCV;grad <- b$GCV1;hess <- b$GCV2
  }
  
  grad <- t(L)%*%grad
  hess <- t(L)%*%hess%*%L

  score.scale <- b$scale.est + score;    
  uconv.ind <- abs(grad) > score.scale*conv.tol
  ## check for all converged too soon, and undo !
  if (!sum(uconv.ind)) uconv.ind <- uconv.ind | TRUE
  kk <- 0 ## counter for QN steps between Newton steps
  for (i in 1:200) {
   
    if (kk==0) { ## time to reset B
      eh <- eigen(hess,symmetric=TRUE)
      d <- eh$values;U <- eh$vectors
      ind <- d < 0
      d[ind] <- -d[ind] ## see Gill Murray and Wright p107/8
      d <- 1/d
      d[d==0] <- min(d)*.Machine$double.eps^.5
      B <- U%*%(d*t(U)) ## Newton based inverse Hessian
    }
     
    kk <- kk + 1
    if (kk > QNsteps) kk <- 0 
 
    ## get the trial step ...
    
    Nstep <- 0 * grad
    Nstep[uconv.ind] <- -drop(B[uconv.ind,uconv.ind]%*%grad[uconv.ind]) # (modified) Newton direction
    
    ms <- max(abs(Nstep))
    if (ms>maxNstep) Nstep <- maxNstep * Nstep/ms

    ## try the step ...
    sc.extra <- 1e-4*sum(grad*Nstep) ## -ve sufficient decrease 
    ii <- 0 ## step halving counter
    step <- Nstep*2
    score1 <- abs(score)*2
    while (score1>score+sc.extra && ii < maxHalf) { ## reject and step halve
      ii <- ii + 1
      step <- step/2
      lsp1 <- lsp + step
  
#      ptm <- proc.time()
      if (kk!=0||ii==1) deriv <- 1 else deriv <- 0
      b1<-gam.fit3(x=X, y=y, sp=L%*%lsp1+lsp0, S=S,rS=rS,off=off, H=H,
          offset = offset,family = family,weights=weights,deriv=deriv,
          control=control,gamma=gamma,scale=scale,
          printWarn=FALSE,mustart=mustart,use.svd=use.svd,scoreType=scoreType,...)
#       ptm <- proc.time()-ptm
#       cat("deriv= ",deriv,"  ",ptm,"\n")
      
      if (scoreType=="REML") {
          score1 <- b1$REML1
      } else if (scoreType=="GACV") {
          score1 <- b1$GACV
      } else if (scoreType=="UBRE") {
          score1 <- b1$UBRE
      } else score1 <- b1$GCV
    } ## accepted step or step failed to lead to decrease

    if (ii < maxHalf) { ## step succeeded 
      mustart <- b1$fitted.values
      if (kk==0) { ## time for a full Newton step ...
#        ptm <- proc.time()
        b<-gam.fit3(x=X, y=y, sp=L%*%lsp1+lsp0, S=S,rS=rS,off=off, H=H,
               offset = offset,family = family,weights=weights,deriv=2,
               control=control,gamma=gamma,scale=scale,
               printWarn=FALSE,mustart=mustart,use.svd=use.svd,scoreType=scoreType,...)
#         ptm <- proc.time()-ptm
#         cat("deriv=2 ",ptm,"\n")

        mustart <- b$fitted.values
        old.score <- score;lsp <- lsp1
        if (scoreType=="REML") {
           score <- b$REML;grad <- b$REML1;hess <- b$REML2 
        } else if (scoreType=="GACV") {
          score <- b$GACV;grad <- b$GACV1;hess <- b$GACV2
        } else if (scoreType=="UBRE") {
          score <- b$UBRE;grad <- b$UBRE1;hess <- b$UBRE2 
        } else { score <- b$GCV;grad <- b$GCV1;hess <- b$GCV2}
        grad <- t(L)%*%grad
        hess <- t(L)%*%hess%*%L
      } else { ## just a BFGS update
        ## first derivatives only.... 
#        ptm <- proc.time()
         if (ii==1) b <- b1 else  
         b<-gam.fit3(x=X, y=y, sp=L%*%lsp1+lsp0, S=S,rS=rS,off=off, H=H,
               offset = offset,family = family,weights=weights,deriv=1,
               control=control,gamma=gamma,scale=scale,
               printWarn=FALSE,mustart=mustart,use.svd=use.svd,scoreType=scoreType,...)
#         ptm <- proc.time()-ptm
#         cat("deriv=1 ",ptm,"\n")

        mustart <- b$fitted.values
        old.score <- score;lsp <- lsp1
        old.grad <- grad
        if (scoreType=="REML") {
          score <- b$REML;grad <- b$REML1 
        } else if (scoreType=="GACV") {
          score <- b$GACV;grad <- b$GACV1
        } else if (scoreType=="UBRE") {
          score <- b$UBRE;grad <- b$UBRE1
        } else { score <- b$GCV;grad <- b$GCV1}
        grad <- t(L)%*%grad
        ## BFGS update of the inverse Hessian...
        yg <- grad-old.grad
        rho <- 1/sum(yg*step)
        B <- B - rho*step%*%(t(yg)%*%B)
        B <- B - rho*(B%*%yg)%*%t(step) + rho*step%*%t(step)
      } ## end of BFGS
    } ## end of successful step updating
    ## test for convergence
    converged <- TRUE
    score.scale <- b$scale.est + score;    
    uconv.ind <- abs(grad) > score.scale*conv.tol
    if (sum(uconv.ind)) converged <- FALSE
    if (abs(old.score-score)>score.scale*conv.tol) { 
      if (converged) uconv.ind <- uconv.ind | TRUE ## otherwise can't progress
      converged <- FALSE      
    }
    if (ii==maxHalf) converged <- TRUE ## step failure
    if (converged) break
  } ## end of iteration loop
  if (ii==maxHalf) ct <- "step failed"
  else if (i==200) ct <- "iteration limit reached" 
  else ct <- "full convergence"
  list(score=score,lsp=lsp,lsp.full=L%*%lsp,grad=grad,hess=hess,iter=i,conv =ct,object=b)
}



gam2derivative <- function(lsp,args,...)
## Performs IRLS GAM fitting for smoothing parameters given in lsp 
## and returns the derivatives of the GCV or UBRE score w.r.t the 
## smoothing parameters for the model.
## args is a list containing the arguments for gam.fit3
## For use as optim() objective gradient
{ if (!is.null(args$L)) {
    lsp <- args$L%*%lsp + args$lsp0
  }
  b<-gam.fit3(x=args$X, y=args$y, sp=lsp, S=args$S,rS=args$rS,off=args$off, H=args$H,
     offset = args$offset,family = args$family,weights=args$w,deriv=1,
     control=args$control,gamma=args$gamma,scale=args$scale,scoreType=args$scoreType,
     use.svd=FALSE,...)
  if (args$scoreType == "deviance") ret <- b$GCV1 else ret <- b$UBRE1
  if (!is.null(args$L)) ret <- t(args$L)%*%ret
  ret
}


gam2objective <- function(lsp,args,...)
## Performs IRLS GAM fitting for smoothing parameters given in lsp 
## and returns the GCV or UBRE score for the model.
## args is a list containing the arguments for gam.fit3
## For use as optim() objective
{ if (!is.null(args$L)) {
    lsp <- args$L%*%lsp + args$lsp0
  }
  b<-gam.fit3(x=args$X, y=args$y, sp=lsp, S=args$S,rS=args$rS,off=args$off, H=args$H,
     offset = args$offset,family = args$family,weights=args$w,deriv=0,
     control=args$control,gamma=args$gamma,scale=args$scale,scoreType=args$scoreType,
     use.svd=FALSE,...)
  if (args$scoreType == "deviance") ret <- b$GCV else ret <- b$UBRE
  attr(ret,"full.fit") <- b
  ret
}



gam4objective <- function(lsp,args,...)
## Performs IRLS GAM fitting for smoothing parameters given in lsp 
## and returns the GCV or UBRE score for the model.
## args is a list containing the arguments for gam.fit3
## For use as nlm() objective
{ if (!is.null(args$L)) {
    lsp <- args$L%*%lsp + args$lsp0
  }
  b<-gam.fit3(x=args$X, y=args$y, sp=lsp, S=args$S,rS=args$rS,off=args$off, H=args$H,
     offset = args$offset,family = args$family,weights=args$w,deriv=1,
     control=args$control,gamma=args$gamma,scale=args$scale,scoreType=args$scoreType,
     use.svd=FALSE,...)
  if (args$scoreType == "deviance") ret <- b$GCV else ret <- b$UBRE
  attr(ret,"full.fit") <- b
  if (args$scoreType == "deviance") at <- b$GCV1 else at <- b$UBRE1
  if (!is.null(args$L)) at <- t(args$L)%*%at
  attr(ret,"gradient") <- at
  ret
}

##
## The following fix up family objects for use with gam.fit3
##


fix.family.link<-function(fam)
# adds d2link the second derivative of the link function w.r.t. mu
# to the family supplied, as well as a 3rd derivative function 
# d3link...
# All d2link and d3link functions have been checked numerically. 
{ if (!inherits(fam,"family")) stop("fam not a family object")
  if (is.null(fam$canonical)) { ## note the canonical link - saves effort in full Newton
    if (fam$family=="gaussian") fam$canonical <- "identity" else
    if (fam$family=="poisson"||fam$family=="quasipoisson") fam$canonical <- "log" else
    if (fam$family=="binomial"||fam$family=="quasibinomial") fam$canonical <- "logit" else
    if (fam$family=="Gamma") fam$canonical <- "inverse" else
    if (fam$family=="inverse.gaussian") fam$canonical <- "1/mu^2" else
    fam$canonical <- "none"
  }
  if (!is.null(fam$d2link)&&!is.null(fam$d3link)&&!is.null(fam$d4link)) return(fam) 
  link <- fam$link
  if (length(link)>1) if (fam$family=="quasi") # then it's a power link
  { lambda <- log(fam$linkfun(exp(1))) ## the power, if > 0
    if (lambda<=0) { fam$d2link <- function(mu) -1/mu^2
      fam$d3link <- function(mu) 2/mu^3
      fam$d4link <- function(mu) -6/mu^4
    }
    else { fam$d2link <- function(mu) lambda*(lambda-1)*mu^(lambda-2)
      fam$d3link <- function(mu) (lambda-2)*(lambda-1)*lambda*mu^(lambda-3)
      fam$d4link <- function(mu) (lambda-3)*(lambda-2)*(lambda-1)*lambda*mu^(lambda-4)
    }
    return(fam)
  } else stop("unrecognized (vector?) link")

  if (link=="identity") {
    fam$d4link <- fam$d3link <- fam$d2link <- 
    function(mu) rep.int(0,length(mu))
    return(fam)
  } 
  if (link == "log") {
    fam$d2link <- function(mu) -1/mu^2
    fam$d3link <- function(mu) 2/mu^3
    fam$d4link <- function(mu) -6/mu^4
    return(fam)
  }
  if (link == "inverse") {
    fam$d2link <- function(mu) 2/mu^3
    fam$d3link <- function(mu) { mu <- mu*mu;-6/(mu*mu)}
    fam$d4link <- function(mu) { mu2 <- mu*mu;24/(mu2*mu2*mu)}
    return(fam)
  }
  if (link == "logit") {
    fam$d2link <- function(mu) 1/(1 - mu)^2 - 1/mu^2
    fam$d3link <- function(mu) 2/(1 - mu)^3 + 2/mu^3
    fam$d4link <- function(mu) 6/(1-mu)^4 - 6/mu^4
    return(fam)
  }
  if (link == "probit") {
    fam$d2link <- function(mu) { 
      eta <- fam$linkfun(mu)
      eta/fam$mu.eta(eta)^2
    }
    fam$d3link <- function(mu) {
      eta <-  fam$linkfun(mu)
      (1 + 2*eta^2)/fam$mu.eta(eta)^3
    }
    fam$d4link <- function(mu) {
       eta <-  fam$linkfun(mu)
       (7*eta + 6*eta^3)/fam$mu.eta(eta)^4
    }
    return(fam)
  }
  if (link == "cloglog") {
    fam$d2link <- function(mu) { l1m <- log(1-mu)
      -1/((1 - mu)^2*l1m) *(1+ 1/l1m)
    }
    fam$d3link <- function(mu) { l1m <- log(1-mu)
       mu3 <- (1-mu)^3
      (-2 - 3*l1m - 2*l1m^2)/mu3/l1m^3
    }
    fam$d4link <- function(mu){
      l1m <- log(1-mu)
      mu4 <- (1-mu)^4
      ( - 12 - 11 * l1m - 6 * l1m^2 - 6/l1m )/mu4  /l1m^3
    }
    return(fam)
  }
  if (link == "sqrt") {
    fam$d2link <- function(mu) -.25 * mu^-1.5
    fam$d3link <- function(mu) .375 * mu^-2.5
    fam$d4link <- function(mu) -0.9375 * mu^-3.5
    return(fam)
  }
  if (link == "cauchit") {
    fam$d2link <- function(mu) { 
     eta <- fam$linkfun(mu)
     2*pi*pi*eta*(1+eta*eta)
    }
    fam$d3link <- function(mu) { 
     eta <- fam$linkfun(mu)
     eta2 <- eta*eta
     2*pi*pi*pi*(1+3*eta2)*(1+eta2)
    }
    fam$d4link <- function(mu) { 
     eta <- fam$linkfun(mu)
     eta2 <- eta*eta
     2*pi^4*(8*eta+12*eta2*eta)*(1+eta2)
    }
    return(fam)
  }
  if (link == "1/mu^2") {
    fam$d2link <- function(mu) 6 * mu^-4
    fam$d3link <- function(mu) -24 * mu^-5
    fam$d4link <- function(mu) 120 * mu^-6
    return(fam)
  }
  if (substr(link,1,3)=="mu^") { ## it's a power link
    ## note that lambda <=0 gives log link so don't end up here
    lambda <- get("lambda",environment(fam$linkfun))
    fam$d2link <- function(mu) (lambda*(lambda-1)) * mu^{lambda-2}
    fam$d3link <- function(mu) (lambda*(lambda-1)*(lambda-2)) * mu^{lambda-3}
    fam$d4link <- function(mu) (lambda*(lambda-1)*(lambda-2)*(lambda-3)) * mu^{lambda-4}
    return(fam)
  }
  stop("link not recognised")
}


fix.family.var<-function(fam)
# adds dvar the derivative of the variance function w.r.t. mu
# to the family supplied, as well as d2var the 2nd derivative of 
# the variance function w.r.t. the mean. (All checked numerically). 
{ if (!inherits(fam,"family")) stop("fam not a family object")
  if (!is.null(fam$dvar)&&!is.null(fam$d2var)&&!is.null(fam$d3var)) return(fam) 
  family <- fam$family
  if (family=="gaussian") {
    fam$d3var <- fam$d2var <- fam$dvar <- function(mu) rep.int(0,length(mu))
    return(fam)
  } 
  if (family=="poisson"||family=="quasipoisson") {
    fam$dvar <- function(mu) rep.int(1,length(mu))
    fam$d3var <- fam$d2var <- function(mu) rep.int(0,length(mu))
    return(fam)
  } 
  if (family=="binomial"||family=="quasibinomial") {
    fam$dvar <- function(mu) 1-2*mu
    fam$d2var <- function(mu) rep.int(-2,length(mu))
    fam$d3var <- function(mu) rep.int(0,length(mu))
    return(fam)
  }
  if (family=="Gamma") {
    fam$dvar <- function(mu) 2*mu
    fam$d2var <- function(mu) rep.int(2,length(mu))
    fam$d3var <- function(mu) rep.int(0,length(mu))
    return(fam)
  }
  if (family=="quasi") {
    fam$dvar <- switch(fam$varfun,
       constant = function(mu) rep.int(0,length(mu)),
       "mu(1-mu)" = function(mu) 1-2*mu,
       mu = function(mu) rep.int(1,length(mu)),
       "mu^2" = function(mu) 2*mu,
       "mu^3" = function(mu) 3*mu^2           
    )
    if (is.null(fam$dvar)) stop("variance function not recognized for quasi")
    fam$d2var <- switch(fam$varfun,
       constant = function(mu) rep.int(0,length(mu)),
       "mu(1-mu)" = function(mu) rep.int(-2,length(mu)),
       mu = function(mu) rep.int(0,length(mu)),
       "mu^2" = function(mu) rep.int(2,length(mu)),
       "mu^3" = function(mu) 6*mu           
    )
    fam$d3var <- switch(fam$varfun,
       constant = function(mu) rep.int(0,length(mu)),
       "mu(1-mu)" = function(mu) rep.int(0,length(mu)),
       mu = function(mu) rep.int(0,length(mu)),
       "mu^2" = function(mu) rep.int(0,length(mu)),
       "mu^3" = function(mu) rep.int(6,length(mu))           
    )
    return(fam)
  }
  if (family=="inverse.gaussian") {
    fam$dvar <- function(mu) 3*mu^2
    fam$d2var <- function(mu) 6*mu
    fam$d3var <- function(mu) rep.int(6,length(mu)) 
    return(fam)
  }
  stop("family not recognised")
}


fix.family.ls<-function(fam)
# adds ls the log saturated likelihood and its derivatives
# w.r.t. the scale parameter to the family object.
{ if (!inherits(fam,"family")) stop("fam not a family object")
  if (!is.null(fam$ls)) return(fam) 
  family <- fam$family
  if (family=="gaussian") {
    fam$ls <- function(y,w,n,scale) c(-sum(w)*log(2*pi*scale)/2,-sum(w)/(2*scale),sum(w)/(2*scale*scale))
    return(fam)
  } 
  if (family=="poisson") {
    fam$ls <- function(y,w,n,scale) {
      res <- rep(0,3)
      res[1] <- sum(dpois(y,y,log=TRUE)*w)
      res
    }
    return(fam)
  } 
  if (family=="binomial") {
    fam$ls <- function(y,w,n,scale) { 
      c(-binomial()$aic(y,n,y,w,0)/2,0,0)
    }
    return(fam)
  }
  if (family=="Gamma") {
    fam$ls <- function(y,w,n,scale) {
      res <- rep(0,3)
      k <- -lgamma(1/scale) - log(scale)/scale - 1/scale
      res[1] <- sum(w*(k-log(y)))
      k <- (digamma(1/scale)+log(scale))/(scale*scale)
      res[2] <- sum(w*k)  
      k <- (-trigamma(1/scale)/(scale) + (1-2*log(scale)-2*digamma(1/scale)))/(scale^3)
      res[3] <- sum(w*k) 
      res
    }
    return(fam)
  }
  if (family=="quasi"||family=="quasipoisson"||family=="quasibinomial") {
    fam$ls <- function(y,w,n,scale) rep(0,3)
    return(fam)
  }
  if (family=="inverse.gaussian") {
    fam$ls <- function(y,w,n,scale) c(-sum(w*log(2*pi*scale*y^3))/2,
     -sum(w)/(2*scale),sum(w)/(2*scale*scale))
    return(fam)
  }
  stop("family not recognised")
}


negbin <- function (theta = stop("'theta' must be specified"), link = "log") { 
## modified from Venables and Ripley's MASS library to work with gam.fit3,
## and to allow a range of `theta' values to be specified
## single `theta' to specify fixed value; 2 theta values (first smaller that second)
## are limits within which to search for theta; otherwise supplied values make up 
## search set.
  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  if (linktemp %in% c("log", "identity", "sqrt")) stats <- make.link(linktemp)
  else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  } else {
    if (inherits(link, "link-glm")) {
       stats <- link
            if (!is.null(stats$name))
                linktemp <- stats$name
        }
        else stop(linktemp, " link not available for negative binomial family; available links are \"identity\", \"log\" and \"sqrt\"")
    }
    env <- new.env(parent = .GlobalEnv)
    assign(".Theta", theta, envir = env)
    variance <- function(mu) mu + mu^2/get(".Theta")
    ## dvaraince/dmu needed as well
    dvar <- function(mu) 1 + 2*mu/get(".Theta")
    ## d2variance/dmu...
    d2var <- function(mu) rep(2/get(".Theta"),length(mu))
    d3var <- function(mu) rep(0,length(mu))
    getTheta <- function() get(".Theta")
    validmu <- function(mu) all(mu > 0)

    dev.resids <- function(y, mu, wt) { Theta <- get(".Theta")
      2 * wt * (y * log(pmax(1, y)/mu) - 
        (y + Theta) * log((y + Theta)/(mu + Theta))) 
    }
    aic <- function(y, n, mu, wt, dev) {
        Theta <- get(".Theta")
        term <- (y + Theta) * log(mu + Theta) - y * log(mu) +
            lgamma(y + 1) - Theta * log(Theta) + lgamma(Theta) -
            lgamma(Theta + y)
        2 * sum(term * wt)
    }
    initialize <- expression({
        if (any(y < 0)) stop("negative values not allowed for the negative binomial family")
        n <- rep(1, nobs)
        mustart <- y + (y == 0)/6
    })
    environment(dvar) <- environment(d2var) <- environment(variance) <- environment(validmu) <- 
                         environment(dev.resids) <- environment(aic) <- environment(getTheta) <- env
    famname <- paste("Negative Binomial(", format(round(theta,3)), ")", sep = "")
    structure(list(family = famname, link = linktemp, linkfun = stats$linkfun,
        linkinv = stats$linkinv, variance = variance,dvar=dvar,d2var=d2var,d3var=d3var, dev.resids = dev.resids,
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize,
        validmu = validmu, valideta = stats$valideta,getTheta = getTheta,canonical="log"), class = "family")
}



totalPenalty <- function(S,H,off,theta,p)
{ if (is.null(H)) St <- matrix(0,p,p)
  else { St <- H; 
    if (ncol(H)!=p||nrow(H)!=p) stop("H has wrong dimension")
  }
  theta <- exp(theta)
  m <- length(theta)
  if (m>0) for (i in 1:m) {
    k0 <- off[i]
    k1 <- k0 + nrow(S[[i]]) - 1
    St[k0:k1,k0:k1] <- St[k0:k1,k0:k1] + S[[i]] * theta[i]
  }
  St
}

mini.roots <- function(S,off,np)
# function to obtain square roots, B[[i]], of S[[i]]'s having as few
# columns as possible. S[[i]]=B[[i]]%*%t(B[[i]]). np is the total number
# of parameters. S is in packed form. 
{ m<-length(S)
  B<-S
  for (i in 1:m)
  { b<-mroot(S[[i]])
    B[[i]]<-matrix(0,np,ncol(b))
    B[[i]][off[i]:(off[i]+nrow(b)-1),]<-b
  }
  B
}


ldTweedie <- function(y,mu=y,p=1.5,phi=1) {
## evaluates log Tweedie density for 1<=p<=2, using series summation of
## Dunn & Smyth (2005) Statistics and Computing 15:267-280.
 
  if (length(p)>1||length(phi)>1) stop("only scalar `p' and `phi' allowed.")
  if (p<1||p>2) stop("p must be in [1,2]")
  ld <- cbind(y,y,y)
  if (p == 2) {
    ld[,1] <- dgamma(y, shape = 1/phi,rate = 1/(phi * mu),log=TRUE)
    ld[,2] <- (digamma(1/phi) + log(phi) - 1 + y/mu - log(y/mu))/(phi*phi)
    ld[,3] <- -2*ld[,2]/phi + (1-trigamma(1/phi)/phi)/(phi^3)
    return(ld)
  }  
  if (p == 1) {
    ## ld[,1] <- dpois(x = y/phi, lambda = mu/phi,log=TRUE)
    bkt <- (y*log(mu/phi) - mu)
    dig <- digamma(y/phi+1)
    trig <- trigamma(y/phi+1)
    ld[,1] <- bkt/phi - lgamma(y/phi+1)
    ld[,2] <- (-bkt - y + dig*y)/(phi*phi)
    ld[,3] <- (2*bkt + 3*y - 2*dig*y - trig *y*y/phi)/(phi^3)
    return(ld) 
  }
  ## .. otherwise need the full series thing....
  ## first deal with the zeros  
  if (length(mu)==1) mu <- rep(mu,length(y))
  ind <- y==0
 
  ld[ind,1] <- -mu[ind]^(2-p)/(phi*(2-p))
  ld[ind,2] <- -ld[ind,1]/phi
  ld[ind,3] <- -2*ld[ind,2]/phi

  if (sum(!ind)==0) return(ld)

  ## now the non-zeros
  y <- y[!ind];mu <- mu[!ind]
  w <- w1 <- w2 <- y*0
  oo <- .C(C_tweedious,w=as.double(w),w1=as.double(w1),w2=as.double(w2),y=as.double(y),
           phi=as.double(phi),p=as.double(p),eps=as.double(.Machine$double.eps),n=as.integer(length(y)))
  theta <- mu^(1-p)
  k.theta <- mu*theta/(2-p)
  theta <- theta/(1-p)
  l.base <-  (y*theta-k.theta)/phi
  ld[!ind,1] <- l.base - log(y) + oo$w
  ld[!ind,2] <- -l.base/phi + oo$w1   
  ld[!ind,3] <- 2*l.base/(phi*phi) + oo$w2
  
  ld
}

Tweedie <- function(p=1,link=power(0)) {
## a restricted Tweedie family
  if (p<1||p>2) stop("`p' can not be less than 1 or greater than 2")
  
  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  okLinks <- c("log", "identity", "sqrt","inverse")
  if (linktemp %in% okLinks)
    stats <- make.link(linktemp) else 
  if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  } else {
    if (inherits(link, "link-glm")) {
       stats <- link
       if (!is.null(stats$name))
          linktemp <- stats$name
        } else {
            stop(gettextf("link \"%s\" not available for poisson family.",
                linktemp, collapse = ""),domain = NA)
        }
    }
    
    variance <- function(mu) mu^p
    dvar <- function(mu) p*mu^(p-1)
    if (p==1) d2var <- function(mu) 0*mu else
      d2var <- function(mu) p*(p-1)*mu^(p-2)
    if (p==1||p==2)  d3var <- function(mu) 0*mu else
      d3var <- function(mu) p*(p-1)*(p-2)*mu^(p-3)
    validmu <- function(mu) all(mu >= 0)

    dev.resids <- function(y, mu, wt) {
        y1 <- y + (y == 0)
        if (p == 1)
            theta <- log(y1/mu)
        else theta <- (y1^(1 - p) - mu^(1 - p))/(1 - p)
        if (p == 2)
            kappa <- log(y1/mu)
        else kappa <- (y^(2 - p) - mu^(2 - p))/(2 - p)
        2 * wt * (y * theta - kappa)
    }
    initialize <- expression({
        n <- rep(1, nobs)
        mustart <- y + 0.1 * (y == 0)
    })
    ls <-  function(y,w,n,scale) {
      power <- p
      colSums(w*ldTweedie(y,y,p=power,phi=scale))
    }

    aic <- function(y, n, mu, wt, dev) {
      power <- p
      scale <- dev/sum(wt)
      -2*sum(ldTweedie(y,mu,p=power,phi=scale)[,1]*wt) + 2
    }
    structure(list(family = paste("Tweedie(",p,")",sep=""), variance = variance, 
              dev.resids = dev.resids,aic = aic, link = linktemp, linkfun = stats$linkfun, linkinv = stats$linkinv,
        mu.eta = stats$mu.eta, initialize = initialize, validmu = validmu,
        valideta = stats$valideta,dvar=dvar,d2var=d2var,d3var=d3var,ls=ls,canonical="none"), class = "family")


}