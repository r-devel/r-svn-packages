## (c) Simon N. Wood (ocat, tw, nb) & Natalya Pya (t.scaled, beta, zip), 
## 2013, 2014. Released under GPL2.

## extended families for mgcv ...

## extended family object for ordered categorical

ocat <- function(theta=NULL,link="identity",R=NULL) {
## extended family object for ordered categorical model.
## one of theta and R must be supplied. length(theta) == R-2.
## weights are ignored.
  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  if (linktemp %in% c("identity")) stats <- make.link(linktemp)
  else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  } else {
    if (inherits(link, "link-glm")) {
      stats <- link
      if (!is.null(stats$name))
            linktemp <- stats$name
    } else stop(linktemp, " link not available for ordered categorical family; available links are \"identity\"")
  }
  if (is.null(theta)&&is.null(R)) stop("Must supply theta or R to ocat")
  if (!is.null(theta)) R <- length(theta) + 2 ## number of catergories
  Theta <-  NULL;n.theta <- R-2
  ## NOTE: data based initialization is in preinitialize...
  if (!is.null(theta)&&sum(theta==0)==0) {
    if (sum(theta<0)) iniTheta <- log(abs(theta)) ## initial theta supplied
    else { 
      iniTheta <- Theta <- log(theta) ## fixed theta supplied
      n.theta <- 0
    }
  } else iniTheta <- rep(-1,length=R-2) ## inital log theta value

  env <- new.env(parent = .GlobalEnv)
  assign(".Theta", iniTheta, envir = env)

  putTheta <-function(theta) assign(".Theta", theta,envir=environment(sys.function()))

  getTheta <-function(trans=FALSE) { 
    theta <- get(".Theta")
    if (trans) { ## transform to (finite) cut points...
      R = length(theta)+2
      alpha <- rep(0,R-1) ## the thresholds
      alpha[1] <- -1
      if (R > 2) { 
        ind <- 2:(R-1)
        alpha[ind] <- alpha[1] + cumsum(exp(theta))
      } 
      theta <- alpha
    }
    theta
  }
  postproc <- expression({
    object$family$family <- 
    paste("Ordered Categorical(",paste(round(object$family$getTheta(TRUE),2),collapse=","),")",sep="")
  })

  validmu <- function(mu) all(is.finite(mu))

  dev.resids <- function(y, mu, wt,theta=NULL) {
    F <- function(x) {
      h <- ind <- x > 0; h[ind] <- 1/(exp(-x[ind]) + 1)
      x <- exp(x[!ind]); h[!ind] <- (x/(1+x))
      h 
    }
    Fdiff <- function(a,b) {
      ## cancellation resistent F(b)-F(a), b>a
      h <- rep(1,length(b)); h[b>0] <- -1; eb <- exp(b*h)
      h <- h*0+1; h[a>0] <- -1; ea <- exp(a*h)
      ind <- b<0;bi <- eb[ind];ai <- ea[ind]
      h[ind] <- bi/(1+bi) - ai/(1+ai)
      ind1 <- a>0;bi <- eb[ind1];ai <- ea[ind1]
      h[ind1] <- (ai-bi)/((ai+1)*(bi+1))
      ind <- !ind & !ind1;bi <- eb[ind];ai <- ea[ind]
      h[ind] <- (1-ai*bi)/((bi+1)*(ai+1))
      h
    }
    if (is.null(theta)) theta <- get(".Theta")
    R = length(theta)+2
    alpha <- rep(0,R+1) ## the thresholds
    alpha[1] <- -Inf;alpha[R+1] <- Inf
    alpha[2] <- -1
    if (R > 2) { 
      ind <- 3:R
      alpha[ind] <- alpha[2] + cumsum(exp(theta))
    } 
    al1 <- alpha[y+1];al0 = alpha[y] ## cut points above and below y
    ## Compute sign for deviance residuals, based on sign of interval
    ## midpoint for each datum minus the fitted value of the latent 
    ## variable. This makes first and last categories like 0s and 1s in
    ## logistic regression....
    s <- sign((al1 + al0)/2-mu) ## sign for deviance residuals
    al1mu <- al1-mu;al0mu <- al0-mu
    f1 <- F(al1mu);f0 <- F(al0mu);
    ##f <- pmax(f1 - f0,.Machine$double.eps)
    f <- Fdiff(al0mu,al1mu)
    ##a1 <- f1^2 - f1;a0 <- f0^2 - f0; a <- a1 -a0
    al1al0 <- (al1-al0)/2;al0al1 <- (al0-al1)/2
    g1 <- F(al1al0);g0 <- F(al0al1)
    ##A <- pmax(g1 - g0,.Machine$double.eps)
    A <- Fdiff(al0al1,al1al0)
    rsd <- 2*(log(A)-log(f))
    attr(rsd,"sign") <- s
    rsd
  } ## end of dev.resids

  Dd <- function(y, mu, theta, wt=NULL, level=0) {
  ## derivatives of the deviance...
    F <- function(x) { ## e^(x)/(1+e^x) without overflow
      h <- ind <- x > 0; h[ind] <- 1/(exp(-x[ind]) + 1)
      x <- exp(x[!ind]); h[!ind] <- (x/(1+x))
      h 
    }
    Fdiff <- function(a,b) {
      ## cancellation resistent F(b)-F(a), b>a 
      h <- rep(1,length(b)); h[b>0] <- -1; eb <- exp(b*h)
      h <- h*0+1; h[a>0] <- -1; ea <- exp(a*h)
      ind <- b<0;bi <- eb[ind];ai <- ea[ind]
      h[ind] <- bi/(1+bi) - ai/(1+ai)
      ind1 <- a>0;bi <- eb[ind1];ai <- ea[ind1]
      h[ind1] <- (ai-bi)/((ai+1)*(bi+1))
      ind <- !ind & !ind1;bi <- eb[ind];ai <- ea[ind]
      h[ind] <- (1-ai*bi)/((bi+1)*(ai+1))
      h
    }
    abcd <- function(x,level=-1) {
      bj <- cj <- dj <- NULL
      ## compute f_j^2 - f_j without cancellation error
      ## x <- 10;F(x)^2-F(x);abcd(x)$aj
      h <- rep(1,length(x)); h[x>0] <- -1; ex <- exp(x*h)
      ex1 <- ex+1;ex1k <- ex1^2
      aj <- -ex/ex1k
      if (level>=0) {
        ## compute f_j - 3 f_j^2 + 2f_j^3 without cancellation error
        ##  x <- 10;F(x)-3*F(x)^2+2*F(x)^3;abcd(x,0)$bj
        ex1k <- ex1k*ex1;ex2 <- ex^2
        bj <- h*(ex-ex^2)/ex1k
        if (level>0) {
          ## compute -f_j + 7 f_j^2 - 12 f_j^3 + 6 f_j^4
          ## x <- 10;-F(x)+7*F(x)^2-12*F(x)^3+6*F(x)^4;abcd(x,1)$cj
          ex1k <- ex1k*ex1;ex3 <- ex2*ex
          cj <- (-ex3 + 4*ex2 - ex)/ex1k
          if (level>1) {    
            ## compute d_j
            ## x <- 10;F(x)-15*F(x)^2+50*F(x)^3-60*F(x)^4+24*F(x)^5;abcd(x,2)$dj
            ex1k <- ex1k*ex1;ex4 <- ex3*ex
            dj <- h * (-ex4 + 11*ex3 - 11*ex2 + ex)/ex1k
          }
        }
      }        
      list(aj=aj,bj=bj,cj=cj,dj=dj)
    }
 
    R = length(theta)+2
    alpha <- rep(0,R+1) ## the thresholds
    alpha[1] <- -Inf;alpha[R+1] <- Inf
    alpha[2] <- -1
    if (R > 2) { 
      ind <- 3:R
      alpha[ind] <- alpha[2] + cumsum(exp(theta))
    } 
    al1 <- alpha[y+1];al0 = alpha[y]
    al1mu <- al1-mu;al0mu <- al0 - mu
    f1 <- F(al1mu);f0 <- F(al0mu);
    ##f <- pmax(f1 - f0,.Machine$double.eps)
    f <- pmax(Fdiff(al0mu,al1mu),.Machine$double.xmin)
    r1 <- abcd(al1mu,level); a1 <- r1$aj
    r0 <- abcd(al0mu,level); a0 <- r0$aj
    ## a1 <- f1^2 - f1;a0 <- f0^2 - f0; 
    a <- a1 - a0
    
    al1al0 <- (al1-al0)/2; al0al1 <- (al0-al1)/2
    g1 <- F(al1al0);g0 <- F(al0al1)
    ##A <- pmax(g1 - g0,.Machine$double.eps)
    A <- Fdiff(al0al1,al1al0)
    if (level>=0) {
      ## b1 <- f1 - 3 * f1^2 + 2 * f1^3;b0 <- f0 - 3 * f0^2 + 2 * f0^3  
      b1 <- r1$bj;b0 <- r0$bj
      b <- b1 - b0
    }
    if (level>0) {
      ##c1 <- -f1 + 7 * f1^2 - 12* f1^3 + 6 * f1^4
      ##c0 <- -f0 + 7 * f0^2 - 12* f0^3 + 6 * f0^4
      c1 <- r1$cj;c0 <- r0$cj
      c <- c1 - c0
      R1 <- abcd(al1al0,level-2) 
      R0 <- abcd(al0al1,level-2)
      ## B <- g1^2 - g1 + g0^2 - g0
      B <- R1$aj + R0$aj
    }
    if (level>1) {
      ##d1 <- f1 - 15 * f1^2 + 50 * f1^3 - 60 * f1^4 + 24 * f1^5
      ##d0 <- f0 - 15 * f0^2 + 50 * f0^3 - 60 * f0^4 + 24 * f0^5
      d1 <- r1$dj;d0 <- r0$dj
      d <- d1 - d0
      ##C <- 2 * g1^3 - 3 * g1^2 + g1 - 2 * g0^3 + 3 * g0^2 - g0
      C <- R1$bj - R0$bj
    }

    oo <- list(D=NULL,Dmu=NULL,Dmu2=NULL,Dth=NULL,Dmuth=NULL,
             Dmu2th=NULL)
    n <- length(y)
    ## deviance...
    oo$D <- 2*(log(A)-log(f))
    if (level >= 0) { ## get derivatives used in coefficient estimation
      oo$Dmu <- -2 * a / f
      a2 <- a^2
      oo$EDmu2 <- oo$Dmu2 <- 2 * (a2/f -  b)/f
    }
    if (R<3) level <- 0 ## no free parameters

    if (level > 0) { ## get first derivative related stuff
      f2 <- f^2;a3 <- a2*a
      oo$Dmu3 <- 2*(- c - 2 * a3/f2 + 3 * a * b/f)/f
      Dmua0 <- 2 * (a0 * a/f -  b0)/f
      Dmua1 <- -2 * (a1 * a /f - b1)/f
      Dmu2a0 <- -2* (c0 + (a0*(2*a2/f - b)- 2*b0*a  )/f)/f
      Dmu2a1 <- 2*(c1  + (2*(a1*a2/f - b1*a) - a1*b)/f)/f
      Da0 <- B/A - 2*a0/f; Da1 <- -B/A + 2 * a1/f 
      ## now transform to derivatives w.r.t. theta...
      oo$Dmu2th <- oo$Dmuth <- oo$Dth <- matrix(0,n,R-2)
      for (k in 1:(R-2)) { 
        etk <- exp(theta[k])
        ind <- y == k+1
        oo$Dth[ind,k] <- Da1[ind]*etk
        oo$Dmuth[ind,k] <- Dmua1[ind]*etk
        oo$Dmu2th[ind,k] <- Dmu2a1[ind]*etk
    
        if (R>k+2) { 
          ind <- y > k+1 & y < R
          oo$Dth[ind,k] <- (Da1[ind]+Da0[ind])*etk
          oo$Dmuth[ind,k] <- (Dmua1[ind]+Dmua0[ind])*etk
          oo$Dmu2th[ind,k] <- (Dmu2a1[ind]+Dmu2a0[ind])*etk
       
        }
        ind <- y == R
        oo$Dth[ind,k] <- Da0[ind]*etk
        oo$Dmuth[ind,k] <- Dmua0[ind]*etk
        oo$Dmu2th[ind,k] <- Dmu2a0[ind]*etk 
      } 
    }  
    if (level >1) { ## and the second derivative components 
      oo$Dmu4 <- 2*((3*b^2 + 4*a*c)/f + a2*(6*a2/f - 12*b)/f2 - d)/f
      Dmu3a0 <-  2 * ((a0*c  + 3*c0*a + 3*b0*b)/f - d0  + 
                      6*a*(a0*a2/f - b0*a - a0*b)/f2 )/f
      Dmu3a1 <- 2 * (d1 - (a1*c + 3*(c1*a + b1*b))/f
                     + 6*a*(b1*a - a1*a2/f  + a1*b)/f2)/f

      Dmua0a0 <- 2*(c0 + (2*a0*(b0 - a0*a/f) - b0*a)/f )/f
      Dmua1a1 <- 2*( (b1*a + 2*a1*(b1 - a1*a/f))/f - c1)/f
      Dmua0a1 <- 2*(a0*(2*a1*a/f - b1) - b0*a1)/f2 

      Dmu2a0a0 <- 2*(d0 + (b0*(2*b0 - b)  + 2*c0*(a0 - a))/f +
                     2*(b0*a2 + a0*(3*a0*a2/f  - 4*b0*a - a0*b))/f2)/f

      Dmu2a1a1 <-  2*( (2*c1*(a + a1) + b1*(2*b1 + b))/f
                + 2*(a1*(3*a1*a2/f  - a1*b) - b1*a*(a + 4*a1))/f2 - d1)/f

      Dmu2a0a1 <- 0 ## (8*a0*b1*a/f^3 + 8*b0*a1*a/f^3 - 12*a0*a1*a/f^4 - 4*b0*b1/f^2 +
                    ## 4*a0*a1*b/f^3 - 2*c0*a1/f^2 - 2*c1*a0/f^2)

      Da0a0 <- 2 * (b0 + a0^2/f)/f + .5 * (C - B^2/A)/A
      Da1a1 <- -2* (b1 - a1^2/f)/f + .5 * (C - B^2/A)/A
      Da0a1 <- -2*a0*a1/f2 - .5 * (C - B^2/A)/A

      ## now transform to derivatives w.r.t. theta...
      n2d <- (R-2)*(R-1)/2
      oo$Dmu3th <- matrix(0,n,R-2)
      oo$Dmu2th2 <- oo$Dmuth2 <- oo$Dth2 <- matrix(0,n,n2d)
      i <- 0
      for (j in 1:(R-2)) for (k in j:(R-2)) { 
        i <- i + 1 ## the second deriv col
        ind <- y >= j ## rest are zero
        ar1.k <- ar.k <- rep(exp(theta[k]),n)
        ar.k[y==R | y <= k] <- 0; ar1.k[y<k+2] <- 0
        ar.j <- ar1.j <- rep(exp(theta[j]),n)
        ar.j[y==R | y <= j] <- 0; ar1.j[y<j+2] <- 0
        ar.kj <- ar1.kj <- rep(0,n)
        if (k==j) {
          ar.kj[y>k&y<R] <- exp(theta[k])
          ar1.kj[y>k+1] <- exp(theta[k])
          oo$Dmu3th[ind,k] <- Dmu3a1[ind]*ar.k[ind]  + Dmu3a0[ind]*ar1.k[ind]
        }
        oo$Dth2[,i] <- Da1a1*ar.k*ar.j + Da0a1*ar.k*ar1.j + Da1 * ar.kj +
                     Da0a0*ar1.k*ar1.j + Da0a1*ar1.k*ar.j + Da0 * ar1.kj
        oo$Dmuth2[,i] <- Dmua1a1*ar.k*ar.j + Dmua0a1*ar.k*ar1.j + Dmua1 * ar.kj +
                     Dmua0a0*ar1.k*ar1.j + Dmua0a1*ar1.k*ar.j + Dmua0 * ar1.kj
        oo$Dmu2th2[,i] <- Dmu2a1a1*ar.k*ar.j + Dmu2a0a1*ar.k*ar1.j + Dmu2a1 * ar.kj +
                     Dmu2a0a0*ar1.k*ar1.j + Dmu2a0a1*ar1.k*ar.j + Dmu2a0 * ar1.kj
      } 
    }
    oo
  } ## end of Dd
 
  aic <- function(y, mu, theta=NULL, wt, dev) {
  
   Fdiff <- function(a,b) {
      ## cancellation resistent F(b)-F(a), b>a 
      h <- rep(1,length(b)); h[b>0] <- -1; eb <- exp(b*h)
      h <- h*0+1; h[a>0] <- -1; ea <- exp(a*h)
      ind <- b<0;bi <- eb[ind];ai <- ea[ind]
      h[ind] <- bi/(1+bi) - ai/(1+ai)
      ind1 <- a>0;bi <- eb[ind1];ai <- ea[ind1]
      h[ind1] <- (ai-bi)/((ai+1)*(bi+1))
      ind <- !ind & !ind1;bi <- eb[ind];ai <- ea[ind]
      h[ind] <- (1-ai*bi)/((bi+1)*(ai+1))
      h
    }
    if (is.null(theta)) theta <- get(".Theta")
    R = length(theta)+2
    alpha <- rep(0,R+1) ## the thresholds
    alpha[1] <- -Inf;alpha[R+1] <- Inf
    alpha[2] <- -1
    if (R > 2) { 
      ind <- 3:R
      alpha[ind] <- alpha[2] + cumsum(exp(theta))
    } 
    al1 <- alpha[y+1];al0 = alpha[y]
    ##f1 <- F(al1-mu);f0 <- F(al0-mu);f <- f1 - f0
    f <- Fdiff(al0-mu,al1-mu)
    -2*sum(log(f))
  } ## end aic

  ls <- function(y,w,n,theta,scale) {
    ## the log saturated likelihood function.
    F <- function(x) {
      h <- ind <- x > 0; h[ind] <- 1/(exp(-x[ind]) + 1)
      x <- exp(x[!ind]); h[!ind] <- (x/(1+x))
      h 
    } 
    Fdiff <- function(a,b) {
      ## cancellation resistent F(b)-F(a), b>a
      h <- rep(1,length(b)); h[b>0] <- -1; eb <- exp(b*h)
      h <- h*0+1; h[a>0] <- -1; ea <- exp(a*h)
      ind <- b<0;bi <- eb[ind];ai <- ea[ind]
      h[ind] <- bi/(1+bi) - ai/(1+ai)
      ind1 <- a>0;bi <- eb[ind1];ai <- ea[ind1]
      h[ind1] <- (ai-bi)/((ai+1)*(bi+1))
      ind <- !ind & !ind1;bi <- eb[ind];ai <- ea[ind]
      h[ind] <- (1-ai*bi)/((bi+1)*(ai+1))
      h
    }
    R = length(theta)+2
    alpha <- rep(0,R+1) ## the thresholds
    alpha[1] <- -Inf;alpha[R+1] <- Inf
    alpha[2] <- -1
    if (R > 2) { 
      ind <- 3:R
      alpha[ind] <- alpha[2] + cumsum(exp(theta))
    } 
    al1 <- alpha[y+1];al0 = alpha[y]
    g1 <- F((al1-al0)/2);g0 <- F((al0-al1)/2)
    ##A <- pmax(g1 - g0,.Machine$double.eps)
    A <- Fdiff((al0-al1)/2,(al1-al0)/2)
    ls <- sum(log(A))
    B <- g1^2 - g1 + g0^2 - g0 
    C <- 2 * g1^3 - 3 * g1^2 + g1 - 2 * g0^3 + 3 * g0^2 - g0
    Da0 <- .5 * B/A ; Da1 <- -0.5 *B/A 
    Da0a0 <- .25 * C/A - .25 * B^2/A^2
    Da1a1 <- .25 * C/A - .25 * B^2/A^2
    Da0a1 <- - .25 * C/A + .25 * B^2/A^2
    i <- 0 
    n2d <- (R-2)*(R-1)/2
    n <- length(y)
    Dth <- matrix(0,n,R-2)
    Dth2 <- matrix(0,n,n2d)
    for (j in 1:(R-2)) for (k in j:(R-2)) { 
      i <- i + 1 ## the second deriv col
      ind <- y >= j ## rest are zero
      ar1.k <- ar.k <- rep(exp(theta[k]),n)
      ar.k[y==R | y <= k] <- 0; ar1.k[y<k+2] <- 0
      ar.j <- ar1.j <- rep(exp(theta[j]),n)
      ar.j[y==R | y <= j] <- 0; ar1.j[y<j+2] <- 0
      ar.kj <- ar1.kj <- rep(0,n)
      if (k==j) {
        ar.kj[y>k&y<R] <- exp(theta[k])
        ar1.kj[y>k+1] <- exp(theta[k])
        Dth[ind,k] <- Da1[ind]*ar.k[ind]  + Da0[ind]*ar1.k[ind]
      }
      Dth2[,i] <- Da1a1*ar.k*ar.j + Da0a1*ar.k*ar1.j + Da1 * ar.kj +
                  Da0a0*ar1.k*ar1.j + Da0a1*ar1.k*ar.j + Da0 * ar1.kj
    } 
    lsth2=colSums(Dth2)
    if (R>2) {
      ls2 <- matrix(0,R-2,R-2);ii <- 0
      for (i in 1:(R-2)) for (j in i:(R-2)) { 
        ii <- ii + 1 
        ls2[i,j] <- ls2[j,i] <- lsth2[ii]
      } 
    }
    list(ls=ls,lsth1=colSums(Dth),lsth2=ls2)
  } ## end of ls
  
  ## initialization is interesting -- needs to be with reference to initial cut-points
  ## so that mu puts each obs in correct category initially...
 
  preinitialize <- expression({
    ocat.ini <- function(R,y) {
    ## initialize theta from raw counts in each category
      if (R<3) return
      p <- cumsum(tabulate(y[is.finite(y)])/length(y[is.finite(y)]))
      eta <- if (p[1]==0) 5 else -1 - log(p[1]/(1-p[1])) ## mean of latent
      theta <- rep(-1,R-1) 
      for (i in 2:(R-1)) theta[i] <- log(p[i]/(1-p[i])) + eta 
      theta <- diff(theta)
      theta[theta <= 0.01] <- 0.01
      theta <- log(theta)
    }
    R3 <- length(G$family$getTheta())+2
    if (R3>2&&G$family$n.theta>0) { 
      Theta <- ocat.ini(R3,G$y)
      G$family$putTheta(Theta)
    } 
  })

  initialize <- expression({ 
    R <- family$n.theta+2
    if (any(y < 1)||any(y>R)) stop("values out of range")
    n <- rep(1, nobs)
    ## get the cut points... 
    alpha <- rep(0,R+1) ## the thresholds
    alpha[1] <- -2;alpha[2] <- -1
    if (R > 2) { 
      ind <- 3:R
      alpha[ind] <- alpha[2] + cumsum(exp(family$getTheta()))
    } 
    alpha[R+1] <- alpha[R] + 1
    mustart <- (alpha[y+1] + alpha[y])/2
  })

  residuals <- function(object,type=c("deviance","working","response")) {
    if (type == "working") { 
      res <- object$residuals 
    } else if (type == "response") {
       theta <- object$family$getTheta()
       mu <- object$linear.predictors
       R = length(theta)+2
       alpha <- rep(0,R+1) ## the thresholds
       alpha[1] <- -Inf;alpha[R+1] <- Inf
       alpha[2] <- -1
       if (R > 2) { 
         ind <- 3:R
         alpha[ind] <- alpha[2] + cumsum(exp(theta))
       } 
       fv <- mu*NA
       for (i in 1:(R+1)) {
         ind <- mu>alpha[i] & mu<=alpha[i+1]
         fv[ind] <- i
       } 
       res <- object$y - fv
    } else if (type == "deviance") { 
      y <- object$y
      mu <- object$fitted.values
      wts <- object$prior.weights
      res <- object$family$dev.resids(y,mu,wts)
      s <- attr(res,"sign")
      if (is.null(s)) s <- sign(y-mu)
      res <- as.numeric(sqrt(pmax(res,0)) * s) 
      
    }
    res
  } ## residuals


  #fv <- function(mu,se=NULL,theta=NULL,predict=FALSE) {
  ## optional function to give fitted values - idea is that 
  ## predict.gam(...,type="response") will use this, as well
  ## as residuals(...,type="response").
  ## predict == TRUE indicates prediction call.
 
  predict <- function(family,se=FALSE,eta=NULL,y=NULL,X=NULL,
                beta=NULL,off=NULL,Vb=NULL,family.data=NULL) {
  ## optional function to give predicted values - idea is that 
  ## predict.gam(...,type="response") will use this, and that
  ## either eta will be provided, or {X, beta, off, Vb}. family.data
  ## contains any family specific extra information. 
    ocat.prob <- function(theta,lp,se=NULL) {
    ## compute probabilities for each class in ocat model
    ## theta is finite cut points, lp is linear predictor, se 
    ## is standard error on lp...
      R <- length(theta) 
      dp <- prob <- matrix(0,length(lp),R+2)
      prob[,R+2] <- 1
      for (i in 1:R) {
        x <- theta[i] - lp
        ind <- x > 0
        prob[ind,i+1] <- 1/(1+exp(-x[ind]))
        ex <- exp(x[!ind])
        prob[!ind,i+1] <- ex/(1+ex)
        dp[,i+1] <- prob[,i+1]*(prob[,i+1]-1)
      }
      prob <- t(diff(t(prob)))
      dp <- t(diff(t(dp))) ## dprob/deta
      if (!is.null(se)) se <- as.numeric(se)*abs(dp)
      list(prob,se)
    } ## ocat.prob

    theta <- family$getTheta(TRUE)
    if (is.null(eta)) { ## return probabilities
      mu <- X%*%beta + off 
      se <- if (se) sqrt(rowSums((X%*%Vb)*X)) else NULL
      ##theta <- cumsum(c(-1,exp(theta)))
      p <- ocat.prob(theta,mu,se)
      if (is.null(se)) return(p) else { ## approx se on prob also returned
        names(p) <- c("fit","se.fit")
        return(p)
      } 
    } else { ## return category implied by eta (i.e mean of latent)
      R = length(theta)+2
      alpha <- rep(0,R) ## the thresholds
      alpha[1] <- -Inf;alpha[R] <- Inf
      fv <- eta*NA
      for (i in 1:(R+1)) {
        ind <- eta>alpha[i] & eta<=alpha[i+1]
        fv[ind] <- i
      } 
      return(fv)
    }
  } ## fv

  rd <- function(mu,wt,scale) {
    ## simulate data given fitted latent variable in mu 
    theta <- get(".Theta") 
    R = length(theta)+2
    alpha <- rep(0,R+1) ## the thresholds
    alpha[1] <- -Inf;alpha[R+1] <- Inf
    alpha[2] <- -1
    if (R > 2) { 
      ind <- 3:R
      alpha[ind] <- alpha[2] + cumsum(exp(theta))
    } 
    ## ... cut points computed, now simulate latent variable, u
    y <- u <- runif(length(mu))
    u <- mu + log(u/(1-u)) 
    ## and allocate categories according to u and cut points...
    for (i in 1:R) {
      y[u > alpha[i]&u <= alpha[i+1]] <- i
    }
    y
  }

  environment(dev.resids) <- environment(aic) <- environment(putTheta) <-
  environment(getTheta) <- environment(rd) <- environment(predict) <- env
  structure(list(family = "Ordered Categorical", link = linktemp, linkfun = stats$linkfun,
        linkinv = stats$linkinv, dev.resids = dev.resids,Dd=Dd,
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize,postproc=postproc,
        preinitialize = preinitialize, ls=ls,rd=rd,residuals=residuals,
        validmu = validmu, valideta = stats$valideta,n.theta=n.theta,
        ini.theta = iniTheta,putTheta=putTheta,predict=predict,step = 1,
        getTheta=getTheta,no.r.sq=TRUE), class = c("extended.family","family"))
} ## end of ocat

#######################
## negative binomial...
#######################

nb <- function (theta = NULL, link = "log") { 
## Extended family object for negative binomial, to allow direct estimation of theta
## as part of REML optimization. Currently the template for extended family objects.
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
  Theta <-  NULL;n.theta <- 1
  if (!is.null(theta)&&theta!=0) {
      if (theta>0) { 
        iniTheta <- Theta <- log(theta) ## fixed theta supplied
        n.theta <- 0 ## signal that there are no theta parameters to estimate
      } else iniTheta <- log(-theta) ## initial theta supplied
  } else iniTheta <- 0 ## inital log theta value
    
    env <- new.env(parent = .GlobalEnv)
    assign(".Theta", iniTheta, envir = env)
    getTheta <- function(trans=FALSE) if (trans) exp(get(".Theta")) else get(".Theta") # get(".Theta")
    putTheta <- function(theta) assign(".Theta", theta,envir=environment(sys.function()))

    variance <- function(mu) mu + mu^2/exp(get(".Theta"))

    validmu <- function(mu) all(mu > 0)

    dev.resids <- function(y, mu, wt,theta=NULL) {
      if (is.null(theta)) theta <- get(".Theta")
      theta <- exp(theta) ## note log theta supplied
      2 * wt * (y * log(pmax(1, y)/mu) - 
        (y + theta) * log((y + theta)/(mu + theta))) 
    }
    
    Dd <- function(y, mu, theta, wt, level=0) {
    ## derivatives of the deviance...
      ltheta <- theta
      theta <- exp(theta)
      yth <- y + theta
      muth <- mu + theta
      r <- list()
      ## get the quantities needed for IRLS. 
      ## Dmu2eta2 is deriv of D w.r.t mu twice and eta twice,
      ## Dmu is deriv w.r.t. mu once, etc...
      r$Dmu <- 2 * wt * (yth/muth - y/mu)
      r$Dmu2 <- -2 * wt * (yth/muth^2 - y/mu^2)
      r$EDmu2 <- 2 * wt * (1/mu - 1/muth) ## exact (or estimated) expected weight
      if (level>0) { ## quantities needed for first derivatives
        r$Dth <- -2 * wt * theta * (log(yth/muth) + (1 - yth/muth) ) 
        r$Dmuth <- 2 * wt * theta * (1 - yth/muth)/muth
        r$Dmu3 <- 4 * wt * (yth/muth^3 - y/mu^3)
        r$Dmu2th <- 2 * wt * theta * (2*yth/muth - 1)/muth^2
      } 
      if (level>1) { ## whole damn lot
        r$Dmu4 <- 2 * wt * (6*y/mu^4 - 6*yth/muth^4)
        r$Dth2 <- -2 * wt * theta * (log(yth/muth) +
                     theta*yth/muth^2 - yth/muth - 2*theta/muth + 1 +
                     theta /yth)
        r$Dmuth2 <- 2 * wt * theta * (2*theta*yth/muth^2 - yth/muth - 2*theta/muth + 1)/muth
        r$Dmu2th2 <- 2 * wt * theta * (- 6*yth*theta/muth^2 + 2*yth/muth + 4*theta/muth - 1) /muth^2
        r$Dmu3th <- 4 * wt * theta * (1 - 3*yth/muth)/muth^3
      }
      r
    }

    aic <- function(y, mu, theta=NULL, wt, dev) {
        if (is.null(theta)) theta <- get(".Theta")
        Theta <- exp(theta)
        term <- (y + Theta) * log(mu + Theta) - y * log(mu) +
            lgamma(y + 1) - Theta * log(Theta) + lgamma(Theta) -
            lgamma(Theta + y)
        2 * sum(term * wt)
    }
    
    ls <- function(y,w,n,theta,scale) {
       ## the log saturated likelihood function.
       Theta <- exp(theta)
       ylogy <- y;ind <- y>0;ylogy[ind] <- y[ind]*log(y[ind])
       term <- (y + Theta) * log(y + Theta) - ylogy +
            lgamma(y + 1) - Theta * log(Theta) + lgamma(Theta) -
            lgamma(Theta + y)
       ls <- -sum(term*w)
       ## first derivative wrt theta...
       yth <- y+Theta
       lyth <- log(yth)
       psi0.yth <- digamma(yth) 
       psi0.th <- digamma(Theta)
       term <- Theta * (lyth - psi0.yth + psi0.th-theta)
       lsth <- -sum(term*w)
       ## second deriv wrt theta...
       psi1.yth <- trigamma(yth) 
       psi1.th <- trigamma(Theta)
       term <- Theta * (lyth - Theta*psi1.yth - psi0.yth + Theta/yth + Theta * psi1.th + psi0.th - theta -1)        
       lsth2 <- -sum(term*w)
       list(ls=ls, ## saturated log likelihood
            lsth1=lsth, ## first deriv vector w.r.t theta - last element relates to scale, if free
            lsth2=lsth2) ## Hessian w.r.t. theta, last row/col relates to scale, if free
    }

    initialize <- expression({
        if (any(y < 0)) stop("negative values not allowed for the negative binomial family")
        n <- rep(1, nobs)
        mustart <- y + (y == 0)/6
    })
  
    postproc <- expression({
      object$family$family <- 
      paste("Negative Binomial(",round(object$family$getTheta(TRUE),3),")",sep="")
    })

    rd <- function(mu,wt,scale) {
      Theta <- exp(get(".Theta"))
      rnbinom(mu,size=Theta,mu=mu)
    }

    qf <- function(p,mu,wt,scale) {
      Theta <- exp(get(".Theta"))
      qnbinom(p,size=Theta,mu=mu)
    }
 

     environment(dev.resids) <- environment(aic) <- environment(getTheta) <- 
     environment(rd)<- environment(qf)<- environment(variance) <- environment(putTheta) <- env
#    famname <- paste("Negative Binomial(", format(round(theta,3)), ")", sep = "")
    structure(list(family = "negative binomial", link = linktemp, linkfun = stats$linkfun,
        linkinv = stats$linkinv, dev.resids = dev.resids,Dd=Dd,variance=variance,
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize,postproc=postproc,ls=ls,
        validmu = validmu, valideta = stats$valideta,n.theta=n.theta, 
        ini.theta = iniTheta,putTheta=putTheta,getTheta=getTheta,rd=rd,qf=qf),
        class = c("extended.family","family"))
} ## nb


## Tweedie....


tw <- function (theta = NULL, link = "log",a=1.01,b=1.99) { 
## Extended family object for Tweedie, to allow direct estimation of p
## as part of REML optimization. 
## p = (a+b*exp(theta))/(1+exp(theta)), i.e. a < p < b
## NOTE: The Tweedie density computation at low phi, low p is susceptible
##       to cancellation error, which seems unavoidable. Furthermore 
##       there are known problems with spurious maxima in the likelihood 
##       w.r.t. p when the data are rounded. 
  
  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  if (linktemp %in% c("log", "identity", "sqrt","inverse")) stats <- make.link(linktemp)
  else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  } else {
    if (inherits(link, "link-glm")) {
       stats <- link
            if (!is.null(stats$name))
                linktemp <- stats$name
        }
        else  stop(gettextf("link \"%s\" not available for Tweedie family.", 
                linktemp, collapse = ""), domain = NA)
  }
  Theta <-  NULL;n.theta <- 1
  if (!is.null(theta)&&theta!=0) {
      if (abs(theta)<=a||abs(theta)>=b) stop("Tweedie p must be in interval (a,b)")
      if (theta>0) { ## fixed theta supplied
        iniTheta <- Theta <- log((theta-a)/(b-theta)) 
        n.theta <- 0 ## so no theta to estimate
      } else iniTheta <- log((-theta-a)/(b+theta)) ## initial theta supplied
  } else iniTheta <- 0 ## inital log theta value
    
  env <- new.env(parent = .GlobalEnv)
  assign(".Theta", iniTheta, envir = env) 
  assign(".a",a, envir = env);assign(".b",b, envir = env)
  getTheta <- function(trans=FALSE) { 
  ## trans transforms to the original scale...
    th <- get(".Theta")
    a <- get(".a");b <- get(".b")
    #a <- 1.01;b <- 1.99
    if (trans) th <- if (th>0) (b+a*exp(-th))/(1+exp(-th)) else (b*exp(th)+a)/(exp(th)+1)
    th
  }
  putTheta <- function(theta) assign(".Theta", theta,envir=environment(sys.function()))

  validmu <- function(mu) all(mu > 0)
  
  variance <- function(mu) { 
    th <- get(".Theta");a <- get(".a");b <- get(".b")
    #a <- 1.01;b <- 1.99
    p <- if (th>0) (b+a*exp(-th))/(1+exp(-th)) else (b*exp(th)+a)/(exp(th)+1)
    mu^p
  }
  
  dev.resids <- function(y, mu, wt,theta=NULL) {
    if (is.null(theta)) theta <- get(".Theta")
    #p <- if (theta>0) (exp(-theta)+2)/(exp(-theta)+1) else (1 + 2*exp(theta))/(1 + exp(theta)) 
    #a <- 1.01;b <- 1.99
    a <- get(".a");b <- get(".b")
    p <- if (theta>0) (b+a*exp(-theta))/(1+exp(-theta)) else (b*exp(theta)+a)/(exp(theta)+1)
    y1 <- y + (y == 0)
    theta <- if (p == 1) log(y1/mu) else (y1^(1 - p) - mu^(1 - p))/(1 - p)
    kappa <- if (p == 2) log(y1/mu) else (y^(2 - p) - mu^(2 - p))/(2 - p)
    2 * (y * theta - kappa) * wt
  }
    
  Dd <- function(y, mu, theta, wt, level=0) {
  ## derivatives of the deviance...
    #a <- 1.01;b <- 1.99
    a <- get(".a");b <- get(".b")
    th <- theta
    p <- if (th>0) (b+a*exp(-th))/(1+exp(-th)) else (b*exp(th)+a)/(exp(th)+1)
    dpth1 <- if (th>0) exp(-th)*(b-a)/(1+exp(-th))^2 else exp(th)*(b-a)/(exp(th)+1)^2
    dpth2 <- if (th>0) ((a-b)*exp(-th)+(b-a)*exp(-2*th))/(exp(-th)+1)^3 else
                   ((a-b)*exp(2*th)+(b-a)*exp(th))/(exp(th)+1)^3
    #p <- if (theta>0) (exp(-theta)+2)/(exp(-theta)+1) else (1 + 2*exp(theta))/(1 + exp(theta))
    #g <- if (theta>0) 1/(1 + exp(-theta)) else exp(theta)/(1+exp(theta)) 
    #dpth1 <- g*(2-p) ## d p / d th
    #dpth2 <- 2*g^2*(p-2) + dpth1  ## d^2 p / d th^2 
    mu1p <- mu^(1-p)
    mup <- mu^p
    r <- list()
    ## get the quantities needed for IRLS. 
    ## Dmu2eta2 is deriv of D w.r.t mu twice and eta twice,
    ## Dmu is deriv w.r.t. mu once, etc...
    ymupi <- y/mup
    r$Dmu <- 2*wt*(mu1p - ymupi)
    r$Dmu2 <- 2*wt*(mu^(-1-p)*p*y + (1-p)/mup)
    r$EDmu2 <- (2*wt)/mup ## expected Dmu2 (weight)
   
    if (level>0) { ## quantities needed for first derivatives
        i1p <- 1/(1-p)
        y1 <- y + (y==0)
        ylogy <- y*log(y1)
        logmu <- log(mu)
        mu2p <- mu * mu1p
        r$Dth <- 2 * wt * ( (y^(2-p)*log(y1) - mu2p*logmu)/(2-p) + 
                            (y*mu1p*logmu - y^(2-p)*log(y1))/(1-p) -
                            (y^(2-p) - mu2p)/(2-p)^2 + 
                            (y^(2-p) - y*mu1p)*i1p^2) *dpth1       

        r$Dmuth <- 2 * wt * logmu * (ymupi - mu1p)*dpth1
        mup1 <-  mu^(-p-1)
        r$Dmu3 <- -2 * wt * mup1*p*(y/mu*(p+1) + 1-p)    
        r$Dmu2th <- 2 * wt  * (mup1*y*(1-p*logmu)-(logmu*(1-p)+1)/mup )*dpth1
      } 
      if (level>1) { ## whole damn lot
        mup2 <- mup1/mu
        r$Dmu4 <- 2 * wt * mup2*p*(p+1)*(y*(p+2)/mu + 1 - p)
        y2plogy <- y^(2-p)*log(y1);y2plog2y <- y2plogy*log(y1)
       
        r$Dth2 <- 2 * wt * (((mu2p*logmu^2-y2plog2y)/(2-p) + (y2plog2y - y*mu1p*logmu^2)/(1-p) +
                             2*(y2plogy-mu2p*logmu)/(2-p)^2 + 2*(y*mu1p*logmu-y2plogy)/(1-p)^2
                             + 2 * (mu2p - y^(2-p))/(2-p)^3+2*(y^(2-p)-y*mu^(1-p))/(1-p)^3)*dpth1^2) +
                              r$Dth*dpth2/dpth1
        
        r$Dmuth2 <- 2 * wt * ((mu1p * logmu^2 - logmu^2*ymupi)*dpth1^2) + r$Dmuth*dpth2/dpth1

        r$Dmu2th2 <- 2 * wt * ( (mup1 * logmu*y*(logmu*p - 2) + logmu/mup*(logmu*(1-p) + 2))*dpth1^2) +
                                r$Dmu2th * dpth2/dpth1

        r$Dmu3th <- 2 * wt * mup1*(y/mu*(logmu*(1+p)*p-p -p-1) +logmu*(1-p)*p + p - 1 + p)*dpth1
      }
      r
    }

    aic <- function(y, mu, theta=NULL, wt, dev) {
        if (is.null(theta)) theta <- get(".Theta")  
        #a <- 1.01;b <- 1.99
        a <- get(".a");b <- get(".b")
        p <- if (theta>0) (b+a*exp(-theta))/(1+exp(-theta)) else (b*exp(theta)+a)/(exp(theta)+1)
        scale <- dev/sum(wt)
        -2 * sum(ldTweedie(y, mu, p = p, phi = scale)[, 1] * 
            wt) + 2
    }

    ls <- function(y, w, n, theta, scale) {
        ## evaluate saturated log likelihood + derivs w.r.t. working params and log(scale)
        #p <- if (theta>0) (exp(-theta)+2)/(exp(-theta)+1) else (1 + 2*exp(theta))/(1 + exp(theta))
        #g <- if (theta>0) 1/(1 + exp(-theta)) else exp(theta)/(1+exp(theta)) 
        #dpth1 <- g*(2-p) ## d p / d th
        #dpth2 <- 2*g^2*(p-2) + dpth1  ## d^2 p / d th^2 
        #a <- 1.01;b <- 1.99;th <- theta
        a <- get(".a");b <- get(".b")
        #p <- if (th>0) (b+a*exp(-th))/(1+exp(-th)) else (b*exp(th)+a)/(exp(th)+1)
        #dpth1 <- if (th>0) exp(-th)*(b-a)/(1+exp(-th))^2 else exp(th)*(b-a)/(exp(th)+1)^2
        #dpth2 <- if (th>0) ((a-b)*exp(-th)+(b-a)*exp(-2*th))/(exp(-th)+1)^3 else
        #           ((a-b)*exp(2*th)+(b-a)*exp(th))/(exp(th)+1)^3
        #LS <- colSums(w * ldTweedie(y, y, p = p, phi = scale))
        LS <- colSums(w * ldTweedie(y, y, rho=log(scale), theta=theta,a=a,b=b))
        #lsth1 <- c(LS[4]*dpth1,LS[2]*scale)
        #d2p <- LS[5]*dpth1^2+LS[4]*dpth2
        #dmix <- LS[6]*dpth1*scale
        #lsth2 <- matrix(c(d2p,dmix,dmix,LS[3]*scale^2+LS[2]*scale),2,2)
        lsth1 <- c(LS[4],LS[2])
        lsth2 <- matrix(c(LS[5],LS[6],LS[6],LS[3]),2,2)
        list(ls=LS[1],lsth1=lsth1,lsth2=lsth2)
    }

 
    initialize <- expression({
        n <- rep(1, nobs)
        mustart <- y + (y == 0)*.1
    })
    postproc <- expression({
      object$family$family <- 
      paste("Tweedie(p=",round(object$family$getTheta(TRUE),3),")",sep="")
    })
  
    rd <- function(mu,wt,scale) {
     th <- get(".Theta") 
     #a <- 1.01;b <- 1.99 
     a <- get(".a");b <- get(".b")
     p <- if (th>0) (b+a*exp(-th))/(1+exp(-th)) else (b*exp(th)+a)/(exp(th)+1)

     if (p == 2) 
            rgamma(mu, shape = 1/scale, scale = mu * scale)
     else
            rTweedie(mu, p = p, phi = scale)
    }

     environment(Dd) <- environment(ls) <-
     environment(dev.resids) <- environment(aic) <- environment(getTheta) <- 
      environment(rd) <- environment(variance) <- environment(putTheta) <- env
#    famname <- paste("Negative Binomial(", format(round(theta,3)), ")", sep = "")
    structure(list(family = "Tweedie", link = linktemp, linkfun = stats$linkfun,
        linkinv = stats$linkinv, dev.resids = dev.resids,Dd=Dd,variance=variance,rd=rd,
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize,postproc=postproc,ls=ls,
        validmu = validmu, valideta = stats$valideta,canonical="none",n.theta=n.theta, 
        ini.theta = iniTheta,putTheta=putTheta,getTheta=getTheta,scale = -1),
        class = c("extended.family","family"))
} ## tw

##################################
## Natalya Pya code from here....
##################################

betar <- function (theta = NULL, link = "logit") { 
## Extended family object for beta regression
## length(theta)=1; log theta supplied
  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  if (linktemp %in% c("logit", "probit", "cloglog", "cauchit", "log")) stats <- make.link(linktemp)
  else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  } else {
    if (inherits(link, "link-glm")) {
       stats <- link
            if (!is.null(stats$name))
                linktemp <- stats$name
        }
        else stop(linktemp, " link not available for beta regression; available links are  \"logit\", \"probit\", \"cloglog\", \"log\" and \"cauchit\"")
    }
   
    Theta <-  NULL; n.theta <- 1
    if (!is.null(theta)&&theta!=0) {
       if (theta>0) {
           iniTheta <- Theta <- log(theta) ## fixed theta supplied
           n.theta <- 0 ## signal that there are no theta parameters to estimate
       } else iniTheta <- log(-theta) ## initial theta supplied
    } else iniTheta <- 0 ##  inital log theta value
    
    env <- new.env(parent = .GlobalEnv)
    assign(".Theta", iniTheta, envir = env)
    getTheta <- function(trans=FALSE) if (trans) exp(get(".Theta")) else get(".Theta") 
    putTheta <- function(theta) assign(".Theta", theta,envir=environment(sys.function()))

    variance <- function(mu) { 
        th <- get(".Theta")
        mu*(1 - mu)/(1+exp(th))
    }
  
    validmu <- function(mu) all(mu > 0 & mu < 1)

    dev.resids <- function(y, mu, wt,theta=NULL) {
    ## '-'2*loglik instead of deviance in REML/ML expression
      if (is.null(theta)) theta <- get(".Theta")
      theta <- exp(theta) ## note log theta supplied
      muth <- mu*theta
      yth <- y*theta
      2* wt * (-lgamma(theta) +lgamma(muth) + lgamma(theta - muth) - muth*log(y/(1-y)) - theta*log(1-y) + log(y*(1-y))) 
    }
    
    Dd <- function(y, mu, theta, wt, level=0) {
    ## derivatives of the -2*loglik...
      ltheta <- theta
      theta <- exp(theta)
      onemu <- 1 - mu;  oney <- 1 - y
      muth <- mu*theta; yth <- y*theta
      onemuth <- onemu*theta  ## (1-mu)*theta
      psi0.th <- digamma(theta)
      psi1.th <- trigamma(theta)
      psi0.muth <- digamma(muth) 
      psi0.onemuth <- digamma(onemuth)
      psi1.muth <- trigamma(muth)
      psi1.onemuth <- trigamma(onemuth)
      psi2.muth <- psigamma(muth,2)
      psi2.onemuth <- psigamma(onemuth,2)
      psi3.muth <- psigamma(muth,3)
      psi3.onemuth <- psigamma(onemuth,3)
      log.yoney <- log(y)-log(oney)
      r <- list()
      ## get the quantities needed for IRLS. 
      ## Dmu2eta2 is deriv of D w.r.t mu twice and eta twice,
      ## Dmu is deriv w.r.t. mu once, etc...
      r$Dmu <- 2 * wt * theta* (psi0.muth - psi0.onemuth - log.yoney)
      r$Dmu2 <- 2 * wt * theta^2*(psi1.muth+psi1.onemuth)
      r$EDmu2 <- r$Dmu2
      if (level>0) { ## quantities needed for first derivatives
        r$Dth <- 2 * wt *theta*(-mu*log.yoney - log(oney)+ mu*psi0.muth+onemu*psi0.onemuth -psi0.th) 
        r$Dmuth <- r$Dmu + 2 * wt * theta^2*(mu*psi1.muth -onemu*psi1.onemuth)
        r$Dmu3 <- 2 * wt *theta^3 * (psi2.muth - psi2.onemuth) 
        r$Dmu2th <- 2* r$Dmu2 + 2 * wt * theta^3* (mu*psi2.muth + onemu*psi2.onemuth)
      } 
      if (level>1) { ## whole lot
        r$Dmu4 <- 2 * wt *theta^4 * (psi3.muth+psi3.onemuth) 
        r$Dth2 <- r$Dth +2 * wt *theta^2* (mu^2*psi1.muth+ onemu^2*psi1.onemuth-psi1.th)
        r$Dmuth2 <- r$Dmuth + 2 * wt *theta^2* (mu^2*theta*psi2.muth+ 2*mu*psi1.muth -
                    theta*onemu^2*psi2.onemuth - 2*onemu*psi1.onemuth)
        r$Dmu2th2 <- 2*r$Dmu2th + 2* wt * theta^3* (mu^2*theta*psi3.muth +3*mu*psi2.muth+ 
                    onemu^2*theta*psi3.onemuth + 3*onemu*psi2.onemuth )
        r$Dmu3th <- 3*r$Dmu3 + 2 * wt *theta^4*(mu*psi3.muth-onemu*psi3.onemuth)
      }
      r
    }

    aic <- function(y, mu, theta=NULL, wt, dev) {
        if (is.null(theta)) theta <- get(".Theta")
        theta <- exp(theta)
        muth <- mu*theta
        term <- -lgamma(theta)+lgamma(muth)+lgamma(theta-muth)-(muth-1)*log(y)-
               (theta-muth-1)*log(1-y) ## `-' log likelihood for each observation
        2 * sum(term * wt)
    }
    
    ls <- function(y,w,n,theta,scale) {
       ## the log saturated likelihood function.
       ## ls is defined as zero for REML/ML expression as deviance is defined as -2*log.lik 
       list(ls=0,## saturated log likelihood
            lsth1=0,  ## first deriv vector w.r.t theta - last element relates to scale
            lsth2=0) ##Hessian w.r.t. theta
     }

   
    ## preinitialization to reset G$y values of <=0 and >=1... 
    preinitialize <- expression({
     ## code to evaluate in estimate.gam...
     ## reset G$y values of <=0 and >= 1 to eps and 1-eps...
     # G$family.data <- list()
      eps <- 1e-7 
      G$y[G$y >= 1-eps] <- 1 - eps
      G$y[G$y<= eps] <- eps
    })

    saturated.ll <- function(y,wt,theta=NULL){
    ## function to find the saturated loglik by Newton method,
    ## searching for the mu (on logit scale) that max loglik given theta and data...

      gbh <- function(y,eta,phi,deriv=FALSE,a=1e-8,b=1-a) {
      ## local function evaluating log likelihood (l), gradient and second deriv 
      ## vectors for beta... a and b are min and max mu values allowed. 
      ## mu = (a + b*exp(eta))/(1+exp(eta))
        ind <- eta>0
        expeta <- mu <- eta;
        expeta[ind] <- exp(-eta[ind]);expeta[!ind] <- exp(eta[!ind])
        mu[ind] <- (a*expeta[ind] + b)/(1+expeta[ind])
        mu[!ind] <- (a + b*expeta[!ind])/(1+expeta[!ind])
        l <- dbeta(y,phi*mu,phi*(1-mu),log=TRUE)
        if (deriv) {
          g <- phi * log(y) - phi * log(1-y) - phi * digamma(mu*phi) + phi * digamma((1-mu)*phi)
          h <- -phi^2*(trigamma(mu*phi)+trigamma((1-mu)*phi))
          dmueta2 <- dmueta1 <-  eta
          dmueta1 <- expeta*(b-a)/(1+expeta)^2 
          dmueta2 <- sign(eta)* ((a-b)*expeta+(b-a)*expeta^2)/(expeta+1)^3
          h <- h * dmueta1^2 + g * dmueta2
          g <- g * dmueta1
        } else g=h=NULL
        list(l=l,g=g,h=h,mu=mu)
      } ## gbh 
      ## now Newton loop...
      eps <- 1e-7
      eta <- y
      a <- eps;b <- 1 - eps
      y[y<eps] <- eps;y[y>1-eps] <- 1-eps
      eta[y<=eps*1.2] <- eps *1.2
      eta[y>=1-eps*1.2] <- 1-eps*1.2
      eta <- log((eta-a)/(b-eta)) 
      mu <- LS <- ii <- 1:length(y)
      for (i in 1:200) {
        ls <- gbh(y,eta,theta,TRUE)
        conv <- abs(ls$g)<mean(abs(ls$l)+.1)*1e-8
        if (sum(conv)>0) { ## some convergences occured
          LS[ii[conv]] <- ls$l[conv] ## store converged
          mu[ii[conv]] <- ls$mu[conv] ## store mu at converged
          ii <- ii[!conv] ## drop indices
          if (length(ii)>0) { ## drop the converged
            y <- y[!conv];eta <- eta[!conv]
            ls$l <- ls$l[!conv];ls$g <- ls$g[!conv];ls$h <- ls$h[!conv]
          } else break ## nothing left to do
        }
        h <- -ls$h
        hmin <- max(h)*1e-4 
        h[h<hmin] <- hmin ## make +ve def
        delta <- ls$g/h   ## step
        ind <- abs(delta)>2
        delta[ind] <- sign(delta[ind])*2 ## step length limit
        ls1 <- gbh(y,eta+delta,theta,FALSE); ## did it work?
        ind <- ls1$l<ls$l ## failure index
        k <- 0
        while (sum(ind)>0&&k<20) { ## step halve only failed steps
          k <- k + 1
          delta[ind] <- delta[ind]/2
          ls1$l[ind] <- gbh(y[ind],eta[ind]+delta[ind],theta,FALSE)$l
          ind <- ls1$l<ls$l
        }
        eta <- eta + delta
      } ## end newton loop
      if (length(ii)>0) { 
        LS[ii] <- ls$l
        warning("saturated likelihood may be inaccurate")
      }

      list(f=sum(wt*LS),term=LS,mu=mu) ## fields f (sat lik) and term (individual datum sat lik) expected
    } ## saturated.ll


    postproc <- expression({
    ## code to evaluate in estimate.gam, to find the saturated
    ## loglik by Newton method
    ## searching for the mu (on logit scale) that max loglik given theta...
      wts <- object$prior.weights
      theta <- object$family$getTheta(trans=TRUE) ## exp theta
      lf <- object$family$saturated.ll(G$y, wts,theta)
      ## storing the saturated loglik for each datum...
      object$family.data <- list(ls = lf$term,mu.ls = lf$mu)   
      l2 <- object$family$dev.resids(G$y,object$fitted.values,wts)
      object$deviance <- 2*lf$f + sum(l2)
      wtdmu <- if (G$intercept) sum(wts * G$y)/sum(wts) 
              else object$family$linkinv(G$offset)
      object$null.deviance <- 2*lf$f + sum(object$family$dev.resids(G$y, wtdmu, wts))
      object$family$family <- 
      paste("Beta regression(",round(theta,3),")",sep="")
    })

    initialize <- expression({
        n <- rep(1, nobs)
        mustart <- y 
    })

    residuals <- function(object,type=c("deviance","working","response","pearson")) {
      if (type == "working") { 
        res <- object$residuals 
      } else if (type == "response") {
        res <- object$y - object$fitted.values
      } else if (type == "deviance") { 
        y <- object$y
        mu <- object$fitted.values
        wts <- object$prior.weights
        sim <- attr(y,"simula")
     #   if (!is.null(sim)) {  ## if response values simulated, Newton search called to get saturated log.lik
           lf <- object$family$saturated.ll(y, wts,object$family$getTheta(TRUE))
           object$family.data$ls <- lf$term  
     #   }
        res <- 2*object$family.data$ls + object$family$dev.resids(y,mu,wts)
        res[res<0] <- 0
        s <- sign(y-mu)
        res <- sqrt(res) * s   
      } else if (type == "pearson") {
        mu <- object$fitted.values
        res <- (object$y - mu)/object$family$variance(mu)^.5
      }
      res
     } ## residuals

    rd <- function(mu,wt,scale) {
     ## simulate data given fitted latent variable in mu 
      Theta <- exp(get(".Theta"))
      r <- rbeta(mu,shape1=Theta*mu,shape2=Theta*(1-mu))
      eps <- 1e-7 ;
      r[r>=1-eps] <- 1 - eps
      r[r<eps] <- eps
      r
    }

    qf <- function(p,mu,wt,scale) {
      Theta <- exp(get(".Theta"))
      q <- qbeta(p,shape1=Theta*mu,shape2=Theta*(1-mu))
      eps <- 1e-7 ;
      q[q>=1-eps] <- 1 - eps
      q[q<eps] <- eps
      q
    }

    environment(dev.resids) <- environment(aic) <- environment(getTheta) <- 
    environment(rd)<- environment(qf) <- environment(variance) <- environment(putTheta) <-
    environment(saturated.ll) <- env

    structure(list(family = "Beta regression", link = linktemp, linkfun = stats$linkfun,
        linkinv = stats$linkinv, dev.resids = dev.resids,Dd=Dd, variance=variance,
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize,ls=ls,
        preinitialize=preinitialize,postproc=postproc, residuals=residuals, saturated.ll=saturated.ll,
        validmu = validmu, valideta = stats$valideta, n.theta=n.theta,  
        ini.theta = iniTheta,putTheta=putTheta,getTheta=getTheta,rd=rd,qf=qf), 
        class = c("extended.family","family"))
} ## betar


  
## scaled t ...

scat <- function (theta = NULL, link = "identity") { 
## Extended family object for scaled t distribution
## length(theta)=2; log theta supplied
  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  if (linktemp %in% c("identity", "log", "inverse")) stats <- make.link(linktemp)
  else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  } else {
    if (inherits(link, "link-glm")) {
       stats <- link
            if (!is.null(stats$name))
                linktemp <- stats$name
        }
        else stop(linktemp, " link not available for scaled t distribution; available links are \"identity\", \"log\",  and \"inverse\"")
    }
    Theta <-  NULL;n.theta <- 2
    if (!is.null(theta)&&sum(theta==0)==0) {
      if (abs(theta[1]<2)) stop("scaled t df must be >2")
      if (sum(theta<0)) { 
        iniTheta <- c(log(abs(theta[1])-2),log(abs(theta[2]))) ## initial theta supplied
      } else { ## fixed theta supplied
        iniTheta <- Theta <- c(log(theta[1]-2),log(theta[2])) 
        n.theta <- 0 ## no thetas to estimate
      }
    } else iniTheta <- c(-2,-1) ## inital log theta value
               
    env <- new.env(parent = .GlobalEnv)
    assign(".Theta", iniTheta, envir = env)
    getTheta <- function(trans=FALSE) { 
    ## trans transforms to the original scale...
      th <- get(".Theta")
      if (trans) { th <- exp(th); th[1] <- th[1] + 2  }
      th
    }
    putTheta <- function(theta) assign(".Theta", theta,envir=environment(sys.function()))

    variance <- function(mu) { 
        th <- get(".Theta")
        nu <- exp(th[1])+2; sig <- exp(th[2])
        sig^2*nu/(nu-2)
    }

   validmu <- function(mu) all(is.finite(mu))

   dev.resids <- function(y, mu, wt,theta=NULL) {
      if (is.null(theta)) theta <- get(".Theta")
      nu <- exp(theta[1])+2; sig <- exp(theta[2])
      wt * (nu + 1)*log(1+(1/nu)*((y-mu)/sig)^2)
    }
    
    Dd <- function(y, mu, theta, wt, level=0) {
    ## derivatives of the deviance...
      ltheta <- theta
    #  theta <- exp(theta); nu <- theta[1]; sig <- theta[2]
      nu <- exp(theta[1])+2; sig <- exp(theta[2])
      nu1 <- nu + 1;  ym <- y - mu; nu2 <- nu - 2;
      a <- 1 + (ym/sig)^2/nu
      oo <- list()
      ## get the quantities needed for IRLS. 
      ## Dmu2eta2 is deriv of D w.r.t mu twice and eta twice,
      ## Dmu is deriv w.r.t. mu once, etc...
      nu1ym <- nu1*ym
      sig2a <- sig^2*a
      nusig2a <- nu*sig2a
      f <- nu1ym/nusig2a
      f1 <- ym/nusig2a
      oo$Dmu <- -2 * wt * f
      oo$Dmu2 <- 2 * wt * nu1*(1/nusig2a- 2*f1^2)  # - 2*ym^2/(nu^2*sig^4*a^2)
      term <- 2*nu1/sig^2/(nu+3)
      n <- length(y) 
      oo$EDmu2 <- rep(term,n)
      if (level>0) { ## quantities needed for first derivatives
        nu1nusig2a <- nu1/nusig2a
        nu2nu <- nu2/nu
        fym <- f*ym; ff1 <- f*f1; f1ym <- f1*ym; fymf1 <- fym*f1
        ymsig2a <- ym/sig2a
        oo$Dmu2th <- oo$Dmuth <- oo$Dth <- matrix(0,n,2)
        oo$Dth[,1] <- 1 * wt * nu2 * (log(a) - fym/nu) 
        oo$Dth[,2] <- -2 * wt * fym    
        oo$Dmuth[,1] <- 2 * wt *(f - ymsig2a - fymf1)*nu2nu
      #  oo$Dmuth[,2] <- 8*wt*nu1*ym*(1/(nu*sig^2*a)- ym^2/(nu^2*sig^4*a^2))
        oo$Dmuth[,2] <- 4* wt* f* (1- f1ym)
      #  oo$Dmu3 <- 8 * wt *nu1*ym * (3/(nu^2*sig^4*a^2) - 4*ym^2/(nu^3*sig^6*a^3)) 
        oo$Dmu3 <- 4 * wt * f * (3/nusig2a - 4*f1^2) 
        oo$Dmu2th[,1] <- 2* wt * (-nu1nusig2a + 1/sig2a + 5*ff1- 2*f1ym/sig2a - 4*fymf1*f1)*nu2nu
      # oo$Dmu2th[,2] <- 8*wt*nu1*(-1/(nu*sig^2*a) + 5*ym^2/(nu^2*sig^4*a^2) - 4*ym^4/(nu^3*sig^6*a^3))
        oo$Dmu2th[,2] <- 4*wt*(-nu1nusig2a + ff1*5 - 4*ff1*f1ym)
      } 
      if (level>1) { ## whole lot
        nu1nu2 <- nu1*nu2; nu1nu <- nu1/nu
        fymf1ym <- fym*f1ym; f1ymf1 <- f1ym*f1
        oo$Dmu4 <- 12 * wt * (-nu1nusig2a/nusig2a + 8*ff1/nusig2a - 8*ff1 *f1^2) 
        n2d <- 3 # number of the 2nd order derivatives
        oo$Dmu3th <- matrix(0,n,2)
        oo$Dmu2th2 <- oo$Dmuth2 <- oo$Dth2 <- matrix(0,n,n2d)
        oo$Dmu3th[,1] <- 4*wt*(-6*f/nusig2a + 3*f1/sig2a + 18*ff1*f1 - 4*f1ymf1/sig2a - 12*nu1ym*f1^4)*nu2nu
        oo$Dmu3th[,2] <- 48*wt* f* (- 1/nusig2a + 3*f1^2 -  2*f1ymf1*f1)
        
        oo$Dth2[,1] <- 1*wt *(nu2*log(a) +nu2nu*ym^2*(-2*nu2-nu1+ 2*nu1*nu2nu - nu1*nu2nu*f1ym)/nusig2a) ## deriv of D w.r.t. theta1 theta1 
  
        oo$Dth2[,2] <- 2*wt*(fym - ym*ymsig2a - fymf1ym)*nu2nu  ## deriv of D wrt theta1 theta2
        oo$Dth2[,3] <- 4 * wt * fym *(1 - f1ym)  ## deriv of D wrt theta2 theta2

        term <- 2*nu2nu - 2*nu1nu*nu2nu -1 + nu1nu
        oo$Dmuth2[,1] <- 2*wt*f1*nu2*(term - 2*nu2nu*f1ym + 4*fym*nu2nu/nu - fym/nu - 2*fymf1ym*nu2nu/nu)
 
        oo$Dmuth2[,2] <- 4*wt* (-f + ymsig2a + 3*fymf1 - ymsig2a*f1ym - 2*fymf1*f1ym)*nu2nu

        oo$Dmuth2[,3] <- 8*wt* f * (-1 + 3*f1ym - 2*f1ym^2)

        oo$Dmu2th2[,1] <- 2*wt*nu2*(-term + 10*nu2nu*f1ym -
               16*fym*nu2nu/nu - 2*f1ym + 5*nu1nu*f1ym - 8*nu2nu*f1ym^2 + 26*fymf1ym*nu2nu/nu - 
               4*nu1nu*f1ym^2 - 12*nu1nu*nu2nu*f1ym^3)/nusig2a
       
        oo$Dmu2th2[,2] <- 4*wt*(nu1nusig2a - 1/sig2a - 11*nu1*f1^2 + 5*f1ym/sig2a + 22*nu1*f1ymf1*f1 - 
                    4*f1ym^2/sig2a - 12*nu1*f1ymf1^2)*nu2nu
        oo$Dmu2th2[,3] <- 8*wt * (nu1nusig2a - 11*nu1*f1^2 + 22*nu1*f1ymf1*f1 - 12*nu1*f1ymf1^2)

       #  term <- nu1nusig2a - 11*nu1*f1^2 + 22*nu1*f1ymf1*f1 - 12*nu1*f1ymf1^2
       #  oo$Dmu2th2[,2] <- 8*wt*(term - 1/sig2a + 5*f1ym/sig2a - 4*f1ym^2/sig2a)*nu2nu
       #  oo$Dmu2th2[,3] <- 16*wt * term
      }
      oo
    }  ## end of Dd

    aic <- function(y, mu, theta=NULL, wt, dev) {
        if (is.null(theta)) theta <- get(".Theta")
        nu <- exp(theta[1])+2; sig <- exp(theta[2])
        term <- -lgamma((nu+1)/2)+ lgamma(nu/2) + log(sig*(pi*nu)^.5) +
           (nu+1)*log(1 + ((y-mu)/sig)^2/nu)/2  ## `-'log likelihood for each observation
        2 * sum(term * wt)
    }
    
    ls <- function(y,w,n,theta,scale) {
       ## the log saturated likelihood function.
     #  Theta <- exp(theta);  nu <- Theta[1]; sig <- Theta[2]
       nu <- exp(theta[1])+2; sig <- exp(theta[2]); nu2 <- nu-2;
       nu2nu <- nu2/nu; nu12 <- (nu+1)/2
       term <- lgamma(nu12) - lgamma(nu/2) - log(sig*(pi*nu)^.5)
       ls <- sum(term*w) 
       ## first derivative wrt theta...
       lsth <- rep(0,2) 
       lsth2 <- matrix(0,2,2)  ## rep(0, 3)
       term <- nu2 * digamma(nu12)/2- nu2 * digamma(nu/2)/2 - 0.5*nu2nu
       lsth[1] <- sum(w*term)
       lsth[2] <- sum(-1*w)
       
       ## second deriv...      
       term <-  nu2^2 * trigamma(nu12)/4 + nu2 * digamma(nu12)/2 -
           nu2^2 * trigamma(nu/2)/4 - nu2 * digamma(nu/2)/2 + 0.5*(nu2nu)^2 - 0.5*nu2nu
       lsth2[1,1] <- sum(term*w)
       lsth2[1,2] <- lsth2[2,1] <- lsth2[2,2] <- 0
       list(ls=ls,## saturated log likelihood
            lsth1=lsth, ## first derivative vector wrt theta
            lsth2=lsth2) ## Hessian wrt theta
    }

    preinitialize <- expression({
      ## initialize theta from raw observations..
       if (G$family$n.theta>0) {
         Theta <- c(-1, log(0.2*var(G$y)^.5))
         G$family$putTheta(Theta)
       } ## otherwise fixed theta supplied
    })

    initialize <- expression({
        if (any(is.na(y))) stop("NA values not allowed for the scaled t family")
        n <- rep(1, nobs)
        mustart <- y + (y == 0)*.1
    })
    postproc <- expression({
      object$family$family <- 
      paste("Scaled t(",paste(round(object$family$getTheta(TRUE),3),collapse=","),")",sep="")
    })
    rd <- function(mu,wt,scale) {
     ## simulate data given fitted latent variable in mu 
      theta <- get(".Theta")
      nu <- exp(theta[1])+2; sig <- exp(theta[2])
      n <- length(mu)
      stats::rt(n=n,df=nu)*sig + mu
    }

    environment(dev.resids) <- environment(aic) <- environment(getTheta) <- 
    environment(rd)<- environment(variance) <- environment(putTheta) <- env

    structure(list(family = "scaled t", link = linktemp, linkfun = stats$linkfun,
        linkinv = stats$linkinv, dev.resids = dev.resids,Dd=Dd,variance=variance,postproc=postproc,
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize,ls=ls, preinitialize=preinitialize,
        validmu = validmu, valideta = stats$valideta,n.theta=n.theta,   
        ini.theta = iniTheta,putTheta=putTheta,getTheta=getTheta, rd=rd),
        class = c("extended.family","family"))
} ## scat



## zero inflated Poisson...


ziP <- function (theta = NULL, link = "log") { 
## Extended family object for zero inflated distribution
## n.theta=2; log theta supplied
  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  if (linktemp %in% c("log")) stats <- make.link(linktemp) ## done only for the "log" link at the moment
  else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  } else {
    if (inherits(link, "link-glm")) {
       stats <- link
            if (!is.null(stats$name))
                linktemp <- stats$name
        }
        else stop(linktemp, " link not available for zero inflated; available link for `lambda' is only  \"log\"")
    }
    Theta <-  NULL;n.theta <- 2
    if (!is.null(theta)) {
      if (theta[2]<0) iniTheta <- c(theta[1],log(-theta[2])) ## initial theta supplied
      else { ## fixed theta supplied
        iniTheta <- Theta <- c(theta[1],log(theta[2]))
        n.theta <- 0 ## no thetas to estimate
      }
    } else iniTheta <- c(1,-1) ## inital theta value
           
    
    env <- new.env(parent = .GlobalEnv)
    assign(".Theta", iniTheta, envir = env)
    getTheta <- function(trans=FALSE) { 
     ## trans transforms to the original scale...
      th <- get(".Theta")
      if (trans) th[2] <- exp(th[2])
      th
    }
    putTheta <- function(theta) assign(".Theta", theta,envir=environment(sys.function()))

   # variance <- function(mu) { 
   #    th <- get(".Theta")
   #    f <- th[1] - exp(th[2])*mu
   #    mu + exp(f)*mu^2
   # }

   validmu <- function(mu) all(mu > 0)

   dev.resids <- function(y, mu, wt,theta=NULL) {
      if (is.null(theta)) theta <- get(".Theta")
      ## note: log theta_2 supplied and working on `lambda' scale rather than `mu'
      th1 <- theta[1]; th2 <- exp(theta[2]); la <- mu;
      if (length(la)==1) la <- rep(mu,length(y))
      yzero <- y == 0
      dr <- rep(0,length(y))
      f <- th1 - th2*la[yzero]; ela <- exp(-la[yzero]);
      h <- ind <- f >0 
      h[ind] <- (1 + ela[ind]*exp(-f[ind]))/(1 + exp(-f[ind]))
      f <- exp(f[!ind])
      h[!ind] <- (f + ela[!ind])/(1+f)
      
      f <- th1 - th2*y[yzero]; ela <- exp(-y[yzero]);
      h0 <- ind <- f >0; 
      h0[ind] <- (1 + ela[ind]*exp(-f[ind]))/(1 + exp(-f[ind]))
      f <- exp(f[!ind])
      h0[!ind] <- (f + ela[!ind])/(1+f)
      dr[yzero] <- log(h0) - log(h);
      dr[!yzero] <- -log(1+exp(th1-th2*y[!yzero])) + log(1+exp(th1-th2*la[!yzero])) - y[!yzero] + 
          la[!yzero] + y[!yzero]*log(y[!yzero]) - y[!yzero]*log(la[!yzero])
      2*dr*wt
    }
    
    Dd <- function(y, mu, theta, wt=NULL, level=0) {
    ## derivatives of the deviance...
      ltheta <- theta
      th1 <- theta[1]; th2 <- exp(theta[2]); la <- mu;
      ## notations for the first terms when y == 0...
      yzero <- y == 0
      f1 <- th1-th2*la[yzero]; 
      a1 <- a2 <- d2 <- t1 <- ind <- f1 > 0
      emf <- exp(-f1[ind]); ef <- exp(f1[!ind]); efa <- exp(f1);
      a1[ind] <- 1/(1+emf); a1[!ind] <- ef/(1+ef)
      b1 <- exp(-la[yzero])
      a2[ind] <- emf/(1+emf)^2; a2[!ind] <- ef/(1+ef)^2;
      d2[ind] <- emf/(1+emf*b1[ind])^2; d2[!ind] <- ef/(b1[!ind]+ef)^2;
      t1[ind] <- (th2 + emf*b1[ind])/(1+emf*b1[ind]); t1[!ind] <- (th2*ef + b1[!ind])/(b1[!ind]+ef);
      
      ## notations for the second terms when y > 0...
      f1.s <- th1-th2*la[!yzero]; 
      a1.s <- a2.s <- ind <- f1.s > 0
      emf.s <- exp(-f1.s[ind]); ef.s <- exp(f1.s[!ind]); efa.s <- exp(f1.s);
      a1.s[ind] <- 1/(1+emf.s); a1.s[!ind] <- ef.s/(1+ef.s)
      a2.s[ind] <- emf.s/(1+emf.s)^2; a2.s[!ind] <- ef.s/(1+ef.s)^2;
      
      oo <- list()
      ## get the quantities needed for IRLS. 
      ## Dmu2eta2 is deriv of D w.r.t mu twice and eta twice,
      ## Dmu is deriv w.r.t. mu once, etc...
      n <- length(y)
      ath <- a1*th2; ath.s <- a1.s*th2
            
      oo$Dmu <- oo$Dmu2 <- rep(0,n)
      oo$Dmu[yzero] <- 2*(-ath + t1)
      oo$Dmu[!yzero] <- 2*(-y[!yzero]/la[!yzero] - ath.s +1)
      oo$Dmu2[yzero] <-  2*(a2*th2^2 - d2*b1*(th2-1)^2)
      oo$Dmu2[!yzero] <-  2*(y[!yzero]/la[!yzero]^2 + a2.s*th2^2)
      f <- th1-th2*la
      a <- aa <- ind <- f > 0
      b <- exp(-la);
      emf <- exp(-f[ind]); ef <- exp(f[!ind]); efa <- exp(f);
      a[ind] <- 1/(1+emf); a[!ind] <- ef/(1+ef)
      c <- 1/(efa + b); s <- 1/(1+efa);
      aa[ind] <- emf/(1+emf)^2; aa[!ind] <- ef/(1+ef)^2;
      oo$EDmu2 <- 2*(aa*th2^2 + s/la- a*b*c*(th2-1)^2)

      if (level>0) { ## quantities needed for first derivatives
        ## notations for the first terms when y == 0...
        a3 <- a23 <- d1 <- t2 <- t3 <-  g1 <- g12 <- ind <- f1 > 0
        emf <- exp(-f1[ind]); ef <- exp(f1[!ind]); efa <- exp(f1);
        c1 <- 1/(efa + b1); s1 <- 1/(1+efa);
      
        a3[ind] <- emf^2/(1+emf)^3; a3[!ind] <- ef/(1+ef)^3; 
        a23[ind] <- emf/(1+emf)^3; a23[!ind] <- ef^2/(1+ef)^3;
        d1[ind] <- 1/(1+emf*b1[ind]); d1[!ind] <- ef/(b1[!ind]+ef);
        t2[ind] <- (th2^2 + emf*b1[ind])/(1+emf*b1[ind]); t2[!ind] <- (th2^2*ef + b1[!ind])/(b1[!ind]+ef);   
        t3[ind] <- (th2^3 + emf*b1[ind])/(1+emf*b1[ind]); t3[!ind] <- (th2^3*ef + b1[!ind])/(b1[!ind]+ef);
        g1 <- a1 - d1; g12 <- a1^2 - d1^2;

        f0 <- th1-th2*y[yzero]; 
        a0 <- a20 <- b0 <- c0 <-  d0 <- g0 <- g02 <- ind <- f0 > 0
        emf0 <- exp(-f0[ind]); ef0 <- exp(f0[!ind]); efa0 <- exp(f0);
        b0 <- exp(-y[yzero])
        c0 <- 1/(efa0 + b0); 
        a0[ind] <- 1/(1+emf0); a0[!ind] <- ef0/(1+ef0) 
        d0[ind] <- 1/(1+ b0[ind]*emf0); d0[!ind] <- ef0/(b0[!ind]+ef0);
        a20[ind] <- emf0/(1+emf0)^2; a20[!ind] <- ef0/(1+ef0)^2;
        g0 <- a0 - d0; g02 <- a0^2 - d0^2
        lath <- la[yzero]*th2; olath <- 1 - lath; tlath <- 2 - lath
 
        ## notations for the second terms when y > 0...
        a3.s <- a23.s <- ind <- f1.s > 0
        a3.s[ind] <- emf.s^2/(1+emf.s)^3; a3.s[!ind] <- ef.s/(1+ef.s)^3; 
        a23.s[ind] <- emf.s/(1+emf.s)^3; a23.s[!ind] <- ef.s^2/(1+ef.s)^3;
       
        f0.s <- th1-th2*y[!yzero]; 
        a0.s <- a20.s <- ind <- f0.s > 0
        emf0 <- exp(-f0.s[ind]); ef0 <- exp(f0.s[!ind]); efa0 <- exp(f0.s);
        a0.s[ind] <- 1/(1+emf0); a0.s[!ind] <- ef0/(1+ef0) 
        a20.s[ind] <- emf0/(1+emf0)^2; a20.s[!ind] <- ef0/(1+ef0)^2;
        lath.s <- la[!yzero]*th2; olath.s <- 1 - lath.s; tlath.s <- 2 - lath.s

        oo$Dmu2th <- oo$Dmuth <- oo$Dth <- matrix(0,n,2)
        b0m <- 1 - b0; b1m <- 1 - b1;
        a2th.s <- a2.s * th2
        d1t1 <- d1 * t1

        oo$Dth[yzero,1] <- 2*(a0*c0*b0m - a1*b1m*c1)
        oo$Dth[!yzero,1] <- 2*(a1.s - a0.s)
        oo$Dth[yzero,2] <- 2*( -y[yzero]*a0*th2*c0*b0m + lath*a1*b1m*c1)
        oo$Dth[!yzero,2] <- 2*(y[!yzero]*a0.s*th2 - lath.s*a1.s)

        oo$Dmuth[yzero,1] <- 2 * (-th2*g1 + th2*a1^2 - d1t1)
        oo$Dmuth[!yzero,1] <- -2*a2th.s
        oo$Dmuth[yzero,2] <- 2 *(-th2*olath*g1 + lath*d1t1 -lath*th2*a1^2)
        oo$Dmuth[!yzero,2] <- -2*th2*(a2.s +a1.s^2 - lath.s*a2.s)

        oo$Dmu3[yzero] <- 2*(-ath*th2^2 + t3+ 3*ath^2*th2 - 2*ath^3 - 3*t1*t2 + 2*t1^3)
        oo$Dmu3[!yzero] <- 2*(-2*y[!yzero]/la[!yzero]^3 - th2^3*(a3.s - a23.s))
        oo$Dmu2th[yzero,1] <- 2*(-2*d1t1*t1 + th2^2*g1 - 3*ath^2 +2*ath^2*a1 + 3*th2^2*d1^2 + 2*th2*b1*d2 + b1*d2)
        oo$Dmu2th[!yzero,1] <- 2 * th2 * ath.s*(1 - 3*a1.s +2*a1.s^2)
        oo$Dmu2th[yzero,2] <- 2* (-lath*d1*t2 + th2^2*tlath*g1 + 2*th2*olath*d1t1 + 2*lath*d1t1*t1 -
            2*lath*ath^2*a1 - ath^2*(2-3*lath))
        oo$Dmu2th[!yzero,2] <- 2*th2^2*(tlath.s*a1.s + a1.s^2*(3*lath.s - 2) - 2*lath.s*a1.s^3)
      } 
      if (level>1) { ## whole damn lot
        d1t2 <- d1*t2              
        term1 <- -d1*th2^4 -b1*c1 + ath*th2^3 + 12*ath^3*th2 - 6*ath^4 
        term2 <- - 7*ath^2*th2^2 + 4*t1*t3 + 3*t2^2 - 12*t1^2*t2 + 6*t1^4
        oo$Dmu4[yzero] <- 2*(term1 + term2)
        oo$Dmu4[!yzero] <- 2*(6*y[!yzero]/la[!yzero]^4 + a2th.s*th2^3 - 6*a1.s^2*th2^4 + 12*a23.s*a1.s*th2^4 
               + 6*ath.s^4) 
        n2d <- 3 # number of 2nd order derivatives
        oo$Dmu3th <- matrix(0,n,2)
        oo$Dmu2th2 <- oo$Dmuth2 <- oo$Dth2 <- matrix(0,n,n2d)
        term1 <- -6*d1t1*t1^2 + 6*th2*d1t1*t1 + th2^3*(d1 - a1+ 7*a1^2 -12*a1^3 + 6*a1^4)
        term2 <- -d1*t3 - 3*th2*d1t2 + 6*d1t1*t2 - 3*th2^2*d1t1
        oo$Dmu3th[yzero,1] <- 2*(term1 + term2)
        oo$Dmu3th[!yzero,1] <- 2*th2^3 *a1.s * (-1 + 7*a1.s - 12*a1.s^2 + 6*a1.s^3)

        oo$Dmu3th[yzero,2] <- 2 * (th2^3*((lath-3)*g1 + 3*a1^2*(3-2*lath) - 6*a1^3*olath) - 3*th2*olath*d1t2 + lath*d1*t3 - 6*lath*d1t1*t2 -3*th2^2*tlath*d1t1 + 6*lath*d1t1*t1^2 + 6*th2*olath*d1t1*t1 + 
             lath*th2^3*a1^2*(6*a1 - 6*a1^2 - 1))
        oo$Dmu3th[!yzero,2] <- 2*ath.s*th2^2*(-3 +lath.s +(9-7*lath.s)*a1.s - 
                6*(1-2*lath.s)*a1.s^2 - 6*lath.s*a1.s^3)

        oo$Dth2[yzero,1] <- 2* (a20*b0m*b0*c0^2 -a20*b0m*d0^2 - a2*b1m*b1*c1^2 + a2*b1m*d1^2)  ## deriv of D w.r.t. theta1 theta1
        oo$Dth2[!yzero,1] <- 2 * (- a20.s +a2.s)
        
        y0th <- y[yzero]*th2; y1th <- y[!yzero]*th2;
        oo$Dth2[yzero,2] <- 2*(y0th*g0 - y0th*g02 - lath*g1 + lath*g12 )## deriv of D wrt theta1 theta2
        oo$Dth2[!yzero,2] <- 2*(y1th*a0.s*(1 - a0.s) - lath.s*a1.s*(1 - a1.s) )    

        oo$Dth2[yzero,3] <- 2*(y0th*(1-y0th)*g0 + y0th^2*g02^2 -
           lath*olath*g1 - lath^2*g12) ## deriv of D wrt theta2 theta2
        oo$Dth2[!yzero,3] <- 2*(y1th*(1-y1th)*a0.s + (y1th*a0.s)^2 - lath.s*a1.s*olath.s - (lath.s*a1.s)^2 )    

        oo$Dmuth2[yzero,1] <- 2* (th2*d1*(1 - 3*d1) - d2*b1 + 2*d1^2*t1 - ath*(1- 3*a1 + 2*a1^2))
        oo$Dmuth2[!yzero,1] <- 2*ath.s*(-1 +3*a1.s - 2*a1.s^2)

        oo$Dmuth2[yzero,2] <- 2* (lath*d1t1*(1-2*d1) -th2*olath*(g1 +d1^2) + lath*th2*(2*a1^3 - g12) +
           th2*(1-2*lath)*a1^2 )
        oo$Dmuth2[!yzero,2] <- 2*ath.s*( -olath.s +(1-3*lath.s)*a1.s + 2*lath.s*a1.s^2)

        oo$Dmuth2[yzero,3] <- 2*th2*(lath*g1 + la[yzero]*olath*d1t1 - olath^2*g1 - 2*lath^2*a1^3 + 
           2*lath*la[yzero]*d1t1*d1 + lath*olath*(2*d1^2 - 3*a1^2))
        oo$Dmuth2[!yzero,3] <- 2*ath.s*(lath.s - olath.s^2 - 3*lath.s*a1.s*olath.s - 2*lath.s^2*a1.s^2)
  
        d1d2b1 <- d1*d2*b1
        term1 <- -th2^2*d1*(1 - 4*d1+ 10*d1^2) + d1t2 + 2*th2*d1t1 - 8*th2*d1d2b1 - 2*d1d2b1 -
            2*d1t1*t1 + 6*d1t1^2 
        term2 <- th2^2*a1*(1 - 7*a1 + 12*a1^2 - 6*a1^3)
        oo$Dmu2th2[yzero,1] <- 2* (term1 + term2)
        oo$Dmu2th2[!yzero,1] <- 2*th2*ath.s*(1 - 7*a1.s + 12*a1.s^2 - 6*a1.s^3)

        oo$Dmu2th2[yzero,2] <- 2 * (lath*d1t2*(2*d1-1) + tlath*th2^2*g1 + 2*th2*olath*d1t1 -
           4*th2*olath*d1^2*t1 + 2*lath*d1*t1^2 - 6*lath*d1^2*t1^2 - 6*lath*th2^2*a1^3*(1-a1) +
           4*lath*th2*d1^2*t1 +lath*th2^2*g12 + th2^2*(4-3*lath)*d1^2 - 6*th2^2*olath*a1^2 +
           2*th2^2*(2-3*lath)*a1^3 )
        oo$Dmu2th2[!yzero,2] <- 2 * th2^2*a1.s*(tlath.s + (7*lath.s-6)*a1.s + 4*(1-3*lath.s)*a1.s^2 
            + 6*lath.s*a1.s^3)   

    #    term1 <- lath*(-th2^2*g1 - 2*lath*d1t2*d1 - 2*th2*d1t1 + 6*lath*d1t1^2 - 2*th2^2*tlath*d1^2) 
    #       + lath*olath*(-d1t2 + 2*d1t1*t1 + 8*th2*d1t1*d1)
        term2 <-  th2^2*tlath^2*g1 + 2*th2*olath^2*d1t1 + 2*lath^2*th2^2*a1^3*(1-3*a1) +
          lath*th2^2*a1^2*(-10*olath*a1 + 7-3*lath) + 2*th2^2*olath^2*(d1^2 - 2*a1^2)
  #      oo$Dmu2th2[yzero,3] <- 2* (term1 + term2 )
        
        oo$Dmu2th2[yzero,3] <- 2* ( lath*(-th2^2*g1 - 2*lath*d1t2*d1 - 2*th2*d1t1 + 6*lath*d1t1^2 - 2*th2^2*tlath*d1^2)  + lath*olath*(-d1t2 + 2*d1t1*t1 + 8*th2*d1t1*d1)  + term2 )

        oo$Dmu2th2[!yzero,3] <- 2*ath.s*th2* (-lath.s + tlath.s^2 + lath.s*a1.s*(7-3*lath.s) +
            2*a1.s^2*lath.s^2 - 10*lath.s*a1.s^2*olath.s - 6*lath.s^2*a1.s^3 - 4*a1.s*olath.s^2) 
      }
      oo
    }  ## end of Dd

    aic <- function(y, mu, theta=NULL, wt, dev) {
        if (is.null(theta)) theta <- get(".Theta")
        th1 <- theta[1]; th2 <- exp(theta[2]); la <- mu;
        yzero <- y == 0
        term <- rep(0,length(y))  ## log likelihood for each observation
        f <- th1 - th2*la[yzero]; ela <- exp(-la[yzero]);
        h <- ind <- f >0 
        h[ind] <- (1 + ela[ind]*exp(-f[ind]))/(1 + exp(-f[ind]))
        f <- exp(f[!ind])
        h[!ind] <- (f + ela[!ind])/(1+f)
        term[yzero] <- log(h)
        f <- exp(th1-th2*la[!yzero])
        term[!yzero] <- -log(1+f) -la[!yzero] +y[!yzero]*log(la[!yzero]) -
                 lfactorial(y[!yzero]) 
        -2 * sum(term * wt)
    }
    
    ls <- function(y,w,n,theta,scale) {
       ## the log saturated likelihood function.
       th1 <- theta[1]; th2 <- exp(theta[2]); 
       yzero <- y == 0
       term <- rep(0,length(y))  ## saturated log likelihood for each observation
       f <- th1 - th2*y[yzero]; ey <- exp(-y[yzero]);
       h <- ind <- f >0 
       h[ind] <- (1 + ey[ind]*exp(-f[ind]))/(1 + exp(-f[ind]))
       f <- exp(f[!ind])
       h[!ind] <- (f + ey[!ind])/(1+f)
       term[yzero] <- log(h)
       f <- exp(th1-th2*y[!yzero])
       term[!yzero] <- -log(1+f) - y[!yzero] +y[!yzero]*log(y[!yzero]) -
                 lfactorial(y[!yzero]) 
       ls <- sum(term*w) 
       ## first derivatives wrt theta...
       lsth <- rep(0,2) 
       f <- th1 - th2*y[yzero]
       a0 <- a20 <- b0 <- c0 <- d0 <- ind <- f > 0
       emf <- exp(-f[ind]); ef <- exp(f[!ind]); efa <- exp(f);
       b0 <- exp(-y[yzero])
       c0 <- 1/(efa + b0); 
       a0[ind] <- 1/(1+emf); a0[!ind] <- ef/(1+ef) 
       a20[ind] <- emf/(1+emf)^2; a20[!ind] <- ef/(1+ef)^2;    
       d0[ind] <- 1/(1+ b0[ind]*emf); d0[!ind] <- ef/(b0[!ind]+ef);   
       
       f0.s <- th1-th2*y[!yzero]; 
       a0.s <- a20.s <- ind <- f0.s > 0
       emf0 <- exp(-f0.s[ind]); ef0 <- exp(f0.s[!ind]); efa0 <- exp(f0.s);
       a0.s[ind] <- 1/(1+emf0); a0.s[!ind] <- ef0/(1+ef0) 
       a20.s[ind] <- emf0/(1+emf0)^2; a20.s[!ind] <- ef0/(1+ef0)^2;

       term <- rep(0,length(y))
       term[yzero] <- a0*(1-b0)*c0;  term[!yzero] <- -a0.s
       lsth[1] <- sum(w*term)
       term[yzero] <- y[yzero]*th2*a0*(b0-1)*c0
       term[!yzero] <- y[!yzero]*th2*a0.s
       lsth[2] <- sum(w*term)
       
       ## second deriv...      
       lsth2 <- matrix(0,2,2)  
       term[yzero] <- tt <- (1-b0)*a20 * (b0*c0^2 - d0^2)
       term[!yzero] <- -a20.s 
       lsth2[1,1] <- sum(term*w)
       term[yzero] <- -y[yzero]*th2*tt
       term[!yzero] <- y[!yzero]*th2*a20.s
       lsth2[1,2] <- lsth2[2,1] <- sum(term*w)
       
       term[yzero] <- y[yzero]*th2*(b0-1)*a20*(d0^2 + a0*b0*c0^2) + (y[yzero]*th2)^2*tt
       term[!yzero] <- y[!yzero]*th2* (a20.s + a0.s^2 - th2*y[!yzero]*a20.s)
       lsth2[2,2] <- sum(term*w)
       list(ls=ls,## saturated log likelihood
         lsth1=lsth,lsth2=lsth2) ## its first and second derivs wrt theta
    }

    initialize <- expression({
        if (any(y < 0)) stop("negative values not allowed for the zero inflated Poisson family")
        n <- rep(1, nobs)
        mustart <- y +exp(-y) # + (y==0)/5
    })
    postproc <- expression({
      object$family$family <- 
      paste("Zero inflated Poisson(",paste(round(object$family$getTheta(TRUE),3),collapse=","),")",sep="")
    })
    fv <- function(mu,theta=NULL) {
    ## optional function to give fitted values - idea is that 
    ## predict.gam(...,type="response") will use this, as well
    ## as residuals(...,type="response")...
      if (is.null(theta)) theta <- get(".Theta")
      th1 <- theta[1]; th2 <- exp(theta[2]); 
      p <- exp(th1- th2*mu)/(1 + exp(th1- th2*mu)) 
      fv <- (1 - p)* mu      
      fv
    } ## fv

    rd <- function(mu,wt,scale) {
    ## simulate data given fitted latent variable in mu 
      theta <- get(".Theta") 
      th1 <- theta[1]; th2 <- exp(theta[2]); 
      n.sim <- length(mu)
      p <- exp(th1- th2*mu)/(1 + exp(th1- th2*mu)) 
      z <- y <- runif(n.sim)
      y <- rep(0,n.sim)
      good <- z >= p
      y[good] <- rpois(sum(good),mu[good])
      y
    }
   
    environment(dev.resids) <- environment(aic) <- environment(getTheta) <- 
    environment(putTheta) <- environment(rd) <- environment(fv) <- env

    structure(list(family = "zero inflated Poisson", link = linktemp, linkfun = stats$linkfun,
        linkinv = stats$linkinv, dev.resids = dev.resids,Dd=Dd, rd=rd,
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize,postproc=postproc,ls=ls,fv=fv,
        validmu = validmu, valideta = stats$valideta,n.theta=n.theta, 
        ini.theta = iniTheta,putTheta=putTheta,getTheta=getTheta),
        class = c("extended.family","family"))
} ## zip
