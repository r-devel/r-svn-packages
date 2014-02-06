## (c) Simon N. Wood 2013. Released under GPL2.

## extended families for mgcv...

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
  Theta <-  NULL
  ## NOTE: following not really right - should initialize from cut points
  if (!is.null(theta)&&sum(theta==0)==0) {
    if (sum(theta<0)) iniTheta <- log(abs(theta)) ## initial theta supplied
    else iniTheta <- Theta <- log(theta) ## fixed theta supplied
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
    if (R3>2) {
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
      res <- sqrt(pmax(res,0)) * s 
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
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize,
        preinitialize = preinitialize, ls=ls,rd=rd,residuals=residuals,
        validmu = validmu, valideta = stats$valideta,n.theta=R-2,
        ini.theta = iniTheta,putTheta=putTheta,predict=predict,
        getTheta=getTheta), class = c("extended.family","family"),step = 1)
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
    Theta <-  NULL
    if (!is.null(theta)&&theta!=0) {
      if (theta>0) iniTheta <- Theta <- log(theta) ## fixed theta supplied
      else iniTheta <- log(-theta) ## initial theta supplied
    } else iniTheta <- 0 ## inital log theta value
    
    env <- new.env(parent = .GlobalEnv)
    assign(".Theta", iniTheta, envir = env)
    getTheta <- function(trans=FALSE) if (trans) exp(.Theta) else .Theta # get(".Theta")
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
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize,ls=ls,
        validmu = validmu, valideta = stats$valideta,n.theta=1, 
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
  Theta <-  NULL
  if (!is.null(theta)&&theta!=0) {
      if (theta>0) iniTheta <- Theta <- log((theta-a)/(b-theta)) ## fixed theta supplied
      else iniTheta <- log((-theta-a)/(b+theta)) ## initial theta supplied
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
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize,ls=ls,
        validmu = validmu, valideta = stats$valideta,canonical="none",n.theta=1, 
        ini.theta = iniTheta,putTheta=putTheta,getTheta=getTheta,scale = -1),
        class = c("extended.family","family"))
} ## tw

