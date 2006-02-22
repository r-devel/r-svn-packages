#### Mallows quasi-likelihood estimator of E. Cantoni and E. Ronchetti (2001)
#### based originally on Eva Cantoni's S-plus code "robGLM"

glmrobMqle <-
    function(X, y, weights = rep(1, nobs), start = NULL,
	     offset = rep(0, nobs), family, weights.on.x = "none",
	     control = glmrobMqle.control(), intercept = TRUE)
{
    ## To DO:
    ## o weights are not really implemented
    ## o offset is not	implemented
    ##
    X <- as.matrix(X)
    xnames <- dimnames(X)[[2]]
    ynames <- names(y)
    conv <- FALSE
    nobs <- NROW(y)
    ncoef <- ncol(X)
    ##	  EMPTY <- ncoef == 0
    if (is.null(weights))
	weights <- rep.int(1, nobs)
    if(any(weights <= 0))
	stop("All weights must be positive")
    if (is.null(offset))
	offset <- rep.int(0, nobs)
    variance <- family$variance
    linkinv <- family$linkinv
    mu.eta <- family$mu.eta
    if (!is.function(variance) || !is.function(linkinv))
	stop("illegal 'family' argument")
    valideta <- family$valideta
    if (is.null(valideta))
	valideta <- function(eta) TRUE
    validmu <- family$validmu
    if (is.null(validmu))
	validmu <- function(mu) TRUE
    ##

## MM: nonsense
##	   if(ncoef == 1)
##	       warning("There is only an intercept: No weights accepted on X")
##	   else {
    w.x <- if(ncoef) {
	switch(weights.on.x,
	       none = rep.int(1, nobs),
	       Hii = glmrobMqleHii(X = X),
	       robCov = glmrobMqleRobDist(X = X, intercept = intercept),
	       stop("Weighting method", sQuote(weights.on.x),
		    " is not implemented"))
    } else ## ncoef == 0
    rep.int(1,nobs)



### Initializations

    tcc <- control$tcc
    eval(family$initialize) ## --> n, mustart, y and weights (=ni)
    ni <- weights
    ##
    if(is.null(start))
	start <- glm.fit(x = X, y = y, weights = weights, offset = offset,
			 family = family)$coefficients
    thetaOld <- theta <- start
    eta <- drop(X %*% theta)
    mu <- linkinv(eta) # mu estimates pi (in [0,1]) at the binomial model
    if (!(validmu(mu) && valideta(eta)))
	stop("Can't find valid starting values: You need help")
    ##
    if(family$family == "poisson") {
	Epsi <- glmrobMqleEpsiPois
	EpsiS <- glmrobMqleEpsiSPois
	Epsi2 <- glmrobMqleEpsi2Pois
    }
    else {
	Epsi <- glmrobMqleEpsiB
	EpsiS <- glmrobMqleEpsiSB
	Epsi2 <- glmrobMqleEpsi2B
    }

### Iterations

    sni <- sqrt(ni)
    if(ncoef) for (nit in 1:control$maxit) {
	Vmu <- variance(mu)
	if (any(is.na(Vmu)))  stop("NAs in V(mu)")
	if (any(Vmu == 0))    stop("0s in V(mu)")
	dmu.deta <- mu.eta(eta)
	if (any(is.na(dmu.deta)))
	    stop("NAs in d(mu)/d(eta)")
	sV <- sqrt(Vmu)
	residP <- (y - mu)* sni/sV
	H <- floor(mu*ni - tcc* sni*sV)
	K <- floor(mu*ni + tcc* sni*sV)
	## Computation of alpha and (7) is done in one apply-loop:
	cpsi <- pmax(-tcc,pmin(residP,tcc)) - Epsi(mu, Vmu, ni, H, K, tcc)
	EEq <- colMeans(cpsi*w.x * sni/sV * dmu.deta * X)
	##
	## 1/n t(X) %*% B %*% X) %*%  delta.coef   = EEq
	DiagB <- EpsiS(mu, Vmu, ni, H, K, tcc) /(sni*sV) * w.x * (ni*dmu.deta)^2
	Dtheta <- solve(crossprod(X, DiagB*X)/nobs, EEq)
	if (any(!is.finite(Dtheta))) {
	    conv <- FALSE
	    warning("Non-finite coefficients at iteration ", nit)
	    break
	}
	theta <- thetaOld + Dtheta
	eta <- drop(X %*% theta) + offset
	mu <- linkinv(eta)
	## Check convergence
	convi <- sqrt(sum((Dtheta)^2)/max(1e-20, sum(thetaOld^2)))
	conv <-	 (convi <= control$acc)
	if(conv)  break
	thetaOld <- theta
    } ## end of iteration
    else { ## ncoef == 0
        conv <- TRUE
        nit <- 0
    }
    if (!conv)
	warning("Algorithm did not converge")
    eps <- 10 * .Machine$double.eps
    if (family$family == "binomial") {
	if (any(mu/weights > 1 - eps) || any(mu/weights < eps))
	    warning("fitted probabilities numerically 0 or 1 occurred")
    }
    if (family$family == "poisson") {
	if (any(mu < eps))
	    warning("fitted rates numerically 0 occurred")
    }

### FIXME: much of the following is "cut & paste" from the iteration part
    Vmu <- variance(mu)
    sV <- sqrt(Vmu)
    residP <- (y - mu)* sni/sV

    ## Estimated asymptotic covariance of the robust estimator

    H <- floor(mu*ni - tcc* sni*sV)
    K <- floor(mu*ni + tcc* sni*sV)
    dmu.deta <- mu.eta(eta)
    w.r <- pmin(1, tcc/abs(residP))
    weights <- w.x*w.r
    if(ncoef) {
	alpha <- colMeans(Epsi(mu, Vmu, ni,H,K, tcc) * w.x *
			  sni/sV * dmu.deta * X)
	DiagA <- Epsi2(mu, Vmu, ni,H,K, tcc) / (ni*Vmu)* w.x^2* (ni*dmu.deta)^2
	matQ  <- crossprod(X, DiagA*X)/nobs - outer(alpha, alpha)

	DiagB <- EpsiS(mu, Vmu, ni,H,K, tcc) / (sni*sV)* w.x * (ni*dmu.deta)^2
	matM <- crossprod(X, DiagB*X)/nobs
	matMinv <- solve(matM)
	asCov <-  matMinv %*% matQ %*% matMinv / nobs
    } else { ## ncoef == 0
	matM <- matQ <- asCov <- matrix(, 0,0)
    }
    ##
    list(coefficients = theta, residuals = residP, fitted.values = mu,
	 w.r = w.r, w.x = w.x, ni = ni, cov = asCov, matM = matM, matQ = matQ,
	 tcc = tcc, family = family, linear.predictors = eta, deviance = NULL,
	 iter = nit, y = y, converged = conv)
}


glmrobMqleHii <- function(X) {
    x <- qr(X)
    Hii <- colSums(qr.qy(x, diag(1, nrow = NROW(y), ncol = x$rank))^2)
    sqrt(1-Hii)
}

glmrobMqleRobDist <- function(X, intercept) {
    if(intercept) {
	Z <- as.matrix(X[, -1])
	Zrc <- cov.rob(Z) ## worse:  method="mcd"
	dist2 <- mahalanobis(Z, center = Zrc$center, cov = Zrc$cov)
	## dist2 <- mahalanobis(Z, center=XX.rcov$center, cov=XX.rcov$cov)
    }
    else {
	Z <- as.matrix(X)
	Zrc <- cov.rob(Z)
	mu <- as.matrix(Zrc$center)
	Mu <- Zrc$cov + mu %*% t(mu)
	dist2 <- mahalanobis(Z, center = rep(0,ncol(Z)), cov = Mu)
    }
    ncoef <- ncol(X)-intercept ## E[chi^2_p] = p
    1/sqrt(pmax(0, (dist2-ncoef)/sqrt(2*ncoef)*8) + 1)
}


## MM: 'acc' seems a misnomer to me, but it's inherited from  MASS::rlm
glmrobMqle.control <-
    function(acc = 1e-04, test.acc = "coef", maxit = 50, tcc = 1.345)
{
    if (!is.numeric(acc) || acc <= 0)
	stop("value of acc must be > 0")
    if (test.acc != "coef")
	stop("Only 'test.acc = \"coef\"' is currently implemented")
    ## if (!(any(test.vec == c("coef", "resid"))))
    ##	  stop("invalid argument for test.acc")
    if (!is.numeric(maxit) || maxit <= 0)
	stop("maximum number of iterations must be > 0")
    if (!is.numeric(tcc) || tcc <= 0)
	stop("value of the tuning constant c (tcc) must be > 0")
    list(acc = acc, test.acc = test.acc, maxit = maxit, tcc = tcc)
}


### ----------------- E[ f(psi ( X ) ) ] -------------------------------

## MM (FIXME): These functions are always called "together" (2 or all 3)
##	       with *identical* argument vectors ---> save computations
##	       [all this was already in Martin's robGLM1 (!), see below
##   FIXME(2): Some of these fail when Huber's "c", 'tcc' is = +Inf
##   -----    --> ../../robGLM1/R/rglm.R

glmrobMqleEpsiPois <- function(mu, Vmu, ni, H, K, tcc) {
    h1 <- mu/sqrt(Vmu)*(dpois(H,mu)- dpois(K,mu))
    tcc*(1 - ppois(K,mu) - ppois(H,mu)) + h1
}

glmrobMqleEpsi2Pois <- function(mu, Vmu, ni, H, K, tcc)
{
    ## Calculation of E(psi^2) for the diagonal elements of A in
    ## matrix Q:
    Epsi2A <- tcc^2*(ppois(H,mu) + 1 - ppois(K,mu))
    Epsi2B <- mu*(dpois(H-1,mu) - dpois(H,mu) -
		  dpois(K-1,mu) + dpois(K,mu))
    Epsi2C <- ppois(K-1,mu)- ppois(H-1,mu)
    (Epsi2A  + Epsi2B + Epsi2C)
}

glmrobMqleEpsiSPois <- function(mu, Vmu, ni, H, K, tcc)
{
    ## Calculation of E(psi*s) for the diagonal elements of B in the
    ## expression matrix M = 1/n t(X) %*% B %*% X:
    EpsiSA <- tcc*(dpois(H,mu) + dpois(K,mu))
    EpsiSB <- mu/sqrt(Vmu)*(dpois(H-1,mu) - dpois(H,mu) -
			    dpois(K-1,mu) + dpois(K,mu))
    EpsiSC <- (ppois(K-1,mu) - ppois(H-1,mu))/sqrt(Vmu)
    (EpsiSA  + EpsiSB + EpsiSC)
}

glmrobMqleEpsiB <- function(mu, Vmu, ni, H, K, tcc)
{
    h1 <- ifelse(ni == 1, (- (H < 0) + (K >= 1) ) * sqrt(Vmu),
		 (pbinom(K-1,pmax(ni-1,1),mu) - pbinom(H-1,pmax(ni-1,1),mu)
		  - pbinom(K,ni,mu) + pbinom(H,ni,mu)) * mu * sqrt(ni/Vmu))
    ## pmax was needed to get numeric returns from pbinom
    ##
    tcc*(1 - pbinom(K,ni,mu) - pbinom(H,ni,mu)) + h1
}

glmrobMqleEpsi2B <- function(mu, Vmu, ni, H, K, tcc)
{
    ## Calculation of E(psi^2) for the diagonal elements of A in
    ## matrix Q:

    Epsi2AC <- tcc^2*(pbinom(H,ni,mu) + 1 - pbinom(K,ni,mu))
    Epsi2B1 <- mu^2*ni/Vmu*(pbinom(K,ni,mu)-pbinom(H,ni,mu))
    Epsi2B2 <- ((mu - 2*mu^2*ni)/Vmu *
		ifelse(ni == 1, (H <= 0)*(K >= 1),
		       (pbinom(K-1,pmax(ni-1,1),mu)-pbinom(H-1,pmax(ni-1,1),mu))))
    Epsi2B3 <- ((ni-1) * mu^2/Vmu*
		ifelse(ni == 2, (H <= 1)*(K >= 2),
		       (pbinom(K-2,pmax(ni-2,1),mu)-pbinom(H-2,pmax(ni-2,1),mu))))
    ##		 pmax was needed to get numeric returns from pbinom
    (Epsi2AC + Epsi2B1 +  Epsi2B2 +  Epsi2B3)
}

glmrobMqleEpsiSB <- function(mu, Vmu, ni, H, K, tcc)
{
    ## Calculation of E(psi*s) for the diagonal elements of B in the
    ## expression matrix M = 1/n t(X) %*% B %*% X:

    EpsiSA <- tcc*mu/Vmu*(pbinom(H,ni,mu) -
			  ifelse(ni == 1, H >= 1, pbinom(H-1,pmax(ni-1,1),mu)))
    EpsiSC <- tcc*mu/Vmu*(pbinom(K,ni,mu) -
			  ifelse(ni == 1, K > 0, pbinom(K-1,pmax(ni-1,1),mu)))
    EpsiSB1 <- mu^2*sqrt(ni)/(Vmu*sqrt(Vmu))*(pbinom(K,ni,mu)-pbinom(H,ni,mu))
    EpsiSB2 <- ((mu/sqrt(ni) - 2*mu^2*sqrt(ni))/(Vmu*sqrt(Vmu)) *
		ifelse(ni == 1, (H <= 0)*(K >= 1),
		       (pbinom(K-1,pmax(ni-1,1),mu)-pbinom(H-1,pmax(ni-1,1),mu))))
    EpsiSB3 <- ((ni-1) * mu^2/(Vmu*sqrt(ni*Vmu))*
		ifelse(ni == 2, (H <= 1)*(K >= 2),
		       (pbinom(K-2,pmax(ni-2,1),mu)-pbinom(H-2,pmax(ni-2,1),mu))))
    ##		 pmax was needed to get numeric returns from pbinom
    (EpsiSA + EpsiSC + EpsiSB1 +  EpsiSB2 +  EpsiSB3)
}

