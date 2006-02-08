glmrobMqle <-
    function(X, y, weights = rep(1, nobs), start = NULL,
             offset = rep(0, nobs), family, weights.on.x = "none",
             control = glmrobMqle.control(), intercept = TRUE)
{
    ## To DO:
    ## o weights are not really implemented
    ## o offset is not  implemented
    ##
    X <- as.matrix(X)
    xnames <- dimnames(X)[[2]]
    ynames <- names(y)
    conv <- FALSE
    nobs <- NROW(y)
    ncoef <- ncol(X)
    ##    EMPTY <- ncoef == 0
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
        stop("illegal `family' argument")
    valideta <- family$valideta
    if (is.null(valideta))
        valideta <- function(eta) TRUE
    validmu <- family$validmu
    if (is.null(validmu))
        validmu <- function(mu) TRUE
    ##
    if(weights.on.x > 0) {
        if(ncoef == 1)
            warning("There is only an intercept: No weights accepted on X")
        else {
            w.x <- switch(weights.on.x,
                          none = rep.int(1, nobs),
                          Hii = glmrobMqleHii(X = X),
                          robCov = glmrobMqleRobDist(X = X, intercept = intercept),
                          stop("This weighting method on X is not implemented"))
        }
    }
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
###
### Iterations
    for (nit in 1:control$maxit) {
        Vmu <- variance(mu)
        if (any(is.na(Vmu)))  stop("NAs in V(mu)")
        if (any(Vmu == 0))    stop("0s in V(mu)")
        dmu.deta <- mu.eta(eta)
        if (any(is.na(dmu.deta)))
            stop("NAs in d(mu)/d(eta)")
        residP <- (y-mu)*sqrt(ni/Vmu)
        H <- floor(mu*ni-tcc*sqrt(ni*Vmu))
        K <- floor(mu*ni+tcc*sqrt(ni*Vmu))
        ## Computation of alpha and (7) is done in one apply-loop:
        cpsi <- pmax(-tcc,pmin(residP,tcc)) - Epsi(mu, Vmu, ni, H, K, tcc)
        EEq <- colMeans(cpsi*w.x * sqrt(ni/Vmu) * dmu.deta*X)
        ##
        ## 1/n t(X) %*% B %*% X) %*%  delta.coef   = EEq
        DiagB <- EpsiS(mu, Vmu, ni, H, K, tcc) /sqrt(ni*Vmu) * w.x * (ni*dmu.deta)^2
        Dtheta  <- solve(crossprod(X, DiagB*X)/nobs, EEq)
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
        conv <-  (convi <= control$acc)
        if(conv)  break
        thetaOld <- theta
    } ## end of iteration
    ##
    if (!conv)
        warning("Algorithm did not converge")
    Vmu <- variance(mu)
    residP <- (y-mu)*sqrt(ni/Vmu)
    eps <- 10 * .Machine$double.eps
    if (family$family == "binomial") {
        if (any(mu/weights > 1 - eps) || any(mu/weights < eps))
            warning("fitted probabilities numerically 0 or 1 occurred")
    }
    if (family$family == "poisson") {
        if (any(mu < eps))
            warning("fitted rates numerically 0 occurred")
    }
###
    ## Estimated asymptotic covariance of the robust estimator
    dmu.deta <- mu.eta(eta)
    w.r <- pmin(1, tcc/abs(residP))
    weights <- w.x*w.r
    H <- floor(mu*ni - tcc*sqrt(ni*Vmu))
    K <- floor(mu*ni + tcc*sqrt(ni*Vmu))
    alpha <- colMeans(Epsi(mu, Vmu, ni, H, K, tcc) * w.x * sqrt(ni/Vmu)
                      * dmu.deta * X)
    DiagA <- Epsi2(mu, Vmu, ni, H, K, tcc) / (ni*Vmu) * w.x^2 * (ni*dmu.deta)^2
    matQ  <- crossprod(X, DiagA*X)/nobs  - outer(alpha, alpha)

    DiagB <- EpsiS(mu, Vmu, ni, H, K, tcc) /sqrt(ni*Vmu) * w.x * (ni*dmu.deta)^2
    matM <- crossprod(X, DiagB*X)/nobs
    matMinv  <- solve(matM)
    asCov <-  matMinv %*% matQ %*% matMinv / nobs
    ##
    list(coefficients = theta, residuals = residP, fitted.values = mu,
         w.r = w.r, w.x = w.x, ni = ni, cov = asCov, matM = matM, matQ = matQ, tcc = tcc,
         family = family, linear.predictors = eta, deviance = NULL,
         iter = nit, y = y, converged = conv)
}

