glmrobMqleHii <- function(X) {
    x <- qr(X)
    Hii <- colSums(qr.qy(x, diag(1, nrow = n, ncol = x$rank))^2)
    sqrt(1-Hii)
}

glmrobMqleEpsiPois <- function(mu, Vmu, ni, H, K, tcc) {
    h1 <- mu/sqrt(Vmu)*(dpois(H,mu)- dpois(K,mu))
    tcc*(1 - ppois(K,mu) - ppois(H,mu)) + h1
}


glmrobMqleEpsi2Pois <- function(mu, Vmu, ni, H, K, tcc) {
    ##
    ## Calculation of E(psi^2) for the diagnonal elements of A in
    ## matrix Q:
    Epsi2A <- tcc^2*(ppois(H,mu) + 1 - ppois(K,mu))
    Epsi2B <- mu*(dpois(H-1,mu) - dpois(H,mu) -
                  dpois(K-1,mu) + dpois(K,mu))
    Epsi2C <- (ppois(K-1,mu)- ppois(H-1,mu))
    (Epsi2A  + Epsi2B + Epsi2C)
}


glmrobMqleEpsiB <- function(mu, Vmu, ni, H, K, tcc) {
    h1 <- ifelse(ni == 1, (- (H < 0) + (K >= 1) ) * sqrt(Vmu),
                 (pbinom(K-1,pmax(ni-1,1),mu) - pbinom(H-1,pmax(ni-1,1),mu)
                  - pbinom(K,ni,mu) + pbinom(H,ni,mu)) * mu * sqrt(ni/Vmu))
    ## pmax was needed to get numeric returns from pbinom
    ##
    tcc*(1 - pbinom(K,ni,mu) - pbinom(H,ni,mu)) + h1
}


glmrobMqleEpsiSPois <- function(mu, Vmu, ni, H, K, tcc) {
    ##
    ## Calculation of E(psi*s) for the diagnonal elements of B in the
    ## expression matrix M = 1/n t(X) %*% B %*% X:
    EpsiSA <- tcc*(dpois(H,mu) + dpois(K,mu))
    EpsiSB <- mu/sqrt(Vmu)*(dpois(H-1,mu) - dpois(H,mu) -
                            dpois(K-1,mu) + dpois(K,mu))
    EpsiSC <- (ppois(K-1,mu)- ppois(H-1,mu))/sqrt(Vmu)
    (EpsiSA  + EpsiSB + EpsiSC)
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


glmrobMqle.control <- function (acc = 1e-04, test.acc = "coef", maxit = 50,
                                 tcc = 1.345)
{
    if (!is.numeric(acc) || acc <= 0)
        stop("value of acc must be > 0")
    if (test.acc != "coef")
        stop("Sorry, only 'coef' is currently implemented")
    ## if (!(any(test.vec == c("coef", "resid"))))
    ##    stop("invalid argument for test.acc")
    if (!is.numeric(maxit) || maxit <= 0)
        stop("maximum number of iterations must be > 0")
    if (!is.numeric(tcc) || tcc <= 0)
        stop("value of the tuning constant c (tcc) must be > 0")
    list(acc = acc, test.acc = test.acc, maxit = maxit, tcc = tcc)
}


glmrobMqleEpsi2B <- function(mu, Vmu, ni, H, K, tcc) {
    ##
    ## Calculation of E(psi^2) for the diagnonal elements of A in
    ## matrix Q:
    Epsi2AC <- tcc^2*(pbinom(H,ni,mu) + 1 - pbinom(K,ni,mu))
    Epsi2B1 <- mu^2*ni/Vmu*(pbinom(K,ni,mu)-pbinom(H,ni,mu))
    Epsi2B2 <- ((mu - 2*mu^2*ni)/Vmu *
                ifelse(ni == 1, (H <= 0)*(K >= 1),
                       (pbinom(K-1,pmax(ni-1,1),mu)-pbinom(H-1,pmax(ni-1,1),mu))))
    Epsi2B3 <- ((ni-1) * mu^2/Vmu*
                ifelse(ni == 2, (H <= 1)*(K >= 2),
                       (pbinom(K-2,pmax(ni-2,1),mu)-pbinom(H-2,pmax(ni-2,1),mu))))
    ##           pmax was needed to get numeric returns from pbinom
    (Epsi2AC + Epsi2B1 +  Epsi2B2 +  Epsi2B3)
}

glmrobMqleEpsiSB <- function(mu, Vmu, ni, H, K, tcc)
{
    ##
    ## Calculation of E(psi*s) for the diagnonal elements of B in the
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

