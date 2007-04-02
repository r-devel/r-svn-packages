#### This is originally from the R package
####
####  rrcov : Scalable Robust Estimators with High Breakdown Point
####
#### by Valentin Todorov

##  I would like to thank Peter Rousseeuw and Katrien van Driessen for
##  providing the initial code of this function.

### This program is free software; you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation; either version 2 of the License, or
### (at your option) any later version.
###
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
### GNU General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program; if not, write to the Free Software
### Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

## hidden in namespace:
quan.f <- function(alpha, n, p) {
    ## Compute h(alpha) := size of subsample, given alpha, (n,p)
    ## Same function for covMcd() and ltsReg()
    n2 <- (n+p+1) %/% 2
    floor(2 * n2 - n + 2 * (n - n2) * alpha)
}

covMcd <- function(x,
           cor = FALSE,
           alpha = 1/2,
           nsamp = 500,
           seed = NULL,
           trace = FALSE,
           use.correction = TRUE,
           control = rrcov.control())
{
    ##   Analyze and validate the input parameters ...

    ## if a control object was supplied, take the option parameters
    ## from it, but if single parameters were passed (not defaults)
    ## they will override the control object.
    ## defCtrl <- rrcov.control()# default control
    # 'control' specified
    if(missing(alpha))  alpha <- control$alpha
    if(missing(nsamp))  nsamp <- control$nsamp
    if(missing(seed))    seed <- control$seed
    if(missing(trace))  trace <- control$trace
    if(missing(use.correction)) use.correction <- control$use.correction

    tolSolve <- control$tolSolve # had 1e-10 hardwired {now defaults to 1e-14}

    if(length(seed) > 0) {
    if(exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE))  {
        seed.keep <- get(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
        on.exit(assign(".Random.seed", seed.keep, envir=.GlobalEnv))
    }
    assign(".Random.seed", seed, envir=.GlobalEnv)
    }

    ##   vt::03.02.2006 - added options "best" and "exact" for nsamp
    ##   nsamp will be further analized in the wrapper .fastmcd()
    if(!missing(nsamp) && is.numeric(nsamp) && nsamp <= 0)
    stop("Invalid number of trials nsamp = ",nsamp, "!")



    if(is.data.frame(x))
    x <- data.matrix(x)
    else if (!is.matrix(x))
    x <- matrix(x, length(x), 1,
            dimnames = list(names(x), deparse(substitute(x))))

    ## drop all rows with missing values (!!) :
    na.x <- !is.finite(x %*% rep(1, ncol(x)))
    ok <- !na.x
    x <- x[ok, , drop = FALSE]
    dx <- dim(x)
    if(!length(dx))
    stop("All observations have missing values!")
    dimn <- dimnames(x)
    n <- dx[1]
    p <- dx[2]
    ## h(alpha) , the size of the subsamples
    quan <- quan.f(alpha, n, p)
    if(n <= p + 1) # ==> floor((n+p+1)/2) > n - 1  -- not Ok
    stop(if (n <= p) # absolute barrier!
         "n <= p -- you can't be serious!"
    else "n == p+1  is too small sample size for MCD")
    ## else
    if(n < 2 * p) { ## p+1 < n < 2p
    warning("n < 2 * p, i.e., possibly too small sample size")
    ## was stop("Need at least 2*(number of variables) observations ")
    }
##     jmin <- (n + p + 1) %/% 2
##     if(alpha < 1/2) ## FIXME? shouldn't we rather test   'alpha < jmin/n' ?
##  stop("The MCD must cover at least", jmin, "observations")
    if(quan > n)
    stop("Sample size n  <  h(alpha; n,p) := size of \"good\" subsample")
    else if(alpha > 1) stop("alpha must be <= 1")

    quantiel <- qchisq(0.975, p)
    ## vt::03.02.2006 - raw.cnp2 and cnp2 are vectors of size 2 and  will
    ##   contain the correction factors (concistency and finite sample)
    ##   for the raw and reweighted estimates respectively. Set them
    ##   initially to 1.  If use.correction is set to FALSE
    ##   (default=TRUE), the finite sample correction factor will not
    ##   be used (neither for the raw estimates nor for the reweighted)
    raw.cnp2 <- cnp2 <- c(1,1)

    ans <- list(method = "Minimum Covariance Determinant Estimator.",
        call = match.call())

    if(alpha == 1) { ## alpha=1: Just compute the classical estimates --------
    mcd <- cov.wt(x)$cov
    loc <- as.vector(colMeans(x))
    obj <- determinant(mcd, log = TRUE)$modulus[1]
    if ( -obj/p > 50 ) {
        ans$cov <- mcd
        dimnames(ans$cov) <- list(dimn[[2]], dimn[[2]])
        if (cor)
        ans$cor <- cov2cor(ans$cov)
        ans$center <- loc
        if(length(dimn[[2]]))
        names(ans$center) <- dimn[[2]]
        ans$n.obs <- n
            ans$singularity <- list(kind = "classical")
            if(trace) cat("classical estimate is singular\n")

        weights <- 1
    }
    else {
        mah <- mahalanobis(x, loc, mcd, tol = tolSolve)
        ## VT:: 01.09.2004 - bug in alpha=1
        weights <- as.numeric(mah < quantiel) # 0/1
        sum.w <- sum(weights)
        ans <- c(ans, cov.wt(x, wt = weights, cor = cor))
        ans$cov <- sum.w/(sum.w - 1) * ans$cov

        ## Consistency factor for reweighted MCD
        if(sum.w != n) {
            ## VT::19.3.2007 - replace this code used several times by a function MCDcons(p, alpha)
            ## - for the reweighted cov use 'sum(weights)/n' instead of alpha
            ##
            ### qdelta.rew <- qchisq(sum.w/n, p)
            ### cdeltainvers.rew <- pgamma(qdelta.rew/2, p/2 + 1) / (sum.w/n)
            ### cnp2[1] <- 1/cdeltainvers.rew
            cnp2[1] <- MCDcons(p, sum.w/n)
            ans$cov <- ans$cov * cnp2[1]
        }
        if( - (determinant(ans$cov, log = TRUE)$modulus[1] - 0)/p > 50) {
                ans$singularity <- list(kind = "reweighted.MCD")
        if(trace) cat("reweighted MCD is singular\n")
        }
        else {
        mah <- mahalanobis(x, ans$center, ans$cov, tol = tolSolve)
        weights <- as.numeric(mah < quantiel) # 0/1
        }
    }

    ans$alpha <- alpha
    ans$quan <- quan
    ans$raw.cov <- mcd
    ans$raw.center <- loc
    if(!is.null(nms <- dimn[[2]])) {
        names(ans$raw.center) <- nms
        dimnames(ans$raw.cov) <- list(nms,nms)
    }
    ans$crit <- exp(obj)
    ans$method <- paste(ans$method,
                "\nThe minimum covariance determinant estimates based on",
                n, "observations \nare equal to the classical estimates.")
    ans$mcd.wt <- rep(NA, length(ok))
    ans$mcd.wt[ok] <- weights
    if(length(dimn[[1]]))
        names(ans$mcd.wt) <- dimn[[1]]
    ans$wt <- NULL
    ans$X <- x
    if(length(dimn[[1]]))
        dimnames(ans$X)[[1]] <- names(ans$mcd.wt)[ok]
    else
        dimnames(ans$X) <- list(seq(along = ok)[ok], NULL)
    if(trace)
        cat(ans$method, "\n")
    ans$raw.cnp2 <- raw.cnp2
    ans$cnp2 <- cnp2
    class(ans) <- "mcd"
    return(ans)
    } ## end {alpha=1} --

    mcd <- .fastmcd(x, quan, nsamp)

    ## Compute the consistency correction factor for the raw MCD
    ##  (see calfa in Croux and Haesbroeck)
    ## VT::19.3.2007 
    ### qalpha <- qchisq(quan/n, p)
    ### calphainvers <- pgamma(qalpha/2, p/2 + 1)/(quan/n)
    ### calpha <- 1/calphainvers
    calpha <- MCDcons(p, quan/n)
    correct <- if(use.correction) MCDcnp2(p, n, alpha) else 1.
    raw.cnp2 <- c(calpha, correct)

    if(p == 1) {
    ## ==> Compute univariate location and scale estimates
    ans$method <- "Univariate location and scale estimation."

    scale <- sqrt(calpha * correct) * as.double(mcd$initcovariance)
    center <- as.double(mcd$initmean)
    if(abs(scale - 0) < 1e-07) {
            ans$singularity <- list(kind = "identicalObs", q = quan)
        ans$raw.cov <- ans$cov <- matrix(0)
        ans$raw.center <- ans$center <- center
        ans$n.obs <- n
        ans$alpha <- alpha
        ans$quan <- quan
        if(!is.null(nms <- dimn[[2]][1])) {
        names(ans$raw.center) <- names(ans$center) <- nms
        dimnames(ans$raw.cov) <- dimnames(ans$cov) <- list(nms,nms)
        }
        ans$crit <- 0
        weights <- as.numeric(abs(x - center) < 1e-07) # 0 / 1
    } ## end { scale ~= 0 }
    else {
        ## Compute the weights for the raw MCD in case p=1
        weights <- as.numeric(((x - center)/scale)^2 < quantiel) # 0/1
        sum.w <- sum(weights)
        ans <- c(ans, cov.wt(x, wt = weights, cor = cor))
        ans$cov <- sum.w/(sum.w - 1) * ans$cov

        ## Apply the correction factor for the reweighted cov
        if(sum.w == n) {
        cdelta.rew <- 1
        correct.rew <- 1
        }
        else {
            ## VT::19.3.2007
            ### qdelta.rew <- qchisq(sum.w/n, p)
            ### cdeltainvers.rew <- pgamma(qdelta.rew/2, p/2 + 1)/(sum.w/n)
            ### cdelta.rew <- 1/cdeltainvers.rew
            cdelta.rew  <- MCDcons(p, sum.w/n)
            correct.rew <- if(use.correction) MCDcnp2.rew(p, n, alpha) else 1.
            cnp2 <- c(cdelta.rew, correct.rew)
        }
        ans$cov <- ans$cov * cdelta.rew * correct.rew
        ans$alpha <- alpha
        ans$quan <- quan
        ans$raw.cov <- as.matrix(scale^2)
        ans$raw.center <- as.vector(center)
        if(!is.null(nms <- dimn[[2]][1])) {
        dimnames(ans$raw.cov) <- list(nms,nms)
        names(ans$raw.center) <- nms
        }
        ans$crit <- 1/(quan - 1) *
                sum(sort((x - as.double(mcd$initmean))^2, partial = quan)[1:quan])
        center <- ans$center
        scale <- as.vector(sqrt(ans$cov))
        weights <- as.numeric(((x - center)/scale)^2 < quantiel)
    } ## end{ scale > 0 }
    } ## end p=1

    else { ## p >= 2 : ----------------------------------------------------------

    ## Apply correction factor to the raw estimates and use them to compute weights
    mcd$initcovariance <- calpha * mcd$initcovariance * correct
    dim(mcd$initcovariance) <- c(p, p)

    ## If not all observations are in general position, i.e. more than
    ## h observations lie on a hyperplane, the program still yields
    ## the MCD location and scatter matrix, the latter being singular
    ## (as it should be), as well as the equation of the hyperplane.
    if(mcd$exactfit != 0) {
        dim(mcd$coeff) <- c(5, p)
        ans$cov <- mcd$initcovariance
        ans$center <- as.vector(mcd$initmean)
        if(!is.null(nms <- dimn[[2]])) {
        dimnames(ans$cov) <- list(nms, nms)
        names(ans$center) <- nms
        }
        ans$n.obs <- n

## no longer relevant:
##      if(mcd$exactfit == -1)
##      stop("The program allows for at most ", mcd$kount, " observations.")
##      if(mcd$exactfit == -2)
##      stop("The program allows for at most ", mcd$kount, " variables.")
            if(!(mcd$exactfit %in% c(1,2)))
                stop("Unexpected 'exactfit' code ", mcd$exactfit, ". Please report!")

            ## new (2007-01) and *instead* of older long 'method' extension;
            ## the old message is stilled *printed* via singularityMessage()
            ##
            ## exactfit is now *passed* to result instead of coded into 'message':
            ans$singularity <-
                list(kind = "on.hyperplane", exactCode = mcd$exactfit,
                     p = p, count = mcd$kount, coeff = mcd$coeff[1,])
        ans$alpha <- alpha
        ans$quan <- quan
        ans$raw.cov <- mcd$initcovariance
        ans$raw.center <- as.vector(mcd$initmean)
        if(!is.null(nms <- dimn[[2]])) {
        names(ans$raw.center) <- nms
        dimnames(ans$raw.cov) <- list(nms,nms)
        }
        ans$crit <- 0
        weights <- mcd$weights
    } ## end exact fit <==>  (mcd$exactfit != 0)

    else { ## exactfit == 0 : have general position ------------------------

            ## FIXME: here we assume that mcd$initcovariance is not singular
            ## ----- but it is for data(mortality, package = "riv") !
        mah <- mahalanobis(x, mcd$initmean, mcd$initcovariance, tol = tolSolve)
        mcd$weights <- weights <- as.numeric(mah < quantiel)
        sum.w <- sum(weights)

        ## Compute and apply the consistency correction factor for
        ## the reweighted cov
        if(sum.w == n) {
        cdelta.rew <- 1
        correct.rew <- 1
        }
        else {
        ## VT::19.3.2007
        ##
        ### qdelta.rew <- qchisq(sum.w/n, p)
        ### cdeltainvers.rew <- pgamma(qdelta.rew/2, p/2 + 1)/(sum.w/n)
        ### cdelta.rew <- 1/cdeltainvers.rew
        cdelta.rew <- MCDcons(p, sum.w/n)
        correct.rew <- if(use.correction) MCDcnp2.rew(p, n, alpha) else 1.
                cnp2 <- c(cdelta.rew, correct.rew)
        }

        ans <- c(ans, cov.wt(x, wt = weights, cor))
        ans$cov <- sum.w/(sum.w - 1) * ans$cov
        ans$cov <- ans$cov * cdelta.rew * correct.rew

        ##vt:: add also the best found subsample to the result list
        ans$best <- sort(as.vector(mcd$best))

        ans$alpha <- alpha
        ans$quan <- quan
        ans$raw.cov <- mcd$initcovariance
        ans$raw.center <- as.vector(mcd$initmean)
        if(!is.null(nms <- dimn[[2]])) {
        names(ans$raw.center) <- nms
        dimnames(ans$raw.cov) <- list(nms,nms)
        }
        ans$raw.weights <- weights
        ans$crit <- mcd$mcdestimate
        ans$raw.mah <- mahalanobis(x, ans$raw.center, ans$raw.cov, tol = tolSolve)

        ## Check if the reweighted scatter matrix is singular.
        if( - (determinant(ans$cov, log = TRUE)$modulus[1] - 0)/p > 50) {
                ans$singularity <- list(kind = "reweighted.MCD")
        if(trace) cat("The reweighted MCD scatter matrix is singular.\n")
        ans$mah <- ans$raw.mah
        }
        else {
        mah <- mahalanobis(x, ans$center, ans$cov, tol = tolSolve)
        ans$mah <- mah
        weights <- as.numeric(mah < quantiel)
        }
    } ## end{ not exact fit }

    } ## end{ p >= 2 }

    ans$mcd.wt <- rep(NA, length(ok))
    ans$mcd.wt[ok] <- weights
    if(length(dimn[[1]]))
    names(ans$mcd.wt) <- dimn[[1]]
    ans$wt <- NULL
    ans$X <- x
    if(length(dimn[[1]]))
    dimnames(ans$X)[[1]] <- names(ans$mcd.wt)[ok]
    else
    dimnames(ans$X) <- list(seq(along = ok)[ok], NULL)
    if(trace)
    cat(ans$method, "\n")
    ans$raw.cnp2 <- raw.cnp2
    ans$cnp2 <- cnp2
    class(ans) <- "mcd"
    return(ans)
}

singularityMsg <- function(singList, n.obs)
{
    stopifnot(is.list(singList))
    switch(singList$kind,
       "classical" = {
           "The classical covariance matrix is singular."
       },
       "reweighted.MCD" = {
           "The reweighted MCD scatter matrix is singular."
       },
       "identicalObs" = {
           sprintf("Initial scale 0 because more than 'quan' (=%d) observations are identical.",
               singList$q)
       },
       "on.hyperplane" = {
           stopifnot(c("p", "count", "coeff") %in% names(singList))

           obsMsg <- function(m, n)
           paste("There are", m,
             "observations (in the entire dataset of",
             n, "obs.) lying on the")
           with(singList,
                    c(switch(exactCode,
                             ## exactfit == 1 :
                             "The covariance matrix of the data is singular.",
                             ## exactfit == 2 :
                             c("The covariance matrix has become singular during",
                               "the iterations of the MCD algorithm.")),

                      if(p == 2) {
                          paste(obsMsg(count, n.obs), "line with equation ",
                                signif(coeff[1], digits= 5), "(x_i1-m_1) +",
                                signif(coeff[2], digits= 5), "(x_i2-m_2) = 0",
                                "with (m_1,m_2) the mean of these observations.")
                      }
                      else if(p == 3) {
                          paste(obsMsg(count, n.obs), "plane with equation ",
                                signif(coeff[1], digits= 5), "(x_i1-m_1) +",
                                signif(coeff[2], digits= 5), "(x_i2-m_2) +",
                                signif(coeff[3], digits= 5), "(x_i3-m_3) = 0",
                                "with (m_1,m_2) the mean of these observations."
                                )
                      }
                      else { ##  p > 3 -----------
                          con <- textConnection("astring", "w")
                          dput(zapsmall(coeff), con)
                          close(con)
                          paste(obsMsg(count, n.obs), "hyperplane with equation ",
                                "a_1*(x_i1 - m_1) + ... + a_p*(x_ip - m_p) = 0",
                                " with (m_1,...,m_p) the mean of these observations",
                                " and coefficients a_i from the vector   a <- ", astring)
                      }))
       },
       ## Otherwise
       stop("illegal 'singularity$kind'")
       ) ## end{switch}
}

print.mcd <- function(x, digits = max(3, getOption("digits") - 3), print.gap = 2, ...)
{
    cat("Minimum Covariance Determinant (MCD) estimator.\n")
    if(!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl)
    }
    cat("-> Method: ", x$method, "\n")
    if(is.list(x$singularity))
        cat(strwrap(singularityMsg(x$singularity, x$n.obs)), sep ="\n")

    ## VT::29.03.2007 - solve a conflict with fastmcd() in package robust - 
    ##      also returning an object of class "mcd"
    xx <- NA
    if(!is.null(x$crit))
        xx <- format(log(x$crit), digits = digits)
    else if (!is.null(x$raw.objective))
        xx <- format(log(x$raw.objective), digits = digits)
        
    cat("\nLog(Det.): ", xx ,"\n")
    cat("Robust Estimate of Location:\n")
    print(x$center, digits = digits, print.gap = print.gap, ...)
    cat("Robust Estimate of Covariance:\n")
    print(x$cov, digits = digits, print.gap = print.gap, ...)
    invisible(x)
}

summary.mcd <- function(object, ...)
{
    class(object) <- c("summary.mcd", class(object))
    object
}

print.summary.mcd <-
    function(x, digits = max(3, getOption("digits") - 3), print.gap = 2, ...)
{
    print.mcd(x, digits = digits, print.gap = print.gap, ...) # see above

    ## hmm, maybe not *such* a good idea :
    if(!is.null(x$cor)) {
      cat("\nRobust Estimate of Correlation: \n")
      dimnames(x$cor) <- dimnames(x$cov)
      print(x$cor, digits = digits, print.gap = print.gap, ...)
    }

    cat("\nEigenvalues:\n")
    print(eigen(x$cov, only.values = TRUE)$values, digits = digits, ...)

    if(!is.null(x$mah)) {
    cat("\nRobust Distances: \n")
    print(as.vector(x$mah), digits = digits, ...)
    }
    invisible(x)
}

## NOTE:  plot.mcd() is in ./covPlot.R !
## ----                    ~~~~~~~~~~~

### --- Namespace hidden (but parsed once and for all) : -------------
MCDcons <- function(p, alpha){
    ## VT::19.3.2007 - replace the code used several times by a function MCDcons(p, alpha)
    ##
    ## Compute the consistency correction factor for the MCD estimate
    ##  (see calfa in Croux and Haesbroeck)
    ##  - alpha = h/n = quan/n
    ##  - use the same function for the reweighted estimates, 
    ##      but instead of 'alpha' call with 'sum(weights)/n'
    
    qalpha <- qchisq(alpha, p)
    calphainvers <- pgamma(qalpha/2, p/2 + 1)/alpha
    1/calphainvers
}

MCDcnp2 <- function(p, n, alpha)
{
    if(p > 2) {
    ##              "alfaq"        "betaq"    "qwaarden"
    coeffqpkwad875 <- matrix(c(-0.455179464070565, 1.11192541278794, 2,
                   -0.294241208320834, 1.09649329149811, 3), ncol = 2)
    coeffqpkwad500 <- matrix(c(-1.42764571687802,  1.26263336932151, 2,
                   -1.06141115981725,  1.28907991440387, 3), ncol = 2)

    y.500 <- log( - coeffqpkwad500[1, ] / p^coeffqpkwad500[2, ] )
    y.875 <- log( - coeffqpkwad875[1, ] / p^coeffqpkwad875[2, ] )

    A.500 <- cbind(1, - log(coeffqpkwad500[3, ] * p^2))
    A.875 <- cbind(1, - log(coeffqpkwad875[3, ] * p^2))
    coeffic.500 <- solve(A.500, y.500)
    coeffic.875 <- solve(A.875, y.875)
    fp.500.n <- 1 - exp(coeffic.500[1]) / n^coeffic.500[2]
    fp.875.n <- 1 - exp(coeffic.875[1]) / n^coeffic.875[2]
    }
    else { ## p <= 2
    if(p == 2) {
        fp.500.n <- 1 - exp( 0.673292623522027) / n^0.691365864961895
        fp.875.n <- 1 - exp( 0.446537815635445) / n^1.06690782995919
    }
    if(p == 1) {
        fp.500.n <- 1 - exp( 0.262024211897096) / n^0.604756680630497
        fp.875.n <- 1 - exp(-0.351584646688712) / n^1.01646567502486
    }
    }

    stopifnot(0.5 <= alpha, alpha <= 1)
    if(alpha <= 0.875)
    fp.alpha.n <- fp.500.n + (fp.875.n - fp.500.n)/0.375 * (alpha - 0.5)
    else ##  0.875 < alpha <= 1
    fp.alpha.n <- fp.875.n + (1 - fp.875.n)/0.125 * (alpha - 0.875)

    return(1/fp.alpha.n)
} ## end{ MCDcnp2 }

MCDcnp2.rew <- function(p, n, alpha)
{
    if(p > 2) {
    ##              "alfaq"        "betaq"    "qwaarden"
    coeffrewqpkwad875 <- matrix(c(-0.544482443573914, 1.25994483222292, 2,
                      -0.343791072183285, 1.25159004257133, 3), ncol = 2)
    coeffrewqpkwad500 <- matrix(c(-1.02842572724793,  1.67659883081926, 2,
                      -0.26800273450853,  1.35968562893582, 3), ncol = 2)

    y.500 <- log( - coeffrewqpkwad500[1, ] / p^ coeffrewqpkwad500[2, ] )
    y.875 <- log( - coeffrewqpkwad875[1, ] / p^ coeffrewqpkwad875[2, ] )

    A.500 <- cbind(1, - log(coeffrewqpkwad500[3, ] * p^2))
    coeffic.500 <- solve(A.500, y.500)
    A.875 <- cbind(1, - log(coeffrewqpkwad875[3, ] * p^2))
    coeffic.875 <- solve(A.875, y.875)
    fp.500.n <- 1 - exp(coeffic.500[1]) / n^ coeffic.500[2]
    fp.875.n <- 1 - exp(coeffic.875[1]) / n^ coeffic.875[2]
    }
    else {
    if(p == 2) {
        fp.500.n <- 1 - exp( 3.11101712909049 ) / n^ 1.91401056721863
        fp.875.n <- 1 - exp( 0.79473550581058 ) / n^ 1.10081930350091
    }
    if(p == 1) {
        fp.500.n <- 1 - exp( 1.11098143415027 ) / n^ 1.5182890270453
        fp.875.n <- 1 - exp( -0.66046776772861) / n^ 0.88939595831888
    }
    }

    stopifnot(0.5 <= alpha, alpha <= 1)
    if(alpha <= 0.875)
    fp.alpha.n <- fp.500.n + (fp.875.n - fp.500.n)/0.375 * (alpha - 0.5)
    else ##  0.875 < alpha <= 1
    fp.alpha.n <- fp.875.n + (1 - fp.875.n)/0.125 * (alpha - 0.875)
    return(1/fp.alpha.n)
} ## end{ MCDcnp2.rew }


.fastmcd <- function(x, quan, nsamp)
{
    dx <- dim(x)
    n <- dx[1]
    p <- dx[2]

    ##   parameters for partitioning {equal to those in Fortran !!}
    kmini <- 5
    nmini <- 300
    km10 <- 10*kmini
    nmaxi <- nmini*kmini

    ##   vt::03.02.2006 - added options "best" and "exact" for nsamp
    if(!missing(nsamp)) {
    if(is.numeric(nsamp) && (nsamp < 0 || nsamp == 0 && p > 1)) {
        warning("Invalid number of trials nsamp= ", nsamp,
                    " ! Using default.\n")
        nsamp <- -1
    } else if(nsamp == "exact" || nsamp == "best") {
        myk <- p + 1 ## was 'p'; but p+1 ("nsel = nvar+1") is correct
        if(n > 2*nmini-1) {
                ## TODO: make 'nmini' user configurable; however that
                ##      currently needs changes to the fortran source
                ##  and re-compilation!
        warning("Options 'best' and 'exact' not allowed for n greater than ",
            2*nmini-1,".\nUsing default.\n")
        nsamp <- -1
        } else {
        nall <- choose(n, myk)
        if(nall > 5000 && nsamp == "best") {
            nsamp <- 5000
            warning("'nsamp = \"best\"' allows maximally 5000 subsets;\n",
                "computing these subsets of size ",
                            myk," out of ",n,"\n")
        } else { ## "exact" or ("best"  &  nall < 5000)
            nsamp <- 0 ## all subsamples
            if(nall > 5000)
            warning("Computing all ", nall, " subsets of size ",myk,
                " out of ",n,
                "\n This may take a very long time!\n",
                immediate. = TRUE)
        }
        }
    }

    if(!is.numeric(nsamp) || nsamp == -1) { # still not defined - set it to the default
        defCtrl <- rrcov.control() # default control
        if(!is.numeric(nsamp))
        warning("Invalid number of trials nsamp= ",nsamp,
            " ! Using default nsamp= ",defCtrl$nsamp,"\n")
        nsamp <- defCtrl$nsamp  # take the default nsamp
    }
    }

    ##   Allocate temporary storage for the Fortran implementation,
    ##   directly in the .Fortran() call.
    ##    (if we used C, we'd rather allocate there, and be quite faster!)

    .Fortran("rffastmcd",
         x = if(is.double(x)) x else as.double(x),
         n =    as.integer(n),
         p =    as.integer(p),   ## = 'nvar'  in Fortran
         nhalff =   as.integer(quan),
         nsamp  =   as.integer(nsamp),# = 'krep'
         initcovariance = double(p * p),
         initmean       = double(p),
         best       = rep.int(as.integer(10000), quan),
         mcdestimate = double(1),    ## = 'det'
         weights   = integer(n),
         exactfit  = integer(1), # output indicator: 0: ok; 1: ..., 2: ..
         coeff     = matrix(double(5 * p), nrow = 5, ncol = p), ## plane
         kount     = integer(1),
         adjustcov = double(p * p), ## << never used -- FIXME
         integer(1),## << 'seed' no longer used -- FIXME
         temp   = integer(n),
         index1 = integer(n),
         index2 = integer(n),
         nmahad = double(n),
         ndist  = double(n),
         am     = double(n),
         am2    = double(n),
         slutn  = double(n),

         med   = double(p),
         mad   = double(p),
         sd    = double(p),
         means = double(p),
         bmeans= double(p),
         w     = double(p),
         fv1   = double(p),
         fv2   = double(p),

         rec   = double(p+1),
         sscp1 = double((p+1)*(p+1)),
         cova1 = double(p * p),
         corr1 = double(p * p),
         cinv1 = double(p * p),
         cova2 = double(p * p),
         cinv2 = double(p * p),
         z     = double(p * p),

         cstock = double(10 * p * p),# (10,nvmax2)
         mstock = double(10 * p),    # (10,nvmax)
         c1stock = double(km10 * p * p), # (km10,nvmax2)
         m1stock = double(km10 * p), # (km10,nvmax)
         dath = double(nmaxi * p),   # (nmaxi,nvmax)

         cutoff = qchisq(0.975, p),
         chimed = qchisq(0.5,   p),

         PACKAGE = "robustbase")[ ## keep the following ones:
         c("initcovariance", "initmean", "best", "mcdestimate",
           "weights", "exactfit", "coeff", "kount", "adjustcov") ]
}
