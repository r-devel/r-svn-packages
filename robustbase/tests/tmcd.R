library(robustbase)
## library(MASS)## only 'Animals' data  plus
##                 MASS::cov.mcd and MASS::mvrnorm {unused}

dodata <- function(nrep = 1, time = nrep >= 3, short = time, full = !short,
		   method = c("FASTMCD", "MASS")) {
    ##@bdescr
    ## Test the function covMcd() on the literature datasets:
    ##
    ## Call covMcd() for all regression datasets available in robustbase and print:
    ##  - execution time (if time)
    ##  - objective function
    ##  - best subsample found (if not short)
    ##  - outliers identified (with cutoff 0.975) (if not short)
    ##  - estimated center and covariance matrix if full)
    ##
    ##@edescr
    ##
    ##@in  nrep              : [integer] number of repetitions to use for estimating the
    ##                                   (average) execution time
    ##@in  time              : [boolean] whether to evaluate the execution time
    ##@in  short             : [boolean] whether to do short output (i.e. only the
    ##                                   objective function value). If short == FALSE,
    ##                                   the best subsample and the identified outliers are
    ##                                   printed. See also the parameter full below
    ##@in  full              : [boolean] whether to print the estimated cente and covariance matrix
    ##@in  method            : [character] select a method: one of (FASTMCD, MASS)

    domcd <- function(x, xname, nrep = 1) {
        n <- dim(x)[1]
        p <- dim(x)[2]
        if(method == "MASS") {
            mcd <- MASS::cov.mcd(x)
            quan <- as.integer(floor((n + p + 1)/2)) #default: floor((n+p+1)/2)
        }
        else {
            mcd <- covMcd(x) # trace = FALSE
            quan <- as.integer(mcd$quan)
        }

        crit <- if(method == "MASS") mcd$crit else log(mcd$crit)

        xres <- sprintf("%*s %3d %3d %3d %12.6f", lname, xname, n, p, quan, crit)
        if(time) {
            xtime <- system.time(dorep(x, nrep, method))[1]/nrep
            xres <- sprintf("%s %10.1f", xres, 1000 * xtime)
        }
        cat(xres, "\n")

        if(!short) {
            cat("Best subsample: \n")
            print(mcd$best)

            ibad <- which(mcd$mcd.wt == 0)
            names(ibad) <- NULL
            nbad <- length(ibad)
            cat("Outliers: ",nbad,"\n")
            if(nbad > 0)
                print(ibad)
            if(full) {
                cat("-------------\n")
                print(mcd)
            }
            cat("--------------------------------------------------------\n")
        }
    }

    lname <- 20
    method <- match.arg(method)

    data(heart)
    data(phosphor)
    data(starsCYG)
    data(stackloss)
    data(coleman)
    data(salinity)
    data(wood)
    data(hbk)

    data(Animals, package = "MASS")
    brain <- Animals[c(1:24, 26:25, 27:28),]
    data(milk)
    data(bushfire)

    ##    data(x1000)
    ##    data(x5000)

    tmp <- sys.call()
    cat("\nCall: ", deparse(substitute(tmp)),"\n")

    cat("Data Set               n   p  Half   LOG(obj)  Time [ms]\n")
    cat("========================================================\n")
    domcd(heart[, 1:2], data(heart), nrep)
    domcd(data.matrix(subset(phosphor, select = -plant)),
          data(phosphor), nrep)
    domcd(starsCYG, data(starsCYG), nrep)
    domcd(stack.x, data(stackloss), nrep)
    domcd(data.matrix(subset(coleman, select = -Y)), data(coleman), nrep)
    domcd(data.matrix(subset(salinity, select = -Y)), data(salinity), nrep)
    domcd(data.matrix(subset(wood, select = -y)), data(wood), nrep)
    domcd(data.matrix(subset(hbk,  select = -Y)), data(hbk), nrep)

    domcd(brain, "Animals", nrep)
    domcd(milk, data(milk), nrep)
    domcd(bushfire, data(bushfire), nrep)
    cat("========================================================\n")
    ##    domcd(x1000$X,data(x1000), nrep)
    ##    domcd(x5000$X,data(x5000), nrep)
}

if(FALSE) { ## the following functions are completely unused here ---------------

#### gendata() ####
## Generates a location contaminated multivariate
## normal sample of n observations in p dimensions
##    (1-eps)*Np(0,Ip) + eps*Np(m,Ip)
## where
##    m = (b,b,...,b)
## Defaults: eps=0 and b=10
##
gendata <- function(n,p,eps = 0,b = 10) {

    if(missing(n) || missing(p))
        stop("Please specify (n,p)")
    if(eps < 0 || eps >= 0.5)
        stop(message = "eps must be in [0,0.5)")
    X <- MASS::mvrnorm(n,rep(0,p),diag(1,nrow = p,ncol = p))
    nbad <- as.integer(eps * n)
    if(nbad > 0) {
        Xbad <- MASS::mvrnorm(nbad,rep(b,p),diag(1,nrow = p,ncol = p))
        xind <- sample(n,nbad)
        X[xind,] <- Xbad
    }
    list(X = X, xind = xind)
}

dogen <- function(nrep = 1, eps = 0.49, method = c("FASTMCD", "MASS")) {

    domcd <- function(x, nrep = 1) {
        xtime <- system.time(dorep(x, nrep, method))[1]/nrep
        cat(sprintf("%6d %3d %10.2f\n", dim(x)[1], dim(x)[2], xtime))
        xtime
    }

    set.seed(1234)

    method <- match.arg(method)

    ap <- c(2, 5, 10, 20, 30)
    an <- c(100, 500, 1000, 10000, 50000)

    tottime <- 0
    cat("     n   p       Time\n")
    cat("=====================\n")
    for(i in 1:length(an)) {
        for(j in 1:length(ap)) {
            n <- an[i]
            p <- ap[j]
            if(5*p <= n) {
                xx <- gendata(n, p, eps)
                X <- xx$X
                tottime <- tottime + domcd(X, nrep)
            }
        }
    }

    cat("=====================\n")
    cat("Total time: ", tottime*nrep, "\n")
}

docheck <- function(n, p, eps) {
    xx <- gendata(n,p,eps)
    mcd <- covMcd(xx$X)
    check(mcd, xx$xind)
}

check <- function(mcd, xind) {
##  check if mcd is robust w.r.t xind, i.e. check how many of xind
##  did not get zero weight
    mymatch <- xind %in% which(mcd$mcd.wt == 0)
    length(xind) - length(which(mymatch))
}

pad.right <- function(z, pads)
{
    ## Pads spaces to right of text
    padding <- paste(rep(" ", pads), collapse = "")
    paste(z, padding, sep = "")
}

}## the above are not used here -------------------------------------------------

dorep <- function(x, nrep = 1, method = c("FASTMCD","MASS")) {

    method <- match.arg(method)
    for(i in 1:nrep)
    if(method == "MASS")
        MASS::cov.mcd(x)
    else
        covMcd(x)
}


## -- now do it:
set.seed(101) # <<-- sub-sampling algorithm now based on R's RNG and seed
dodata()
dodata(nrep = 12)
dodata(nrep = 12, method = "MASS")

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
