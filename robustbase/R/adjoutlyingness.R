#### MC-adjusted Outlyingness
#### ------------------------
###
### Original code from  the web site from the Antwerpen Statistics Group :
###  http://www.agoras.ua.ac.be/Robustn.htm
### which has a "MC" section and for the software links to
### ftp://ftp.win.ua.ac.be/pub/software/agoras/newfiles/mc.tar.gz
### and that contains  mcrsoft/adjoutlyingness.R

## MM [ FIXME ]:
## -----------

## 1)   Use  *transposed*  B[] and A[] matrices

## 2)   use  IQR() instead of   quantile(., .75) - quantile(., .25)

##-->  but only *after* testing original code
##     ^^^^^^^^^^^^^^^^^^^^^^^^

adjOutlyingness <- function(x, ndir=250, clower=3, cupper=4)
## Skewness-Adjusted Outlyingness
{
    x <- data.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    stopifnot(n > 0, p > 0, is.numeric(x))
    if (p < n) {
        i <- 1
        it <- 0
        B <- matrix(0,nrow = ndir,ncol = p)
        A <- matrix(0,nrow = ndir,ncol = p)
        while (i <= ndir  &&  (it <- it+1) < 1000) {
            xx <- sample(n)
            P <- x[xx[1:p],]
            if (qr(P, 0.000000000001)$rank == p) {
                B[i,] <- solve(P,matrix(1,nrow = p,ncol = 1))
                i <- i+1
            }
        }
        Bnorm <- matrix(0,nrow = ndir,ncol = 1)
        for (i in 1:ndir) {
            Bnorm[i] <- sum(abs(B[i,])^2)^(1/2)
        }
        Bnormr <- Bnorm[Bnorm > 0.000000000001]
        B      <-     B[Bnorm > 0.000000000001,]
        for (i in 1:ndir) {
            A[i,] <- Bnormr[i]^(-1)*B[i,]
        }
    }
    else {
        stop('More dimensions than observations, currently not implemented')
    }
    Y <- x %*% t(A)
    YZ <- matrix(0,nrow = n,ncol = ndir)
browser()
    tmc <- mc(Y)
    ##     ===
    tme <- apply(Y,MARGIN = 2,median)
    tp3 <- apply(Y,MARGIN = 2,quantile,0.75)
    tp1 <- apply(Y,MARGIN = 2,quantile,0.25)
    tiq <- tp3-tp1
    tup <- (tp3+1.5*exp(cupper*tmc*(tmc >= 0)-clower*tmc*(tmc < 0))*tiq)-tme
    tlo <- tme-(tp1-1.5*exp(-clower*tmc*(tmc >= 0)+cupper*tmc*(tmc < 0))*tiq)
    for (i in 1:ndir) {
        tup[i] <- -min(-Y[Y[,i] < (tup[i]+tme[i]),i])-tme[i]
        tlo[i] <- tme[i]-min(Y[Y[,i] > (tme[i]-tlo[i])])
    }
    for (j in 1:n) {
        YZ[j,] <- abs(Y[j,]-tme)/((Y[j,] >= tme)*tup+(Y[j,] < tme)*tlo)
    }
    adjout <- apply(YZ,MARGIN = 1,max)
    mcadjout <- mc(adjout)
    ##          ===
    if (mcadjout > 0) {
        cutoff <- quantile(adjout,0.75)+1.5*exp(cupper*mcadjout)*(quantile(adjout,0.75)-quantile(adjout,0.25))
    }
    else {
        cutoff <- quantile(adjout,0.75)+1.5*(quantile(adjout,0.75)-quantile(adjout,0.25))
    }
    flag <- (adjout <= cutoff)
    res <- list(adjout = adjout,cutoff = cutoff,flag = flag)
    res
}
