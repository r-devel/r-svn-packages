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

adjOutlyingness <- function(x, ndir=250, clower=3, cupper=4,
                            alpha.cutoff = 0.75, coef = 1.5)
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
            if (qr(P, 1e-12)$rank == p) {
                B[i,] <- solve(P,matrix(1,nrow = p,ncol = 1))
                i <- i+1
            }
        }
        Bnorm <- matrix(0,nrow = ndir,ncol = 1)
        for (i in 1:ndir) {
            Bnorm[i] <- sum(abs(B[i,])^2)^(1/2)
        }
        Bnormr <- Bnorm[Bnorm > 1e-12]
        B      <-     B[Bnorm > 1e-12,]
 ## FIXME: no loop :
        for (i in 1:ndir) {
            A[i,] <- Bnormr[i]^(-1)*B[i,]
        }
    }
    else {
        stop('More dimensions than observations: not yet implemented')
        ## MM: In LIBRA(matlab) they have it implemented:
        ##    seed=0;
        ##    nrich1=n*(n-1)/2;
        ##    ndirect=min(250,nrich1);
        ##    true = (ndirect == nrich1);
        ##    B=extradir(x,ndir,seed,true); %n*ri
        ##      ======== % Calculates ndirect directions through
        ##               % two random choosen data points from data
        ##    for i=1:size(B,1)
        ##        Bnorm(i)=norm(B(i,:),2);
        ##    end
        ##    Bnormr=Bnorm(Bnorm > 1.e-12); %ndirect*1
        ##    B=B(Bnorm > 1.e-12,:);       %ndirect*n
        ##    A=diag(1./Bnormr)*B;         %ndirect*n

    }
    Y <- x %*% t(A)
    YZ <- matrix(0, nrow = n, ncol = ndir)
    tmc <- apply(Y, MARGIN = 2, mc)
    ##                          ==
    ## MM: Quart <- apply(Y, quantile, c(1:3)/4, names = FALSE)
    tme <- apply(Y, MARGIN = 2, median)
    tp3 <- apply(Y, MARGIN = 2, quantile,0.75)
    tp1 <- apply(Y, MARGIN = 2, quantile,0.25)
    tiq <- tp3-tp1
    tup <- -tme +(tp3 + coef*tiq*exp( cupper*tmc*(tmc >= 0)-clower*tmc*(tmc < 0)))
    tlo <- tme - (tp1 - coef*tiq*exp(-clower*tmc*(tmc >= 0)+cupper*tmc*(tmc < 0)))
    for (i in 1:ndir) {
        tup[i] <- -min(-Y[Y[,i] < (tup[i]+tme[i]),i])-tme[i]
        tlo[i] <- tme[i]-min(Y[Y[,i] > (tme[i]-tlo[i])])
    }
## FIXME (without loop!):
    for (j in 1:n) {
        YZ[j,] <- abs(Y[j,]-tme)/((Y[j,] >= tme)*tup+(Y[j,] < tme)*tlo)
    }
    adjout <- apply(YZ, MARGIN = 1, max)
    Qadj <- quantile(adjout, probs = c(1 - alpha.cutoff, alpha.cutoff))

    mcadjout <- mc(adjout)
    ##          ===
    if (mcadjout > 0) {
        cutoff <- Qadj[2]+ coef*exp(cupper*mcadjout)*(Qadj[2]-Qadj[1])
    }
    else {
        cutoff <- Qadj[2]+ coef*(Qadj[2]-Qadj[1])
    }

    list(adjout = adjout, cutoff = cutoff,
         flag = (adjout <= cutoff))
}
