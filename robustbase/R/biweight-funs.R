#### These Chi() and Psi() functions are currently used by lmrob() functions
#### FIXME: integrate these with the psi-rho objects --> ./psi-rho-funs.R

## Chi'() is just a scaled version of psi(); i.e. chi() = rho() ??

## With current scale (new for psi()):
##	 i)  Chi'(x, c) == (6/c^2) Psi(x,c)
## ==>	 ii) Chi''(x,c) == (6/c^2) Psi'(x,c)

lmrob.Chi <- function(x, cc, deriv = 0)
{
    x <- x / cc
    x2 <- x*x
    out <- x2 > 1
    switch(deriv + 1,
       {  ## deriv = 0
	   r <- x2*(3 + x2*(-3 + x2))
	   r[out] <- 1
       },
       {  ## deriv = 1
	   r <- 6/cc * x * (1-x2)^2
	   r[out] <- 0
       },
       {  ## deriv = 2
	   r <- 6/(cc^2) * (1 - x2) * (1 - 5*x2)
	   r[out] <- 0
       },
       stop("deriv must be in {0,1,2}"))
    r
}

lmrob.Psi <- function(x, cc, deriv = 0)
{
    ## This version of psi() is scaled such that psi'(0) = 1
    x2 <- (x / cc)^2
    r <- if(deriv == 0) x*(1 - x2)^2 else (1 - x2)*(1 - 5*x2)
    r[ x2 > 1 ] <- 0
    r
}
