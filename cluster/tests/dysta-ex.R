#### Interface to intermediate level Fortran Routines
#### Purpose : explore different versions of  dysta(), dysta3(), .. in
####	       Fortran code.

library(cluster)

if(FALSE) # we no longer C-export dysta() and dysta3
dysta <- function(x, kind = c("euclidean","manhattan", "SqEuclidean"),
                  dystaK = "dysta")
{
    ## Purpose:
    ## -------------------------------------------------------------------------
    ## Arguments:
    ## -------------------------------------------------------------------------
    ## Author: Martin Maechler, Date:  3 Sep 2002, 08:21

    kind <- match.arg(kind)
    ndyst <- which(kind == eval(formals()$kind))# 1, 2, or 3
    n <- nrow(x <- as.matrix(x))
    p <- ncol(x)
    storage.mode(x) <- "double"
    hasNA <- apply(is.na(x), 2, any) # == apply(x, 2, anyNA)
    if(any(hasNA)) {
	ina <- is.na(x)
	x[ina] <- valmd <- -1.1*max(abs(range(x, na.rm = TRUE)))
	valmd <- rep(valmd, p)
    } else valmd <- 0.

    dys <- double(1 + n*(n-1)/2)
    jtmd <- as.integer(ifelse(hasNA, -1, 1))

    r <-
        if(dystaK == "dysta3") {
            .C(cluster:::dysta3,
               n,
               jp = p,
               x,
               dys= dys,
               ndyst= ndyst,
               jtmd=  jtmd,
               valmd= valmd,
	       jhalt= integer(1))[c("dys", "jhalt")]
        } else {
            .C(cluster:::dysta,
                     n,
                     jp = p,
                     x,
                     dys= dys,
                     ndyst= ndyst,
                     jtmd=  jtmd,
                     valmd= valmd,
		     jhalt= integer(1))[c("dys", "jhalt")]
        }
    if(r$jhalt) {
	cat("'jhalt' was ", r$jhalt,
	    " -- some dissimilarities will be missing.\n")
	r$dys[r$dys == -1.] <- NA
    }
    r$dys
}

(x <- cbind(c(0:6,NA), c(1,2,NA,7,NA,8:9,8)))
if(FALSE) {
dysta(x)
(d1 <- dysta(x, kind = "m"))
(d3 <- dysta(x, kind = "m",dystaK = "dysta3"))
identical(sort(d1), sort(d3)) # TRUE
}
cbind(# d1=d1[-1], d3=d3[-length(d3)],
      dist=dist(x,"manhattan"), daisy= daisy(x,"manhattan"))

if(FALSE)
identical(d3[-length(d3)],
	  c(dist(x,"manhattan")))# !
identical(c(daisy(x,"manhattan")), c(dist(x,"manhattan")))
identical(c(daisy(x,"euclidean")), c(dist(x,"euclidean")))
if(FALSE)
identical(dysta(x, dystaK="dysta3")[-length(d3)],
	  c(dist(x,"euclidean")))# !
