#### Interface to intermediate level Fortran Routine
####
## Purpose: explore different versions of  dysta(), dysta2(), .. in
##          Fortran code.
dysta <- function(x, kind = c("euclidean","manhattan"))
{
  ## Purpose:
  ## -------------------------------------------------------------------------
  ## Arguments:
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date:  3 Sep 2002, 08:21

    kind <- arg.match(kind)
    n <- nrow(x <- as.matrix(x))
    p <- ncol(x)
    storage.mode(x) <- "double"
    hasNA <- apply(x, function(cl) any(is.na(cl)), 2)
    if(any(hasNA)) {

    }
    .Fortran("dysta",
             n,
             jp = p,
             x,
             dys= double(1 + n*(n-1)/2),
             ndyst=as.integer(kind == "euclidean"),
             jtmd=as.integer(hasNA)
             valmd=double()
             jhalt=integer(1)
             )

         }
