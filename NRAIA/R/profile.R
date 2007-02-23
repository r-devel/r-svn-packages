## functions to work with profile.nls objects
as.data.frame.profile.nls <-
    function(x, row.names = NULL, optional = FALSE, ...)
{
    df <- do.call("rbind", lapply(x, "[[", "par.vals"))
    tau <- lapply(x, "[[", "tau")
    pnames <- factor(rep(names(x), sapply(tau, length)), levels = names(x))
    pars <- df[cbind(seq_len(nrow(df)),
                     match(as.character(pnames), colnames(df)))]
    df <- as.data.frame(df)
    df$.par <- pars
    df$.pnm <- pnames
    df$.tau <- unlist(tau)
    df
}
    
plot.profile.nls <-
    function (x, levels, conf = c(99, 95, 90, 80, 50)/100, 
              absVal = TRUE, ...) 
{
    dfres <- attr(x, "summary")$df[2]
    confstr <- NULL
    if (missing(levels)) {
        levels <- sqrt(qf(pmax(0, pmin(1, conf)), 1, dfres))
        confstr <- paste(format(100 * conf), "%", sep = "")
    }
    if (any(levels <= 0)) {
        levels <- levels[levels > 0]
        warning("levels truncated to positive values only")
    }
    if (is.null(confstr)) {
        confstr <- paste(format(100 * pf(levels^2, 1, dfres)), 
            "%", sep = "")
    }
    levels <- sort(levels)
    spl <- lapply(x, function(x)
                  splines::interpSpline(x$par.vals[, attr(x, "parameters")$par],
                                        x$tau))
    bspl <- lapply(spl, splines::backSpline)
    tau <- c(-rev(levels), 0, levels)
    df <- data.frame(tau = rep.int(tau, length(x)),
                     pval = unlist(lapply(bspl,
                     function(sp) predict(sp, x= tau)$y)),
                     pnm = gl(length(x), length(tau), labels = names(x)))
    ylab <- expression(tau)
    if (absVal) {
        df$tau <- abs(df$tau)
        ylab <- expression("|" * tau * "|")
    }
    xyplot(tau ~ pval | pnm, df,
           scales = list(x = list(relation = 'free'), y = list(rot = 0)),
           ylab = ylab, xlab = "", panel = function(x, y, ...)
       {
           pfun <- function(x) predict(spl[[panel.number()]], x = x)$y
           panel.grid(h = -1, v = -1)
           lsegments(x, y, x, 0, ...)
           if (absVal) {
               lsegments(x, y, rev(x), y)
               pfun <- function(x) abs(predict(spl[[panel.number()]], x = x)$y)
           } else {
               panel.abline(h = 0, ...)
           }
           panel.curve(pfun, ...)
       }, ...)
}
