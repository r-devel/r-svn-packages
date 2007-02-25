## functions and methods for profile.nls objects

as.data.frame.profile.nls <-
    function(x, row.names = NULL, optional = FALSE, ...)
{
    fr <- do.call("rbind", lapply(x, "[[", "par.vals"))
    tau <- lapply(x, "[[", "tau")
    pnames <- factor(rep(names(x), sapply(tau, length)), levels = names(x))
    pars <- fr[cbind(seq_len(nrow(fr)),
                     match(as.character(pnames), colnames(fr)))]
    fr <- as.data.frame(fr)
    fr$.tau <- unlist(tau)
    fr$.par <- pars
    fr$.pnm <- pnames
    fr
}

## A lattice-based plot method for profile.nls objects
## FIXME use pmax.int and pmin.int after 2.5.0 is released
plot.profile.nls <-
    function (x, levels = qf(pmax(0, pmin(1, conf)), 1, df[2]),
              conf = c(50, 80, 90, 95, 99)/100, 
              absVal = TRUE, ...) 
{
    df <- attr(x, "summary")$df
    levels <- sort(levels[is.finite(levels) && levels > 0])
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

dp <- function(x = NULL,
               varname = NULL, limits, at = NULL, lab = NULL,
               draw = TRUE,

               varname.col = add.text$col,
               varname.cex = add.text$cex,
               varname.lineheight = add.text$lineheight,
               varname.font = add.text$font,
               varname.fontfamily = add.text$fontfamily,
               varname.fontface = add.text$fontface,
               
               axis.text.col = axis.text$col,
               axis.text.alpha = axis.text$alpha,
               axis.text.cex = axis.text$cex,
               axis.text.font = axis.text$font,
               axis.text.fontfamily = axis.text$fontfamily,
               axis.text.fontface = axis.text$fontface,
               
               axis.line.col = axis.line$col,
               axis.line.alpha = axis.line$alpha,
               axis.line.lty = axis.line$lty,
               axis.line.lwd = axis.line$lwd,
               ...)
{
    add.text <- trellis.par.get("add.text")
    axis.line <- trellis.par.get("axis.line")
    axis.text <- trellis.par.get("axis.text")

    if (!is.null(varname))
        grid.text(varname,
                  gp =
                  gpar(col = varname.col,
                       cex = varname.cex,
                       lineheight = varname.lineheight,
                       fontface = lattice:::chooseFace(varname.fontface,
                       varname.font),
                       fontfamily = varname.fontfamily))

    if (FALSE) ## plot axes
    ## if (draw)    
    {
        rot <- c(90, 0)
        if (is.null(at))
        {
            at <- 
                if (is.character(limits)) seq_along(limits)
                else pretty(limits)
        }
        if (is.null(lab))
        {
            lab <- 
                if (is.character(limits)) limits
                else {
                    rot <- 0
                    format(at, trim = TRUE)
                }
        }
        for (side in c("left", "top", "right", "bottom"))
            panel.axis(side = side,
                       at = at,
                       labels = lab,
                       tick = TRUE,
                       half = TRUE,

                       tck = 1, ## from scales ?
                       rot = rot, 

                       text.col = axis.text.col,
                       text.alpha = axis.text.alpha,
                       text.cex = axis.text.cex,
                       text.font = axis.text.font,
                       text.fontfamily = axis.text.fontfamily,
                       text.fontface = axis.text.fontface,

                       line.col = axis.line.col,
                       line.alpha = axis.line.alpha,
                       line.lty = axis.line.lty,
                       line.lwd = axis.line.lwd)
    }
}


splom.profile.nls <-
    function (x, levels = qf(pmax(0, pmin(1, conf)), df[1], df[2]),
              conf = c(50, 80, 90, 95, 99)/100, ...)
{
    df <- attr(x, "summary")$df
    levels <- sort(levels[is.finite(levels) && levels > 0])
    mlev <- max(levels)
    pfr <- do.call("expand.grid", lapply(x, function(el) c(-mlev, mlev)))
    spl <- lapply(x, function(x)
                  splines::interpSpline(x$par.vals[, attr(x, "parameters")$par],
                                        x$tau))
    bspl <- lapply(spl, splines::backSpline)
    fr <- as.data.frame(x)
    nms <- names(spl)
    for (nm in nms) fr[[nm]] <- predict(spl[[nm]], fr[[nm]])$y
    lp <- function(x, y, groups, subscripts, ...)
    {
        browser()
        i <- eval.parent(expression(i))
        j <- eval.parent(expression(j))
        fri <- subset(fr, .pnm == nms[i])
        sij <- interpSpline(fri[ , i], fri[ , j])
        psij <- predict(sij)
        panel.lines(psij$y, psij$x, ...)
        frj <- subset(fr, .pnm == nms[j])
        sji <- interpSpline(frj[ , j], frj[ , i])
        psji <- predict(sji)
        panel.lines(psji$x, psji$y, ...)
    }
    up <- function(...) {}
    
    splom(~ pfr, lower.panel = lp, upper.panel = up, diag.panel = dp, ...)
}

