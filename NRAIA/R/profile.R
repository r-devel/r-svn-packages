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
    function (x, levels = sqrt(qf(pmax(0, pmin(1, conf)), 1, df[2])),
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
    fr <- data.frame(tau = rep.int(tau, length(x)),
                     pval = unlist(lapply(bspl,
                     function(sp) predict(sp, x= tau)$y)),
                     pnm = gl(length(x), length(tau), labels = names(x)))
    ylab <- expression(tau)
    if (absVal) {
        fr$tau <- abs(fr$tau)
        ylab <- expression("|" * tau * "|")
    }
    xyplot(tau ~ pval | pnm, fr,
           scales = list(x = list(relation = 'free')),
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
    j <- eval.parent(expression(j))
    n.var <- eval.parent(expression(n.var))
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

    if (draw)    
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
        sides <- c("right", "bottom")
        if (j == 1) sides <- "right"
        if (j == n.var) sides <- "bottom"
        for (side in sides)
            panel.axis(side = side,
                       at = at,
                       labels = lab,
                       tick = TRUE,
                       check.overlap = TRUE,
                       half = !(j %in% c(1, n.var)),

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

## convert the x-cosine and y-cosine to an average and difference,
## ensuring that the difference is positive by flipping signs if
## necessary
ad <- function(xc, yc)
{
    a <- (xc + yc)/2
    d <- (xc - yc)
    cbind(ifelse(d > 0, a, -a), abs(d))
}

## convert d versus a (as an xyVector) and level to a matrix of taui and tauj
tauij <- function(xy, lev) lev * cos(xy$x + outer(xy$y/2, c(-1, 1)))

## safe arc-cosine
sacos <- function(x) acos(pmax(-1, pmin(1, x)))

## extract only the y component from a prediction
predy <- function(sp, vv) predict(sp, vv)$y

cont <- function(sij, sji, levels, nseg = 101)
{
    ada <- array(0, c(length(levels), 2, 4))
    ada[, , 1] <- ad(0, sacos(predy(sij,  levels)/levels))
    ada[, , 2] <- ad(sacos(predy(sji, levels)/levels), 0)
    ada[, , 3] <- ad(pi, sacos(predy(sij, -levels)/levels))
    ada[, , 4] <- ad(sacos(predy(sji, -levels)/levels), pi)
    pts <- array(0, c(length(levels), nseg + 1, 2))
    for (i in seq_along(levels))
        pts[i, ,] <- tauij(predict(periodicSpline(ada[i, 1, ], ada[i, 2, ]),
                                   nseg = nseg), levels[i])
    levs <- c(-rev(levels), 0, levels)
    list(tki = predict(sij, levs), tkj = predict(sji, levs), pts = pts)
}

splom.profile.nls <-
    function (x, data,  ## unused - for compatibility with generic only
              levels = sqrt(df[1] * qf(pmax(0, pmin(1, conf)), df[1], df[2])),
              conf = c(50, 80, 90, 95, 99)/100, ...)
{
    df <- attr(x, "summary")$df
    levels <- sort(levels[is.finite(levels) && levels > 0])
    mlev <- max(levels)
    pfr <- do.call("expand.grid", lapply(x, function(el) c(-mlev, mlev)))
    spl <- lapply(x, function(x)
                  interpSpline(x$par.vals[, attr(x, "parameters")$par], x$tau))
    bsp <- lapply(spl, backSpline)
    fr <- as.data.frame(x)
    nms <- names(spl)
    for (nm in nms) fr[[nm]] <- predy(spl[[nm]], fr[[nm]])
    lp <- function(x, y, groups, subscripts, ...) {
        i <- eval.parent(expression(i))
        j <- eval.parent(expression(j))
        fri <- subset(fr, .pnm == nms[i])
        sij <- interpSpline(fri[ , i], fri[ , j])
        frj <- subset(fr, .pnm == nms[j])
        sji <- interpSpline(frj[ , j], frj[ , i])
        psij <- predict(sij)
        ll <- cont(sij, sji, levels)
        dd <- sapply(current.panel.limits(), diff)/50
        ## now do the actual plotting
        panel.grid(h = -1, v = -1)
        llines(psij$y, psij$x, ...)
        llines(predict(sji), ...)
        with(ll$tki, lsegments(y - dd[1], x, y + dd[1], x, ...))
        with(ll$tkj, lsegments(x, y - dd[2], x, y + dd[2], ...))
        for (k in seq_along(levels)) llines(ll$pts[k, , ], ...)
    }
    up <- function(x, y, groups, subscripts, ...) {
        i <- eval.parent(expression(j)) ## plots are transposed
        j <- eval.parent(expression(i))
        fri <- subset(fr, .pnm == nms[i])
        sij <- interpSpline(fri[ , i], fri[ , j])
        frj <- subset(fr, .pnm == nms[j])
        sji <- interpSpline(frj[ , j], frj[ , i])
        psij <- predict(sij)
        psji <- predict(sji)
        ll <- cont(sij, sji, levels)
        pts <- ll$pts
        ## do the actual plotting
        pushViewport(viewport(xscale = range(predy(bsp[[i]], x)),
                              yscale = range(predy(bsp[[j]], y))))
        limits <- current.panel.limits()
        panel.grid(h = -1, v = -1)
        llines(predy(bsp[[i]], psij$y), predy(bsp[[j]], psij$x), ...)
        llines(predy(bsp[[i]], psji$x), predy(bsp[[j]], psji$y), ...)
        for (k in seq_along(levels))
            llines(predy(bsp[[i]], pts[k, , 1]),
                   predy(bsp[[j]], pts[k, , 2]), ...)
        popViewport(1)
    }        

    splom(~ pfr, lower.panel = lp, upper.panel = up, diag.panel = dp, ...)
}

