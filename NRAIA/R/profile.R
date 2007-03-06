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

## extract only the y component from a prediction
predy <- function(sp, vv) predict(sp, vv)$y

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
                     pval = unlist(lapply(bspl, predy, tau)),
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
           pfun <- function(x) predy(spl[[panel.number()]], x)
           panel.grid(h = -1, v = -1)
           lsegments(x, y, x, 0, ...)
           if (absVal) {
               lsegments(x, y, rev(x), y)
               pfun <- function(x) abs(predy(spl[[panel.number()]], x))
           } else {
               panel.abline(h = 0, ...)
           }
           panel.curve(pfun, ...)
       }, ...)
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
sacos <- function(x) acos(pmax(-0.999, pmin(0.999, x)))

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
    function (x, data, ## unused - only for compatibility with generic
              levels = sqrt(df[1] * qf(pmax(0, pmin(1, conf)), df[1], df[2])),
              conf = c(50, 80, 90, 95, 99)/100, ...)
{
    df <- attr(x, "summary")$df
    levels <- sort(levels[is.finite(levels) && levels > 0])
    mlev <- max(levels)
    spl <- lapply(x, function(x)
                  interpSpline(x$par.vals[, attr(x, "parameters")$par], x$tau))
    bsp <- lapply(spl, backSpline)
    pfr <- do.call("cbind", lapply(bsp, predy, c(-mlev, mlev)))
    fr <- as.data.frame(x)
    nms <- names(spl)
    for (nm in nms) fr[[nm]] <- predy(spl[[nm]], fr[[nm]])
    np <- length(nms)
    ## create a list of lists with the names of the parameters
    traces <- lapply(x, function(el) lapply(x, function(el1) list()))
    for (j in seq_along(nms)[-1]) {
        for (i in seq_len(j - 1)) {
            fri <- subset(fr, .pnm == nms[i])
            sij <- interpSpline(fri[ , i], fri[ , j])
            frj <- subset(fr, .pnm == nms[j])
            sji <- interpSpline(frj[ , j], frj[ , i])
            ll <- cont(sij, sji, levels)
            traces[[j]][[i]] <- list(sij = sij, sji = sji, ll = ll)
        }
    }
    ## panel function for lower triangle
    lp <- function(x, y, groups, subscripts, ...) {
        tr <- traces[[eval.parent(expression(j))]][[
            eval.parent(expression(i))]]
        pushViewport(viewport(xscale = c(-1.07, 1.07) * mlev,
                              yscale = c(-1.07, 1.07) * mlev))
        dd <- sapply(current.panel.limits(), diff)/50
        psij <- predict(tr$sij)
        ll <- tr$ll
        ## now do the actual plotting
        panel.grid(h = -1, v = -1)
        llines(psij$y, psij$x, ...)
        llines(predict(tr$sji), ...)
        with(ll$tki, lsegments(y - dd[1], x, y + dd[1], x, ...))
        with(ll$tkj, lsegments(x, y - dd[2], x, y + dd[2], ...))
        for (k in seq_along(levels)) llines(ll$pts[k, , ], ...)
        popViewport(1)
    }
    ## panel function for upper triangle
    up <- function(x, y, groups, subscripts, ...) {
        ## panels are transposed so reverse i and j
        i <- eval.parent(expression(j))
        j <- eval.parent(expression(i))
        tr <- traces[[j]][[i]]
        ll <- tr$ll
        pts <- ll$pts
        limits <- current.panel.limits()
        psij <- predict(tr$sij)
        psji <- predict(tr$sji)        
        ## do the actual plotting
        panel.grid(h = -1, v = -1)
        llines(predy(bsp[[i]], psij$y), predy(bsp[[j]], psij$x), ...)
        llines(predy(bsp[[i]], psji$x), predy(bsp[[j]], psji$y), ...)
        for (k in seq_along(levels))
            llines(predy(bsp[[i]], pts[k, , 1]),
                   predy(bsp[[j]], pts[k, , 2]), ...)
    }        
    dp <- function(x = NULL,            # diagonal panel
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
            at <- pretty(limits)
            sides <- c("left", "top")
            if (j == 1) sides <- "top"
            if (j == n.var) sides <- "left"
            for (side in sides)
                panel.axis(side = side,
                           at = at,
                           labels = format(at, trim = TRUE),
                           tick = TRUE,
                           check.overlap = TRUE,
                           half = side == "top" && j > 1,
                           
                           tck = 1, rot = 0, 
                           
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
            lims <- c(-1.07, 1.07) * mlev
            pushViewport(viewport(xscale = lims, yscale = lims))
            side <- ifelse(j == 1, "right", "bottom")
            which.half <- ifelse(j == 1, "lower", "upper")
            at <- pretty(lims)
            panel.axis(side = side, at = at, labels = format(at, trim = TRUE),
                       tick = TRUE, half = TRUE, which.half = which.half,
                       tck = 1, rot = 0,

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
            popViewport(1)
        }
    }

    splom(~ pfr, lower.panel = lp, upper.panel = up, diag.panel = dp, ...)
}

