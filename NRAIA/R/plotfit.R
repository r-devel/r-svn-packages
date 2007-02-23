plotfit <- function(fm, ...) UseMethod("plotfit")

plotfit.nls <- function(fm, ...)
{
    mm <- fm$m
    cc <- fm$call
    pnms <- names(mm$getPars())
    form <- cc$formula
    rhsnms <- all.vars(form[[3]])
    vnms <- rhsnms[!(rhsnms %in% pnms)]
    if (length(vnms) > 1)
        stop("plotfit not yet implements for >1 covariate")
    predfun <- function(x) {
        ll <- list(x)
        names(ll) <- vnms
        predict(fm, ll)
    }
    xyplot(eval(substitute(y ~ x, list(y = form[[2]],
                                       x = as.name(vnms)))),
           mm$getEnv(),
           panel = function(x, y, ...) {
               panel.grid(h = -1, v = -1)
               panel.points(x, y, ...)
               panel.curve(predfun, ...)
           }, ...)
}

    
