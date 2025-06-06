library(nlme)

## from example(nlme.nlsList)
fm1 <- nlsList(SSasymp, Loblolly)
fm1

stopifnot(
    all.equal(as.matrix(coef(fm1)),
              array(c(94.1282, 94.9406, 89.8849, 110.699, 111.003, 109.986, 101.056,
                      127.134, 101.087, 95.6669, 95.5563, 113.514, 105.718, 99.1719,
                      -8.25075, -7.75749, -8.75902, -8.16943, -8.46261, -8.55854, -8.44363,
                      -7.67936, -8.50234, -9.07824, -9.66503, -7.59562, -8.90608, -9.91665,
                      -3.21758, -3.22932, -3.08622, -3.39034, -3.39757, -3.36252, -3.23282,
                      -3.57533, -3.21402, -3.11638, -3.09227, -3.35282, -3.22296, -3.08484),
                    dim = c(14L, 3L),
                    dimnames = list(c("329", "327", "325", "307", "331", "311", "315",
                                      "321", "319", "301", "323", "309", "303", "305"),
                                    c("Asym", "R0", "lrc"))),
              tol = 1e-5)
    ,
    all.equal(pooledSD(fm1), structure(0.70039649, df = 42), tol = 1e-5)
    ,
    84 == sum(lengths(lapply(fm1, fitted))) # total deg.freedom
)

## confint():
if (getRversion() >= "4.4.0" || # confint.nls() taken from MASS
    (requireNamespace("MASS") && packageVersion("MASS") < "7.3-60.1")) {
    system.time(cnL1 <- confint(fm1)) # 0.4 sec
    stopifnot(exprs = {
        is.list(cnL1)
        identical(names(cnL1), names(fm1))
        vapply(cnL1, is.matrix, NA)
        vapply(cnL1, dim, c(NA,0L)) == 3:2
        identical(unname(sapply(cnL1, dim)), matrix(3:2, 2, length(fm1)))
        sapply(cnL1, is.finite)
    })
}

## build a random-effects model from the stratified nls fits
fm2 <- nlme(fm1, random = Asym ~ 1)
fm2
stopifnot(
    all.equal(fixef(fm2),
              c(Asym = 101.4483, R0 = -8.6274937, lrc = -3.2337304), tol = 4e-7)
    ,
    all.equal(logLik(fm2),
              structure(-114.743, class = "logLik", nall = 84, nobs = 84, df = 5),
              tol = 4e-6)
    ,
    all.equal(sigma(fm2), 0.71886244, tol = 1e-6)
    )

pm2.0 <- predict(fm2, Loblolly, level=0)## failed in nlme 3.1-123
stopifnot(all.equal(head(pm2.0),
		    c(3.64694, 11.0597, 27.2258, 40.5006, 51.4012, 60.3522),
		    tol = 1e-5)) # 4e-7 {64b nb-mm4}

## same as fm2 but fitted via nlme's formula interface
fm3 <- nlme(height ~ SSasymp(age, Asym, R0, lrc),
            data = Loblolly,
            fixed = Asym + R0 + lrc ~ 1,
            random = Asym ~ 1, groups = ~Seed)
## explicitly specifying 'groups' failed in 3.1-153 with
## Error in nlme::nlsList(model = height ~ SSasymp(age, Asym, R0, lrc), data = Loblolly,  : 
##   unused argument (groups = ~Seed)
fm2$origCall <- NULL
stopifnot(all.equal(fm2, fm3))


## reproduce a simple random-intercept lme() using nlme()
BW1 <- subset(BodyWeight, Diet == 1, -Diet)
fm1BW.lme <- lme(weight ~ 1, data = BW1, random = ~ 1 | Rat, method = "ML")
fm1BW.nlme <- nlme(weight ~ b0, fixed = b0 ~ 1, random = b0 ~ 1 | Rat,
                   data = BW1, start = fixef(fm1BW.lme), method = "ML")
stopifnot(exprs = {
    all.equal(logLik(fm1BW.lme), logLik(fm1BW.nlme))
    all.equal(fixef(fm1BW.lme), fixef(fm1BW.nlme), check.names = FALSE)
    all.equal(ranef(fm1BW.lme), ranef(fm1BW.nlme), check.attributes = FALSE)
    all.equal(as.numeric(VarCorr(fm1BW.lme)), as.numeric(VarCorr(fm1BW.nlme)))
})

## now the same with an additional spatial correlation structure (PR#18192)
fm2BW.lme <- update(fm1BW.lme, correlation = corExp(form = ~ Time))
fm2BW.nlme <- try(update(fm1BW.nlme, correlation = corExp(form = ~ Time)))
## nlme() <= 3.1-166 would sometimes crash with varying memory errors, e.g.:
## segfault, possibly also reporting many "*** recursive gc invocation",
## or crash with "corrupted size vs. prev_size" or "double free or corruption (!prev)"
if (!inherits(fm2BW.nlme, "try-error")) # may not have converged (seen on Apple M1)
stopifnot(exprs = {
    all.equal(logLik(fm2BW.lme), logLik(fm2BW.nlme))
    all.equal(fixef(fm2BW.lme), fixef(fm2BW.nlme), check.names = FALSE)
    all.equal(ranef(fm2BW.lme), ranef(fm2BW.nlme), check.attributes = FALSE,
              tolerance = 1e-5) # seen 4e-05
    all.equal(as.numeric(VarCorr(fm1BW.lme)), as.numeric(VarCorr(fm1BW.nlme)))
})
