#### Tests for gam() and its methods, with  NA  data,
#### notably with  na.action = na.exclude {which some say is *the* reasonable na.action}
library(mgcv)

set.seed(11)
n <- 500
## NA's and outlier
iOut <- 200:205
## non-outlier indices -- add NAs only to those
iIn <- which(is.na(match(1:n, iOut)))
d.4 <- within(gamSim(1, n = n), {
    y[iOut] <- y[iOut] + 20
    y[c(2, sample(iIn, n/20))] <- NA
})
for(xn in paste("x",0:3, sep=""))
    d.4[sample(iIn, 5), xn] <- NA
head(d.4, 16)

m1 <- gam(y ~ s(x0) + s(x1) + s(x2) + s(x3), data = d.4,
          na.action = "na.exclude")
## Large residuals at exactly the "known" locations [because of "na.exclude"]:
stopifnot(identical(iOut, which(unname(residuals(m1)) > 10)))
gam.check(m1) ## gave error in earlier versions
