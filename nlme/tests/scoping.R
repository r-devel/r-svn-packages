library("nlme")

## PR#18157
mygnls <- function (mydata)
    gnls(weight ~ SSlogis(Time, Asym, xmid, scal), data = mydata)
fm1 <- mygnls(Soybean) # failed in 3.1-153 with
## Error in stats::nls(formula = weight ~ SSlogis(Time, Asym, xmid, scal),  : 
##   object 'mydata' not found
