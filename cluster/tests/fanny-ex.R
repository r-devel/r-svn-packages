
library(cluster)

###--- An extension of  example(fanny) : -------------------

## generate 25 objects, divided into two clusters, and 3 objects lying
## between those clusters.
set.seed(171)
x <- rbind(cbind(rnorm(10, 0, 0.5), rnorm(10, 0, 0.5)),
           cbind(rnorm(15, 5, 0.5), rnorm(15, 5, 0.5)),
           cbind(rnorm( 3,3.5,0.5), rnorm( 3,3.5,0.5)))

data(ruspini)

.proctime00 <- proc.time()

(fannyx <- fanny(x, 2))
summary(fannyx)
str(fannyx)
summary(fanny(x,3))# one extra

summary(fanny(ruspini, 3), digits = 9)
summary(fanny(ruspini, 4), digits = 9)# `correct' #{clusters}
summary(fanny(ruspini, 5), digits = 9)

## Last Line:
cat('Time elapsed: ', proc.time() - .proctime00,'\n')
