### Clusplot() RESULT checking ...


## plotting votes.diss(dissimilarity) in a bivariate plot and
## partitioning into 2 clusters
data(votes.repub)
votes.diss <- daisy(votes.repub)
votes.clus <- pam(votes.diss, 2, diss = TRUE)$clustering
print.default(clusplot(votes.diss, votes.clus, diss = TRUE, shade = TRUE))

## plotting iris (dataframe) in a 2-dimensional plot and partitioning
## into 3 clusters.
data(iris)
iris.x <- iris[, 1:4]
print.default(clusplot(iris.x, pam(iris.x, 3)$clustering, diss = FALSE))

.Random.seed <- c(0,rep(7654,3))

## generate 25 objects, divided into 2 clusters.
x <- rbind(cbind(rnorm(10,0,0.5), rnorm(10,0,0.5)),
           cbind(rnorm(15,5,0.5), rnorm(15,5,0.5)))
print.default(clusplot(pam(x, 2)))


