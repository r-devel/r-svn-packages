## Originally inspired by  Kasper Fischer-Rasmussen 's  clara_gower.html  [html from Rmd]
'
  Gower’s distance with CLARA
'
library(cluster)
packageDescription("cluster")

## carefully getting  150 + 200 + 150 = 500 obs. from the 3  xclara clusters :
str(dd <- xclara[c(1:150, 1001:1200, 2101:2250), ])
dim(dd) # 500 2

set.seed(47)
cl_manhat <- clara(dd, 3, metric = "manhattan", rngR=TRUE, pamLike=TRUE, samples = 500)
cl_gower  <- clara(dd, 3, metric = "gower",     rngR=TRUE, pamLike=TRUE, samples = 500)

stopifnot( cl_manhat$cluster == cl_gower $cluster )## they *are* the same

## There are no distinction between the clusters depending on the usage of
## Manhattan vs. Gower’s distance. Thus the changes survived hte initial
## testing.

"
 Using “cluster”’s built in tools to compute Gower’s distance.
"

cl_gower_full <- clara(dd, k = 3, metric = "gower", rngR = TRUE, pamLike = TRUE, samples = 500, sampsize = nrow(dd))
dist_cl_full <- as.matrix(cl_gower_full$diss)
i_full <- rownames(dist_cl_full)
d_full <- data.frame(CLARA = as.vector(cl_gower_full$diss),
                     DAISY = as.vector(daisy(dd[i_full, ], metric = "gower")))

## MM: instead of all this, just
all.equal(d_full$CLARA,
          d_full$DAISY, tol=0) # ] "Mean relative difference: 2.17e-16"

stopifnot( all.equal(d_full$CLARA,
                     d_full$DAISY, tol = 1e-15) ) ## equal up to  15 digits!

## We can see that the distance measurements are exactly identical when the
## whole data is used in the clustering. This is because the Gower distance
## scales the distances measurements with the range of each feature. Due to
## the subsampling, approximate ranges are calculated based on each
## subsample explaining the deviations.


## MM: compare -- with pam():
dGow <- daisy(dd, metric="gower")
cl_full <- clara(dd, k = 3, metric = "gower", rngR = TRUE, pamLike = TRUE, samples = 1, sampsize = nrow(dd))

all.equal(c(dGow) , c(cl_full$diss), tol=0) # "Mean relative difference: 2.171402e-16"

pam_3 <- pam(dGow, k = 3, variant = "faster")
stopifnot(exprs = {
    identical(pam_3  $ clustering,
              cl_full$ clustering)
    all.equal(c(dGow) , c(cl_full$diss), tol = 1e-15)
})

