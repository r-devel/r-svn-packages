library(cluster)

data(animals)
(mani <- mona(animals, trace=TRUE))

str(mani)

set.seed(1)
ani.N1 <- animals; ani.N1[cbind(sample.int(20, 10), sample.int(6, 10, replace=TRUE))] <- NA
(maniN <- mona(ani.N1, trace=TRUE))

for(seed in c(2:20)) {
    set.seed(seed); cat("seed = ", seed,"\n")
    ani.N2 <- animals
    ani.N2[cbind(sample.int(20, 9),
                 sample.int( 6, 9, replace=TRUE))] <- NA
    try(print(maniN2 <- mona(ani.N2, trace=TRUE)))
}

## Check all "too many NA" and other illegal cases
ani.NA   <- animals; ani.NA[4,] <- NA
aniNA    <- within(animals, { end[2:9] <- NA })
aniN2    <- animals; aniN2[cbind(1:6, c(3, 1, 4:6, 2))] <- NA
ani.non2 <- within(animals, end[7] <- 3 )
ani.idNA <- within(animals, end[!is.na(end)] <- 1 )
## use tools::assertError() {once you don't use *.Rout.save anymore}
try( mona(ani.NA)   )
try( mona(aniNA)    )
try( mona(aniN2)    )
try( mona(ani.non2) )
try( mona(ani.idNA) )

if(require(MASS)) {

    if(R.version$major != "1" || as.numeric(R.version$minor) >= 7)
        RNGversion("1.6")
    set.seed(253)
    n <- 512; p <- 3
    Sig <- diag(p); Sig[] <- 0.8 ^ abs(col(Sig) - row(Sig))
    x3 <- mvrnorm(n, rep(0,p), Sig) >= 0
    x <- cbind(x3, rbinom(n, size=1, prob = 1/2))

    print(sapply(as.data.frame(x), table))

    mx <- mona(x)
    str(mx)
    print(lapply(mx[c(1,3,4)], table))
}

if(FALSE)
    mona(cbind(1:0), trace=2)
## Loop npass = 1: (ka,kb) = (1,13)
##   for(j ..) -> jma=1, jtel(.,z) = (8, -1777201152) --> splitting: (nel; jres, ka, km) = (1; -6, 1, 9)
##  1  2  2  2  2 ..... [infinite loop]
