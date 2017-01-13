library(cluster)

data(animals)
(mani <- mona(animals))

str(mani)

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
