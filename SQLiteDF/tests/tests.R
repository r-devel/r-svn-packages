library(datasets)

# test if SQLiteDF successfully loaded. it should create (if it does not exist)
# or open workspace.db in the current directory (assuming workspace.db is 
# a valid workspace file)
library(SQLiteDF)
stopifnot(file.exists("workspace.db"))

# test creating unnamed sdfs
u1.sdf <- sqlite.data.frame(iris)
stopifnot(file.exists("data1.db"))
u2.sdf <- sqlite.data.frame(attenu)
stopifnot(file.exists("data2.db"))

# test classes
stopifnot(class(u1.sdf) == "sqlite.data.frame",
          class(u2.sdf) == "sqlite.data.frame")

# test row.names, head, tail
stopifnot(all(as.character(row.names(u1.sdf))==row.names(iris)),
          all(as.character(row.names(u2.sdf))==row.names(attenu)),
          all(head(u1.sdf) == head(iris)),
          all(head(u2.sdf) == head(attenu)),
          all(tail(u1.sdf) == tail(iris)),
          all(tail(u2.sdf) == tail(attenu)))


compareSdfToDf <- function(sdf, df, with.names=TRUE) {
    ncols <- ncol(df)
    nrows <- nrow(df)
    for (i in 1:nrows) { for (j in 1:ncols) {
        if (df[i,j] != sdf[i,j]) stop("Not equal on ", i, ",", j, "\n")
    }}
    for (j in 1:ncols) { 
        sv <- sdf[[j]]; # test sqlite.vector
        stopifnot(class(sv) == "sqlite.vector", has.typeSvec(sv, class(df[[j]])[1]))
        if (length(sv) != nrows) stop("Unexpected # of rows for col", j, "\n")
        for (i in 1:nrows)
            if (sv[i] != df[i,j]) stop("Not equal on", i, "on col", j, "\n")
    }
    if (with.names) for (j in names(df)) { 
        sv <- sdf[[j]]; # test sqlite.vector
        stopifnot(class(sv) == "sqlite.vector", has.typeSvec(sv, class(df[[j]])[1]))
        if (length(sv) != nrows) stop("Unexpected # of rows for col", j, "\n")
        for (i in 1:nrows)
            if (sv[i] != df[i,j]) stop("Not equal on", i, "on col", j, "\n")
    }
}

# test creating named sdfs
iris.sdf <- sqlite.data.frame(iris, "iris")
stopifnot(file.exists("iris.db"))
iris2.sdf <- sqlite.data.frame(iris, "iris")
stopifnot(file.exists("iris1.db"))

# test lsSdf(), w/c queries workspace.db
lsSdf() 
stopifnot((sort(lsSdf()) == c("data1", "data2", "iris", "iris1")))

# test attach/detach sdf
detachSdf("data1")
stopifnot(sort(lsSdf()) == c("data2", "iris", "iris1"))
attachSdf("data1.db")
stopifnot(sort(lsSdf()) == sort(c("data2", "iris", "iris1", "data1")))

# test sdf methods
stopifnot(length(iris.sdf) == length(iris),
          length(row.names(iris.sdf)) == length(row.names(iris)),
          all(dim(iris.sdf) == dim(iris)),
          all(names(iris.sdf) == names(iris)))

# test sdf indexers
compareSdfToDf(iris.sdf, iris)

# test operators
stopifnot(all((iris.sdf[,1] + iris.sdf[,2]) == (iris[,1] + iris[,2])))
stopifnot(all((iris.sdf[,1]*10 %/% iris.sdf[,2]) == (iris[,1]*10 %/% iris[,2])))
stopifnot(all(signif(log(iris.sdf[,1],17)) == signif(log(iris[,1],17))))
stopifnot(all(round(log(iris.sdf[,1]),5) == round(log(iris[,1]),5)))
stopifnot(all(with(iris.sdf, Sepal.Length*Sepal.Width - Petal.Length/Petal.Width) == with(iris, Sepal.Length*Sepal.Width - Petal.Length/Petal.Width)))

stopifnot(sapply(iris.sdf[,1:4],sum) == sapply(iris[,1:4],sum))

if (require(RSQLite)) {
    dr <- SQLite()
    con <- dbConnect(dr, dbname="example.db")
    i1 <- sdfImportDBI(con, "select * from iris")
    compareSdfToDf(i1[,1:4], iris[,1:4], with.names=FALSE)
    stopifnot(all(as.character(i1[,5]) == as.character(iris[,5])))

    i2 <- sdfImportDBI(con, "select * from iris", 30)  # test rbindSdf
    compareSdfToDf(i2[,1:4], iris[,1:4], with.names=FALSE)
    stopifnot(all(as.character(i2[,5]) == as.character(iris[,5])))
}

# test summary
#for (j in 1:5) stopifnot(all(summary(iris.sdf[,j]) == summary(iris[,j])))
#iris.sdf.summary <- summary(iris.sdf)
#iris.summary <- summary(iris)
#stopifnot(all(iris.sdf.summary[,1:4] == iris.summary[,1:4]))
#stopifnot(all(na.exclude(iris.sdf.summary[,5])==na.exclude(iris.summary[,5])))
