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

# test creating named sdfs
iris.sdf <- sqlite.data.frame(iris, "iris")
stopifnot(file.exists("iris.db"))
iris2.sdf <- sqlite.data.frame(iris, "iris")
stopifnot(file.exists("iris1.db"))

# test lsSdf(), w/c queries workspace.db
lsSdf() 
stopifnot((lsSdf() == c("data1", "data2", "iris", "iris1")))

# test attach/detach sdf
detachSdf("data1")
stopifnot(lsSdf() == c("data2", "iris", "iris1"))
attachSdf("data1.db")
stopifnot(lsSdf() == c("data2", "iris", "iris1", "data1"))

# test sdf methods
stopifnot(length(iris.sdf) == length(iris),
          length(row.names(iris.sdf)) == length(row.names(iris)),
          all(dim(iris.sdf) == dim(iris)),
          all(names(iris.sdf) == names(iris)))

# test sdf indexers
ncols <- length(iris)
nrows <- length(row.names(iris))
for (i in 1:nrows) { for (j in 1:ncols) {
    if (iris[i,j] != iris.sdf[i,j]) stop("Not equal on ", i, ",", j, "\n")
}}
for (j in 1:ncols) { 
    sv <- iris.sdf[[j]]; # test sqlite.vector
    if (length(sv) != nrows) stop("Unexpected # of rows for col", j, "\n")
    for (i in 1:nrows)
        if (sv[i] != iris[i,j]) stop("Not equal on", i, "on col", j, "\n")
}
for (j in names(iris)) { 
    sv <- iris.sdf[[j]]; # test sqlite.vector
    if (length(sv) != nrows) stop("Unexpected # of rows for col", j, "\n")
    for (i in 1:nrows)
        if (sv[i] != iris[i,j]) stop("Not equal on", i, "on col", j, "\n")
}
#iris.sdf.project <- iris[,1:4]
