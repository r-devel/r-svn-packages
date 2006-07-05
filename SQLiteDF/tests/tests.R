library(datasets)

# test if SQLiteDF successfully loaded. it should create (if it does not exist)
# or open workspace.db in the current directory (assuming workspace.db is 
# a valid workspace file)
library(SQLiteDF)
file.exists("workspace.db")

# test creating unnamed sdfs
u1.sdf <- sqlite.data.frame(iris)
file.exists("data1.db")
u2.sdf <- sqlite.data.frame(attenu)
file.exists("data2.db")

# test creating named sdfs
iris.sdf <- sqlite.data.frame(iris, "iris")
file.exists("iris.db")
iris2.sdf <- sqlite.data.frame(iris, "iris")
file.exists("iris1.db")

# test lsSdf(), w/c queries workspace.db
lsSdf() 
all(lsSdf() == c("data1", "data2", "iris", "iris1"))

# test attach/detach sdf
detachSdf("data1")
lsSdf() 
all(lsSdf() == c("data2", "iris", "iris1"))
attachSdf("data1.db")
lsSdf() 
all(lsSdf() == c("data2", "iris", "iris1", "data1"))

# test sdf methods
length(iris.sdf) == length(iris)
length(row.names(iris.sdf)) == length(row.names(iris))
all(dim(iris.sdf) == dim(iris))
all(names(iris.sdf) == names(iris))

# test sdf indexers
ncols <- length(iris)
nrows <- length(row.names(iris))
for (i in 1:nrows) { for (j in 1:ncols) {
    if (iris[i,j] != iris.sdf[i,j]) cat("Not equal on", i, ",", j, "\n")
}}
for (j in 1:ncols) { 
    sv <- iris.sdf[[j]]; # test sqlite.vector
    if (length(sv) != nrows) cat("Unexpected # of rows for col", j, "\n")
    for (i in 1:nrows)
        if (sv[i] != iris[i,j]) cat("Not equal on", i, "on col", j, "\n")
}
for (j in names(iris)) { 
    sv <- iris.sdf[[j]]; # test sqlite.vector
    if (length(sv) != nrows) cat("Unexpected # of rows for col", j, "\n")
    for (i in 1:nrows)
        if (sv[i] != iris[i,j]) cat("Not equal on", i, "on col", j, "\n")
}
