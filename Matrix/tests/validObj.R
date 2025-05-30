library(Matrix)

### Do all kinds of object creation and coercion

source(system.file("test-tools.R", package = "Matrix"))

## the empty ones:
chk.matrix(new("dgeMatrix"))

## "dge"
assertError( new("dgeMatrix", Dim = c(2,2), x= 1:4) )# double 'Dim'
assertError( new("dgeMatrix", Dim = as.integer(c(2,2)), x= 1:4) )# int 'x'
assertError( new("dgeMatrix", Dim = 2:2, x=as.double(1:4)) )# length(Dim) !=2
assertError( new("dgeMatrix", Dim = as.integer(c(2,2)), x= as.double(1:5)))

chk.matrix(m1 <- Matrix(1:6, ncol=2))
chk.matrix(m2 <- Matrix(1:7 +0, ncol=3)) # a (desired) warning
stopifnot(unique(is(m1)) == c("dgeMatrix", "ddenseMatrix", "generalMatrix",
	    "dMatrix", "denseMatrix", "Matrix", "compMatrix"),
	  dim(t(m1)) == 2:3,
	  identical(m1, t(t(m1))))
c.nam <- paste("C",1:2, sep='')
dimnames(m1) <- list(NULL, c.nam)
stopifnot(colnames(m1) == c.nam,
	  identical(dimnames(t(m1)), list(c.nam, NULL)),
	  identical(m1, t(t(m1))))

## an example of *named* dimnames
(t34N <- as(unclass(table(x = gl(3,4), y=gl(4,3))), "dgeMatrix"))
stopifnot(identical(dimnames(t34N),
		    dimnames(as(t34N, "matrix"))),
          identical(t34N, t(t(t34N))))

## "dpo"
chk.matrix(cm <- crossprod(m1))
chk.matrix(cp <- as(cm, "dppMatrix"))# 'dpp' + factors
chk.matrix(cs <- as(cm, "dsyMatrix"))# 'dsy' + factors
chk.matrix(dcm <- as(cm, "dgeMatrix"))#'dge'
chk.matrix(mcm <- as(cm, "dMatrix")) # 'dsy' + factors -- buglet? rather == cm?
chk.matrix(mc. <- as(cm, "Matrix"))
stopifnot(identical(mc., mcm),
          identical4(2*cm, cm + cp, cp + cs, mcm * 2))

chk.matrix(eq <- cm == cs)
stopifnot(all(eq@x),
	  identical3(eq, cs == cp, cm == cp),
	  as.logical(!(cs < cp)),
	  identical4(!(cs < cp), !(cp > cs), cp <= cs, cs >= cp))

## Coercion to 'dpo' should give an error if result would be invalid
M <- Matrix(diag(4) - 1)
if(FALSE)## FIXME?: dsy -> dpo works here{chol() check too expensive ?}
assertError(as(M, "dpoMatrix"))
M. <- as(M, "dgeMatrix")
M.[1,2] <- 10 # -> not even symmetric anymore
assertError(as(M., "dpoMatrix"))


## Cholesky
chk.matrix(ch <- chol(cm))
chk.matrix(ch2 <- chol(as(cm, "dsyMatrix")))
chk.matrix(ch3 <- chol(as(cm, "dgeMatrix")))
stopifnot(is.all.equal3(as(ch, "matrix"), as(ch2, "matrix"), as(ch3, "matrix")))
### Very basic	triangular matrix stuff

assertError( new("dtrMatrix", Dim = c(2,2), x= 1:4) )# double 'Dim'
assertError( new("dtrMatrix", Dim = as.integer(c(2,2)), x= 1:4) )# int 'x'
## This caused a segfault (before revision r1172 in ../src/dtrMatrix.c):
assertError( new("dtrMatrix", Dim = 2:2, x=as.double(1:4)) )# length(Dim) !=2
assertError( new("dtrMatrix", Dim = as.integer(c(2,2)), x= as.double(1:5)))

tr22 <- new("dtrMatrix", Dim = as.integer(c(2,2)), x=as.double(1:4))
tt22 <- t(tr22)
(tPt <- tr22 + tt22)
stopifnot(identical(10 * tPt, tPt * 10),
	  (t.22 <- (tr22 / .5)* .5)@x == c(1,0,3,4),
	  TRUE) ## not yet: class(t.22) == "dtrMatrix")

## non-square triagonal Matrices --- are forbidden ---
assertError(new("dtrMatrix", Dim = 2:3,
		x=as.double(1:6), uplo="L", diag="U"))

n <- 3:3
assertError(new("dtCMatrix", Dim = c(n,n), diag = "U"))
stopifnot(validObject(T <- new("dtTMatrix", Dim = c(n,n), diag = "U")),
	  identical(as.mat(T), diag(n)),
	  validObject(M <- new("dtCMatrix", Dim = c(n,n), diag = "U",
			       p = rep.int(0:0, n+1)))
	  )

set.seed(3) ; (p9 <- as(sample(9), "pMatrix"))
## Check that the correct error message is triggered
ind.try <- try(p9[1,1] <- 1, silent = TRUE)
stopifnot(grep("replacing.*sensible", ind.try[1]) == 1,
	  is.logical(p9[1,]),
	  isTRUE(p9[-c(1:6, 8:9), 1]),
	  identical(t(p9), solve(p9)),
	  identical(p9[TRUE,], as(p9, "ngTMatrix")),
	  identical(as(diag(9), "pMatrix"), as(1:9, "pMatrix"))
	  )
