#include "cs_utils.h"

/* Borrowed from one of Tim Davis' examples in the CSparse Demo directory */
/* 1 if A is square & upper tri., -1 if square & lower tri., 0 otherwise */
static int is_sym (cs *A)
{
    int is_upper, is_lower, j, p, n = A->n, m = A->m, *Ap = A->p, *Ai = A->i ;
    if (m != n) return (0) ;
    is_upper = 1 ;
    is_lower = 1 ;
    for (j = 0 ; j < n ; j++)
    {
	for (p = Ap [j] ; p < Ap [j+1] ; p++)
	{
	    if (Ai [p] > j) is_upper = 0 ;
	    if (Ai [p] < j) is_lower = 0 ;
	}
    }
    return (is_upper ? 1 : (is_lower ? -1 : 0)) ;
}

/**
 * Create a cs object with the contents of x.  Typically called via  AS_CSP()
 *
 * @param ans pointer to a cs struct.  This is allocated in the caller
 * so it is easier to keep track of where it should be freed - in many
 * applications the memory can be allocated with alloca and
 * automatically freed on exit from the caller.
 * @param x pointer to an object that inherits from CsparseMatrix
 *
 * @return pointer to a cs object that contains pointers
 * to the slots of x.
 */
cs *Matrix_as_cs(cs *ans, SEXP x)
{
    char *valid[] = {"dgCMatrix", "dtCMatrix", ""};/* had also "dsCMatrix", but that
						    * only stores one triangle */
    int *dims, ctype = Matrix_check_class(class_P(x), valid);
    SEXP islot;

    if (ctype < 0) error("invalid class of object to Matrix_as_cs");
				/* dimensions and nzmax */
    dims = INTEGER(GET_SLOT(x, Matrix_DimSym));
    ans->m = dims[0]; ans->n = dims[1];
    islot = GET_SLOT(x, Matrix_iSym);
    ans->nz = -1;		/* indicates compressed column storage */
    ans->nzmax = LENGTH(islot);
    ans->i = INTEGER(islot);
    ans->p = INTEGER(GET_SLOT(x, Matrix_pSym));
    ans->x = REAL(GET_SLOT(x, Matrix_xSym));

    return ans;
}

/**
 * Copy the contents of a to an appropriate CsparseMatrix object and,
 * optionally, free a or free both a and the pointers to its contents.
 *
 * @param a matrix to be converted
 * @param cl the name of the S4 class of the object to be generated
 * @param dofree 0 - don't free a; > 0 cs_free a; < 0 Free a
 *
 * @return SEXP containing a copy of a
 */
SEXP Matrix_cs_to_SEXP(cs *a, char *cl, int dofree)
{
    SEXP ans;
    char *valid[] = {"dgCMatrix", "dsCMatrix", "dtCMatrix", ""};
    int *dims, ctype = Matrix_check_class(cl, valid), nz;

    if (ctype < 0)
	error(_("invalid class of object to Matrix_cs_to_SEXP"));
    ans = PROTECT(NEW_OBJECT(MAKE_CLASS(cl)));
				/* allocate and copy common slots */
    dims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2));
    dims[0] = a->m; dims[1] = a->n;
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_pSym, INTSXP, a->n + 1)),
	   a->p, a->n + 1);
    nz = a->p[a->n];
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_iSym, INTSXP, nz)), a->i, nz);
    Memcpy(REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, nz)), a->x, nz);
    if (ctype > 0) { /* dsC or dtC */
	int uplo = is_sym(a);
	if (!uplo)
	    error(_("cs matrix not compatible with class '%s'"), valid[ctype]);
	if (ctype == 2) /* dtC* */
	    SET_SLOT(ans, Matrix_diagSym, mkString("N"));
	SET_SLOT(ans, Matrix_uploSym, mkString(uplo < 0 ? "L" : "U"));
    }
    if (dofree > 0) cs_spfree(a);
    if (dofree < 0) Free(a);
    UNPROTECT(1);
    return ans;
}

#if 0 				/* unused ------------------------------------*/
/**
 * Populate a css object with the contents of x.
 *
 * @param ans pointer to a csn object
 * @param x pointer to an object of class css_LU or css_QR.
 *
 * @return pointer to a cs object that contains pointers
 * to the slots of x.
 */
css *Matrix_as_css(css *ans, SEXP x)
{
    char *cl = class_P(x);
	*valid[] = {"css_LU", "css_QR", ""};
    int *nz = INTEGER(GET_SLOT(x, install("nz"))),
	ctype = Matrix_check_class(cl, valid);

    if (ctype < 0) error("invalid class of object to Matrix_as_css");
    ans->q = INTEGER(GET_SLOT(x, install("Q")));
    ans->m2 = nz[0]; ans->lnz = nz[1]; ans->unz = nz[2];
    switch(ctype) {
    case 0:			/* css_LU */
	ans->pinv = (int *) NULL;
	ans->parent = (int *) NULL;
	ans->cp = (int *) NULL;
	break;
    case 1:			/* css_QR */
	ans->pinv = INTEGER(GET_SLOT(x, install("Pinv")));
	ans->parent = INTEGER(GET_SLOT(x, install("parent")));
	ans->cp = INTEGER(GET_SLOT(x, install("cp")));
	break;
    default:
	error("invalid class of object to Matrix_as_css");
    }
    return ans;
}

/**
 * Populate a csn object with the contents of x.
 *
 * @param ans pointer to a csn object
 * @param x pointer to an object of class csn_LU or csn_QR.
 *
 * @return pointer to a cs object that contains pointers
 * to the slots of x.
 */
csn *Matrix_as_csn(csn *ans, SEXP x)
{
    char *valid[] = {"csn_LU", "csn_QR", ""};
    int ctype = Matrix_check_class(class_P(x), valid);

    if (ctype < 0) error("invalid class of object to Matrix_as_csn");
    ans->U = Matrix_as_cs(GET_SLOT(x, install("U")));
    ans->L = Matrix_as_cs(GET_SLOT(x, install("L")));
    switch(ctype) {
    case 0:
	ans->B = (double*) NULL;
	ans->pinv = INTEGER(GET_SLOT(x, install("Pinv")));
	break;
    case 1:
	ans->B = REAL(GET_SLOT(x, install("beta")));
	ans->pinv = (int*) NULL;
	break;
    default:
	error("invalid class of object to Matrix_as_csn");
    }
    return ans;
}

/**
 * Copy the contents of S to a css_LU or css_QR object and,
 * optionally, free S or free both S and the pointers to its contents.
 *
 * @param a css object to be converted
 * @param cl the name of the S4 class of the object to be generated
 * @param dofree 0 - don't free a; > 0 cs_free a; < 0 Free a
 * @param m number of rows in original matrix
 * @param n number of columns in original matrix
 *
 * @return SEXP containing a copy of S
 */
SEXP Matrix_css_to_SEXP(css *S, char *cl, int dofree, int m, int n)
{
    SEXP ans;
    char *valid[] = {"css_LU", "css_QR", ""};
    int *nz, ctype = Matrix_check_class(cl, valid);

    if (ctype < 0)
	error("Inappropriate class '%s' for Matrix_css_to_SEXP", cl);
    ans = PROTECT(NEW_OBJECT(MAKE_CLASS(cl)));
				/* allocate and copy common slots */
    Memcpy(INTEGER(ALLOC_SLOT(ans, install("Q"), INTSXP, n)), S->q, n);
    nz = INTEGER(ALLOC_SLOT(ans, install("nz"), INTSXP, 3));
    nz[0] = S->m2; nz[1] = S->lnz; nz[2] = S->unz;
    switch(ctype) {
    case 0:
	break;
    case 1:
	Memcpy(INTEGER(ALLOC_SLOT(ans, install("Pinv"), INTSXP, m)),
	       S->pinv, m);
	Memcpy(INTEGER(ALLOC_SLOT(ans, install("parent"), INTSXP, n)),
	       S->parent, n);
	Memcpy(INTEGER(ALLOC_SLOT(ans, install("cp"), INTSXP, n)),
	       S->cp, n);
	break;
    default:
	error("Inappropriate class '%s' for Matrix_css_to_SEXP", cl);
    }
    if (dofree > 0) cs_sfree(S);
    if (dofree < 0) Free(S);
    UNPROTECT(1);
    return ans;
}

/**
 * Copy the contents of N to a csn_LU or csn_QR object and,
 * optionally, free N or free both N and the pointers to its contents.
 *
 * @param a csn object to be converted
 * @param cl the name of the S4 class of the object to be generated
 * @param dofree 0 - don't free a; > 0 cs_free a; < 0 Free a
 *
 * @return SEXP containing a copy of S
 */
SEXP Matrix_csn_to_SEXP(csn *N, char *cl, int dofree)
{
    SEXP ans;
    char *valid[] = {"csn_LU", "csn_QR", ""};
    int ctype = Matrix_check_class(cl, valid), n = (N->U)->n;

    if (ctype < 0)
	error("Inappropriate class '%s' for Matrix_csn_to_SEXP", cl);
    ans = PROTECT(NEW_OBJECT(MAKE_CLASS(cl)));
				/* allocate and copy common slots */
    /* FIXME: Use the triangular matrix classes for csn_LU */
    SET_SLOT(ans, install("L"),	/* these are free'd later if requested */
	     Matrix_cs_to_SEXP(N->L, "dgCMatrix", 0));
    SET_SLOT(ans, install("U"),
	     Matrix_cs_to_SEXP(N->U, "dgCMatrix", 0));
    switch(ctype) {
    case 0:
	Memcpy(INTEGER(ALLOC_SLOT(ans, install("Pinv"), INTSXP, n)),
	       N->pinv, n);
	break;
    case 1:
	Memcpy(REAL(ALLOC_SLOT(ans, install("beta"), REALSXP, n)),
	       N->B, n);
	break;
    default:
	error("Inappropriate class '%s' for Matrix_csn_to_SEXP", cl);
    }
    if (dofree > 0) cs_nfree(N);
    if (dofree < 0) {
	Free(N->L); Free(N->U); Free(N);
    }
    UNPROTECT(1);
    return ans;
}

#endif	/* unused */
