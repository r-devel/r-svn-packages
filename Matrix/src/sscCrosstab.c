#include "sscCrosstab.h"

SEXP sscCrosstab(SEXP flist, SEXP upper)
{
    int **fpt, *Ap, *Gp, *TTi, *Ti, *Tj, *dims, i,
	ncol = 0,
	nfac = length(flist),
	nfc2 = (nfac * (nfac - 1))/2, /* nfac choose 2 */
	nobs = length(VECTOR_ELT(flist, 0)),
	ntrpl, nz, pos, up = asLogical(upper);
    double *TTx, *Tx;
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("sscCrosstab")));

    if (!isNewList(flist) || nfac < 1)
	error("flist must be a non-empty list");
    SET_SLOT(val, Matrix_GpSym, allocVector(INTSXP, nfac + 1));
    Gp = INTEGER(GET_SLOT(val, Matrix_GpSym));
    fpt = (int **) R_alloc(nfac, sizeof(int *));
    for (i = 0; i < nfac; i++) {
	SEXP el = VECTOR_ELT(flist, i);
	if (!inherits(el, "factor"))
	    error("flist must be a non-empty list of factors");
	if (length(el) != nobs)
	    error("All elements of flist must have the same length");
	Gp[i] = ncol;
	ncol += length(getAttrib(el, R_LevelsSymbol));
	fpt[i] = INTEGER(el);
    }
    Gp[nfac] = ncol;
    SET_SLOT(val, Matrix_uploSym, ScalarString(mkChar(up ? "U" : "L")));
    SET_SLOT(val, Matrix_DimSym, allocVector(INTSXP, 2));
    dims = INTEGER(GET_SLOT(val, Matrix_DimSym));
    dims[0] = dims[1] = ncol;
    ntrpl = nfc2 * nobs + ncol;
    Ti = Calloc(ntrpl, int); Tj = Calloc(ntrpl, int); TTi = Calloc(ntrpl, int);
    Tx = Calloc(ntrpl, double); TTx = Calloc(ntrpl, double);
				/* Generate the triplet form of the result */
    for (i = 0; i < ncol; i++) {
	Ti[i] = Tj[i] = i;	/* The diagonals - these will store counts */
	Tx[i] = 0.0;
    }
    pos = ncol;
    for (i = 0; i < nobs ; i++) {
	int j, jcol, k;
	for (j = 0; j < nfac; j++) {
	    jcol = Gp[j] + fpt[j][i] - 1;
	    Tx[jcol] += 1.;	/* increment diagonal count */
	    for (k = j + 1; k < nfac; k++) { /* off-diagonals */
		int irow = Gp[k] + fpt[k][i] - 1;
		if (up) {
		    Ti[pos] = jcol; Tj[pos] = irow;
		} else {
		    Tj[pos] = jcol; Ti[pos] = irow;
		}
		Tx[pos] = 1.;
		pos++;
	    }
	}
    }
    SET_SLOT(val, Matrix_pSym, allocVector(INTSXP, ncol + 1));
    Ap = INTEGER(GET_SLOT(val, Matrix_pSym));
    triplet_to_col(ncol, ncol, ntrpl, Ti, Tj, Tx, Ap, TTi, TTx);
    nz = Ap[ncol];		/* non-zeros in Z'Z crosstab */
    SET_SLOT(val, Matrix_iSym, allocVector(INTSXP, nz));
    SET_SLOT(val, Matrix_xSym, allocVector(REALSXP, nz));
    Memcpy(INTEGER(GET_SLOT(val, Matrix_iSym)), TTi, nz);
    Memcpy(REAL(GET_SLOT(val, Matrix_xSym)), TTx, nz);
    Free(Ti); Free(Tj); Free(Tx);
    Free(TTi); Free(TTx);
    
    UNPROTECT(1);
    return val;
}

extern void ssclme_fill_LIp(int n, const int Parent[], int LIp[]);

SEXP sscCrosstab_L_LI_sizes(SEXP ctab, SEXP permexp)
{
    SEXP ans = PROTECT(allocVector(INTSXP, 4));
    int *Ai = INTEGER(GET_SLOT(ctab, Matrix_iSym)),
	*Ap = INTEGER(GET_SLOT(ctab, Matrix_pSym)),
	*aa = INTEGER(ans),
	*perm = INTEGER(permexp),
	n = INTEGER(GET_SLOT(ctab, Matrix_DimSym))[1],
	*Lp = Calloc(n + 1, int),
	*Parent = Calloc(n, int),
	*Lnz = Calloc(n, int),
	*Flag = Calloc(n, int);

    ldl_symbolic(n, Ap, Ai, Lp, Parent, Lnz, Flag,
		 (int *) NULL, (int *) NULL); /* P & Pinv */
    aa[0] = Lp[n];
    ssclme_fill_LIp(n, Parent, Lp);
    aa[1] = Lp[n];
    ssc_symbolic_permute(n, 1, perm, Ap, Ai);
    ldl_symbolic(n, Ap, Ai, Lp, Parent, Lnz, Flag,
		 (int *) NULL, (int *) NULL); /* P & Pinv */
    aa[2] = Lp[n];
    ssclme_fill_LIp(n, Parent, Lp);
    aa[3] = Lp[n];
    Free(Flag); Free(Lnz); Free(Parent); Free(Lp); 
    UNPROTECT(1);
    return ans;
}

/** 
 * Check for a nested series of grouping factors in the sparse,
 * symmetric representation of the pairwise cross-tabulations.
 * 
 * @param n size of pairwise cross-tabulation matrix
 * @param nf number of groups of columns in pairwise cross-tabulation
 * @param upper non-zero if the upper triangle is stored
 * @param Ap array of pointers to columns
 * @param Ai row indices
 * @param Gp array of pointers to groups
 * 
 * @return 0 for non-nested groups, 1 for nested groups
 */
static
int ctab_isNested(int n, int nf, int upper,
		  const int Ap[], const int Ai[], const int Gp[])
{
    if (nf > 1) {  /* single factor always nested */
	int  i;
	if (upper) {
	    int *nnz = (int *) R_alloc(n, sizeof(int)), nz = Ap[n];
				/* count number of nonzeros in each row */
	    for (i = 0; i < n; i++) nnz[i] = 0;
	    for (i = 0; i < nz; i++) nnz[Ai[i]]++;
	    for (i = 0; i < nf; i++) {
		int j, p2 = Gp[i+1], target = nf - i;
		for (j = Gp[i]; j < p2; j++) {
		    if (nnz[j] != target) return 0;
		}
	    }
	} else {		/* lower triangle - the easy case */
	    for (i = 0; i < nf; i++) {
		int j, p2 = Gp[i+1], target = nf - i;
		for (j = Gp[i]; j < p2; j++) {
		    if ((Ap[j+1] - Ap[j]) != target)
			return 0;
		}
	    }
	}
    }
    return 1;
}

/** 
  * Determine the inverse of the permutation of the rows of the sparse,
  * column-oriented matrix represented by m, n, Ap and Ai that will
  * concentrate non-zeros in the lower rows.
  * 
  * @param m number of rows in the matrix
  * @param n number of columns in the matrix
  * @param Ap column pointers in Ai (modified)
  * @param Ai row indices (modified)
  * @param rperm on return contains the permutation of the rows
  * @param cperm if non-null, on return contains the permutation of the columns
  */
static
void pair_perm(int m, int n, int Ap[], int Ai[], int rperm[], int cperm[])
{
    int *cc = Calloc(n, int),	/* column counts */
	*cm = Calloc(n, int),	/* column removed */
	*sm = Calloc(m, int),	/* sum of column removals for this row */
	ii, j, pc,
	*rc = Calloc(m, int);	/* row counts */
    
    for (j = 0; j < n; j++) {	/* initialize col counts */
	cc[j] = Ap[j+1] - Ap[j];
	cm[j] = 0;
    }
    pc = 0;
    for (ii = m - 1; 0 <= ii; ii--) { /* fill rperm from RHS */
	int maxrc, rr;
	int i, mincc, p1, p3;
	
	mincc = m + 1;		/* find minimum positive cc */
	for (j = 0; j < n; j++) {
	    if (0 < cc[j] && cc[j] < mincc) mincc = cc[j];
	    if (mincc < 2) break;
	}
	
	for (i = 0; i < m; i++) {sm[i] = rc[i] = 0;}
	for (j = 0; j < n; j++) { /* row counts for cols where cc = mincc */
	    if (cc[j] == mincc) {
		int p2 = Ap[j+1];
		for (i = Ap[j]; i < p2; i++) {
		    rc[Ai[i]]++;
		    sm[Ai[i]] += cm[j];
		}
	    }
	}
	maxrc = -1;		/* find rows with rc[i] == max(rc) */
	for (i = 0; i < m; i++) {
	    int ic = rc[i];
				/* Choose first on row count.  Ties go
	                         * to smaller sum of moved.  Ties
	                         * there go to the last one. */
	    if (ic > maxrc || (ic == maxrc && sm[i] >= sm[rr])) {
		maxrc = ic;
		rr = i;
	    }
	}

	rperm[rr] = ii;

	p1 = p3 = 0;		/* update cc, Ap and Ai */
	for (j = 0; j < n; j++) {
	    int p2 = Ap[j+1];
	    for (i = Ap[j]; i < p2; i++) {
		if (Ai[i] == rr) {
		    cc[j]--; cm[j]++; /* move from count to removed */
		    if (cperm && cc[j] < 1) cperm[j] = pc++;
		} else {
		    if (i != p1) Ai[p1] = Ai[i];
		    p1++;
		}
	    }
	    Ap[j] = p3;
	    p3 = p1;		/* save current pos for next iteration */
	}
	Ap[n] = p3;
    }
    Free(cc); Free(cm); Free(rc); 
}

SEXP sscCrosstab_groupedPerm(SEXP ctab)
{
    SEXP
	GpSlot = GET_SLOT(ctab, Matrix_GpSym),
	iSlot = GET_SLOT(ctab, Matrix_iSym),
	pSlot = GET_SLOT(ctab, Matrix_pSym);
    int *Ai = INTEGER(iSlot),
	*Ap = INTEGER(pSlot),
	*Gp = INTEGER(GpSlot),
	i,
	n = length(pSlot) - 1,	/* number of columns */
	nf = length(GpSlot) - 1, /* number of factors */
	up;
    SEXP ans = PROTECT(allocVector(INTSXP, n));

    up = toupper(*CHAR(STRING_ELT(GET_SLOT(ctab, Matrix_uploSym), 0))) != 'L';
    if (nf > 1 && up) {			/* transpose */
	int nz = length(iSlot);
	int *ai = Calloc(nz, int),
	    *ap = Calloc(n + 1, int);
	double *ax = Calloc(nz, double);

	csc_components_transpose(n, n, nz, Ap, Ai,
				 REAL(GET_SLOT(ctab, Matrix_xSym)),
				 ap, ai, ax);
	Ap = ap;
	Ai = ai;
	Free(ax);		/* don't need values, only positions */
    }
    for (i = 0; i < n; i++) {
	INTEGER(ans)[i] = i;    /* initialize permutation to identity */
    }
    for (i = 1; i < nf; i++) {
	col_metis_order(Gp[i - 1], Gp[i], Gp[i+1], Ap, Ai, INTEGER(ans));
    }
    if (0) {			/* skip this part */
	if (!ctab_isNested(n, nf, 0, Ap, Ai, Gp)) {
	    int *np = Calloc(n + 1, int), /* column pointers */
		*ni = Calloc(length(iSlot) - n, int); /* row indices */
	    
	    np[0] = 0;
	    for (i = 1; i < nf; i++) { /* adjacent pairs of grouping factors */
		int j, k, p0 = 0, p1 = Gp[i-1], p2 = Gp[i], p3 = Gp[i+1];
		
		for (j = p1; j < p2; j++) { /* for this set of columns */
		    int lk = Ap[j+1];
		    for (k = Ap[j]; k < lk; k++) {
			int ii = Ai[k];
			if (p2 <= ii && ii < p3) { /* check the row */
			    ni[p0++] = ii - p2;
			}
		    }
		    np[j + 1 - p1] = p0;
		}
		pair_perm(p3 - p2, p2 - p1, np, ni,
			  INTEGER(ans) + p2,
			  (i == 1) ? INTEGER(ans) + p1 : (int *) 0);
		for (j = p2; j < p3; j++) INTEGER(ans)[j] += p2;
	    }
	    Free(np); Free(ni);
	}
    }
    if (nf > 1 && up) {Free(Ap); Free(Ai);}
    UNPROTECT(1);
    return ans;
}
