/* Produced by
 * $Id$
 */
/* spannel.f -- translated by f2c (version 20010821).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)


-- be brave, try without >>>
 * #include "f2c.h" /* <<<<< ------*/
#include <math.h>
#ifndef max
# define	max(a, b) 		((a) < (b) ? (b) : (a))
#endif
#ifndef min
# define	min(a, b)		((a) > (b) ? (b) : (a))
#endif
#ifndef abs
# define	abs(x)			((x) >= 0 ? (x) : -(x))
#endif


/* Table of constant values */

static int c__0 = 0;

/* Compute the SPANNing ELlipsoid -- 
 --------------- for clusplot.default(*, span = TRUE) 

 Subroutine */ int spannel_(int *ncas, int *ndep, double *dat,

	double *dstopt, double *cov, double *varsum, double *
	varss, double *prob, double *work, double *eps, int *
	maxit, int *ierr)
{
    /* System generated locals */
    int dat_dim1, dat_offset, cov_dim1, cov_offset, i__1, i__2, i__3;
    double d__1;

    /* Builtin functions */
    double sqrt(double);

    /* Local variables */
    double scal, dmax__, aver, dist;
    int loop, i__, j, k;
    double p, deter;
    extern /* Subroutine */ int sweep_(double *, int *, int *,

	    int *, double *);
    double tempo;


/*    ncas = number of objects. 
    ndep = number of variables. 
    maxit= maximal # iterations (and returns #{iter.}) 
    dstopt = squared distances 
 Local Variables 
 (in clusplot's call,)  dat[i,0] are all  == 1 
 Scale Data dat[i,j] to mean = 0 and var{1/n} = 1  -- for j= 1:ndep (not j=0!) 
     Parameter adjustments */
    --prob;
    --dstopt;
    --varss;
    --varsum;
    cov_dim1 = *ndep - 0 + 1;
    cov_offset = 0 + cov_dim1 * 0;
    cov -= cov_offset;
    dat_dim1 = *ncas;
    dat_offset = 1 + dat_dim1 * 0;
    dat -= dat_offset;

    /* Function Body */
    for (j = 1; j <= *ndep; ++j) { /* f2c-clean: s {i__1} {*ndep} */
	varsum[j] = 0.;
	varss[j] = 0.;
    }
    for (i__ = 1; i__ <= *ncas; ++i__) { /* f2c-clean: s {i__1} {*ncas} */
	for (j = 1; j <= *ndep; ++j) { /* f2c-clean: s {i__2} {*ndep} */
	    varsum[j] += dat[i__ + j * dat_dim1];
/* Computing 2nd power */
	    d__1 = dat[i__ + j * dat_dim1];
	    varss[j] += d__1 * d__1;
	}
    }
    for (j = 1; j <= *ndep; ++j) { /* f2c-clean: s {i__2} {*ndep} */
	aver = varsum[j] / *ncas;
	scal = sqrt(varss[j] / *ncas - aver * aver);
	for (i__ = 1; i__ <= *ncas; ++i__) { /* f2c-clean: s {i__1} {*ncas} */
	    dat[i__ + j * dat_dim1] = (dat[i__ + j * dat_dim1] - aver) / scal;
	}
    }
    p = 1.f / (double) (*ncas);
    for (i__ = 1; i__ <= *ncas; ++i__) { /* f2c-clean: s {i__1} {*ncas} */
	prob[i__] = p;
    }
    *ierr = 0;
    p = (double) (*ndep);
/* ---- Repeat { ... up to `maxit' times ] */
    loop = 0;
L160:
    ++loop;
/*     Cov[,] = weighted covariance of dat[,]  {weights = prob[]} */
    for (j = 0; j <= *ndep; ++j) { /* f2c-clean: s {i__1} {*ndep} */
	for (k = 0; k <= j; ++k) { /* f2c-clean: s {i__2} {j} */
	    cov[k + j * cov_dim1] = 0.f;
	}
    }
    for (i__ = 1; i__ <= *ncas; ++i__) { /* f2c-clean: s {i__2} {*ncas} */
	for (j = 0; j <= *ndep; ++j) { /* f2c-clean: s {i__1} {*ndep} */
	    work[j] = dat[i__ + j * dat_dim1];
	    tempo = prob[i__] * work[j];
	    for (k = 0; k <= j; ++k) { /* f2c-clean: s {i__3} {j} */
		cov[k + j * cov_dim1] += tempo * work[k];
	    }
	}
    }
    for (j = 0; j <= *ndep; ++j) { /* f2c-clean: s {i__2} {*ndep} */
	for (k = 0; k <= j; ++k) { /* f2c-clean: s {i__3} {j} */
	    cov[j + k * cov_dim1] = cov[k + j * cov_dim1];
	}
    }
    deter = 1.f;
    for (i__ = 0; i__ <= *ndep; ++i__) { /* f2c-clean: s {i__3} {*ndep} */
	if (deter <= 0.f) {
	    *ierr = 2;
	    return 0;
	}
	sweep_(&cov[cov_offset], ndep, &c__0, &i__, &deter);
    }
    dmax__ = 0.f;
    for (i__ = 1; i__ <= *ncas; ++i__) { /* f2c-clean: s {i__3} {*ncas} */
	dist = -1.f;
	for (j = 0; j <= *ndep; ++j) { /* f2c-clean: s {i__2} {*ndep} */
	    work[j] = 0.;
	    for (k = 0; k <= *ndep; ++k) { /* f2c-clean: s {i__1} {*ndep} */
		work[j] -= cov[j + k * cov_dim1] * dat[i__ + k * dat_dim1];
	    }
/*           work(j) = - sum_{k=0}^p  dat(i,k) * cov(k,j) { = cov(j,k) }, 
     i.e., work_j = - X[i,] %*% COV[,j] */
	    dist += work[j] * dat[i__ + j * dat_dim1];
	}
	dstopt[i__] = dist;
/*        Dist{opt}_i = -1 - t(X[i,]) %*% COV %*% X[i,] */
	if (dist > dmax__) {
	    dmax__ = dist;
	}
    }
/*     dmax = max{ dstopt[i] } */
    if (dmax__ > p + *eps) {
/*     not yet converged */
	for (i__ = 1; i__ <= *ncas; ++i__) { /* f2c-clean: s {i__3} {*ncas} */
	    prob[i__] = prob[i__] * dstopt[i__] / p;
	}
	if (loop < *maxit) {
	    goto L160;
	}
    }
    *maxit = loop;
    return 0;
} /* spannel_ 

 Subroutine */ int sweep_(double *cov, int *nord, int *ixlo,

	int *nel, double *deter)
{
    /* System generated locals */
    int cov_dim1, cov_offset, i__1, i__2;

    /* Local variables */
    double temp;
    int i__, j;



    /* Parameter adjustments */
    cov_dim1 = *nord - 0 + 1;
    cov_offset = 0 + cov_dim1 * 0;
    cov -= cov_offset;

    /* Function Body */
    temp = cov[*nel + *nel * cov_dim1];
    *deter *= temp;
    if (*nord <= 1) {
	cov[cov_dim1 + 1] = 1.f / temp;
	return 0;
    }
/* else : nord > 1 */
    for (i__ = *ixlo; i__ <= *nord; ++i__) { /* f2c-clean: s {i__1} {*nord} */
	if (i__ == *nel) {
	    goto L30;
	}
	for (j = *ixlo; j <= i__; ++j) { /* f2c-clean: s {i__2} {i__} */
	    if (j == *nel) {
		goto L20;
	    }
	    cov[j + i__ * cov_dim1] = cov[i__ + j * cov_dim1] - cov[i__ + *
		    nel * cov_dim1] * cov[*nel + j * cov_dim1] / temp;
	    cov[i__ + j * cov_dim1] = cov[j + i__ * cov_dim1];
L20:
	    ;
	}
L30:
	;
    }
    cov[*nel + *nel * cov_dim1] = 1.f;
    for (i__ = *ixlo; i__ <= *nord; ++i__) { /* f2c-clean: s {i__1} {*nord} */
	cov[*nel + i__ * cov_dim1] = -cov[i__ + *nel * cov_dim1] / temp;
	cov[i__ + *nel * cov_dim1] = cov[*nel + i__ * cov_dim1];
    }
    return 0;
} /* sweep_ */

