/* Compute the SPANNing ELlipsoid --
 --------------- for clusplot.default(*, span = TRUE)
*/

/* Original spannel.f -- translated by f2c (version 20010821).
 * and f2c-clean,v 1.10 2002/03/28 16:37:27 maechler
 */

#include <math.h>

#include "cluster.h"

void sweep_(double *, int *, int *, int *, double *);

/* Table of constant values */

static int c__0 = 0;

void spannel_(int *ncas, /* = number of objects */
	      int *ndep, /* = number of variables */
	      double *dat,/* [ncas, 0:ndep] */
	      double *dstopt, /* = squared distances [1:ncas] */
	      double *cov,/* matrix [0:ndep, 0:ndep] */
	      double *varsum,	/* [1:ndep] */
	      double *varss,	/* [1:ndep] */
	      double *prob, 	/* [1:ncas] */
	      double *work, 	/* [0:ndep] */
	      double *eps,
	      int *maxit, /* = maximal # iterations (and returns #{iter.})*/
	      int *ierr)
{
    /* System generated locals */
    int dat_dim1, dat_offset, cov_dim1, cov_offset, i__1, i__2, i__3;
    double d__1;

    /* Builtin functions */
    double sqrt(double);

    /* Local variables */
    double scal, dmax__, aver, dist;
    int loop, i__, j, k;
    double p, deter, tempo;

    /* Parameter adjustments */
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
    for (i__ = 1; i__ <= *ncas; ++i__) {
	for (j = 1; j <= *ndep; ++j) {
	    varsum[j] += dat[i__ + j * dat_dim1];
/* Computing 2nd power */
	    d__1 = dat[i__ + j * dat_dim1];
	    varss[j] += d__1 * d__1;
	}
    }
    for (j = 1; j <= *ndep; ++j) {
	aver = varsum[j] / *ncas;
	scal = sqrt(varss[j] / *ncas - aver * aver);
	for (i__ = 1; i__ <= *ncas; ++i__) {
	    dat[i__ + j * dat_dim1] = (dat[i__ + j * dat_dim1] - aver) / scal;
	}
    }
    p = 1.f / (double) (*ncas);
    for (i__ = 1; i__ <= *ncas; ++i__)
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
	    return;
	}
	sweep_(&cov[cov_offset], ndep, &c__0, &i__, &deter);
    }
    dmax__ = 0.f;
    for (i__ = 1; i__ <= *ncas; ++i__) {
	dist = -1.f;
	for (j = 0; j <= *ndep; ++j) {
	    work[j] = 0.;
	    for (k = 0; k <= *ndep; ++k)
		work[j] -= cov[j + k * cov_dim1] * dat[i__ + k * dat_dim1];

/*          work(j) = - sum_{k=0}^p  dat(i,k) * cov(k,j) { = cov(j,k) },
     i.e.,  work_j = - X[i,] %*% COV[,j] */
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

