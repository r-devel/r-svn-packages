/* Compute the SPANNing ELlipsoid --
 --------------- for clusplot.default(*, span = TRUE)
*/

/* Original spannel.f -- translated by f2c (version 20010821).
 * and f2c-clean,v 1.10 2002/03/28 16:37:27 maechler
 */

#include <math.h>

#include "cluster.h"

void sweep(double *, int *, int *, int *, double *);

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
    int dat_dim1, dat_offset, cov_dim1;

    /* Local variables */
    int loop, i, j, k;
    double scal, dmax, aver, dist, p, deter, tempo;

#define COV(i,j) cov[i + j * cov_dim1]
    cov_dim1 = *ndep + 1;

    /* Parameter adjustments */
    --prob;
    --dstopt;
    --varss;
    --varsum;

    dat_dim1 = *ncas;
    dat_offset = 1 + dat_dim1 * 0;
    dat -= dat_offset;

    /* Function Body */

/* (in clusplot's call,)  dat[i,0] are all  == 1
 *
 * Scale Data dat[i,j] to mean = 0 and var{1/n} = 1  -- for j= 1:ndep (not j=0!)
 */
    for (j = 1; j <= *ndep; ++j) {
	varsum[j] = 0.;
	varss[j] = 0.;
    }
    for (i = 1; i <= *ncas; ++i) {
	for (j = 1; j <= *ndep; ++j) {
	    p = dat[i + j * dat_dim1];
	    varsum[j] += p;
	    varss [j] += p * p;
	}
    }
    for (j = 1; j <= *ndep; ++j) {
	aver = varsum[j] / *ncas;
	scal = sqrt(varss[j] / *ncas - aver * aver);
	for (i = 1; i <= *ncas; ++i) {
	    dat[i + j * dat_dim1] = (dat[i + j * dat_dim1] - aver) / scal;
	}
    }
    p = 1. / (double) (*ncas);
    for (i = 1; i <= *ncas; ++i)
	prob[i] = p;
    *ierr = 0;
    p = (double) (*ndep);
/* ---- Repeat { ... up to `maxit' times ] */
    loop = 0;
L160:
    ++loop;
/*     Cov[,] = weighted covariance of dat[,]  {weights = prob[]} */
    for (j = 0; j <= *ndep; ++j) {
	for (k = 0; k <= j; ++k)
	    COV(k,j) = 0.;
    }
    for (i = 1; i <= *ncas; ++i) {
	for (j = 0; j <= *ndep; ++j) {
	    work[j] = dat[i + j * dat_dim1];
	    tempo = prob[i] * work[j];
	    for (k = 0; k <= j; ++k)
		COV(k,j) += tempo * work[k];
	}
    }
    for (j = 0; j <= *ndep; ++j)
	for (k = 0; k <= j; ++k)
	    COV(j,k) = COV(k,j);

    deter = 1.;
    for (i = 0; i <= *ndep; ++i) {
	if (deter <= 0.) { *ierr = 2; return; }
	sweep(cov, ndep, &c__0, &i, &deter);
    }
    dmax = 0.;
    for (i = 1; i <= *ncas; ++i) {
	dist = -1.;
	for (j = 0; j <= *ndep; ++j) {
	    work[j] = 0.;
	    for (k = 0; k <= *ndep; ++k)
		work[j] -= COV(j,k) * dat[i + k * dat_dim1];

/*          work(j) = - sum_{k=0}^p  dat(i,k) * cov(k,j) { = cov(j,k) },
     i.e.,  work_j = - X[i,] %*% COV[,j] */
	    dist += work[j] * dat[i + j * dat_dim1];
	}
	dstopt[i] = dist;
/*        Dist{opt}_i = -1 - t(X[i,]) %*% COV %*% X[i,] */
	if (dmax < dist)
	    dmax = dist;
    }/*     dmax = max{ dstopt[i] } */

    if (dmax > p + *eps) { /*     not yet converged */
	for (i = 1; i <= *ncas; ++i)
	    prob[i] *= (dstopt[i] / p);
	if (loop < *maxit) {
	    goto L160;
	}
    }
    *maxit = loop;
    return;
} /* spannel_ */

/* This is currently also called from R : ../tests/sweep-ex.R
 * ==> keep pointers !*/
void sweep(double *cov, int *nord, int *ixlo, int *nel, double *deter)
{
    /* System generated locals */
    int cov_dim1;

    /* Local variables */
    double temp;
    int i, j;
    cov_dim1 = *nord + 1;

    temp = COV(*nel,*nel);
    *deter *= temp;
    if (*nord <= 1) {
	COV(1,1) = 1. / temp;
    }
    else { /* nord > 1 */
	for (i = *ixlo; i <= *nord; ++i) {
	    if (i != *nel) {
		for (j = *ixlo; j <= i; ++j) {
		    if (j != *nel) {
			COV(j,i) = COV(i,j) - COV(i,*nel) * COV(*nel,j) / temp;
			COV(i,j) = COV(j,i);
		    }
		}
	    }
	}
	COV(*nel,*nel) = 1.;
	for (i = *ixlo; i <= *nord; ++i) {
	    COV(*nel,i) = -COV(i,*nel) / temp;
	    COV(i,*nel) = COV(*nel,i);
	}
    }
    return;
} /* sweep */

