/* FANNY : program for Fuzzy cluster ANalysis */

/* was $Id$
 * fanny.f -- translated by f2c (version 20020621).
 * and treated by  f2c-clean v 1.10, and manually by Martin Maechler
 */

#include <Rmath.h>
#include "cluster.h"

/* dysta3_() is in cluster.h ! */

static void
caddy(int *, double *, int *,
      int *, int *, int *, double *);
static void
fygur(int *, int *, int *,
      int *, int *, int *, int *, int *, double *,
      double *, double *, double *, double *, double *, double *);
static void
fuzzy(int *, int *, double *, double *, double *, double *,
      double *, double *, double *, double *, int *,
      double *, double *, double *, int *);


void fanny(int *nn, int *jpp, int *kk,
	   double *x, double *dss, int *jdyss, double *valmd,
	   int *jtmd, int *ndyst, int *nsend, int *nelem,
	   int *negbr, double *syl, double *p, double *dp,
	   double *pt, int *nfuzz, double *esp, double *ef,
	   double *dvec, double *ttsyl, double *eda, double *edb,
	   double *obj, int *ncluv, double *sylinf, double *r,
	   double *tol, int *maxit)
{
/* Arguments
   nn	= number of objects
   jpp	= number of variables for clustering
   kk	= number of clusters
*/

    /* Local variables */
    int l, nhalf, jhalt, ktrue;
    double s;


    if (*jdyss != 1) { /* compute dissimilarities from data */
	jhalt = 0;
	dysta3_(nn, jpp, x, dss, ndyst, jtmd, valmd, &jhalt);
	if (jhalt != 0) {
	    *jdyss = -1;
	    return;
	}
    }

    nhalf = *nn * (*nn - 1) / 2;

    fuzzy(nn, &nhalf, p, dp, pt, dss, esp,
	  ef, eda, edb, kk, obj, r, tol, maxit);

    caddy(nn, p, kk, &ktrue, nfuzz, ncluv, pt);

    /*  Compute "silhouette": */
    if (2 <= ktrue && ktrue < *nn) {
	/* s := max( dss[i,j] ) */
	for(l = 0, s = 0.; l < nhalf; l++)
	    if (s < dss[l])
		s = dss[l];
	fygur(&ktrue, nn, kk, &nhalf, ncluv, nsend, nelem,
	      negbr, syl, dvec, pt, ttsyl, dss, &s, sylinf);
    }
    return;
} /* fanny */


int dysta3_(int *nn, int *jpp, double *x, double *dys,
	    int *ndyst, int *jtmd, double *valmd, int *jhalt)
{
    /* System generated locals */
    int x_dim1, x_offset;


    /* Local variables */
    int j, k, l, nlk, npres;
    double clk, rpres;

    /* Parameter adjustments */
    --dys;
    --valmd;
    --jtmd;
    x_dim1 = *nn;
    x_offset = 1 + x_dim1 * 1;
    x -= x_offset;

    /* Function Body */
    nlk = 0;
    for (l = 1; l <= (*nn - 1); ++l) {
	for (k = l + 1; k <= *nn; ++k) {
	    clk = 0.;
	    ++nlk;
	    npres = 0;
	    for (j = 1; j <= *jpp; ++j) {
		if (jtmd[j] < 0) {
		    if (x[l + j * x_dim1] == valmd[j] ||
			x[k + j * x_dim1] == valmd[j])

			continue; /* next j */
		}
		++npres;
		if (*ndyst == 1)
		    clk += (x[l + j * x_dim1] - x[k + j * x_dim1]) *
			   (x[l + j * x_dim1] - x[k + j * x_dim1]);
		else
		    clk += fabs(x[l + j * x_dim1] - x[k + j * x_dim1]);
	    }
	    if (npres == 0) {
		*jhalt = 1;
		dys[nlk] = -1.;
	    } else {
		double p_r = (*jpp) / (double) npres;
		if (*ndyst == 1)
		    dys[nlk] = sqrt(clk * p_r);
		else
		    dys[nlk] = clk * p_r;
	    }
	}
    }
    return 0;
} /* dysta3_ */


static
void fuzzy(int *nn, int *hh, double *p,
	   double *dp, double *pt, double *dss, double *esp,
	   double *ef, double *eda, double *edb, int *k,
	   double *obj,
	   double *r,  /* the exponent, > 1. -- was fixed to 2 originally */
	   double *tol,/* the precision for the iterations */
	   int *nit)   /* the maximal number of iterations -- was fixed to 500*/
{
    /* Local variables */
    double p0, dt, zk, xx, ddd, crt, reen, cryt, rvers;
    int j, l, m, nd, jm, it, lx, ndk;

    /* System generated locals */
    int p_dim1, p_offset, dp_dim1, dp_offset;
    double d__1;
    /* Parameter adjustments */
    --dss;
    --ef;
    --esp;
    --pt;
    dp_dim1 = *nn;
    dp_offset = 1 + dp_dim1 * 1;
    dp -= dp_offset;
    p_dim1 = *nn;
    p_offset = 1 + p_dim1 * 1;
    p -= p_offset;


    /* Function Body */
    rvers = 1. / *r;
    reen = 1. / (*r - 1.);

/*     initial fuzzy clustering */

    p0 = 0.1 / (*k - 1);
    for (m = 1; m <= *nn; ++m) {
	for (l = 1; l <= *k; ++l) {
	    dp[m + l * dp_dim1] = 0.;
	    p[m + l * p_dim1] = p0;
	}
    }
    ndk = *nn / *k;
    nd = ndk;
    l = 1;
    for (m = 1; m <= *nn; ++m) {
	p[m + l * p_dim1] = 0.9;
	if (m >= nd) {
	    nd += ndk;
	    ++l;
	    if (l == *k) {
		nd = *nn;
	    }
	}
	for (lx = 1; lx <= *k; ++lx) {
	    p[m + lx * p_dim1] = pow(p[m + lx * p_dim1], *r);
	}
    }

/*     initial criterion value */

    cryt = 0.;
    for (l = 1; l <= *k; ++l) {
	esp[l] = 0.;
	ef[l] = 0.;
	for (m = 1; m <= *nn; ++m) {
	    esp[l] += p[m + l * p_dim1];
	    for (j = 1; j <= *nn; ++j) {
		if (j != m) {
		    jm = imin2(m,j);
		    jm = (jm - 1) * *nn - jm * (jm + 1) / 2 + imax2(m,j);
		    dp[m + l * dp_dim1] += p[j + l * p_dim1] * dss[jm];
		    ef[l] += p[j + l * p_dim1] * p[m + l * p_dim1] * dss[jm];
		}
	    }
	}
	cryt += ef[l] / (esp[l] * 2.);
    }
    crt = cryt;

/*     start of iterations */

    it = 1;
    m = 0;

    do {
	/* the new membership coefficients of the objects are calculated,
	   and the resulting value of the criterion is computed. */
	++m;
	dt = 0.;
	for (l = 1; l <= *k; ++l) {
	    pt[l] = pow(esp[l] * 2. * esp[l] /
			(esp[l] * 2. * dp[m + l * dp_dim1] - ef[l]),
			reen);
	    dt += pt[l];
	}
	xx = 0.;
	for (l = 1; l <= *k; ++l) {
	    pt[l] /= dt;
	    if (pt[l] < 0.)
		xx += pt[l];
	}
	/* now: sum_l (pt[l]) == 1;  xx := sum_{pt[l] < 0} pt[l] */
	for (l = 1; l <= *k; ++l) {
	    pt[l] = (pt[l] > 0.) ? pow(pt[l] / (1 - xx), *r) : 0.;
	    esp[l] += pt[l] - p[m + l * p_dim1];
	    for (j = 1; j <= *nn; ++j) {
		if (j != m) {
		    jm = imin2(m,j);
		    jm = (jm - 1) * *nn - jm * (jm + 1) / 2 + imax2(m,j);
		    ddd = (pt[l] - p[m + l * p_dim1]) * dss[jm];
		    dp[j + l * dp_dim1] += ddd;
		    ef[l] += p[j + l * p_dim1] * 2. * ddd;
		}
	    }
	    p[m + l * p_dim1] = pt[l];
	}

	if (m >= *nn) {
	    cryt = 0.;
	    *eda = 0.;
	    for (l = 1; l <= *k; ++l) {
		*eda += esp[l] / (double) *nn;
		cryt += ef[l] / (esp[l] * 2.);
	    }

	    /* Convergence check */
	    if (fabs(cryt - crt) <= *tol * cryt)
		break;
	    if (it >= *nit) { /* non-convergence in max_it iterations */
		*nit = -1;
		break;
	    }
	    m = 0;
	    ++it;
	    crt = cryt;
	}

    } while(1);


    /* non-fuzzyness index of libert is computed */

    obj[0] = (double) it;
    obj[1] = cryt;
    zk = (double) (*k);
    *edb = (zk * *eda - 1.) / (zk - 1.);
    for (m = 1; m <= *nn; ++m)
	for (l = 1; l <= *k; ++l)
	    p[m + l * p_dim1] = pow(p[m + l * p_dim1], rvers);

    return;
} /* fuzzy */


static
void caddy(int *nn, double *p, int *k, int *ktrue,
	   int *nfuzz, int *ncluv, double *rdraw)
{
    /* System generated locals */
    int p_dim1, p_offset;

    /* Local variables */
    Rboolean stay;
    int l, m, ktry, kleft, kwalk, nbest, lfuzz;
    double pbest;


    /* Parameter adjustments */
    --ncluv;
    --rdraw;
    --nfuzz;
    p_dim1 = *nn;
    p_offset = 1 + p_dim1 * 1;
    p -= p_offset;

    /* Function Body */
    pbest = p[p_dim1 + 1];
    nbest = 1;
    for (l = 2; l <= *k; ++l) {
	if (pbest < p[l * p_dim1 + 1]) {
	    pbest = p[l * p_dim1 + 1];
	    nbest = l;
	}
    }
    nfuzz[1] = nbest;
    ncluv[1] = 1;
    *ktrue = 1;
    for (m = 2; m <= *nn; ++m) {
	pbest = p[m + p_dim1];
	nbest = 1;
	for (l = 2; l <= *k; ++l) {
	    if (pbest < p[m + l * p_dim1]) {
		pbest = p[m + l * p_dim1];
		nbest = l;
	    }
	}
	stay = FALSE;
	for (ktry = 1; ktry <= *ktrue; ++ktry) {
	    if (nfuzz[ktry] == nbest) {
		stay = TRUE;
		ncluv[m] = ktry;
	    }
	}
	if (! stay) {
	    (*ktrue)++;
	    nfuzz[*ktrue] = nbest;
	    ncluv[m] = *ktrue;
	}
    }
    if (*ktrue < *k) {
	for (kwalk = *ktrue + 1; kwalk <= *k; ++kwalk) {
	    for (kleft = 1; kleft <= *k; ++kleft) {
		stay = FALSE;
		for (ktry = 1; ktry <= (kwalk - 1); ++ktry) {
		    if (nfuzz[ktry] == kleft)
			stay = TRUE;
		}
		if (! stay) {
		    nfuzz[kwalk] = kleft;
		    break;
		}
	    }
	}
    }
    for (m = 1; m <= *nn; ++m) {
	for (l = 1; l <= *k; ++l) {
	    lfuzz = nfuzz[l];
	    rdraw[l] = p[m + lfuzz * p_dim1];
	}
	for (l = 1; l <= *k; ++l) {
	    p[m + l * p_dim1] = rdraw[l];
	}
    }
    return;
} /* caddy */

/* -----------------------------------------------------------

     Compute Silhouette Information :

 TODO  cleanup: this is almost identical to black() in	pam.c()
   -- only difference : different  dys() indexing !
*/
static
void fygur(int *ktrue, int *nn, int *kk, int *hh,
	   int *ncluv, int *nsend, int *nelem, int *negbr,
	   double *syl, double *srank, double *avsyl, double *ttsyl,
	   double *dss, double *s, double *sylinf)
{
    /* System generated locals */
    int sylinf_dim1, sylinf_offset;

    /* Local variables */
    int j, l, nj, nl, nbb, mjl, njl, ntt, lang = -1 /*Wall*/;
    int nclu, lplac, numcl, nsylr;
    double db, btt, dysa, dysb, symax;

    /* Parameter adjustments */
    sylinf_dim1 = *nn;
    sylinf_offset = 1 + sylinf_dim1 * 1;
    sylinf -= sylinf_offset;
    --srank;
    --syl;
    --negbr;
    --nelem;
    --nsend;
    --ncluv;
    --avsyl;
    --dss;

    /* Function Body */
    nsylr = 0;
    *ttsyl = 0.;
    for (numcl = 1; numcl <= *ktrue; ++numcl) {
	ntt = 0;
	for (j = 1; j <= *nn; ++j) {
	    if (ncluv[j] == numcl) {
		++ntt;
		nelem[ntt] = j;
	    }
	}
	for (j = 1; j <= ntt; ++j) {
	    nj = nelem[j];
	    dysb = *s * 1.1 + 1.;
	    negbr[j] = -1;
	    for (nclu = 1; nclu <= *ktrue; ++nclu) {
		if (nclu != numcl) {
		    nbb = 0;
		    db = 0.;
		    for (l = 1; l <= *nn; ++l) {
			if (ncluv[l] == nclu) {
			    ++nbb;
			    if (l < nj) {
				mjl = *nn * (l - 1) + nj - l * (l + 1) / 2;
				db += dss[mjl];
			    } else if (l > nj) {
				mjl = *nn * (nj - 1) + l - nj * (nj + 1) / 2;
				db += dss[mjl];
			    } /* else dss(.)=0 ; nothing to add */
			}
		    }
		    btt = (double) nbb;
		    db /= btt;
		    if (dysb > db) {
			dysb = db;
			negbr[j] = nclu;
		    }
		}
	    }
	    if (ntt > 1) {
		dysa = 0.;
		for (l = 1; l <= ntt; ++l) {
		    nl = nelem[l];
		    if (nj < nl) {
			njl = *nn * (nj - 1) + nl - nj * (nj + 1) / 2;
			dysa += dss[njl];
		    } else if (nj > nl) {
			njl = *nn * (nl - 1) + nj - nl * (nl + 1) / 2;
			dysa += dss[njl];
		    }/* else dss(.)=0 ; nothing to add */
		}
		dysa /= ntt - 1;
		if (dysa > 0.) {
		    if (dysb > 0.) {
			if (dysb > dysa) {
			    syl[j] = 1. - dysa / dysb;
			} else if (dysb < dysa) {
			    syl[j] = dysb / dysa - 1.;
			} else { /* dysb == dysa: */
			    syl[j] = 0.;
			}
			if (syl[j] <= -1.)
			    syl[j] = -1.;
			else if (syl[j] >= 1.)
			    syl[j] = 1.;
		    } else {
			syl[j] = -1.;
		    }
		} else if (dysb > 0.) {
		    syl[j] = 1.;
		} else {
		    syl[j] = 0.;
		}
	    } else { /* ntt == 1: */
		syl[j] = 0.;
	    }
	}
	avsyl[numcl] = 0.;
	for (j = 1; j <= ntt; ++j) {
	    symax = -2.;
	    for (l = 1; l <= ntt; ++l) {
		if (symax < syl[l]) {
		    symax = syl[l];
		    lang = l;
		}
	    }
	    nsend[j] = lang;
	    srank[j] = syl[lang];
	    avsyl[numcl] += srank[j];
	    syl[lang] = -3.;
	}
	*ttsyl += avsyl[numcl];
	avsyl[numcl] /= (double) ntt;
	if (ntt < 2) {
	    ++nsylr;
	    sylinf[nsylr + sylinf_dim1] = (double) numcl;
	    sylinf[nsylr + (sylinf_dim1 << 1)] = (double) negbr[1];
	    sylinf[nsylr + sylinf_dim1 * 3] = 0.;
	    sylinf[nsylr + (sylinf_dim1 << 2)] = (double) nelem[1];
	}
	else {
	    for (l = 1; l <= ntt; ++l) {
		++nsylr;
		lplac = nsend[l];
		sylinf[nsylr + sylinf_dim1] = (double) numcl;
		sylinf[nsylr + (sylinf_dim1 << 1)] = (double) negbr[lplac];
		sylinf[nsylr + sylinf_dim1 * 3] = srank[l];
		sylinf[nsylr + (sylinf_dim1 << 2)] = (double) nelem[lplac];
	    }
	}
    }
    *ttsyl /= *nn;
    return;
} /* fygur */

