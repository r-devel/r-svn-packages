/* -*- mode: c; kept-new-versions: 25; kept-old-versions: 20 -*- */

/*
 * PAM := Partitioning Around Medoids
 *
 * $Id$
 * original Id: pam.f,v 1.16 2003/06/03 13:40:56 maechler translated by
 * f2c (version 20031025) and run through f2c-clean,v 1.10 2002/03/28
 */

#include <R_ext/Print.h>/* for diagnostics */

#include "cluster.h"

void pam(int *nn, int *jpp, int *kk,
	 double *x, double *dys, int *jdyss, double *valmd, int *jtmd,
	 int *ndyst, int *nsend, int *nrepr, int *nelem,
	 double *radus, double *damer, double *ttd, double *separ,
	 double *ttsyl, int *med, double *obj, int *ncluv,
	 double *clusinf, double *sylinf, int *nisol)
{
    /* System generated locals */
    int x_dim1, x_offset, clusinf_dim1, clusinf_offset, sylinf_dim1,
	sylinf_offset;

    /* Local variables */
    static int k, l;
    static double s, sky;
    static int nhalf, jhalt;

/*     carries out a clustering using the k-medoid approach.


 jdyss = 0 : compute distances from x
	= 1 : distances provided  in x
     Parameter adjustments */
    sylinf_dim1 = *nn;
    sylinf_offset = 1 + sylinf_dim1;
    sylinf -= sylinf_offset;
    --ncluv;
    --separ;
    --ttd;
    --damer;
    --radus;
    --nelem;
    --nrepr;
    --nsend;
    --dys;
    --jtmd;
    --valmd;
    x_dim1 = *nn;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --nisol;
    clusinf_dim1 = *kk;
    clusinf_offset = 1 + clusinf_dim1;
    clusinf -= clusinf_offset;
    --med;
    --obj;

    /* Function Body */
    if (*jdyss != 1) {
	jhalt = 0;
	F77_CALL(dysta)(nn, jpp, &x[x_offset], &dys[1], ndyst,
			&jtmd[1], &valmd[1], &jhalt);
	if (jhalt != 0) {
	    *jdyss = -1;
	    return;
	}
    }
/*     nhalf := #{distances}+1 = length(dys) */
    nhalf = *nn * (*nn - 1) / 2 + 1;
/*     s := max( dys[.] ), the largest distance */
    s = 0.;
    for (l = 2; l <= nhalf; ++l)
	if (s < dys[l])
	    s = dys[l];

/*     Build + Swap : */
    bswap(kk, nn, &nrepr[1], &radus[1], &damer[1], &ttd[1], &nhalf, &dys[1],
	   &sky, &s, &obj[1]);

/*     Compute STATs : */
    cstat(kk, nn, &nsend[1], &nrepr[1], &radus[1], &damer[1], &ttd[1],
	  &separ[1], &sky, &s, &nhalf, &dys[1], &ncluv[1], &nelem[1], &med[1],
	  &nisol[1]);
    for (k = 1; k <= *kk; ++k) {
	clusinf[k + clusinf_dim1] = (double) nrepr[k];
	clusinf[k + (clusinf_dim1 << 1)] = radus[k];
	clusinf[k + clusinf_dim1 * 3] = ttd[k];
	clusinf[k + (clusinf_dim1 << 2)] = damer[k];
	clusinf[k + clusinf_dim1 * 5] = separ[k];
    }
    if (1 < *kk && *kk < *nn) {
/*	 Compute Silhouette info : */
	dark(kk, nn, &nhalf, &ncluv[1], &nsend[1], &nelem[1], &nrepr[1], &
		radus[1], &damer[1], &ttd[1], ttsyl, &dys[1], &s, &sylinf[
		sylinf_offset]);
    }
} /* pam */

/* NOTE: dysta() is *kept* in Fortran for now --> ./dysta.f */

/* -----------------------------------------------------------

     bswap(): the clustering algorithm in 2 parts:  I. build,	II. swap
*/
void bswap(int *kk, int *nn, int *nrepr,
	   double *dysma, double *dysmb, double *beter, int *hh,
	   double *dys, double *sky, double *s, double *obj)
{
    int i__, j, ij, k, kj, kbest, nbest, njn, nmax;
    double ammax, small, cmd, dz, dzsky;

/* nrepr[]: here is boolean (0/1): 1 = "is representative object"  */

     /* Parameter adjustments */
    --beter;
    --dysmb;
    --dysma;
    --nrepr;
    --dys;
    --obj;

    /* -Wall: */
    nbest = -1;
    kbest = -1;

/*     first algorithm: build. */

    for (i__ = 1; i__ <= *nn; ++i__) {
	nrepr[i__] = 0;
	dysma[i__] = *s * 1.1f + 1.;
    }
    for (k = 1; k <= *kk; ++k) {
	for (i__ = 1; i__ <= *nn; ++i__) {
	    if (nrepr[i__] == 0) {
		beter[i__] = 0.;
		for (j = 1; j <= *nn; ++j) {
		    ij = F77_CALL(meet)(&i__, &j);
		    cmd = dysma[j] - dys[ij];
		    if (cmd > 0.) {
			beter[i__] += cmd;
		    }
		}
	    }
	}
	ammax = 0.;
	for (i__ = 1; i__ <= *nn; ++i__) {
	    if (nrepr[i__] == 0 && ammax <= beter[i__]) {
/*		    does < (instead of <= ) work too? -- NO! */
		ammax = beter[i__];
		nmax = i__;
	    }
	}
	nrepr[nmax] = 1;/* = .true. : *is* a representative */
	for (j = 1; j <= *nn; ++j) {
	    njn = F77_CALL(meet)(&nmax, &j);
	    if (dysma[j] > dys[njn]) {
		dysma[j] = dys[njn];
	    }
	}
    }
    *sky = 0.;
    for (j = 1; j <= *nn; ++j) {
	*sky += dysma[j];
    }
    obj[1] = *sky / *nn;
    if (*kk > 1) {

/*     second algorithm: swap.

 --   Loop : */
L60:
	for (j = 1; j <= *nn; ++j) {
	    dysma[j] = *s * 1.1f + 1.;
	    dysmb[j] = *s * 1.1f + 1.;
	    for (i__ = 1; i__ <= *nn; ++i__) {
		if (nrepr[i__] == 1) {
		    ij = F77_CALL(meet)(&i__, &j);
		    if (dys[ij] < dysma[j]) {
			dysmb[j] = dysma[j];
			dysma[j] = dys[ij];
		    } else {
			if (dys[ij] < dysmb[j]) {
			    dysmb[j] = dys[ij];
			}
		    }
		}
	    }
	}
	dzsky = 1.;
	for (k = 1; k <= *nn; ++k) {
	    if (nrepr[k] == 0) {
		for (i__ = 1; i__ <= *nn; ++i__) {
		    if (nrepr[i__] == 1) {
			dz = 0.;
			for (j = 1; j <= *nn; ++j) {
			    ij = F77_CALL(meet)(&i__, &j);
			    kj = F77_CALL(meet)(&k, &j);
			    if (dys[ij] == dysma[j]) {
				small = dysmb[j];
				if (small > dys[kj]) {
				    small = dys[kj];
				}
				dz = dz - dysma[j] + small;
			    } else {
				if (dys[kj] < dysma[j]) {
				    dz = dz - dysma[j] + dys[kj];
				}
			    }
			}
			if (dz < dzsky) {
			    dzsky = dz;
			    kbest = k;
			    nbest = i__;
			}
		    }
		}
	    }
	}
	if (dzsky < 0.) {
	    nrepr[kbest] = 1;
	    nrepr[nbest] = 0;
	    *sky += dzsky;
	    goto L60;
	}
    }
    obj[2] = *sky / *nn;
} /* bswap */

/* -----------------------------------------------------------
 cstat(): Compute STATistics (numerical output) concerning each partition
*/
void cstat(int *kk, int *nn, int *nsend, int *nrepr,
	   double *radus, double *damer, double *ttd,
	   double *separ, double *z__, double *s, int *hh,
	   double *dys, int *ncluv, int *nelem, int *med, int *nisol)
{
    /*logical*/int kand;
    int j, k, m, ja, jb, jk, jndz;
    int ksmal, numcl, mevj, njaj, nel, njm, nvn, ntt, nvna, numl, nplac;
    double aja, ajb, dam, dsmal, sep, rnn, rtt, ttt;

    /* Parameter adjustments */
    --nisol;
    --med;
    --nelem;
    --ncluv;
    --separ;
    --ttd;
    --damer;
    --radus;
    --nrepr;
    --nsend;
    --dys;

    /* Function Body */
    for (j = 1; j <= *nn; ++j) {
	if (nrepr[j] == 0) {
	    dsmal = *s * 1.1f + 1.;
	    for (k = 1; k <= *nn; ++k) {
		if (nrepr[k] == 1) {
		    njaj = F77_CALL(meet)(&k, &j);
		    if (dys[njaj] < dsmal) {
			dsmal = dys[njaj];
			ksmal = k;
		    }
		}
	    }
	    nsend[j] = ksmal;
	} else {
	    nsend[j] = j;
	}
    }
    jk = 1;
    nplac = nsend[1];
    for (j = 1; j <= *nn; ++j) {
	ncluv[j] = 0;
	if (nsend[j] == nplac) {
	    ncluv[j] = 1;
	}
    }
    for (ja = 2; ja <= *nn; ++ja) {
	nplac = nsend[ja];
	if (ncluv[nplac] == 0) {
	    ++jk;
	    for (j = 2; j <= *nn; ++j) {
		if (nsend[j] == nplac) {
		    ncluv[j] = jk;
		}
	    }
	    if (jk == *kk) {
		goto L148;
	    }
	}
    }

/*     analysis of the clustering. */

L148:
    for (numcl = 1; numcl <= *kk; ++numcl) {
	ntt = 0;
	radus[numcl] = -1.;
	ttt = 0.;
	for (j = 1; j <= *nn; ++j) {
	    if (ncluv[j] == numcl) {
		++ntt;
		m = nsend[j];
		nelem[ntt] = j;
		njm = F77_CALL(meet)(&j, &m);
		ttt += dys[njm];
		if (dys[njm] > radus[numcl]) {
		    radus[numcl] = dys[njm];
		}
	    }
	}
	rtt = (double) ntt;
	ttd[numcl] = ttt / rtt;
	med[numcl] = m;
    }
    rnn = (double) (*nn);
    if (*kk == 1) {
	damer[1] = *s;
	nrepr[1] = *nn;
	return;
    }
/*  ELSE	  kk > 1 :

     numl = number of l-clusters. */

    numl = 0;
    for (k = 1; k <= *kk; ++k) {
	/*
	  identification of cluster k:
	  nel  = number of objects
	  nelem= vector of objects */

	nel = 0;
	for (j = 1; j <= *nn; ++j) {
	    if (ncluv[j] == k) {
		++nel;
		nelem[nel] = j;
	    }
	}
	nrepr[k] = nel;
	if (nel == 1) {
	    nvn = nelem[1];
	    damer[k] = 0.;
	    separ[k] = *s * 1.1f + 1.;
	    for (j = 1; j <= *nn; ++j) {
		if (j != nvn) {
		    mevj = F77_CALL(meet)(&nvn, &j);
		    if (separ[k] > dys[mevj]) {
			separ[k] = dys[mevj];
		    }
		}
	    }

/* Is cluster k
	1) an L-cluster	 or
	2) an L*-cluster ? */
	    if (separ[k] == 0.) {
		++numl;
	    }
	} else {
/*	       nel != 1 : */
	    dam = -1.;
	    sep = *s * 1.1f + 1.;
	    kand = (1);
	    for (ja = 1; ja <= nel; ++ja) {
		nvna = nelem[ja];
		aja = -1.;
		ajb = *s * 1.1f + 1.;
		for (jb = 1; jb <= *nn; ++jb) {
		    jndz = F77_CALL(meet)(&nvna, &jb);
		    if (ncluv[jb] == k) {
			if (dys[jndz] > aja) {
			    aja = dys[jndz];
			}
		    } else {
			if (dys[jndz] < ajb) {
			    ajb = dys[jndz];
			}
		    }
		}
		if (kand && aja >= ajb) {
		    kand = (0);
		}
		if (dam < aja) {
		    dam = aja;
		}
		if (sep > ajb) {
		    sep = ajb;
		}
	    }
	    separ[k] = sep;
	    damer[k] = dam;
	    if (kand) {
		++numl;
		if (dam >= sep) {
/*		  L-cluster */
		    nisol[k] = 1;
		} else {
/*		  L*-cluster */
		    nisol[k] = 2;
		}
		goto L40;
	    }
	}
	nisol[k] = 0;
L40:
	;
    }
} /* cstat */

/* -----------------------------------------------------------
     Compute Silhouette Information :
 */
void dark(int *kk, int *nn, int *hh, int *ncluv,
	  int *nsend, int *nelem, int *negbr,
	  double *syl, double *srank, double *avsyl, double *ttsyl,
	  double *dys, double *s, double *sylinf)
{
    /* System generated locals */
    int sylinf_dim1, sylinf_offset;

    /* Local variables */
    int j, l, lang, lplac, nclu, nj, nl, nbb, mjl, njl, ntt, numcl, nsylr;
    double db, btt, rtt, dysa, dysb, symax;

/* Parameter adjustments */
    sylinf_dim1 = *nn;
    sylinf_offset = 1 + sylinf_dim1;
    sylinf -= sylinf_offset;
    --avsyl;
    --srank;
    --syl;
    --negbr;
    --nelem;
    --nsend;
    --ncluv;
    --dys;

    /* Function Body */
    nsylr = 0;
    *ttsyl = 0.;
    for (numcl = 1; numcl <= *kk; ++numcl) {
	ntt = 0;
	for (j = 1; j <= *nn; ++j) {
	    if (ncluv[j] == numcl) {
		++ntt;
		nelem[ntt] = j;
	    }
	}
	for (j = 1; j <= ntt; ++j) {
	    nj = nelem[j];
	    dysb = *s * 1.1f + 1.;
	    negbr[j] = -1;
	    for (nclu = 1; nclu <= *kk; ++nclu) {
		if (nclu != numcl) {
		    nbb = 0;
		    db = 0.;
		    for (l = 1; l <= *nn; ++l) {
			if (ncluv[l] == nclu) {
			    ++nbb;
			    if (l != nj) {
				mjl = F77_CALL(meet)(&nj, &l);
				db += dys[mjl];
			    }
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
		    if (nj != nl) {
			njl = F77_CALL(meet)(&nj, &nl);
			dysa += dys[njl];
		    }
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
			if (syl[j] <= -1.) {
			    syl[j] = -1.;
			}
			if (syl[j] >= 1.) {
			    syl[j] = 1.;
			}
		    } else {
			syl[j] = -1.;
		    }
		} else if (dysb > 0.) {
		    syl[j] = 1.;
		} else {
		    syl[j] = 0.;
		}
	    } else {
/*     ntt == 1: */
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
	rtt = (double) ntt;
	avsyl[numcl] /= rtt;
	if (ntt < 2) {
	    ++nsylr;
	    sylinf[nsylr + sylinf_dim1] = (double) numcl;
	    sylinf[nsylr + (sylinf_dim1 << 1)] = (double) negbr[1];
	    sylinf[nsylr + sylinf_dim1 * 3] = 0.;
	    sylinf[nsylr + (sylinf_dim1 << 2)] = (double) nelem[1];
	} else {
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
} /* dark */

