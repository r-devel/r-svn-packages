/* Produced by
 * $Id$
 */
/* pam.f -- translated by f2c (version 20031025).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip


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


/* $Id$

 PAM := Partitioning Around Medoids

 Subroutine */ int pam_(int *nn, int *jpp, int *kk, double *
	x, double *dys, int *jdyss, double *valmd, int *jtmd,

	int *ndyst, int *nsend, int *nrepr, int *nelem,

	double *radus, double *damer, double *ttd, double *
	separ, double *ttsyl, int *med, double *obj, int *
	ncluv, double *clusinf, double *sylinf, int *nisol)
{
    /* System generated locals */
    int x_dim1, x_offset, clusinf_dim1, clusinf_offset, sylinf_dim1,

	    sylinf_offset, i__1;

    /* Local variables */
    static int k, l;
    static double s, sky;
    extern /* Subroutine */ int dark_(int *, int *, int *,

	    int *, int *, int *, int *, double *,

	    double *, double *, double *, double *,

	    double *, double *);
    static int nhalf, jhalt;
    extern /* Subroutine */ int bswap_(int *, int *, int *,

	    double *, double *, double *, int *, double *,
	     double *, double *, double *), cstat_(int *,

	    int *, int *, int *, double *, double *,

	    double *, double *, double *, double *, int *,
	     double *, int *, int *, int *, int *),

	    dysta_(int *, int *, double *, double *, int *
	    , int *, double *, int *);


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
	dysta_(nn, jpp, &x[x_offset], &dys[1], ndyst, &jtmd[1], &valmd[1], &
		jhalt);
	if (jhalt != 0) {
	    *jdyss = -1;
	    return 0;
	}
    }
/*     nhalf := #{distances}+1 = length(dys) */
    nhalf = *nn * (*nn - 1) / 2 + 1;
/*     s := max( dys[.] ), the largest distance */
    s = 0.f;
    for (l = 2; l <= nhalf; ++l) { /* f2c-clean: s {i__1} {nhalf} */
	if (s < dys[l]) {
	    s = dys[l];
	}
    }
/*     Build + Swap : */
    bswap_(kk, nn, &nrepr[1], &radus[1], &damer[1], &ttd[1], &nhalf, &dys[1],

	    &sky, &s, &obj[1]);
/*     Compute STATs : */
    cstat_(kk, nn, &nsend[1], &nrepr[1], &radus[1], &damer[1], &ttd[1], &
	    separ[1], &sky, &s, &nhalf, &dys[1], &ncluv[1], &nelem[1], &med[1]
	    , &nisol[1]);
    for (k = 1; k <= *kk; ++k) { /* f2c-clean: s {i__1} {*kk} */
	clusinf[k + clusinf_dim1] = (double) nrepr[k];
	clusinf[k + (clusinf_dim1 << 1)] = radus[k];
	clusinf[k + clusinf_dim1 * 3] = ttd[k];
	clusinf[k + (clusinf_dim1 << 2)] = damer[k];
	clusinf[k + clusinf_dim1 * 5] = separ[k];
    }
    if (1 < *kk && *kk < *nn) {
/* 	 Compute Silhouette info : */
	dark_(kk, nn, &nhalf, &ncluv[1], &nsend[1], &nelem[1], &nrepr[1], &
		radus[1], &damer[1], &ttd[1], ttsyl, &dys[1], &s, &sylinf[
		sylinf_offset]);
    }
    return 0;
} /* pam_

     -----------------------------------------------------------
     Compute Distances from X matrix {also for agnes() and diana()}:

 Subroutine */ int dysta_(int *nn, int *jpp, double *x,

	double *dys, int *ndyst, int *jtmd, double *valmd,

	int *jhalt)
{
    /* System generated locals */
    int x_dim1, x_offset, i__1, i__2, i__3;
    double d__1;

    /* Builtin functions */
    double sqrt(double);

    /* Local variables */
    static int j, k, l;
    static double pp, clk;
    static int nlk, npres, lsubt;
    static double rpres;

/* ndyst = 1 : euclidean
 "else"    : manhattan
 VARs
     Parameter adjustments */
    --dys;
    --valmd;
    --jtmd;
    x_dim1 = *nn;
    x_offset = 1 + x_dim1;
    x -= x_offset;

    /* Function Body */
    nlk = 1;
    dys[1] = 0.f;
    pp = (double) (*jpp);
    for (l = 2; l <= *nn; ++l) { /* f2c-clean: s {i__1} {*nn} */
	lsubt = l - 1;
	for (k = 1; k <= lsubt; ++k) { /* f2c-clean: s {i__2} {lsubt} */
	    clk = 0.f;
	    ++nlk;
	    npres = 0;
	    for (j = 1; j <= *jpp; ++j) { /* f2c-clean: s {i__3} {*jpp} */
		if (jtmd[j] < 0) {
		    if (x[l + j * x_dim1] == valmd[j]) {
			goto L30;
		    }
		    if (x[k + j * x_dim1] == valmd[j]) {
			goto L30;
		    }
		}
		++npres;
		if (*ndyst == 1) {
		    clk += (x[l + j * x_dim1] - x[k + j * x_dim1]) * (x[l + j

			    * x_dim1] - x[k + j * x_dim1]);
		} else {
		    clk += (d__1 = x[l + j * x_dim1] - x[k + j * x_dim1], abs(
			    d__1));
		}
L30:
		;
	    }
	    rpres = (double) npres;
	    if (npres == 0) {
		*jhalt = 1;
		dys[nlk] = -1.f;
	    } else {
		if (*ndyst == 1) {
		    dys[nlk] = sqrt(clk * (pp / rpres));
		} else {
		    dys[nlk] = clk * (pp / rpres);
		}
	    }
	}
    }
    return 0;
} /* dysta_

     -----------------------------------------------------------

     bswap(): the clustering algorithm in 2 parts:  I. build,	II. swap

 Subroutine */ int bswap_(int *kk, int *nn, int *nrepr,

	double *dysma, double *dysmb, double *beter, int *hh,

	double *dys, double *sky, double *s, double *obj)
{
    /* System generated locals */
    int i__1, i__2, i__3;

    /* Local variables */
    static int i__, j, k, ij, kj;
    static double dz, cmd;
    static int njn;
    extern int meet_(int *, int *);
    static int nmax;
    static double ammax;
    static int kbest;
    static double small;
    static int nbest;
    static double dzsky;

/*     nrepr[]: here is boolean (0/1): 1 = "is representative object"
     function called

     -Wall:
     Parameter adjustments */
    --beter;
    --dysmb;
    --dysma;
    --nrepr;
    --dys;
    --obj;

    /* Function Body */
    nbest = -1;
    kbest = -1;

/*     first algorithm: build. */

    for (i__ = 1; i__ <= *nn; ++i__) { /* f2c-clean: s {i__1} {*nn} */
	nrepr[i__] = 0;
	dysma[i__] = *s * 1.1f + 1.f;
    }
    for (k = 1; k <= *kk; ++k) { /* f2c-clean: s {i__1} {*kk} */
	for (i__ = 1; i__ <= *nn; ++i__) { /* f2c-clean: s {i__2} {*nn} */
	    if (nrepr[i__] == 0) {
		beter[i__] = 0.f;
		for (j = 1; j <= *nn; ++j) { /* f2c-clean: s {i__3} {*nn} */
		    ij = meet_(&i__, &j);
		    cmd = dysma[j] - dys[ij];
		    if (cmd > 0.f) {
			beter[i__] += cmd;
		    }
		}
	    }
	}
	ammax = 0.f;
	for (i__ = 1; i__ <= *nn; ++i__) { /* f2c-clean: s {i__2} {*nn} */
	    if (nrepr[i__] == 0 && ammax <= beter[i__]) {
/* 		    does < (instead of <= ) work too? -- NO! */
		ammax = beter[i__];
		nmax = i__;
	    }
	}
	nrepr[nmax] = 1;
/* = .true. : *is* a representative */
	for (j = 1; j <= *nn; ++j) { /* f2c-clean: s {i__2} {*nn} */
	    njn = meet_(&nmax, &j);
	    if (dysma[j] > dys[njn]) {
		dysma[j] = dys[njn];
	    }
	}
    }
    *sky = 0.f;
    for (j = 1; j <= *nn; ++j) { /* f2c-clean: s {i__1} {*nn} */
	*sky += dysma[j];
    }
    obj[1] = *sky / *nn;
    if (*kk > 1) {

/*     second algorithm: swap.

 --   Loop : */
L60:
	for (j = 1; j <= *nn; ++j) { /* f2c-clean: s {i__1} {*nn} */
	    dysma[j] = *s * 1.1f + 1.f;
	    dysmb[j] = *s * 1.1f + 1.f;
	    for (i__ = 1; i__ <= *nn; ++i__) { /* f2c-clean: s {i__2} {*nn} */
		if (nrepr[i__] == 1) {
		    ij = meet_(&i__, &j);
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
	dzsky = 1.f;
	for (k = 1; k <= *nn; ++k) { /* f2c-clean: s {i__1} {*nn} */
	    if (nrepr[k] == 0) {
		for (i__ = 1; i__ <= *nn; ++i__) { /* f2c-clean: s {i__2} {*nn} */
		    if (nrepr[i__] == 1) {
			dz = 0.f;
			for (j = 1; j <= *nn; ++j) { /* f2c-clean: s {i__3} {*nn} */
			    ij = meet_(&i__, &j);
			    kj = meet_(&k, &j);
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
	if (dzsky < 0.f) {
	    nrepr[kbest] = 1;
	    nrepr[nbest] = 0;
	    *sky += dzsky;
	    goto L60;
	}
    }
    obj[2] = *sky / *nn;
    return 0;
} /* bswap_

     -----------------------------------------------------------

 cstat():
 Compute STATistics (numerical output) concerning each partition

 Subroutine */ int cstat_(int *kk, int *nn, int *nsend, int

	*nrepr, double *radus, double *damer, double *ttd,

	double *separ, double *z__, double *s, int *hh,

	double *dys, int *ncluv, int *nelem, int *med,

	int *nisol)
{
    /* System generated locals */
    int i__1, i__2, i__3;

    /* Local variables */
    static int j, k, m, ja, jb, jk;
    static double aja, ajb, dam;
    static int nel, njm;
    static double sep, rnn;
    static int nvn, ntt;
    static double rtt, ttt;
    static /*logical*/int kand;
    static int njaj;
    extern int meet_(int *, int *);
    static int mevj, nvna, jndz, numl, nplac;
    static double dsmal;
    static int ksmal, numcl;

/* function called

     Parameter adjustments */
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
    for (j = 1; j <= *nn; ++j) { /* f2c-clean: s {i__1} {*nn} */
	if (nrepr[j] == 0) {
	    dsmal = *s * 1.1f + 1.f;
	    for (k = 1; k <= *nn; ++k) { /* f2c-clean: s {i__2} {*nn} */
		if (nrepr[k] == 1) {
		    njaj = meet_(&k, &j);
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
    for (j = 1; j <= *nn; ++j) { /* f2c-clean: s {i__1} {*nn} */
	ncluv[j] = 0;
	if (nsend[j] == nplac) {
	    ncluv[j] = 1;
	}
    }
    for (ja = 2; ja <= *nn; ++ja) { /* f2c-clean: s {i__1} {*nn} */
	nplac = nsend[ja];
	if (ncluv[nplac] == 0) {
	    ++jk;
	    for (j = 2; j <= *nn; ++j) { /* f2c-clean: s {i__2} {*nn} */
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
    for (numcl = 1; numcl <= *kk; ++numcl) { /* f2c-clean: s {i__1} {*kk} */
	ntt = 0;
	radus[numcl] = -1.f;
	ttt = 0.f;
	for (j = 1; j <= *nn; ++j) { /* f2c-clean: s {i__2} {*nn} */
	    if (ncluv[j] == numcl) {
		++ntt;
		m = nsend[j];
		nelem[ntt] = j;
		njm = meet_(&j, &m);
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
	return 0;
    }
/*  ELSE	  kk > 1 :

     numl = number of l-clusters. */

    numl = 0;
    for (k = 1; k <= *kk; ++k) { /* f2c-clean: s {i__1} {*kk}

     identification of cluster k:
     nel  = number of objects
     nelem= vector of objects */

	nel = 0;
	for (j = 1; j <= *nn; ++j) { /* f2c-clean: s {i__2} {*nn} */
	    if (ncluv[j] == k) {
		++nel;
		nelem[nel] = j;
	    }
	}
	nrepr[k] = nel;
	if (nel == 1) {
	    nvn = nelem[1];
	    damer[k] = 0.f;
	    separ[k] = *s * 1.1f + 1.f;
	    for (j = 1; j <= *nn; ++j) { /* f2c-clean: s {i__2} {*nn} */
		if (j != nvn) {
		    mevj = meet_(&nvn, &j);
		    if (separ[k] > dys[mevj]) {
			separ[k] = dys[mevj];
		    }
		}
	    }

/* Is cluster k
 	1) an L-cluster	 or
 	2) an L*-cluster ? */
	    if (separ[k] == 0.f) {
		++numl;
	    }
	} else {
/* 	       nel != 1 : */
	    dam = -1.f;
	    sep = *s * 1.1f + 1.f;
	    kand = (1);
	    for (ja = 1; ja <= nel; ++ja) { /* f2c-clean: s {i__2} {nel} */
		nvna = nelem[ja];
		aja = -1.f;
		ajb = *s * 1.1f + 1.f;
		for (jb = 1; jb <= *nn; ++jb) { /* f2c-clean: s {i__3} {*nn} */
		    jndz = meet_(&nvna, &jb);
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
/* 		  L-cluster */
		    nisol[k] = 1;
		} else {
/* 		  L*-cluster */
		    nisol[k] = 2;
		}
		goto L40;
	    }
	}
	nisol[k] = 0;
L40:
	;
    }
    return 0;
} /* cstat_

     -----------------------------------------------------------

     Compute Silhouette Information :

 Subroutine */ int dark_(int *kk, int *nn, int *hh, int *
	ncluv, int *nsend, int *nelem, int *negbr, double *
	syl, double *srank, double *avsyl, double *ttsyl,

	double *dys, double *s, double *sylinf)
{
    /* System generated locals */
    int sylinf_dim1, sylinf_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static int j, l;
    static double db;
    static int nj, nl, nbb, mjl, njl;
    static double btt;
    static int ntt;
    static double rtt;
    static int lang;
    extern int meet_(int *, int *);
    static double dysa;
    static int nclu;
    static double dysb;
    static int lplac, numcl;
    static double symax;
    static int nsylr;

/*     function called

     Parameter adjustments */
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
    *ttsyl = 0.f;
    for (numcl = 1; numcl <= *kk; ++numcl) { /* f2c-clean: s {i__1} {*kk} */
	ntt = 0;
	for (j = 1; j <= *nn; ++j) { /* f2c-clean: s {i__2} {*nn} */
	    if (ncluv[j] == numcl) {
		++ntt;
		nelem[ntt] = j;
	    }
	}
	for (j = 1; j <= ntt; ++j) { /* f2c-clean: s {i__2} {ntt} */
	    nj = nelem[j];
	    dysb = *s * 1.1f + 1.f;
	    negbr[j] = -1;
	    for (nclu = 1; nclu <= *kk; ++nclu) { /* f2c-clean: s {i__3} {*kk} */
		if (nclu != numcl) {
		    nbb = 0;
		    db = 0.f;
		    for (l = 1; l <= *nn; ++l) { /* f2c-clean: s {i__4} {*nn} */
			if (ncluv[l] == nclu) {
			    ++nbb;
			    if (l != nj) {
				mjl = meet_(&nj, &l);
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
		dysa = 0.f;
		for (l = 1; l <= ntt; ++l) { /* f2c-clean: s {i__3} {ntt} */
		    nl = nelem[l];
		    if (nj != nl) {
			njl = meet_(&nj, &nl);
			dysa += dys[njl];
		    }
		}
		dysa /= ntt - 1;
		if (dysa > 0.f) {
		    if (dysb > 0.f) {
			if (dysb > dysa) {
			    syl[j] = 1.f - dysa / dysb;
			} else if (dysb < dysa) {
			    syl[j] = dysb / dysa - 1.f;
			} else {
/*     dysb == dysa: */
			    syl[j] = 0.f;
			}
			if (syl[j] <= -1.f) {
			    syl[j] = -1.f;
			}
			if (syl[j] >= 1.f) {
			    syl[j] = 1.f;
			}
		    } else {
			syl[j] = -1.f;
		    }
		} else if (dysb > 0.f) {
		    syl[j] = 1.f;
		} else {
		    syl[j] = 0.f;
		}
	    } else {
/*     ntt == 1: */
		syl[j] = 0.f;
	    }
	}
	avsyl[numcl] = 0.f;
	for (j = 1; j <= ntt; ++j) { /* f2c-clean: s {i__2} {ntt} */
	    symax = -2.f;
	    for (l = 1; l <= ntt; ++l) { /* f2c-clean: s {i__3} {ntt} */
		if (symax < syl[l]) {
		    symax = syl[l];
		    lang = l;
		}
	    }
	    nsend[j] = lang;
	    srank[j] = syl[lang];
	    avsyl[numcl] += srank[j];
	    syl[lang] = -3.f;
	}
	*ttsyl += avsyl[numcl];
	rtt = (double) ntt;
	avsyl[numcl] /= rtt;
	if (ntt < 2) {
	    ++nsylr;
	    sylinf[nsylr + sylinf_dim1] = (double) numcl;
	    sylinf[nsylr + (sylinf_dim1 << 1)] = (double) negbr[1];
	    sylinf[nsylr + sylinf_dim1 * 3] = 0.f;
	    sylinf[nsylr + (sylinf_dim1 << 2)] = (double) nelem[1];
	} else {
	    for (l = 1; l <= ntt; ++l) { /* f2c-clean: s {i__2} {ntt} */
		++nsylr;
		lplac = nsend[l];
		sylinf[nsylr + sylinf_dim1] = (double) numcl;
		sylinf[nsylr + (sylinf_dim1 << 1)] = (double) negbr[lplac]
			;
		sylinf[nsylr + sylinf_dim1 * 3] = srank[l];
		sylinf[nsylr + (sylinf_dim1 << 2)] = (double) nelem[lplac]
			;
	    }
	}
    }
    *ttsyl /= *nn;
    return 0;
} /* dark_ */

