/* pam.f -- translated by f2c (version 20031025).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* $Id$ */

/* PAM := Partitioning Around Medoids */

/* Subroutine */ int pam_(integer *nn, integer *jpp, integer *kk, doublereal *
	x, doublereal *dys, integer *jdyss, doublereal *valmd, integer *jtmd,
	integer *ndyst, integer *nsend, integer *nrepr, integer *nelem,
	doublereal *radus, doublereal *damer, doublereal *ttd, doublereal *
	separ, doublereal *ttsyl, integer *med, doublereal *obj, integer *
	ncluv, doublereal *clusinf, doublereal *sylinf, integer *nisol)
{
    /* System generated locals */
    integer x_dim1, x_offset, clusinf_dim1, clusinf_offset, sylinf_dim1,
	    sylinf_offset, i__1;

    /* Local variables */
    static integer k, l;
    static doublereal s, sky;
    extern /* Subroutine */ int dark_(integer *, integer *, integer *,
	    integer *, integer *, integer *, integer *, doublereal *,
	    doublereal *, doublereal *, doublereal *, doublereal *,
	    doublereal *, doublereal *);
    static integer nhalf, jhalt;
    extern /* Subroutine */ int bswap_(integer *, integer *, integer *,
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *, doublereal *), cstat_(integer *,
	    integer *, integer *, integer *, doublereal *, doublereal *,
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, integer *, integer *, integer *, integer *),
	    dysta_(integer *, integer *, doublereal *, doublereal *, integer *
	    , integer *, doublereal *, integer *);


/*     carries out a clustering using the k-medoid approach. */


/* jdyss = 0 : compute distances from x */
/* 	= 1 : distances provided  in x */
    /* Parameter adjustments */
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
    i__1 = nhalf;
    for (l = 2; l <= i__1; ++l) {
	if (s < dys[l]) {
	    s = dys[l];
	}
/* L10: */
    }
/*     Build + Swap : */
    bswap_(kk, nn, &nrepr[1], &radus[1], &damer[1], &ttd[1], &nhalf, &dys[1],
	    &sky, &s, &obj[1]);
/*     Compute STATs : */
    cstat_(kk, nn, &nsend[1], &nrepr[1], &radus[1], &damer[1], &ttd[1], &
	    separ[1], &sky, &s, &nhalf, &dys[1], &ncluv[1], &nelem[1], &med[1]
	    , &nisol[1]);
    i__1 = *kk;
    for (k = 1; k <= i__1; ++k) {
	clusinf[k + clusinf_dim1] = (doublereal) nrepr[k];
	clusinf[k + (clusinf_dim1 << 1)] = radus[k];
	clusinf[k + clusinf_dim1 * 3] = ttd[k];
	clusinf[k + (clusinf_dim1 << 2)] = damer[k];
	clusinf[k + clusinf_dim1 * 5] = separ[k];
/* L135: */
    }
    if (1 < *kk && *kk < *nn) {
/* 	 Compute Silhouette info : */
	dark_(kk, nn, &nhalf, &ncluv[1], &nsend[1], &nelem[1], &nrepr[1], &
		radus[1], &damer[1], &ttd[1], ttsyl, &dys[1], &s, &sylinf[
		sylinf_offset]);
    }
    return 0;
} /* pam_ */

/*     ----------------------------------------------------------- */
/*     Compute Distances from X matrix {also for agnes() and diana()}: */

/* Subroutine */ int dysta_(integer *nn, integer *jpp, doublereal *x,
	doublereal *dys, integer *ndyst, integer *jtmd, doublereal *valmd,
	integer *jhalt)
{
    /* System generated locals */
    integer x_dim1, x_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer j, k, l;
    static doublereal pp, clk;
    static integer nlk, npres, lsubt;
    static doublereal rpres;

/* ndyst = 1 : euclidean */
/* "else"    : manhattan */
/* VARs */
    /* Parameter adjustments */
    --dys;
    --valmd;
    --jtmd;
    x_dim1 = *nn;
    x_offset = 1 + x_dim1;
    x -= x_offset;

    /* Function Body */
    nlk = 1;
    dys[1] = 0.f;
    pp = (doublereal) (*jpp);
    i__1 = *nn;
    for (l = 2; l <= i__1; ++l) {
	lsubt = l - 1;
	i__2 = lsubt;
	for (k = 1; k <= i__2; ++k) {
	    clk = 0.f;
	    ++nlk;
	    npres = 0;
	    i__3 = *jpp;
	    for (j = 1; j <= i__3; ++j) {
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
	    rpres = (doublereal) npres;
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
/* L20: */
	}
/* L100: */
    }
    return 0;
} /* dysta_ */

/*     ----------------------------------------------------------- */

/*     bswap(): the clustering algorithm in 2 parts:  I. build,	II. swap */

/* Subroutine */ int bswap_(integer *kk, integer *nn, integer *nrepr,
	doublereal *dysma, doublereal *dysmb, doublereal *beter, integer *hh,
	doublereal *dys, doublereal *sky, doublereal *s, doublereal *obj)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, ij, kj;
    static doublereal dz, cmd;
    static integer njn;
    extern integer meet_(integer *, integer *);
    static integer nmax;
    static doublereal ammax;
    static integer kbest;
    static doublereal small;
    static integer nbest;
    static doublereal dzsky;

/*     nrepr[]: here is boolean (0/1): 1 = "is representative object" */
/*     function called */

/*     -Wall: */
    /* Parameter adjustments */
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

    i__1 = *nn;
    for (i__ = 1; i__ <= i__1; ++i__) {
	nrepr[i__] = 0;
	dysma[i__] = *s * 1.1f + 1.f;
/* L17: */
    }
    i__1 = *kk;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *nn;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (nrepr[i__] == 0) {
		beter[i__] = 0.f;
		i__3 = *nn;
		for (j = 1; j <= i__3; ++j) {
		    ij = meet_(&i__, &j);
		    cmd = dysma[j] - dys[ij];
		    if (cmd > 0.f) {
			beter[i__] += cmd;
		    }
/* L21: */
		}
	    }
/* L22: */
	}
	ammax = 0.f;
	i__2 = *nn;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (nrepr[i__] == 0 && ammax <= beter[i__]) {
/* 		    does .lt. (instead of .le.) work too? -- NO! */
		ammax = beter[i__];
		nmax = i__;
	    }
/* L31: */
	}
	nrepr[nmax] = 1;
/* = .true. : *is* a representative */
	i__2 = *nn;
	for (j = 1; j <= i__2; ++j) {
	    njn = meet_(&nmax, &j);
	    if (dysma[j] > dys[njn]) {
		dysma[j] = dys[njn];
	    }
/* L41: */
	}
/* L20: */
    }
    *sky = 0.f;
    i__1 = *nn;
    for (j = 1; j <= i__1; ++j) {
	*sky += dysma[j];
/* L51: */
    }
    obj[1] = *sky / *nn;
    if (*kk > 1) {

/*     second algorithm: swap. */

/* --   Loop : */
L60:
	i__1 = *nn;
	for (j = 1; j <= i__1; ++j) {
	    dysma[j] = *s * 1.1f + 1.f;
	    dysmb[j] = *s * 1.1f + 1.f;
	    i__2 = *nn;
	    for (i__ = 1; i__ <= i__2; ++i__) {
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
/* L62: */
	    }
/* L63: */
	}
	dzsky = 1.f;
	i__1 = *nn;
	for (k = 1; k <= i__1; ++k) {
	    if (nrepr[k] == 0) {
		i__2 = *nn;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    if (nrepr[i__] == 1) {
			dz = 0.f;
			i__3 = *nn;
			for (j = 1; j <= i__3; ++j) {
			    ij = meet_(&i__, &j);
			    kj = meet_(&k, &j);
			    if (dys[ij] == dysma[j]) {
				small = dysmb[j];
				if (small > dys[kj]) {
				    small = dys[kj];
				}
				dz = dz - dysma[j] + small;
			    } else {
/* L70: */
				if (dys[kj] < dysma[j]) {
				    dz = dz - dysma[j] + dys[kj];
				}
			    }
/* L71: */
			}
			if (dz < dzsky) {
			    dzsky = dz;
			    kbest = k;
			    nbest = i__;
			}
		    }
/* L72: */
		}
	    }
/* L73: */
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
} /* bswap_ */

/*     ----------------------------------------------------------- */

/* cstat(): */
/* Compute STATistics (numerical output) concerning each partition */

/* Subroutine */ int cstat_(integer *kk, integer *nn, integer *nsend, integer
	*nrepr, doublereal *radus, doublereal *damer, doublereal *ttd,
	doublereal *separ, doublereal *z__, doublereal *s, integer *hh,
	doublereal *dys, integer *ncluv, integer *nelem, integer *med,
	integer *nisol)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer j, k, m, ja, jb, jk;
    static doublereal aja, ajb, dam;
    static integer nel, njm;
    static doublereal sep, rnn;
    static integer nvn, ntt;
    static doublereal rtt, ttt;
    static logical kand;
    static integer njaj;
    extern integer meet_(integer *, integer *);
    static integer mevj, nvna, jndz, numl, nplac;
    static doublereal dsmal;
    static integer ksmal, numcl;

/* function called */

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
    i__1 = *nn;
    for (j = 1; j <= i__1; ++j) {
	if (nrepr[j] == 0) {
	    dsmal = *s * 1.1f + 1.f;
	    i__2 = *nn;
	    for (k = 1; k <= i__2; ++k) {
		if (nrepr[k] == 1) {
		    njaj = meet_(&k, &j);
		    if (dys[njaj] < dsmal) {
			dsmal = dys[njaj];
			ksmal = k;
		    }
		}
/* L110: */
	    }
	    nsend[j] = ksmal;
	} else {
	    nsend[j] = j;
	}
/* L130: */
    }
    jk = 1;
    nplac = nsend[1];
    i__1 = *nn;
    for (j = 1; j <= i__1; ++j) {
	ncluv[j] = 0;
	if (nsend[j] == nplac) {
	    ncluv[j] = 1;
	}
/* L135: */
    }
    i__1 = *nn;
    for (ja = 2; ja <= i__1; ++ja) {
	nplac = nsend[ja];
	if (ncluv[nplac] == 0) {
	    ++jk;
	    i__2 = *nn;
	    for (j = 2; j <= i__2; ++j) {
		if (nsend[j] == nplac) {
		    ncluv[j] = jk;
		}
/* L140: */
	    }
	    if (jk == *kk) {
		goto L148;
	    }
	}
/* L145: */
    }

/*     analysis of the clustering. */

L148:
    i__1 = *kk;
    for (numcl = 1; numcl <= i__1; ++numcl) {
	ntt = 0;
	radus[numcl] = -1.f;
	ttt = 0.f;
	i__2 = *nn;
	for (j = 1; j <= i__2; ++j) {
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
/* L150: */
	}
	rtt = (doublereal) ntt;
	ttd[numcl] = ttt / rtt;
	med[numcl] = m;
/* L160: */
    }
/* L230: */
    rnn = (doublereal) (*nn);
    if (*kk == 1) {
	damer[1] = *s;
	nrepr[1] = *nn;
	return 0;
    }
/*  ELSE	  kk > 1 : */

/*     numl = number of l-clusters. */

/* L240: */
    numl = 0;
    i__1 = *kk;
    for (k = 1; k <= i__1; ++k) {

/*     identification of cluster k: */
/*     nel  = number of objects */
/*     nelem= vector of objects */

	nel = 0;
	i__2 = *nn;
	for (j = 1; j <= i__2; ++j) {
	    if (ncluv[j] == k) {
		++nel;
		nelem[nel] = j;
	    }
/* L23: */
	}
	nrepr[k] = nel;
	if (nel == 1) {
	    nvn = nelem[1];
	    damer[k] = 0.f;
	    separ[k] = *s * 1.1f + 1.f;
	    i__2 = *nn;
	    for (j = 1; j <= i__2; ++j) {
		if (j != nvn) {
		    mevj = meet_(&nvn, &j);
		    if (separ[k] > dys[mevj]) {
			separ[k] = dys[mevj];
		    }
		}
/* L250: */
	    }

/* Is cluster k */
/* 	1) an L-cluster	 or */
/* 	2) an L*-cluster ? */
	    if (separ[k] == 0.f) {
		++numl;
	    }
	} else {
/* 	       nel != 1 : */
	    dam = -1.f;
	    sep = *s * 1.1f + 1.f;
	    kand = TRUE_;
	    i__2 = nel;
	    for (ja = 1; ja <= i__2; ++ja) {
		nvna = nelem[ja];
		aja = -1.f;
		ajb = *s * 1.1f + 1.f;
		i__3 = *nn;
		for (jb = 1; jb <= i__3; ++jb) {
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
/* L25: */
		}
		if (kand && aja >= ajb) {
		    kand = FALSE_;
		}
		if (dam < aja) {
		    dam = aja;
		}
		if (sep > ajb) {
		    sep = ajb;
		}
/* L26: */
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
/* L27: */
		    nisol[k] = 2;
		}
		goto L40;
	    }
	}
	nisol[k] = 0;
L40:
	;
    }
/* L300: */
    return 0;
} /* cstat_ */

/*     ----------------------------------------------------------- */

/*     Compute Silhouette Information : */

/* Subroutine */ int dark_(integer *kk, integer *nn, integer *hh, integer *
	ncluv, integer *nsend, integer *nelem, integer *negbr, doublereal *
	syl, doublereal *srank, doublereal *avsyl, doublereal *ttsyl,
	doublereal *dys, doublereal *s, doublereal *sylinf)
{
    /* System generated locals */
    integer sylinf_dim1, sylinf_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer j, l;
    static doublereal db;
    static integer nj, nl, nbb, mjl, njl;
    static doublereal btt;
    static integer ntt;
    static doublereal rtt;
    static integer lang;
    extern integer meet_(integer *, integer *);
    static doublereal dysa;
    static integer nclu;
    static doublereal dysb;
    static integer lplac, numcl;
    static doublereal symax;
    static integer nsylr;

/*     function called */

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
    *ttsyl = 0.f;
    i__1 = *kk;
    for (numcl = 1; numcl <= i__1; ++numcl) {
	ntt = 0;
	i__2 = *nn;
	for (j = 1; j <= i__2; ++j) {
	    if (ncluv[j] == numcl) {
		++ntt;
		nelem[ntt] = j;
	    }
/* L30: */
	}
	i__2 = ntt;
	for (j = 1; j <= i__2; ++j) {
	    nj = nelem[j];
	    dysb = *s * 1.1f + 1.f;
	    negbr[j] = -1;
	    i__3 = *kk;
	    for (nclu = 1; nclu <= i__3; ++nclu) {
		if (nclu != numcl) {
		    nbb = 0;
		    db = 0.f;
		    i__4 = *nn;
		    for (l = 1; l <= i__4; ++l) {
			if (ncluv[l] == nclu) {
			    ++nbb;
			    if (l != nj) {
				mjl = meet_(&nj, &l);
				db += dys[mjl];
			    }
			}
/* L43: */
		    }
		    btt = (doublereal) nbb;
		    db /= btt;
		    if (dysb > db) {
			dysb = db;
			negbr[j] = nclu;
		    }
		}
/* L41: */
	    }
	    if (ntt > 1) {
		dysa = 0.f;
		i__3 = ntt;
		for (l = 1; l <= i__3; ++l) {
		    nl = nelem[l];
		    if (nj != nl) {
			njl = meet_(&nj, &nl);
			dysa += dys[njl];
		    }
/* L45: */
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
/* L40: */
	}
	avsyl[numcl] = 0.f;
	i__2 = ntt;
	for (j = 1; j <= i__2; ++j) {
	    symax = -2.f;
	    i__3 = ntt;
	    for (l = 1; l <= i__3; ++l) {
		if (symax < syl[l]) {
		    symax = syl[l];
		    lang = l;
		}
/* L70: */
	    }
	    nsend[j] = lang;
	    srank[j] = syl[lang];
	    avsyl[numcl] += srank[j];
	    syl[lang] = -3.f;
/* L60: */
	}
	*ttsyl += avsyl[numcl];
	rtt = (doublereal) ntt;
	avsyl[numcl] /= rtt;
	if (ntt < 2) {
	    ++nsylr;
	    sylinf[nsylr + sylinf_dim1] = (doublereal) numcl;
	    sylinf[nsylr + (sylinf_dim1 << 1)] = (doublereal) negbr[1];
	    sylinf[nsylr + sylinf_dim1 * 3] = 0.f;
	    sylinf[nsylr + (sylinf_dim1 << 2)] = (doublereal) nelem[1];
	} else {
	    i__2 = ntt;
	    for (l = 1; l <= i__2; ++l) {
		++nsylr;
		lplac = nsend[l];
		sylinf[nsylr + sylinf_dim1] = (doublereal) numcl;
		sylinf[nsylr + (sylinf_dim1 << 1)] = (doublereal) negbr[lplac]
			;
		sylinf[nsylr + sylinf_dim1 * 3] = srank[l];
		sylinf[nsylr + (sylinf_dim1 << 2)] = (doublereal) nelem[lplac]
			;
/* L80: */
	    }
	}
/* L100: */
    }
    *ttsyl /= *nn;
    return 0;
} /* dark_ */

