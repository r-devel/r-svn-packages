/* Produced by
 * $Id$
 */
/* clara.f -- translated by f2c (version 20010821).
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


/* $Id$
     Clustering LARge Applications
     ~          ~~~   ~
     Clustering program based upon the k-medoid approach,
     and suitable for data sets of at least 100 objects.
     (for smaller data sets, please use program pam.)

 Subroutine */ int clara_(int *nn, int *jpp, int *kk,

	double *x, int *nran, int *nsam, double *dys, int

	*mdata, double *valmd, int *jtmd, int *ndyst, int *
	nrepr, int *nsel, int *nbest, int *nr, int *nrx,

	double *radus, double *ttd, double *ratt, double *
	ttbes, double *rdbes, double *rabes, int *mtt, double

	*obj, double *avsyl, double *ttsyl, double *sylinf,

	int *jstop, double *tmp1, double *tmp2, double *tmp3,

	int *ntmp1, int *ntmp2, int *ntmp3, int *ntmp4,

	int *ntmp5, int *ntmp6)
{
    /* System generated locals */
    int sylinf_dim1, sylinf_offset, i__1, i__2, i__3;

    /* Local variables */
    static int kall;
    static /*logical*/int nafs;
    static int nadv, jran, kran, kans, nneq, less, nsub, nrun, j, l;
    extern /* Subroutine */ int black_(int *, int *, int *,

	    int *, int *, double *, double *, double *,

	    double *, double *, double *, int *, int *,

	    int *, int *, double *, double *);
    static double s;
    static int nhalf;
    static double z__;
    extern /* Subroutine */ int selec_(int *, int *, int *,

	    int *, double *, int *, int *, int *,

	    double *, int *, int *, double *, double *,

	    int *, /*logical*/int *, double *, double *, double *,

	    int *, int *, int *, int *, int *, int *,

	    double *, double *);
    static int nsamb;
    extern /* Subroutine */ int randm_(int *, double *);
    static int jhalt, nadvp, nexap, nexbp, nunfs;
    extern /* Subroutine */ int resul_(int *, int *, int *,

	    int *, int *, int *, double *, double *,

	    int *, int *), bswap2_(int *, int *, int *,

	    double *, double *, double *, double *,

	    double *, double *), dysta2_(int *, int *,

	    int *, double *, int *, double *, int *,

	    int *, double *, int *);
    static int jk, jn, js;
    static double zb, sx;
    static int nad, jjb;
    static double zba;
    static int jkk;
    static double ran;
    static int kkm, kkp, jsm, nsm;
    static double rnn;
    static int ntt;

/*     nn   = number of objects
     jpp  = number of variables
     kk   = number of clusters
     nran = number of random samples drawn (= `samples' in S)
     nsam = number of objects drawn from data set (= `sampsize' in R)

     mdata= {0,1};  1 : min(x) is missing value (NA);  0: no NA
     ndyst= {1,2};  1 : euclidean;  2 : manhattan
     valmd[j] = "missing value code" (instead of NA) for x[,j]
     jtmd [j] = {-1,1}:  -1 : x[,j] has NA ;  1 : no NAs in x[,j]
 Var
     Parameter adjustments */
    --jtmd;
    --valmd;
    --x;
    --avsyl;
    --mtt;
    --rabes;
    --rdbes;
    --ttbes;
    --ratt;
    --ttd;
    --radus;
    --nrx;
    --nr;
    --ntmp6;
    --ntmp5;
    --ntmp4;
    --ntmp3;
    --ntmp2;
    --ntmp1;
    --tmp3;
    --tmp2;
    --tmp1;
    sylinf_dim1 = *nsam;
    sylinf_offset = 1 + sylinf_dim1 * 1;
    sylinf -= sylinf_offset;
    --nbest;
    --nsel;
    --nrepr;
    --dys;

    /* Function Body */
    *jstop = 0;
    rnn = (double) (*nn);
    if (*nn == *nsam) {
	nneq = 1;
    } else {
	nneq = 0;
    }
/*     nhalf := size of distance array dys() */
    nhalf = *nsam * (*nsam - 1) / 2 + 1;
    nsamb = *nsam << 1;
    if (*nn < nsamb) {
	less = *nn - *nsam;
    } else {
	less = *nsam;
    }
    nunfs = 0;
    kall = 0;
/*     nrun : this is the ``random seed'' of the very simple randm() below */
    nrun = 0;

/*     Loop :  random subsamples are drawn and partitioned into kk clusters */
    for (jran = 1; jran <= *nran; ++jran) { /* f2c-clean: s {i__1} {*nran} */
	jhalt = 0;
	if (nneq == 0) {
	    goto L140;
	}
	if (nneq == 2) {
	    goto L400;
	}
/*     else nneq == 1 when above nn == nsam : */
	nneq = 2;
	for (j = 1; j <= *nsam; ++j) { /* f2c-clean: s {i__2} {*nsam} */
	    nsel[j] = j;
	}
	goto L320;
/*     nneq = 0 : */
L140:
	ntt = 0;
	if (jran != 1 && nunfs != jran && *nn >= nsamb) {
	    for (jk = 1; jk <= *kk; ++jk) { /* f2c-clean: s {i__2} {*kk} */
		nsel[jk] = nrx[jk];
	    }
	    kkm = *kk - 1;
	    for (jk = 1; jk <= kkm; ++jk) { /* f2c-clean: s {i__2} {kkm} */
		nsm = nsel[jk];
		kkp = jk + 1;
		jsm = jk;
		for (jkk = kkp; jkk <= *kk; ++jkk) { /* f2c-clean: s {i__3} {*kk} */
		    if (nsel[jkk] >= nsm) {
			goto L160;
		    }
		    nsm = nsel[jkk];
		    jsm = jkk;
L160:
		    ;
		}
		nsel[jsm] = nsel[jk];
		nsel[jk] = nsm;
	    }
	    ntt = *kk;
	} else {
L180:
	    randm_(&nrun, &ran);
	    kran = (int) (rnn * ran + 1.f);
	    if (kran > *nn) {
		kran = *nn;
	    }
	    if (jran != 1) {
		for (jk = 1; jk <= *kk; ++jk) { /* f2c-clean: s {i__2} {*kk} */
		    if (kran == nrx[jk]) {
			goto L180;
		    }
		}
	    }
	    ++ntt;
	    nsel[ntt] = kran;
	    if (less == ntt) {
		goto L290;
	    }
	}
/*     Loop */
L210:
	randm_(&nrun, &ran);
	kran = (int) (rnn * ran + 1.f);
	if (kran > *nn) {
	    kran = *nn;
	}
	if (jran == 1 || *nn >= nsamb) {
	    goto L230;
	}
	for (jk = 1; jk <= *kk; ++jk) { /* f2c-clean: s {i__2} {*kk} */
	    if (kran == nrx[jk]) {
		goto L210;
	    }
	}
L230:
	for (kans = 1; kans <= ntt; ++kans) { /* f2c-clean: s {i__2} {ntt} */
	    if (nsel[kans] < kran) {
		goto L260;
	    }
	    if (nsel[kans] == kran) {
		goto L210;
	    }
	    goto L270;
L260:
	    ;
	}
	++ntt;
	nsel[ntt] = kran;
	goto L290;
L270:
	for (nad = kans; nad <= ntt; ++nad) { /* f2c-clean: s {i__2} {ntt} */
	    nadv = ntt - nad + kans;
	    nadvp = nadv + 1;
	    nsel[nadvp] = nsel[nadv];
	}
	++ntt;
	nsel[kans] = kran;
L290:
	if (ntt < less) {
	    goto L210;
	}
	if (*nn >= nsamb) {
	    goto L320;
	}
	nexap = 1;
	nexbp = 1;
/*     do 305  jn=1, nn */
	jn = 0;
L300:
	++jn;
	if (nsel[nexap] == jn) {
	    ++nexap;
	} else {
	    nrepr[nexbp] = jn;
	    ++nexbp;
	}
	if (jn < *nn) {
	    goto L300;
	}
/*     305  continue */
	for (nsub = 1; nsub <= *nsam; ++nsub) { /* f2c-clean: s {i__2} {*nsam} */
	    nsel[nsub] = nrepr[nsub];
	}
L320:
	dysta2_(nsam, jpp, &nsel[1], &x[1], nn, &dys[1], ndyst, &jtmd[1], &
		valmd[1], &jhalt);
	if (jhalt == 1) {
	    goto L400;
	}
	kall = 1;
	s = 0.f;
	l = 1;
L340:
	++l;
	if (dys[l] > s) {
	    s = dys[l];
	}
	if (l < nhalf) {
	    goto L340;
	}
	bswap2_(kk, nsam, &nrepr[1], &dys[1], &z__, &s, &tmp1[1], &tmp2[1], &
		tmp3[1]);
	selec_(kk, nn, jpp, ndyst, &zb, nsam, mdata, &jtmd[1], &valmd[1], &
		nrepr[1], &nsel[1], &dys[1], &x[1], &nr[1], &nafs, &ttd[1], &
		radus[1], &ratt[1], &ntmp1[1], &ntmp2[1], &ntmp3[1], &ntmp4[1]
		, &ntmp5[1], &ntmp6[1], &tmp1[1], &tmp2[1]);
	if (nafs) {
	    ++nunfs;
	    goto L400;
	}
	if (jran != 1) {
	    if (zb >= zba) {
		goto L400;
	    }
	}
	zba = zb;
	for (jjb = 1; jjb <= *kk; ++jjb) { /* f2c-clean: s {i__2} {*kk} */
	    ttbes[jjb] = ttd[jjb];
	    rdbes[jjb] = radus[jjb];
	    rabes[jjb] = ratt[jjb];
	}
	for (jk = 1; jk <= *kk; ++jk) { /* f2c-clean: s {i__2} {*kk} */
	    nrx[jk] = nr[jk];
	}
	for (js = 1; js <= *nsam; ++js) { /* f2c-clean: s {i__2} {*nsam} */
	    nbest[js] = nsel[js];
	}
	sx = s;
L400:
	;
    }
/* --- end random sampling loop */
    if (nunfs >= *nran) {
	*jstop = 1;
	return 0;
    }

/*     for the best subsample, the objects of the entire data set
     are assigned to their clusters */

    if (kall != 1) {
	*jstop = 2;
	return 0;
    }
    *obj = zba / rnn;
    dysta2_(nsam, jpp, &nbest[1], &x[1], nn, &dys[1], ndyst, &jtmd[1], &valmd[
	    1], &jhalt);
    resul_(kk, nn, jpp, ndyst, mdata, &jtmd[1], &valmd[1], &x[1], &nrx[1], &
	    mtt[1]);
    if (*kk > 1) {
	black_(kk, jpp, nn, nsam, &nbest[1], &dys[1], &sx, &x[1], &avsyl[1],

		ttsyl, &sylinf[sylinf_offset], &ntmp1[1], &ntmp2[1], &ntmp3[1]
		, &ntmp4[1], &tmp1[1], &tmp2[1]);
    }
    return 0;
} /* clara_

     --- of clara() ---------------------------------------------------


 Subroutine */ int dysta2_(int *nsam, int *jpp, int *nsel,

	double *x, int *nn, double *dys, int *ndyst, int *
	jtmd, double *valmd, int *jhalt)
{
    /* System generated locals */
    int i__1, i__2, i__3;
    double d__1;

    /* Builtin functions */
    double sqrt(double);

    /* Local variables */
    static int lsel, ksel, j, k, l, npres, lsubt;
    static double rpres;
    static int kj, lj;
    static double pp, clk;
    static int nlk;



    /* Parameter adjustments */
    --dys;
    --nsel;
    --valmd;
    --jtmd;
    --x;

    /* Function Body */
    pp = (double) (*jpp);
    nlk = 1;
    dys[1] = 0.f;
    for (l = 2; l <= *nsam; ++l) { /* f2c-clean: s {i__1} {*nsam} */
	lsubt = l - 1;
	lsel = nsel[l];
	for (k = 1; k <= lsubt; ++k) { /* f2c-clean: s {i__2} {lsubt} */
	    ksel = nsel[k];
	    clk = 0.f;
	    ++nlk;
	    npres = 0;
	    for (j = 1; j <= *jpp; ++j) { /* f2c-clean: s {i__3} {*jpp} */
		lj = (lsel - 1) * *jpp + j;
		kj = (ksel - 1) * *jpp + j;
		if (jtmd[j] < 0) {
/* in the following line, x(-1) ==> seg.fault {BDR to R-core, Sat, 3 Aug 2002} */
		    if (x[lj] == valmd[j]) {
			goto L30;
		    }
		    if (x[kj] == valmd[j]) {
			goto L30;
		    }
		}
		++npres;
		if (*ndyst == 1) {
		    clk += (x[lj] - x[kj]) * (x[lj] - x[kj]);
		} else {
		    clk += fabs(x[lj] - x[kj]);
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
} /* dysta2_

     -----------------------------------------------------------

 Subroutine */ int randm_(int *nrun, double *ran)
{
    static int k;
    static double ry;

/*   we programmed this generator ourselves because we wanted it
   to be machine independent. it should run on most computers
   because the largest int used is less than 2**30 . the period
   is 2**16=65536, which is good enough for our purposes. */
    *nrun = *nrun * 5761 + 999;
    k = *nrun / 65536;
    *nrun -= k << 16;
    ry = (double) (*nrun);
    *ran = ry / 65536.f;
    return 0;
} /* randm_

     -----------------------------------------------------------

 Subroutine */ int bswap2_(int *kk, int *nsam, int *nrepr,

	double *dys, double *sky, double *s, double *dysma,

	double *dysmb, double *beter)
{
    /* System generated locals */
    int i__1, i__2, i__3;

    /* Local variables */
    static int njaj;
    extern int meet_(int *, int *);
    static double rsam;
    static int nmax;
    static double asky;
    static int j, k;
    static double ammax;
    static int kbest;
    static double small;
    static int nbest;
    static double dzsky;
    static int ja;
    static double dz, cmd;
    static int nkj, njn, nny;


/* Var

     first algorithm: build.

     Parameter adjustments */
    --beter;
    --dysmb;
    --dysma;
    --dys;
    --nrepr;

    /* Function Body */
    nny = 0;
    for (j = 1; j <= *nsam; ++j) { /* f2c-clean: s {i__1} {*nsam} */
	nrepr[j] = 0;
	dysma[j] = *s * 1.1f + 1.f;
    }
/* -- LOOP --------------------- */
L20:
    for (ja = 1; ja <= *nsam; ++ja) { /* f2c-clean: s {i__1} {*nsam} */
	if (nrepr[ja] != 0) {
	    goto L22;
	}
	beter[ja] = 0.f;
	for (j = 1; j <= *nsam; ++j) { /* f2c-clean: s {i__2} {*nsam} */
	    njaj = meet_(&ja, &j);
	    cmd = dysma[j] - dys[njaj];
	    if (cmd > 0.f) {
		beter[ja] += cmd;
	    }
	}
L22:
	;
    }
    ammax = 0.f;
    for (ja = 1; ja <= *nsam; ++ja) { /* f2c-clean: s {i__1} {*nsam} */
	if (nrepr[ja] != 0) {
	    goto L31;
	}
	if (beter[ja] < ammax) {
	    goto L31;
	}
	ammax = beter[ja];
	nmax = ja;
L31:
	;
    }
    nrepr[nmax] = 1;
    ++nny;
    for (j = 1; j <= *nsam; ++j) { /* f2c-clean: s {i__1} {*nsam} */
	njn = meet_(&nmax, &j);
	if (dys[njn] < dysma[j]) {
	    dysma[j] = dys[njn];
	}
    }
    if (nny != *kk) {
	goto L20;
    }
/* -- */
    *sky = 0.f;
    for (j = 1; j <= *nsam; ++j) { /* f2c-clean: s {i__1} {*nsam} */
	*sky += dysma[j];
    }
    if (*kk == 1) {
	return 0;
    }
    rsam = (double) (*nsam);
    asky = *sky / rsam;

/*     second algorithm: swap.

 -- LOOP */
L60:
    for (j = 1; j <= *nsam; ++j) { /* f2c-clean: s {i__1} {*nsam} */
	dysma[j] = *s * 1.1f + 1.f;
	dysmb[j] = *s * 1.1f + 1.f;
	for (ja = 1; ja <= *nsam; ++ja) { /* f2c-clean: s {i__2} {*nsam} */
	    if (nrepr[ja] == 0) {
		goto L62;
	    }
	    njaj = meet_(&ja, &j);
	    if (dys[njaj] >= dysma[j]) {
		goto L61;
	    }
	    dysmb[j] = dysma[j];
	    dysma[j] = dys[njaj];
	    goto L62;
L61:
	    if (dys[njaj] >= dysmb[j]) {
		goto L62;
	    }
	    dysmb[j] = dys[njaj];
L62:
	    ;
	}
    }
    dzsky = 1.f;
    for (k = 1; k <= *nsam; ++k) { /* f2c-clean: s {i__1} {*nsam} */
	if (nrepr[k] == 1) {
	    goto L73;
	}
	for (ja = 1; ja <= *nsam; ++ja) { /* f2c-clean: s {i__2} {*nsam} */
	    if (nrepr[ja] == 0) {
		goto L72;
	    }
	    dz = 0.f;
	    for (j = 1; j <= *nsam; ++j) { /* f2c-clean: s {i__3} {*nsam} */
		njaj = meet_(&ja, &j);
		nkj = meet_(&k, &j);
		if (dys[njaj] != dysma[j]) {
		    goto L70;
		}
		small = dysmb[j];
		if (dys[njaj] < small) {
		    small = dys[nkj];
		}
		dz = dz - dysma[j] + small;
		goto L71;
L70:
		if (dys[nkj] < dysma[j]) {
		    dz = dz - dysma[j] + dys[nkj];
		}
L71:
		;
	    }
	    if (dz >= dzsky) {
		goto L72;
	    }
	    dzsky = dz;
	    kbest = k;
	    nbest = ja;
L72:
	    ;
	}
L73:
	;
    }
    if (dzsky >= 0.f) {
	return 0;
    }
    nrepr[kbest] = 1;
    nrepr[nbest] = 0;
    *sky += dzsky;
    goto L60;
} /* bswap2_

     --- of bswap2() --------------------------------------------------
     selec() : called once [per random sample] from clara()

 Subroutine */ int selec_(int *kk, int *nn, int *jpp, int *
	ndyst, double *zb, int *nsam, int *mdata, int *jtmd,

	double *valmd, int *nrepr, int *nsel, double *dys,

	double *x, int *nr, /*logical*/int *nafs, double *ttd,

	double *radus, double *ratt, int *nrnew, int *nsnew,

	int *npnew, int *ns, int *np, int *new__, double *
	ttnew, double *rdnew)
{
    /* System generated locals */
    int i__1, i__2;
    double d__1;

    /* Builtin functions */
    double sqrt(double);

    /* Local variables */
    static int npab;
    extern int meet_(int *, int *);
    static int newf, jnew, nrjk;
    static double dsum, pres;
    static int j, jkabc;
    static double dnull;
    static int nstrt, ka, kb, na, nb, jk, jn, jp;
    static double pp, abc;
    static int npa, npb, njk;
    static double tra, rns;


/* Var

     nafs = .true. if a distance cannot be calculated
     Parameter adjustments */
    --ratt;
    --radus;
    --ttd;
    --nr;
    --x;
    --valmd;
    --jtmd;
    --rdnew;
    --ttnew;
    --new__;
    --np;
    --ns;
    --npnew;
    --nsnew;
    --nrnew;
    --dys;
    --nsel;
    --nrepr;

    /* Function Body */
    *nafs = (0);

/*     identification of representative objects, and initializations */

    jk = 0;
    for (j = 1; j <= *nsam; ++j) { /* f2c-clean: s {i__1} {*nsam} */
	if (nrepr[j] == 0) {
	    goto L10;
	}
	++jk;
	nr[jk] = nsel[j];
	ns[jk] = 0;
	ttd[jk] = 0.f;
	radus[jk] = -1.f;
	np[jk] = j;
L10:
	;
    }

/*     assignment of the objects of the entire data set to a cluster,
     computation of some statistics, determination of the
     new ordering of the clusters */

    *zb = 0.f;
    pp = (double) (*jpp);
    newf = 0;
    jn = 0;
L15:
    ++jn;
    if (*mdata == 0) {
	for (jk = 1; jk <= *kk; ++jk) { /* f2c-clean: s {i__1} {*kk} */
	    dsum = 0.f;
	    nrjk = nr[jk];
	    if (nrjk != jn) {
		for (jp = 1; jp <= *jpp; ++jp) { /* f2c-clean: s {i__2} {*jpp} */
		    na = (nrjk - 1) * *jpp + jp;
		    nb = (jn - 1) * *jpp + jp;
		    tra = fabs(x[na] - x[nb]);
		    if (*ndyst == 1) {
			tra *= tra;
		    }
		    dsum += tra;
		}
		if (jk != 1) {
		    if (dsum >= dnull) {
			goto L30;
		    }
		}
	    }
	    dnull = dsum;
	    jkabc = jk;
L30:
	    ;
	}
    } else {
	pres = 0.f;
	for (jk = 1; jk <= *kk; ++jk) { /* f2c-clean: s {i__1} {*kk} */
	    dsum = 0.f;
	    nrjk = nr[jk];
	    if (nrjk != jn) {
		abc = 0.f;
		for (jp = 1; jp <= *jpp; ++jp) { /* f2c-clean: s {i__2} {*jpp} */
		    na = (nrjk - 1) * *jpp + jp;
		    nb = (jn - 1) * *jpp + jp;
		    if (jtmd[jp] < 0) {
			if (x[na] == valmd[jp]) {
			    goto L50;
			}
			if (x[nb] == valmd[jp]) {
			    goto L50;
			}
		    }
		    abc += 1.f;
		    tra = fabs(x[na] - x[nb]);
		    if (*ndyst == 1) {
			tra *= tra;
		    }
		    dsum += tra;
L50:
		    ;
		}
		if (abc < .5f) {
		    goto L70;
		}
		dsum = dsum * abc / pp;
	    }
	    if (pres <= .5f) {
		pres = 1.f;
	    } else {
		if (dsum >= dnull) {
		    goto L70;
		}
	    }
	    dnull = dsum;
	    jkabc = jk;
L70:
	    ;
	}
	if (pres > .5f) {
	    goto L80;
	}
	*nafs = (1);
	return 0;
    }
L80:
    if (*ndyst == 1) {
	dnull = sqrt(dnull);
    }
    *zb += dnull;
    ttd[jkabc] += dnull;
    if (dnull > radus[jkabc]) {
	radus[jkabc] = dnull;
    }
    ++ns[jkabc];
    if (newf < *kk) {
	if (newf != 0) {
	    for (jnew = 1; jnew <= newf; ++jnew) { /* f2c-clean: s {i__1} {newf} */
		if (jkabc == new__[jnew]) {
		    goto L90;
		}
	    }
	}
	++newf;
	new__[newf] = jkabc;
    }
L90:
    if (jn < *nn) {
	goto L15;
    }

/*     a permutation is carried out on vectors nr,ns,np,ttd,radus
     using the information in vector new. */

    for (jk = 1; jk <= *kk; ++jk) { /* f2c-clean: s {i__1} {*kk} */
	njk = new__[jk];
	nrnew[jk] = nr[njk];
	nsnew[jk] = ns[njk];
	npnew[jk] = np[njk];
	ttnew[jk] = ttd[njk];
	rdnew[jk] = radus[njk];
    }
    for (jk = 1; jk <= *kk; ++jk) { /* f2c-clean: s {i__1} {*kk} */
	nr[jk] = nrnew[jk];
	ns[jk] = nsnew[jk];
	np[jk] = npnew[jk];
	ttd[jk] = ttnew[jk];
	radus[jk] = rdnew[jk];
    }
    for (j = 1; j <= *kk; ++j) { /* f2c-clean: s {i__1} {*kk} */
	rns = (double) ns[j];
	ttd[j] /= rns;
    }
    if (*kk == 1) {
	goto L150;
    }

/*     computation of minimal distance of medoid ka to any
     other medoid for comparison with the radius of cluster ka. */

    for (ka = 1; ka <= *kk; ++ka) { /* f2c-clean: s {i__1} {*kk} */
	nstrt = 0;
	npa = np[ka];
	for (kb = 1; kb <= *kk; ++kb) { /* f2c-clean: s {i__2} {*kk} */
	    if (kb == ka) {
		goto L110;
	    }
	    npb = np[kb];
	    npab = meet_(&npa, &npb);
	    if (nstrt == 0) {
		nstrt = 1;
	    } else {
		if (dys[npab] >= ratt[ka]) {
		    goto L110;
		}
	    }
	    ratt[ka] = dys[npab];
	    if (ratt[ka] != 0.f) {
		goto L110;
	    }
	    ratt[ka] = -1.f;
L110:
	    ;
	}
	if (ratt[ka] > -.5f) {
	    ratt[ka] = radus[ka] / ratt[ka];
	}
    }
L150:
    return 0;
} /* selec_

     -----------------------------------------------------------

 Subroutine */ int resul_(int *kk, int *nn, int *jpp, int *
	ndyst, int *mdata, int *jtmd, double *valmd, double *
	x, int *nrx, int *mtt)
{
    /* System generated locals */
    int i__1, i__2;
    double d__1;

    /* Builtin functions */
    double sqrt(double);

    /* Local variables */
    static int njnb, nxja, nrjk;
    static double dsum;
    static int j, nrjka;
    static double dnull;
    static int jksky, ja, ka, na, nb, jk, jn;
    static double pp, abc;
    static int jna;
    static double tra;

/* Var
     Parameter adjustments */
    --mtt;
    --nrx;
    --x;
    --valmd;
    --jtmd;

    /* Function Body */
    pp = (double) (*jpp);

/*     clustering vector is incorporated into x, and printed. */

    jn = 0;
L100:
    ++jn;
    njnb = (jn - 1) * *jpp;
    for (jk = 1; jk <= *kk; ++jk) { /* f2c-clean: s {i__1} {*kk} */
	if (nrx[jk] == jn) {
	    goto L220;
	}
    }
    jna = (jn - 1) * *jpp + 1;
    if (*mdata != 0) {
	goto L170;
    }
    for (jk = 1; jk <= *kk; ++jk) { /* f2c-clean: s {i__1} {*kk} */
	dsum = 0.f;
	nrjk = (nrx[jk] - 1) * *jpp;
	for (j = 1; j <= *jpp; ++j) { /* f2c-clean: s {i__2} {*jpp} */
	    na = nrjk + j;
	    nb = njnb + j;
	    tra = fabs(x[na] - x[nb]);
	    if (*ndyst == 1) {
		tra *= tra;
	    }
	    dsum += tra;
	}
	if (*ndyst == 1) {
	    dsum = sqrt(dsum);
	}
	if (jk == 1) {
	    dnull = dsum + .1f;
	}
	if (dsum >= dnull) {
	    goto L160;
	}
	dnull = dsum;
	jksky = jk;
L160:
	;
    }
    goto L200;
L170:
    for (jk = 1; jk <= *kk; ++jk) { /* f2c-clean: s {i__1} {*kk} */
	dsum = 0.f;
	nrjk = (nrx[jk] - 1) * *jpp;
	abc = 0.f;
	for (j = 1; j <= *jpp; ++j) { /* f2c-clean: s {i__2} {*jpp} */
	    na = nrjk + j;
	    nb = njnb + j;
	    if (jtmd[j] >= 0) {
		goto L185;
	    }
	    if (x[na] == valmd[j]) {
		goto L180;
	    }
	    if (x[nb] == valmd[j]) {
		goto L180;
	    }
L185:
	    abc += 1.f;
	    tra = fabs(x[na] - x[nb]);
	    if (*ndyst == 1) {
		tra *= tra;
	    }
	    dsum += tra;
L180:
	    ;
	}
	if (*ndyst == 1) {
	    dsum = sqrt(dsum);
	}
	dsum = dsum * abc / pp;
	if (jk == 1) {
	    dnull = dsum + .1f;
	}
	if (dsum >= dnull) {
	    goto L190;
	}
	dnull = dsum;
	jksky = jk;
L190:
	;
    }
L200:
    x[jna] = (double) jksky;
L220:
    if (jn < *nn) {
	goto L100;
    }
    for (jk = 1; jk <= *kk; ++jk) { /* f2c-clean: s {i__1} {*kk} */
	nrjk = nrx[jk];
	nrjka = (nrjk - 1) * *jpp + 1;
	x[nrjka] = (double) jk;
    }
    for (ka = 1; ka <= *kk; ++ka) { /* f2c-clean: s {i__1} {*kk} */
	mtt[ka] = 0;
	j = 0;
L325:
	++j;
	ja = (j - 1) * *jpp + 1;
	nxja = (int) (x[ja] + .1f);
	if (nxja == ka) {
	    ++mtt[ka];
	}
	if (j < *nn) {
	    goto L325;
	}
    }
    return 0;
} /* resul_

     -----------------------------------------------------------

 Subroutine */ int black_(int *kk, int *jpp, int *nn, int *
	nsam, int *nbest, double *dys, double *sx, double *x,

	double *avsyl, double *ttsyl, double *sylinf, int *
	ncluv, int *nsend, int *nelem, int *negbr, double *
	syl, double *srank)
{
    /* System generated locals */
    int sylinf_dim1, sylinf_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static int lang;
    extern int meet_(int *, int *);
    static double dysa;
    static int nclu;
    static double dysb, rsam;
    static int j, l, ncase, lplac, numcl;
    static double symax;
    static int nsylr;
    static double db;
    static int nj, nl, nbb, jna;
    static double att, btt;
    static int ntt;
    static double rtt;


/* Silhouettes computation and "drawing"  --> syl() and sylinf()


 Var

     construction of clustering vector (ncluv)
     of selected sample (nbest).

     Parameter adjustments */
    --avsyl;
    --x;
    --srank;
    --syl;
    --negbr;
    --nelem;
    --nsend;
    --ncluv;
    sylinf_dim1 = *nsam;
    sylinf_offset = 1 + sylinf_dim1 * 1;
    sylinf -= sylinf_offset;
    --dys;
    --nbest;

    /* Function Body */
    for (l = 1; l <= *nsam; ++l) { /* f2c-clean: s {i__1} {*nsam} */
	ncase = nbest[l];
	jna = (ncase - 1) * *jpp + 1;
	ncluv[l] = (int) (x[jna] + .1f);
    }

/*     drawing of the silhouettes */

    nsylr = 0;
    *ttsyl = 0.f;
    for (numcl = 1; numcl <= *kk; ++numcl) { /* f2c-clean: s {i__1} {*kk} */
	ntt = 0;
	for (j = 1; j <= *nsam; ++j) { /* f2c-clean: s {i__2} {*nsam} */
	    if (ncluv[j] != numcl) {
		goto L30;
	    }
	    ++ntt;
	    nelem[ntt] = j;
L30:
	    ;
	}
	for (j = 1; j <= ntt; ++j) { /* f2c-clean: s {i__2} {ntt} */
	    nj = nelem[j];
	    dysb = *sx * 1.1f + 1.f;
	    negbr[j] = -1;
	    for (nclu = 1; nclu <= *kk; ++nclu) { /* f2c-clean: s {i__3} {*kk} */
		if (nclu == numcl) {
		    goto L41;
		}
		nbb = 0;
		db = 0.f;
		for (l = 1; l <= *nsam; ++l) { /* f2c-clean: s {i__4} {*nsam} */
		    if (ncluv[l] != nclu) {
			goto L43;
		    }
		    ++nbb;
		    db += dys[meet_(&nj, &l)];
L43:
		    ;
		}
		btt = (double) nbb;
		db /= btt;
		if (db >= dysb) {
		    goto L41;
		}
		dysb = db;
		negbr[j] = nclu;
L41:
		;
	    }
	    if (ntt == 1) {
		goto L50;
	    }
	    dysa = 0.f;
	    for (l = 1; l <= ntt; ++l) { /* f2c-clean: s {i__3} {ntt} */
		nl = nelem[l];
		dysa += dys[meet_(&nj, &nl)];
	    }
	    att = (double) (ntt - 1);
	    dysa /= att;
	    if (dysa > 0.f) {
		goto L51;
	    }
	    if (dysb > 0.f) {
		goto L52;
	    }
L50:
	    syl[j] = 0.f;
	    goto L40;
L52:
	    syl[j] = 1.f;
	    goto L40;
L51:
	    if (dysb > 0.f) {
		if (dysb > dysa) {
		    syl[j] = 1.f - dysa / dysb;
		}
		if (dysb < dysa) {
		    syl[j] = dysb / dysa - 1.f;
		}
		if (dysb == dysa) {
		    syl[j] = 0.f;
		}
	    } else {
		syl[j] = -1.f;
	    }
	    if (syl[j] < -1.f) {
		syl[j] = -1.f;
	    }
	    if (syl[j] > 1.f) {
		syl[j] = 1.f;
	    }
L40:
	    ;
	}
	avsyl[numcl] = 0.f;
	for (j = 1; j <= ntt; ++j) { /* f2c-clean: s {i__2} {ntt} */
	    symax = -2.f;
	    for (l = 1; l <= ntt; ++l) { /* f2c-clean: s {i__3} {ntt} */
		if (syl[l] > symax) {
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
	if (ntt >= 2) {
	    goto L75;
	}
	ncase = nelem[1];
	++nsylr;
	sylinf[nsylr + sylinf_dim1] = (double) numcl;
	sylinf[nsylr + (sylinf_dim1 << 1)] = (double) negbr[1];
	sylinf[nsylr + sylinf_dim1 * 3] = 0.f;
	sylinf[nsylr + (sylinf_dim1 << 2)] = (double) nbest[ncase];
	goto L100;
L75:
	for (l = 1; l <= ntt; ++l) { /* f2c-clean: s {i__2} {ntt} */
	    lplac = nsend[l];
	    ncase = nelem[lplac];
	    ++nsylr;
	    sylinf[nsylr + sylinf_dim1] = (double) numcl;
	    sylinf[nsylr + (sylinf_dim1 << 1)] = (double) negbr[lplac];
	    sylinf[nsylr + sylinf_dim1 * 3] = srank[l];
	    sylinf[nsylr + (sylinf_dim1 << 2)] = (double) nbest[ncase];
	}
L100:
	;
    }
    rsam = (double) (*nsam);
    *ttsyl /= rsam;
    return 0;
} /* black_ */

