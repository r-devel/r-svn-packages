/* clara.f -- translated by f2c (version 20010821).
 * and run through  f2c-clean,v 1.10 2002/03/28 16:37:27 maechler
 */

#include <math.h>
#include "cluster.h"

/* $Id$
     Clustering LARge Applications
     ~		~~~   ~
     Clustering program based upon the k-medoid approach,
     and suitable for data sets of at least 100 objects.
     (for smaller data sets, please use program pam.)
*/
void clara_(int *nn, /* = number of objects */
	    int *jpp,/* = number of variables */
	    int *kk, /* = number of clusters */
	    double *x,	/* */
	    int *nran,	/* = #{random samples} drawn (= `samples' in R)*/
	    int *nsam,	/* = #{objects} drawn from data set (`sampsize' in R) */
	    double *dys,/* */
	    int *mdata,/*= {0,1}; 1: min(x) is missing value (NA);  0: no NA */
	    double *valmd,/*[j]= missing value code (instead of NA) for x[,j]*/
	    int *jtmd,/*[j] = {-1,1};  -1: x[,j] has NA;  1: no NAs in x[,j] */
	    int *ndyst, /* = {1,2};  1 : euclidean;  2 : manhattan*/
	    int *nrepr,	/* */
	    int *nsel, int *nbest, int *nr, int *nrx,
	    double *radus, double *ttd, double *ratt,
	    double *ttbes, double *rdbes, double *rabes,
	    int *mtt, double *obj, double *avsyl, double *ttsyl,
	    double *sylinf, int *jstop,
	    double *tmp1, double *tmp2, double *tmp3,
	    int *ntmp1, int *ntmp2, int *ntmp3, int *ntmp4,
	    int *ntmp5, int *ntmp6)
{

    /* Local variables */

    /*logical*/int nafs;
    int j, jjb, jk, jkk, jn, js, jsm, jhalt;
    int kkm, kkp, nsm, ntt;
    int kall, nadv, jran, kran, kans, nneq, less, nsub, nrun, l;
    int nhalf, nsamb, nad, nadvp, nexap, nexbp, nunfs;
    double s, z__, zb, zba, ran, rnn, sx;

/* Parameter adjustments */
    --rabes; --rdbes; --ttbes;
    --ratt;  --radus; --ttd;

    --nrx; --nr;

    --nbest; --nsel;

    --nrepr; --dys;


    *jstop = 0;
    rnn = (double) (*nn);
    /* nhalf := size of distance array dys[] */
    nhalf = *nsam * (*nsam - 1) / 2 + 1;

    if (*nn == *nsam)
	nneq = 1;
    else
	nneq = 0;

    nsamb = *nsam << 1;
    if (*nn < nsamb) {
	less = *nn - *nsam;
    } else {
	less = *nsam;
    }
    nunfs = 0;
    kall = 0;
    /* nrun : this is the ``random seed'' of the very simple randm() below */
    nrun = 0;

/* __LOOP__ :  random subsamples are drawn and partitioned into kk clusters */

    for (jran = 1; jran <= *nran; ++jran) {
	jhalt = 0;
	if (nneq == 0) {
	    goto L140;
	}
	if (nneq == 2)  continue;/* random sample*/

	/* else nneq == 1 when above nn == nsam : */
	nneq = 2;
	for (j = 1; j <= *nsam; ++j) {
	    nsel[j] = j;
	}
	goto L320;

/*     nneq = 0 : */
L140:
	ntt = 0;
	if (jran != 1 && nunfs != jran && *nn >= nsamb) {
	    for (jk = 1; jk <= *kk; ++jk) {
		nsel[jk] = nrx[jk];
	    }
	    kkm = *kk - 1;
	    for (jk = 1; jk <= kkm; ++jk) {
		nsm = nsel[jk];
		kkp = jk + 1;
		jsm = jk;
		for (jkk = kkp; jkk <= *kk; ++jkk) {
		    if (nsel[jkk] < nsm) {
			nsm = nsel[jkk];
			jsm = jkk;
		    }
		}
		nsel[jsm] = nsel[jk];
		nsel[jk] = nsm;
	    }
	    ntt = *kk;
	}
	else {
L180:
	    randm(&nrun, &ran);
	    kran = (int) (rnn * ran + 1.);
	    if (kran > *nn) {
		kran = *nn;
	    }
	    if (jran != 1) {
		for (jk = 1; jk <= *kk; ++jk) {
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

/* Loop -- */
L210:
	randm(&nrun, &ran);
	kran = (int) (rnn * ran + 1.);
	if (kran > *nn) {
	    kran = *nn;
	}
	if (jran != 1 && *nn < nsamb) {
	    for (jk = 1; jk <= *kk; ++jk) {
		if (kran == nrx[jk])
		    goto L210;
	    }
	}

	for (kans = 1; kans <= ntt; ++kans) {
	    if (nsel[kans] >= kran) {
		if (nsel[kans] == kran)
		    goto L210;
		else
		    goto L270;
	    }
	}
	++ntt;
	nsel[ntt] = kran;
	goto L290;
L270:
	for (nad = kans; nad <= ntt; ++nad) {
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
/* end Loop 210 */

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
	for (nsub = 1; nsub <= *nsam; ++nsub) {
	    nsel[nsub] = nrepr[nsub];
	}
L320:
	dysta2(nsam, jpp, &nsel[1], x, nn, &dys[1], ndyst, jtmd,
	       valmd, &jhalt);
	if (jhalt == 1)
	    continue;/* random sample*/

	kall = 1;
	s = 0.;
	l = 1;
	do {
	    ++l;
	    if (s < dys[l])
		s = dys[l];
	} while (l < nhalf);

	bswap2(kk, nsam, &nrepr[1], &dys[1], &z__, &s, tmp1, tmp2, tmp3);
	selec(kk, nn, jpp, ndyst, &zb, nsam, mdata, jtmd,
	      valmd, &nrepr[1], &nsel[1], &dys[1], x, &nr[1],
	      &nafs, &ttd[1], &radus[1], &ratt[1],
	      ntmp1, ntmp2, ntmp3, ntmp4, ntmp5, ntmp6, tmp1, tmp2);

	if (nafs) {
	    ++nunfs;
	}
	else if (!((jran != 1) && (zb >= zba))) {

	    zba = zb;
	    for (jjb = 1; jjb <= *kk; ++jjb) {
		ttbes[jjb] = ttd  [jjb];
		rdbes[jjb] = radus[jjb];
		rabes[jjb] = ratt [jjb];
	    }
	    for (jk = 1; jk <= *kk; ++jk)
		nrx[jk] = nr[jk];
	    for (js = 1; js <= *nsam; ++js)
		nbest[js] = nsel[js];
	    sx = s;
	}
    }
/* --- end random sampling loop */

    if (nunfs >= *nran) { *jstop = 1; return; }

/*     for the best subsample, the objects of the entire data set
     are assigned to their clusters */

    if (kall != 1) { *jstop = 2; return; }

    *obj = zba / rnn;
    dysta2(nsam, jpp, &nbest[1], x, nn, &dys[1], ndyst, jtmd,
	   valmd, &jhalt);
    resul(kk, nn, jpp, ndyst, mdata, jtmd, valmd, x, &nrx[1], mtt);
    if (*kk > 1) {
	black(kk, jpp, nn, nsam, &nbest[1], &dys[1], &sx, x, avsyl,
	      ttsyl, sylinf, ntmp1, ntmp2, ntmp3, ntmp4, tmp1, tmp2);
    }
    return;
} /* End clara() ---------------------------------------------------*/


void dysta2(int *nsam, int *jpp, int *nsel,
	    double *x, int *nn, double *dys, int *ndyst, int *jtmd,
	    double *valmd, int *jhalt)
{

    /* Local variables */
    int j, k, kj, l, lj, ksel, lsel, nlk, npres, lsubt;
    double pp, clk, rpres;

    /* Parameter adjustments */
    --dys;
    --nsel;
    --valmd;
    --jtmd;
    --x;

    /* Function Body */
    pp = (double) (*jpp);
    nlk = 1;
    dys[1] = 0.;
    for (l = 2; l <= *nsam; ++l) {
	lsubt = l - 1;
	lsel = nsel[l];
	for (k = 1; k <= lsubt; ++k) {
	    ksel = nsel[k];
	    clk = 0.;
	    ++nlk;
	    npres = 0;
	    for (j = 1; j <= *jpp; ++j) {
		lj = (lsel - 1) * *jpp + j;
		kj = (ksel - 1) * *jpp + j;
		if (jtmd[j] < 0) {
/* in the following line, x[-1] ==> seg.fault {BDR to R-core, Sat, 3 Aug 2002} */
		    if (x[lj] == valmd[j]) {
			continue /* next j */;
		    }
		    if (x[kj] == valmd[j]) {
			continue /* next j */;
		    }
		}
		++npres;
		if (*ndyst == 1)
		    clk += (x[lj] - x[kj]) * (x[lj] - x[kj]);
		else
		    clk += fabs(x[lj] - x[kj]);
	    }
	    rpres = (double) npres;
	    if (npres == 0) {
		*jhalt = 1;
		dys[nlk] = -1.;
	    } else {
		if (*ndyst == 1) {
		    dys[nlk] = sqrt(clk * (pp / rpres));
		} else {
		    dys[nlk] = clk * (pp / rpres);
		}
	    }
	}
    }
    return;
} /* End dysta2() -----------------------------------------------------------*/

void randm(int *nrun, double *ran)
{
    int k;
    double ry;

/*   we programmed this generator ourselves because we wanted it
   to be machine independent. it should run on most computers
   because the largest int used is less than 2**30 . the period
   is 2**16=65536, which is good enough for our purposes. */
    *nrun = *nrun * 5761 + 999;
    k = *nrun / 65536;
    *nrun -= k << 16;
    ry = (double) (*nrun);
    *ran = ry / 65536.f;
    return;
} /* randm() */

void bswap2(int *kk, int *nsam, int *nrepr,
	    double *dys, double *sky, double *s, double *dysma,
	    double *dysmb, double *beter)
{

    int j, ja, k, kbest, nbest;
    int nkj, njn, njaj, nny, nmax;

    double ammax, rsam, small, asky, dzsky, dz, cmd;


/* ====== first algorithm: build. ====== */

/* Parameter adjustments */
    --beter;
    --dysmb;
    --dysma;
    --dys;
    --nrepr;

    /* Function Body */
    nny = 0;
    for (j = 1; j <= *nsam; ++j) {
	nrepr[j] = 0;
	dysma[j] = *s * 1.1 + 1.;
    }
    /* --- Loop --- */
    do {
	for (ja = 1; ja <= *nsam; ++ja) {
	    if (nrepr[ja] == 0) {
		beter[ja] = 0.;
		for (j = 1; j <= *nsam; ++j) {
		    njaj = meet_(&ja, &j);
		    cmd = dysma[j] - dys[njaj];
		    if (cmd > 0.)
			beter[ja] += cmd;
		}
	    }
	}
	ammax = 0.;
	for (ja = 1; ja <= *nsam; ++ja) {
	    if (nrepr[ja] == 0 && ammax <= beter[ja]) {
		ammax = beter[ja];
		nmax = ja;
	    }
	}
	nrepr[nmax] = 1;
	++nny;
	for (j = 1; j <= *nsam; ++j) {
	    njn = meet_(&nmax, &j);
	    if (dys[njn] < dysma[j]) {
		dysma[j] = dys[njn];
	    }
	}
    } while (nny != *kk);
    /* --- */

    *sky = 0.;
    for (j = 1; j <= *nsam; ++j)
	*sky += dysma[j];

    if (*kk == 1)
	return;

    rsam = (double) (*nsam);
    asky = *sky / rsam;

/*------- second algorithm: SWAP. -------- */

/* Big LOOP : */
L60:

    for (j = 1; j <= *nsam; ++j) {
	dysma[j] = *s * 1.1 + 1.;
	dysmb[j] = *s * 1.1 + 1.;
	for (ja = 1; ja <= *nsam; ++ja) {
	    if (nrepr[ja] != 0) {
		njaj = meet_(&ja, &j);
		if (dys[njaj] < dysma[j]) {
		    dysmb[j] = dysma[j];
		    dysma[j] = dys[njaj];
		} else if (dys[njaj] < dysmb[j])
		    dysmb[j] = dys[njaj];
	    }
	}
    }

    dzsky = 1.;
    for (k = 1; k <= *nsam; ++k) {
	if (nrepr[k] != 1) {
	    for (ja = 1; ja <= *nsam; ++ja) {
		if (nrepr[ja] != 0) {
		    dz = 0.;
		    for (j = 1; j <= *nsam; ++j) {
			njaj = meet_(&ja, &j);
			nkj = meet_(&k, &j);
			if (dys[njaj] == dysma[j]) {
			    small = dysmb[j];
			    if (small > dys[njaj])
				small = dys[nkj];
			    dz = dz - dysma[j] + small;
			}
			else if (dys[nkj] < dysma[j])
			    dz = dz - dysma[j] + dys[nkj];
		    }
		    if (dz < dzsky) {
			dzsky = dz;
			kbest = k;
			nbest = ja;
		    }
		}
	    }
	}
    }
    if (dzsky >= 0.)
	return;

    nrepr[kbest] = 1;
    nrepr[nbest] = 0;
    *sky += dzsky;
    goto L60;

} /* End of bswap2() -------------------------------------------------- */

/* selec() : called once [per random sample] from clara() */
void selec(int *kk, int *nn, int *jpp, int *ndyst,
	   double *zb, int *nsam, int *mdata, int *jtmd,
	   double *valmd, int *nrepr, int *nsel, double *dys,
	   double *x, int *nr, /*logical*/int *nafs, double *ttd,
	   double *radus, double *ratt, int *nrnew, int *nsnew,
	   int *npnew, int *ns, int *np, int *new,
	   double *ttnew, double *rdnew)
{
/* nafs = .true. if a distance cannot be calculated */

    /* Local variables */
    int j, jk, jn, jp, jkabc, jnew, ka, kb;
    int newf, nrjk,  npab, nstrt, na, nb, npa, npb, njk;

    double abc, dnull, dsum, pres, pp, tra, rns;


/* Parameter adjustments */
    --ratt;
    --radus;
    --ttd;
    --nr;
    --x;
    --valmd;
    --jtmd;
    --rdnew;
    --ttnew;
    --new;
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

/* identification of representative objects, and initializations */

    jk = 0;
    for (j = 1; j <= *nsam; ++j) {
	if (nrepr[j] != 0) {
	    ++jk;
	    nr	 [jk] = nsel[j];
	    ns	 [jk] = 0;
	    ttd	 [jk] = 0.;
	    radus[jk] = -1.;
	    np	 [jk] = j;
	}
    }

/* assignment of the objects of the entire data set to a cluster,
     computation of some statistics, determination of the
     new ordering of the clusters */

    *zb = 0.;
    pp = (double) (*jpp);
    newf = 0;
    jn = 0;
/* L15:
 *     ++jn;
 */
    do { jn++;
	if (*mdata == 0) {
	    for (jk = 1; jk <= *kk; ++jk) {
		dsum = 0.;
		nrjk = nr[jk];
		if (nrjk != jn) {
		    for (jp = 1; jp <= *jpp; ++jp) {
			na = (nrjk - 1) * *jpp + jp;
			nb = (jn - 1) * *jpp + jp;
			tra = fabs(x[na] - x[nb]);
			if (*ndyst == 1) {
			    tra *= tra;
			}
			dsum += tra;
		    }
		    if (jk != 1 && dsum >= dnull)
			continue /* next jk */;
		}
		dnull = dsum;
		jkabc = jk;
	    }
	}
	else {
	    pres = 0.;
	    for (jk = 1; jk <= *kk; ++jk) {
		dsum = 0.;
		nrjk = nr[jk];
		if (nrjk != jn) {
		    abc = 0.;
		    for (jp = 1; jp <= *jpp; ++jp) {
			na = (nrjk - 1) * *jpp + jp;
			nb = (jn   - 1) * *jpp + jp;
			if (jtmd[jp] < 0) {
			    if (x[na] == valmd[jp]) {
				continue /* next jp */;
			    }
			    if (x[nb] == valmd[jp]) {
				continue /* next jp */;
			    }
			}
			abc += 1.;
			tra = fabs(x[na] - x[nb]);
			if (*ndyst == 1) {
			    tra *= tra;
			}
			dsum += tra;
		    }
		    if (abc < 0.5) {
			continue /* next jk */;
		    }
		    dsum = dsum * abc / pp;
		}
		if (pres <= 0.5)
		    pres = 1.;
		else if (dsum >= dnull)
		    continue /* next jk */;

		dnull = dsum;
		jkabc = jk;
	    }/* for(jk ..) */

	    if (pres <= 0.5) {
		*nafs = (1);
		return;
	    }
	} /* if (*mdata..) else */

	if (*ndyst == 1)
	    dnull = sqrt(dnull);

	*zb += dnull;
	ttd[jkabc] += dnull;
	if (radus[jkabc] < dnull)
	    radus[jkabc] = dnull;

	++ns[jkabc];
	if (newf < *kk) {
	    if (newf != 0) {
		for (jnew = 1; jnew <= newf; ++jnew) {
		    if (jkabc == new[jnew]) {
			goto L90;/* next jn */
		    }
		}
	    }
	    ++newf;
	    new[newf] = jkabc;
	}
    L90:
	;
    } while(jn < *nn);


/*     a permutation is carried out on vectors nr,ns,np,ttd,radus
     using the information in vector new. */

    for (jk = 1; jk <= *kk; ++jk) {
	njk = new[jk];
	nrnew[jk] = nr[njk];
	nsnew[jk] = ns[njk];
	npnew[jk] = np[njk];
	ttnew[jk] = ttd[njk];
	rdnew[jk] = radus[njk];
    }
    for (jk = 1; jk <= *kk; ++jk) {
	nr[jk] = nrnew[jk];
	ns[jk] = nsnew[jk];
	np[jk] = npnew[jk];
	ttd[jk] = ttnew[jk];
	radus[jk] = rdnew[jk];
    }
    for (j = 1; j <= *kk; ++j) {
	rns = (double) ns[j];
	ttd[j] /= rns;
    }

    if (*kk > 1) {

	/* computation of minimal distance of medoid ka to any
	   other medoid for comparison with the radius of cluster ka. */

	for (ka = 1; ka <= *kk; ++ka) {
	    nstrt = 0;
	    npa = np[ka];
	    for (kb = 1; kb <= *kk; ++kb) {
		if (kb == ka)
		    continue /* next kb */;

		npb = np[kb];
		npab = meet_(&npa, &npb);
		if (nstrt == 0)
		    nstrt = 1;
		else if (dys[npab] >= ratt[ka])
		    continue /* next kb */;

		ratt[ka] = dys[npab];
		if (ratt[ka] != 0.)
		    continue /* next kb */;

		ratt[ka] = -1.;
	    }
	    if (ratt[ka] > -0.5)
		ratt[ka] = radus[ka] / ratt[ka];
	}
    }
    return;
} /* End selec() -----------------------------------------------------------*/

void resul(int *kk, int *nn, int *jpp, int *ndyst,
	   int *mdata, int *jtmd, double *valmd,
	   double *x, int *nrx, int *mtt)
{
    /* Local variables */
    int j, jksky, ja, jk, jn, jna, ka, na, nb, njnb, nrjka, nxja, nrjk;
    double pp, abc, dnull, dsum, tra;

/* Parameter adjustments */
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
    for (jk = 1; jk <= *kk; ++jk) {
	if (nrx[jk] == jn) {
	    goto L220;
	}
    }
    jna = (jn - 1) * *jpp + 1;
    if (*mdata != 0) {
	goto L170;
    }
    for (jk = 1; jk <= *kk; ++jk) {
	dsum = 0.;
	nrjk = (nrx[jk] - 1) * *jpp;
	for (j = 1; j <= *jpp; ++j) {
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
	if (dnull > dsum) {
	    dnull = dsum;
	    jksky = jk;
	}
    }
    goto L200;
L170:
    for (jk = 1; jk <= *kk; ++jk) {
	dsum = 0.;
	nrjk = (nrx[jk] - 1) * *jpp;
	abc = 0.;
	for (j = 1; j <= *jpp; ++j) {
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
	    abc += 1.;
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
    for (jk = 1; jk <= *kk; ++jk) {
	nrjk = nrx[jk];
	nrjka = (nrjk - 1) * *jpp + 1;
	x[nrjka] = (double) jk;
    }
    for (ka = 1; ka <= *kk; ++ka) {
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
    return;
} /* end resul() -----------------------------------------------------------*/


void black(int *kk, int *jpp, int *nn, int *nsam, int *nbest,
	   double *dys, double *sx, double *x,
	   double *avsyl, double *ttsyl, double *sylinf,
	   int *ncluv, int *nsend, int *nelem, int *negbr,
	   double *syl, double *srank)
{
/* Silhouettes computation and "drawing"  --> syl() and sylinf() */

    /* System generated locals */
    int sylinf_dim1, sylinf_offset;

    /* Local variables */
    int lang;

    double att, btt, rtt, db, dysa, dysb, rsam, symax;
    int j, jna, l, lplac, nj, nl, nbb, ncase, nclu, numcl, nsylr, ntt;

/* Parameter adjustments */
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

/*
     construction of clustering vector (ncluv)
     of selected sample (nbest).
*/

    /* Function Body */
    for (l = 1; l <= *nsam; ++l) {
	ncase = nbest[l];
	jna = (ncase - 1) * *jpp + 1;
	ncluv[l] = (int) (x[jna] + .1f);
    }

/*     drawing of the silhouettes */

    nsylr = 0;
    *ttsyl = 0.;
    for (numcl = 1; numcl <= *kk; ++numcl) {
	ntt = 0;
	for (j = 1; j <= *nsam; ++j) {
	    if (ncluv[j] == numcl) {
		++ntt;
		nelem[ntt] = j;
	    }
	}
	for (j = 1; j <= ntt; ++j) {
	    nj = nelem[j];
	    dysb = *sx * 1.1 + 1.;
	    negbr[j] = -1;

	    for (nclu = 1; nclu <= *kk; ++nclu) {
		if (nclu != numcl) {
		    nbb = 0;
		    db = 0.;
		    for (l = 1; l <= *nsam; ++l) {
			if (ncluv[l] == nclu) {
			    ++nbb;
			    db += dys[meet_(&nj, &l)];
			}
		    }
		    btt = (double) nbb;
		    db /= btt;
		    if (db < dysb) {
			dysb = db;
			negbr[j] = nclu;
		    }
		}
	    }

	    if (ntt == 1) {
		syl[j] = 0.;	    continue /* j */;
	    }
	    dysa = 0.;
	    for (l = 1; l <= ntt; ++l) {
		nl = nelem[l];
		dysa += dys[meet_(&nj, &nl)];
	    }
	    att = (double) (ntt - 1);
	    dysa /= att;
	    if (dysa <= 0.) {
		if (dysb > 0.)
		    syl[j] = 1.;
		else
		    syl[j] = 0.;

		continue /* j */;
	    }

	    if (dysb > 0.) {
		if (dysb > dysa)
		    syl[j] = 1. - dysa / dysb;
		else if (dysb < dysa)
		    syl[j] = dysb / dysa - 1.;
		else /* (dysb == dysa) */
		    syl[j] = 0.;
	    }
	    else {
		syl[j] = -1.;
	    }
	    if (syl[j] < -1.)
		syl[j] = -1.;

	    if (syl[j] > 1.)
		syl[j] = 1.;
	} /* for(j ..) */

	avsyl[numcl] = 0.;
	for (j = 1; j <= ntt; ++j) {
	    symax = -2.;
	    for (l = 1; l <= ntt; ++l) {
		if (syl[l] > symax) {
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

	if (ntt >= 2) {
	    for (l = 1; l <= ntt; ++l) {
		lplac = nsend[l];
		ncase = nelem[lplac];
		++nsylr;
		sylinf[nsylr + sylinf_dim1] = (double) numcl;
		sylinf[nsylr + (sylinf_dim1 << 1)] = (double) negbr[lplac];
		sylinf[nsylr + sylinf_dim1 * 3] = srank[l];
		sylinf[nsylr + (sylinf_dim1 << 2)] = (double) nbest[ncase];
	    }
	}
	else {
	    ncase = nelem[1];
	    ++nsylr;
	    sylinf[nsylr + sylinf_dim1] = (double) numcl;
	    sylinf[nsylr + (sylinf_dim1 << 1)] = (double) negbr[1];
	    sylinf[nsylr + sylinf_dim1 * 3] = 0.;
	    sylinf[nsylr + (sylinf_dim1 << 2)] = (double) nbest[ncase];
	}

    }
    rsam = (double) (*nsam);
    *ttsyl /= rsam;
    return;
} /* black */
