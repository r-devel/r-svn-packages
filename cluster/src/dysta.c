/* Produced by f2c and Martin's  f2c-clean,v 1.11 2012/05/04 --
   and some pretty editing
*/

#include <math.h>

/* Dysta() :

 Compute Distances from X matrix {also for agnes() and diana()}:
 -----------------------------------------------------------

 was part of  pam.f (then called both from Fortran & C)
 "keep in sync" with  daisy.c {now both moved to C}
 ((*no* longer called from ../tests/dysta-ex.R ...))
*/
int dysta(int nn, int p,
	  double *x,  // [nn x p] matrix
	  double *dys,// length 1 + nn*(nn-1)/2 -- with first entry := 0.  (FIXME: drop that)
	  int ndyst, /* ndyst = 1 : euclidean ; "else"(2) : manhattan */
	  int *jtmd, double *valmd)
{
    int x_dim1 = nn;
    int nlk = 0, jhalt = 0;
/*    dys[nlk] = 0.; */
/*  -----------== is used potentially for  d[i,i] == dys[0] == 0  (FIXME) */
    double pp = (double) p;
    for (int l = 1; l < nn; ++l) {
	for (int k = 0; k < l; ++k) {
	    double clk = 0.;
	    int npres = 0;
	    for (int j = 0; j < p; ++j) {
		if (jtmd[j] < 0) { /* some  x(*,j) are missing (NA) */
		    if((x[l + j * x_dim1] == valmd[j]) ||
		       (x[k + j * x_dim1] == valmd[j])) continue; // next j
		}
		++npres;
		if (ndyst == 1) {
		    clk += (x[l + j * x_dim1] - x[k + j * x_dim1]) *
			   (x[l + j * x_dim1] - x[k + j * x_dim1]);
		} else {
		    clk += fabs(x[l + j * x_dim1] - x[k + j * x_dim1]);
		}
	    }
	    if (npres == 0) {
		jhalt = 1;
		dys[nlk] = -1.;
	    } else if(ndyst == 1) {
		dys[nlk] = sqrt(clk * (pp / (double) npres));
	    } else {
		dys[nlk] = clk * (pp / (double) npres);
	    }
	    ++nlk;
	} // for( k )
    }
    return jhalt;
} /* dysta_ */

