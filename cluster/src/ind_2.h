/* inlined, to be included in pam.c, clara.c, and twins.c */

static R_INLINE
int ind_2(int l, int j)
{
/* Utility, originally FORTRAN,	 called "meet"; called from CLARA, PAM & TWINS.
 * Original code had MEET(), MEET2(), and MEET3() in the 3 source files.

 * ind_2(l,j) returns the *index* of dys() where diss. d(l,j) is stored:
 *	  d(l,j) == dys[ind_2(l,j)]
 *
 * MM: changed to work with 0-origin matrices  dys[], but l,j are >= 1
 */
#ifdef was_orig
    if(l > j)
	return (l-2)*(l-1)/2 + j;
    else if(l == j)
	return 0;/* and the first element, dys[0] is := 0.  permanently! */
    else /* l < j */
	return (j-2)*(j-1)/2 + l;
#else
    /* from Li Long -- optimizes particularly well on Itanium2 */

    /* max_m check by MMächler: max_m  is the largest integer m
     * ----- such that (m-2)*(m-1) < MAX_INT : */
# define max_m 46342

    int result = 0;
    if (l != j) {
	int m = l>j ? l : j;
	int n = l>j ? j : l;
	result = (m <= max_m)
	    ? (m-2)*(m-1)/2 + n
	    : (int) (((double) m-2)*(m-1)/2 + n);
    }
    return result;
#endif
}

/* A version of  dys[ ind_2(l, j) ]   which does *not* address "dys[0]" == 0. : */
static R_INLINE
double dys_2(const double dys[], int l, int j)
{
    double res = 0.;
    if (l != j) {
	int m = l>j ? l : j; // m := max(l,j) >= 2
	    if (l>j)  l = j; // l := min(l,j) >= 1
	res = dys[(m <= max_m)
		        ? (m-2)*(m-1)/2 + l
	: (int) (((double) m-2)*(m-1)/2 + l)];
    }
    return res;
}
#undef max_m
