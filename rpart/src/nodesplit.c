/*
 * nodesplit -- Split the node in two, and keep a count as we do of how
 *  many splits are determined by each surrogate variable.
 *
 * me	  : pointer to the node of the tree
 * nodenum: the node number of the current node, used to update "which"
 * n1, n2 : starting and ending indices for the observation numbers
 * nnleft, nnright: at the end, how many obs were sent to the left and
 *            to the right.  Beware-  for certain settings of the
 *            usesurrogate option, some obs go neither left nor right
 */
#include "rpart.h"
#include "node.h"
#include "rpartproto.h"
#include <stdio.h>

void
nodesplit(pNode me, int nodenum, int n1, int n2, int *nnleft, int *nnright)
{
    pSplit tsplit = me->primary;
    int someleft = 0;
    int leftson = 2 * nodenum;  /* the label that will go with the left son */
    int rightson = leftson + 1;

    /*
     * Walk through the variables (primary, then surrogate 1, then surr 2...)
     *   and reassign "rp.which"
     */
    int pvar = tsplit->var_num;     /* primary variable */

    int nleft = 0, nright = 0;
    if (rp.numcat[pvar] > 0) {  /* categorical primary variable */
	int *index = tsplit->csplit;
	for (int i = n1; i < n2; i++) {
	    int j = rp.sorts[pvar][i];
	    if (j < 0)
		someleft++;     /* missing value */
	    else
		switch (index[(int) rp.xdata[pvar][j] - 1]) {
		case LEFT:
		    rp.which[j] = leftson;
		    nleft++;
		    break;
		case RIGHT:
		    rp.which[j] = leftson + 1;
		    nright++;
		    break;
		}
	}
    } else {
	double psplit = tsplit->spoint;        /* value of split point */
	int extra = tsplit->csplit[0];
	for (int i = n1; i < n2; i++) {
	    int j = rp.sorts[pvar][i];
	    if (j < 0)
		someleft++;
	    else {
		int k = (rp.xdata[pvar][j] < psplit) ? extra : -extra;
		if (k == LEFT) {
		    rp.which[j] = leftson;
		    nleft++;
		} else {
		    rp.which[j] = leftson + 1;
		    nright++;
		}
	    }
	}
    }

    /*
     * Now the surrogates
     *   Usually, there aren't a lot of observations that need to
     *   be split.  So it is more efficient to make one 1:n pass,
     *   with multiple runs through the surrogate list.
     */
    if (someleft > 0 && rp.usesurrogate > 0) {
	for (int i = n1; i < n2; i++) {
	    int j = rp.sorts[pvar][i];
	    if (j >= 0)
		continue;       /* already split */

	    j = -(j + 1);       /* obs number - search for surrogates */
	    for (tsplit = me->surrogate; tsplit; tsplit = tsplit->nextsplit) {
		int var = tsplit->var_num;
		if (!R_FINITE(rp.xdata[var][j]))
		    continue;
		/* surrogate not missing - process it */

		if (rp.numcat[var] > 0) {       /* categorical surrogate */
		    int *index = tsplit->csplit;
		    int k = (int) rp.xdata[var][j]; /* the value of the surrogate  */
		    /*
		     * The need for the if stmt below may not be obvious.
		     * The surrogate's value must not be missing, AND there
		     * must have been at least 1 person with both this
		     * level of the surrogate and a primary split value
		     * somewhere in the node.  If everyone in this node
		     * with level k of the surrogate also had a missing
		     * value of the primary variable, then index[k-1] will
		     * be zero.
		     */
		    if (index[k - 1]) {
			tsplit->count++;
			if (index[k - 1] == LEFT) {
			    rp.which[j] = leftson;
			    nleft++;
			} else {
			    rp.which[j] = leftson + 1;
			    nright++;
			}
			someleft--;
			break;
		    }
		} else {
		    double psplit = tsplit->spoint;  /* continuous surrogate */
		    int extra = tsplit->csplit[0];
		    tsplit->count++;
		    int k = (rp.xdata[var][j] < psplit) ? extra : -extra;
		    if (k == LEFT) {
			rp.which[j] = leftson;
			nleft++;
		    } else {
			rp.which[j] = rightson;
			nright++;
		    }
		    someleft--;
		    break;
		}
	    }
	}
    }
    if (someleft > 0 && rp.usesurrogate == 2) {
	/* all surrogates missing, use the default */
	int i = me->lastsurrogate;
	if (i) {           /* 50-50 splits are possible - there is no
			    * "default" */
	    int lastisleft;
	    if (i < 0) {
		lastisleft = leftson;
		nleft += someleft;
	    } else {
		lastisleft = rightson;
		nright += someleft;
	    }

	    for (int i = n1; i < n2; i++) {
		int j = rp.sorts[pvar][i];
		/*
		 * only those who weren't split by the primary (j < 0) and
		 * weren't split by a surrogate (rp.which == nodenum) need to be
		 * assigned
		 */
		if (j < 0) {
		    j = -(j + 1);
		    if (rp.which[j] == nodenum)
			rp.which[j] = lastisleft;
		}
	    }
	}
    }
    /*
     * Last part of the work is to update the rp.sorts matrix
     *
     * Say that n1=5, n2=12, 4 go left, 3 go right, and one obs
     *   stays home, and the data looks like this:
     *
     *   sorts[var][5] = 17    rp.which[17]= 17 = 2*nodenum +1 = rightson
     *       "            4    rp.which[4] = 16 = 2*nodenum    = leftson
     *                   21          21   16
     *                    6           6   17
     *      "		  7           7   16
     *                   30          30   17
     *                    8           8   16
     *	sorts[var][12]=	-11    rp.which[11] = 8 = nodenum  (X = missing)
     *
     *  Now, every one of the rows of the sorts contains these same
     *    elements -- 4,6,7,8,11,17,21,30 -- in some order, rp.which
     *    represents both how they are sorted and any missings via
     *    negative numbers.
     *  We need to reorder this as {{goes left}, {goes right}, {stays
     *    here}}, preserving order within the first two groups.  The
     *    order within the "stay here" group doesn't matter since they
     *    won't be looked at again for splitting.
     *  So the result in this case should be
     *       4, 21, 7, 8,   17, 6, 30,  -11
     *  The algorithm is the opposite of a merge sort.
     *
     *  Footnote: if no surrogate variables were used, then one could
     *   skip this process for the primary split variable, as that
     *   portion of "sorts" would remain unchanged.  It's not worth
     *   the bother of checking, however.
     */
    for (int k = 0; k < rp.nvar; k++) {
	int *sindex = rp.sorts[k];   /* point to variable k */
	int i1 = n1,i2 = i1 + nleft, i3 = i2 + nright;
	for (int i = n1; i < n2; i++) {
	    int j = sindex[i];
	    if (j < 0)
		j = -(j + 1);
	    if (rp.which[j] == leftson)
		sindex[i1++] = sindex[i];
	    else {
		if (rp.which[j] == rightson)
		    rp.tempvec[i2++] = sindex[i];
		else
		    rp.tempvec[i3++] = sindex[i];       /* went nowhere */
	    }
	}
	for (int i = n1 + nleft; i < n2; i++)
	    sindex[i] = rp.tempvec[i];
    }

    *nnleft = nleft;
    *nnright = nright;
}
