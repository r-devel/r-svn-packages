/* Code to implement kd tree based nearest neighbour routines in C.
   Based on approach of Press et al Numerical Recipes 3rd ed. 21.2, but 
   re-implemented in vanilla C. Design is based on principles from Press
   et al. but due to licensing, this is a re-write. So results need not
   correspond to Press et al. code in detail (e.g. exactly which point 
   is in which box of tree, and exact indexing details). Efficiency 
   should be same, however.

   R CMD SHLIB kd-tree.c 
   to build.
   dyn.load("kd-tree.so")
   to load into R.
*/

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <math.h>
#include <stdlib.h>

typedef struct { /* defines structure for kd-tree box */
  double *lo,*hi;    /* box defining co-ordinates */
  int parent,child1,child2, /* indices of parent and 2 offspring */
      p0,p1;         /* indices of first and last point in box */
} box_type; 



typedef struct {
  box_type *box;
  int *ind, /* index of points in coordinate matrix tree relates to */
      *rind, /* where is ith row of X in ind? */
      n_box, /* number of boxes */
      d; /* dimension */
  double huge; /* number indicating an open boundary */
} kdtree_type;

void k_order(int *k,int *ind,double *x,int *n) {
/* ind is of length n.
   x is of length max(ind)+1, at least.
   On exit x[ind[0..k-2]] <= x[ind[k-1]] <= x[ind[k..n-1]],
   by permuting ind, but not x. i.e. we split x into two 
   groups... those less than or equal to the kth largest, and 
   those greater than or equal to the kth largest.

   This works by picking an arbitrary pivot value from x, and arranging 
   ind so that x[ind] is partitioned about that value. Now we know 
   immediately which of the two partitions must contain the kth 
   largest value, so we repeat the process on this one. Eventually 
   the partition that must contain the kth largest value has only 1 
   or 2 elements in it... 

   R test code...
   dyn.load("kd-tree.so")
   n <- 20
   x <- runif(n)
   k <- 1
   oo <- .C("k_order",k=as.integer(k),ind=as.integer(0:19),x=as.double(x),n=as.integer(n))
   ind <- oo$ind+1
   sort(x)[k]
   x[ind[k]]
   ## following is TRUE if it's worked...
   sum(!c(x[ind[1:(k-1)]]<=x[ind[k]],x[ind[(k+1):n]]>=x[ind[k]]))==0

*/ 

  int l,r,m,ip,ri,li,dum;
  double xp;
  l = 0;   /* leftmost point of current partition */
  r = *n-1; /* rightmost point of current partitions */
  while (1) {
    if (r > l+1) { /* partition large enough to need work still */
      m = (l+r) / 2; /* pick a point from partition midpoint (by location not value) 
                        (Press et al say to do this to avoid poor behaviour on already 
                         sorted x).*/ 
      dum = ind[l+1];ind[l+1] = ind[m];ind[m] = dum; /* swap points m and l+1 */
      /* now re-arrange so that x[ind[l]] < x[ind[l+1]] < x[ind[r]]... */
      if (x[ind[l]] > x[ind[r]]) { /* swap r and l */
        dum = ind[r];ind[r] = ind[l];ind[l] = dum;
      }
      if (x[ind[l]] > x[ind[l+1]]) { /* swap l and l+1 */
        dum = ind[l];ind[l] = ind[l+1];ind[l+1] = dum;
      } 
      else if (x[ind[l+1]] > x[ind[r]]) { /* swap l+1 and r */
        dum = ind[l+1];ind[l+1] = ind[r];ind[r] = dum;
      }
      ip = ind[l+1]; /* index of pivot */
      xp = x[ip];    /* pivot value */
      /* so pivot is xp = x[ind[l+1]]. start proccess of shuffling array into
         two partitions containing all the values less than xp, and all 
         those larger than xp... */ 
      ri = r-1;  /* start searching down partition from here for wrongly located 
                    values (pos r above pivot already) */
      li = l+2; /* start searching up from here (pos l is already below pivot, l+1 is pivot)*/
      while (1) {
        while(x[ind[li]]<xp) li++; /* move up until value on wrong side found */
        while(x[ind[ri]]>xp) ri--; /* move down until value on wrong side found */
        if (ri<=li) break; /* partitions correct now */
        dum = ind[ri];ind[ri] = ind[li];ind[li] = dum; /* swap ri and li (to correct sides) */
      } /* end of partitioning loop */
      /* now put pivot into right place (and value in that place in correct partition) */
      ind[l+1] = ind[ri];
      ind[ri] = ip;
      /* Now select the partition in which kth largest must lie, by setting new end points */
      if (ri > *k) r = ri - 1; else l = li;
    } else { /* the partition can only contain 1 or 2 points */
      if (r == l+1 && x[ind[r]] < x[ind[l]]) { /* contains two points, but in wrong order */
         dum = ind[r];ind[r] = ind[l];ind[l] = dum; /* so swap indices */
      } 
      return; /* x[ind[k]] is kth largest value in x */ 
    }
  } /* end while(1) - main loop */
}

void free_kdtree(kdtree_type kd) {
/* free a kdtree */
  free(kd.ind);free(kd.rind);
  free(kd.box[0].lo); /* storage for box coordinates */
  free(kd.box);
}

void kd_tree(double *X,int *n, int *d,kdtree_type *kd) {
/* Create a kd tree for the points in n by d matrix X.
   X is in column order. Each row is one point. 
   At end of process... 
   * box[i] contains points indexed by 
     ind[box[i].p0..box[i].p1] 
   * box[i] has one parent and 2 children, unless it contains only one 
     or 2 points, in which case it has no children.

   
*/
  int *ind,*rind,*p,i,m,todo[50],todo_d[50],item,bi,nb,np,k,dim,b,p0,p1;
  box_type *box;
  double huge=1e100,*pd,*x,*dum1,*dum2,*dum3;
  /* create index for points... */
  ind = (int *)calloc((size_t) *n,sizeof(int)); 
  for (i=0,p=ind;i < *n;i++,p++) *p = i; 
  /* Find the number of boxes in the tree */
  m=2;while (m < *n) m *= 2;
  nb = *n * 2 - m / 2 - 1;
  if (nb > m-1) nb = m - 1; 
  /* Now make an array of boxes (all cleared to zero)... */
  box = (box_type *)calloc((size_t)nb,sizeof(box_type));
  /* allocate storage for box defining coordinates... */ 
  pd = (double *)calloc((size_t)nb * (2 * *d),sizeof(double));
  for (i=0;i<nb;i++) {
    box[i].lo = pd;pd += *d;
    box[i].hi = pd;pd += *d;
  } /* note box[0].lo is now the thing to free when done (not p!) */

  /* set up box[0]... */
  for (i=0;i< *d;i++) { /* arbitrary bounding box */
    box[0].lo[i] = -huge;box[0].hi[i] = huge;
  }
  box[0].p1 = *n-1; /* last index item in this box (.p0 is first) */
  todo[0]=0;   /* put box[0] on todo list for processing */ 
  todo_d[0]=0; /* which dimension to start with */
  item=0;      /* which item of todo list to do next (item+1 is number of items on list) */
  /* each box contains 1 or 2 points, or it gets split into
     2 smaller boxes. If the smaller boxes contain 3 or more 
     points then they are added to the todo list. The todo 
     list is always worked on from the end (i.e by processing 
     box todo[item]) */
  bi=0; /* box index, counting number of boxes created so far */
  while (item >= 0) {   /* todo list still has items */
    b = todo[item];     /* current box */
    dim = todo_d[item]; /* dimension on which to split box */ 
    p0 = box[b].p0;p1=box[b].p1; 
    np = p1-p0+1;      /* number of points in box k */
    x = X + dim * *n;  /* array of co-ordinates for current dimension to sort on */
    k = (np-1)/2;          /* split the box around kth value in box */ 
    /* next line re-orders the point index for this box only.
       after reordering the index is split into two parts, indexing
       points below and obove the kth largest value */  
 
    k_order(&k,ind+p0,x,&np); 
  
    /*... so the box is now split at a plane/line through x[ind[p0+k-1]] */
    item--; /* basically done that item */
    /* create the offspring boxes... */
    
    bi++; /* lower box first */ 
    
    if (bi>nb-1) Rprintf("too many boxes!!");
    box[b].child1=bi;/* record box relationships */
    /* copy box coordinates... */
    for (dum1=box[bi].lo,dum2=dum1 + *d,dum3=box[b].lo;dum1<dum2;dum1++,dum3++) *dum1 = *dum3;
    for (dum1=box[bi].hi,dum2=dum1 + *d,dum3=box[b].hi;dum1<dum2;dum1++,dum3++) *dum1 = *dum3;
    box[bi].hi[dim] = x[ind[p0+k]]; /* split location */
    box[bi].parent=b; /* record box relationships */
    box[bi].p0=box[b].p0;
    box[bi].p1=box[b].p0+k;
    if (k>1) { /* more than two points , so more work needed */ 
      item++;
      todo[item] = bi;
      todo_d[item] = dim+1;
      if (todo_d[item] == *d) todo_d[item] = 0;
    }
    bi++; /* now the higher box */ 
    if (bi>nb-1) Rprintf("too many boxes!!");
    box[b].child2=bi;/* record box relationships */
    /* copy box coordinates... */
    for (dum1=box[bi].lo,dum2=dum1 + *d,dum3=box[b].lo;dum1<dum2;dum1++,dum3++) *dum1 = *dum3;
    for (dum1=box[bi].hi,dum2=dum1 + *d,dum3=box[b].hi;dum1<dum2;dum1++,dum3++) *dum1 = *dum3;
    box[bi].lo[dim] = x[ind[p0+k]]; /* split location */
    box[bi].parent=b; 
    box[bi].p1=box[b].p1;
    box[bi].p0=box[b].p0+k+1;
    if (np-k>3) { /* more than two points , so more work needed */ 
      item++;
      todo[item] = bi;
      todo_d[item] = dim+1;
      if (todo_d[item] == *d) todo_d[item] = 0;
    }
  }  
  if (bi!=nb-1) Rprintf("bi not equal to nb-1 %d %d\n",bi,nb-1);
  rind = (int *)calloc((size_t) *n,sizeof(int));
  /* now create index of where ith row of X is in ind */
  for (i=0;i<*n;i++) rind[ind[i]]=i; 
  /* now put tree into kd object */
  kd->box = box;kd->ind = ind;kd->rind = rind;kd->n_box = nb;kd->huge = huge;
  kd->d = *d;
} /* end of kd_tree */


void Rkdtree(double *X,int *n, int *d,double *lo,double *hi,int *ind, int *rind) { 
/* Routine to export kdtree data to R 
     m <- 2;
     while (m<n) m <- m*2
     nb <- min(m-1,2*n-m/2-1)
     hi and lo are nb by k matrices
     ind and rind are n vectors
     X is an n by d matrix
*/
  kdtree_type kd;
  int i,j;
  kd_tree(X,n,d,&kd);
  for (i=0;i<*n;i++) {ind[i] = kd.ind[i];rind[i]=kd.rind[i];}
  for (j=0;j<*d;j++) for (i=0;i<kd.n_box;i++,lo++,hi++) { 
    *lo = kd.box[i].lo[j];
    *hi = kd.box[i].hi[j];
  }
  free_kdtree(kd);
}


inline void update_heap(double *h,int *ind,int n) {
/* h contains n numbers, such that h[i] > h[2*i+1] and
   h[i] > h[2*i+2] (each applying whenever elements exist).
   The exception is that h[0], may not obey these conditions. 
   This function re-arranges h so that it does. It also 
   applies the same re-arrangement to ind. Figure 8.3.1
   of Press et al (2007) shows what's going on.
*/
  double h0;
  int i,i0,ind0;
  h0 = h[0]; /* h0 should be largest element, in properly ordered heap */
  ind0 = ind[0];  /* index vector to re-shuffle exactly as h vector */
  i0 = 0; /* current position of h0 */
  i = 1; /* index for first child node of i0 */
  while (i < n) { /* work through to end of heap */
    if (i < n-1&&h[i]<h[i+1]) i++; /* i indexes largest of i0's child nodes */ 
    if (h0 > h[i]) break; /* h0 should be at h[i0] */
    /* since h0 <= h[i], move h[i] 'up' heap into h[i0], and move i0, the nominal 
       position for h0 'down' heap to i */
    h[i0] = h[i];
    ind[i0] = ind[i];
    i0 = i;
    i = 2*i+1; /* now move on to first child of h[i]... */
  }
  h[i0] = h0; /* put h0 into location it should occupy in heap */
  ind[i0] = ind0;
}

inline double box_dist(box_type *box,double *x,int d) {
/* find distance from d dimensional box to point x */
  double d2 = 0.0,z,*bl,*bh,*xd;
  for (xd=x+d,bl=box->lo,bh=box->hi; x < xd;x++,bl++,bh++) {
    if (*x < *bl) { z = *x - *bl;d2 += z*z;}
    if (*x > *bh) { z = *x - *bh;d2 += z*z;}
  } 
  return(sqrt(d2));
}

inline int which_box(kdtree_type *kd,int j) {
/* Finds smallest box in kd tree containing jth point 
   from point set used to create tree */ 
  int i,bi,b1;
  i = kd->rind[j]; /* where jth point is in kd->ind */
  bi=0;
  
  while (kd->box[bi].child1) { /* still haven't reached smallest */
   
    b1 = kd->box[bi].child1;   /* index of first child */
    if (kd->box[b1].p1>=i) bi = b1; /* point is in child1 */
    else bi = kd->box[bi].child2; /* mkd->box[bi].child1ust be in child2 */
  }
  return(bi); /* index of smallest box containing jth point */
}


int xbox(kdtree_type *kd,double *x) {
/* which box of the kd tree is point x located in? 
   For maximal efficiency use the fact that nested boxes are
   split along one dimension, and that the split dimensions are
   cycled through in the same order, while descending the tree. 
*/
  int bi,d,b1;
  box_type *box;
  bi=0; /* root of the tree - the big box */
  box = kd->box;
  d=0;  /* dimension for first split */
  while (box[bi].child1) { /* still not reached the outermost twig - smallest box*/
    b1 = box[bi].child1;
    /* note that points on boundary are in lower box (child1) */
    if (x[d] <= box[b1].hi[d]) bi = b1; else
    bi = box[bi].child2;
    d ++; if (d == kd->d) d=0;
  }
  return(bi);
}


inline double ijdist(int i, int j, double *X,int n,int d) {
/* return Euclidian distance between ith and jth rows of n by d 
   matrix X */
  double *pi,*pj,*pil,dist=0.0,x;
  for (pi=X+i,pil=pi+n*d,pj=X+j;pi<pil;pi+=n,pj+=n) {x = *pi - *pj;dist += x*x;} 
  return(sqrt(dist));
} 


void p_area(double *a,double *X,kdtree_type kd,int n,int d) {
/* Associates the volume of its kd box with each point. If the point 
   shares a box then the volume is split. If the box has an open boundary 
   then that boundary is shrunk so that the point is enclosed in it.
   Results returned in a.
*/
  double *wa,*lo,*hi,*x0,*x1,av_vol=0.0,min_w,x;
  int np,bi,i,j,k,ok=1,*count,check;

  wa = (double *)calloc((size_t)d,sizeof(double));
  lo = (double *)calloc((size_t)d,sizeof(double));
  hi = (double *)calloc((size_t)d,sizeof(double));
  x0 = (double *)calloc((size_t)d,sizeof(double));
  x1 = (double *)calloc((size_t)d,sizeof(double));
  count = (int *)calloc((size_t)d,sizeof(int));

  /* get average box widths, for fallback purposes */
  for (bi=0;bi<kd.n_box;bi++) {
    for (j=0;j<d;j++) {
      if (kd.box[bi].lo[j]!= - kd.huge && kd.box[bi].hi[j]!= kd.huge) {
        count[j]++; wa[j] += kd.box[bi].hi[j] - kd.box[bi].lo[j];
      }
    }
  }
  for (j=0;j<d;j++) wa[j] /= count[j];

  for (i=0;i<n;i++) {
    bi = which_box(&kd,i); /* locate smallest box containing point i */
    for (j=0;j<d;j++) { /* make copies of box boundaries to work on */
      lo[j] = kd.box[bi].lo[j];
      if (lo[j]==-kd.huge) ok = 0;
      hi[j] = kd.box[bi].hi[j];
      if (hi[j]==kd.huge) ok = 0;
    }
    np = kd.box[bi].p1-kd.box[bi].p0+1; /* number of points in box */
    if (!ok) { /* box is not finite */
      check = 0;
      k = kd.ind[kd.box[bi].p0]; /* row of X that first point is in */
      if (k==i) check=1;
      for (j=0;j<d;j++) x0[j] = X[k + j * n];
      if (np>1) { /* there is a second point to consider */
        k = kd.ind[kd.box[bi].p1];
        if (k==i) check=1;
        for (j=0;j<d;j++) x1[j] = X[k + j * n];
      }
      if (!check) Rprintf("indexing error in p_area!\n");
      /* now work through trying to shrink the limits */
      /* first shrink limits to points unless that collapses a dimension */
      ok=1; min_w = -1.0;
      for (j=0;j<d;j++) {
        if (lo[j] == -kd.huge) { /* attempt to shrink boundary to (lowest) point */
           x = x0[j]; if (np>1 && x1[j]<x) x = x1[j];  
           if (x < hi[j]) lo[j] = x; /* sorted! */
           else ok=0; /* not sorted! */
        } 
        if (hi[j] == kd.huge) { /* attempt to shrink boundary to (highest) point */
           x = x0[j]; if (np>1 && x1[j]>x) x = x1[j];  
           if (x > lo[j]) hi[j] = x; /* sorted! */
           else ok=0; /* not sorted! */
        } 
        if (lo[j] != -kd.huge && hi[j] != kd.huge) {
           x = hi[j]-lo[j];
           if (min_w < 0 || x < min_w) min_w = x;
        }
      } /* end of first pass through limits */
      if (!ok) { /* then there are unfixed limits left to deal with */
        for (j=0;j<d;j++) {
          if (lo[j] == -kd.huge) { /* attempt to shrink boundary to (lowest) point */
            x = x0[j]; if (np>1 && x1[j]<x) x = x1[j]; 
            if (min_w>0) x -= min_w; else x -= wa[j];
            lo[j] = x; /* sorted! */
          } 
          if (hi[j] == kd.huge) { /* attempt to shrink boundary to (highest) point */
            x = x0[j]; if (np>1 && x1[j]>x) x = x1[j];  
            if (min_w>0) x += min_w; else x += wa[j];
            hi[j] = x; /* sorted! */
          }
        } 
      } /* all limits now reset */
    } /* box is now finite */
    /* compute box volume */
    for (x=1.0,j=0;j<d;j++) x*= hi[j]-lo[j];
    x /= np; /* share it out */
    a[i] = x; /* store it */ 
  }
  free(count);
  free(x0);free(x1);free(lo);free(hi);free(wa);
}


void k_nn_work(kdtree_type kd,double *X,double *dist,int *ni,int *n,int *d,int *k) {
/* Given a kd tree, this routine does the actual work of finding the nearest neighbours.
*/
  int i,j,bi,*ik,bii,todo[100],item,pcount,*ind;
  box_type *box;
  double *dk,huge,*p,*p1,*p2,dij,*x;
 
  huge = kd.huge;
  ind = kd.ind;
  box = kd.box;

  dk = (double *)calloc((size_t)*k,sizeof(double)); /* distance k-array */
  ik = (int *)calloc((size_t)*k,sizeof(int)); /* corresponding index array */
  x = (double *)calloc((size_t)*n,sizeof(double)); /* array for current point */    
  pcount=0;

  for (i=0;i < *n;i++) { /* work through all the points in X */
    for (p=X+i,p1=x,p2=p1 + *d;p1<p2;p1++, p+= *n) *p1 = *p; /* copy ith point (ith row of X) to x */
    for (p=dk,p1=dk + *k;p<p1;p++) *p = huge; /* initialize distances to huge */
    /* here I have followed Press et al. in descending tree to smallest 
       box and then re-ascending to find box with enough points. This is
       probably more efficient than checking box size all the way down 
       if n >> k .... */
   
    bi = which_box(&kd,i); /* bi is smallest box containing ith point */
 

 /*   for (j=0;j<*d;j++) 
    if (x[j]<kd.box[bi].lo[j]||x[j]>kd.box[bi].hi[j]) { 
      Rprintf("%d  ",i);
      for (j=0;j<*d;j++) Rprintf("%g  ",x[j]);
      for (j=0;j<*d;j++) Rprintf("%g  ",kd.box[bi].lo[j]);
      for (j=0;j<*d;j++) Rprintf("%g  ",kd.box[bi].hi[j]);
      Rprintf("\n");
    }

    Rprintf("%d  ",bi);*/

    while (box[bi].p1-box[bi].p0 < *k) bi = box[bi].parent; /* note k does not include self */    
    
    /*  Rprintf("Initial box %d contains %d need %d\n",bi,kd.box[bi].p1-kd.box[bi].p0+1,*k); */
    /* now find k nearest points in the box and put in dk... */    
 
    for (j=box[bi].p0;j<=box[bi].p1;j++) 
    if (ind[j]!=i) { /* avoid self! */
      pcount++;
      dij = ijdist(i,ind[j],X,*n,*d); /* distance between points i and j */ 
      if (dij<dk[0]) { /* distance smaller than top of heap */
        dk[0] = dij;       /* so replace top of distance heap */
        ik[0] = ind[j]; /* and put index on index heap */
        if (*k>1) update_heap(dk,ik,*k); /* update heap so it still obeys heap ordering */  
      } 
    } /* finished initialising heap (dk, ik) */
    
    /* Now search the rest of the tree. Basic idea is that if a box is further from the 
       ith point than dk[0] (the largest of the current neighbour distances), then we 
       can ignore all the points it contains (and hence its descendents) */ 
    todo[0] = 0; /* index of root box... first to check */
    item=0;
    bii = bi; /* index of initializing box */ 
    while (item>=0) { /* items on the todo list */
      if (todo[item]==bii) { /* this is the initializing box - already dealt with */
        item--;
      } else {
        bi = todo[item]; /* box to deal with now */
        item--;
        if (box_dist(box+bi,x,*d)<dk[0]) { /* box edge is closer than some of existing points 
                                              -- need to check further */
          if (box[bi].child1) { /* box has children --- add to todo list */
            item++;
            todo[item] = box[bi].child1;
            item++;
            todo[item] = box[bi].child2;
          } else { /* at smallest box end of tree */
            for (j=box[bi].p0;j<=box[bi].p1;j++) {
              pcount++;
              dij = ijdist(i,ind[j],X,*n,*d); /* distance between points i and j */ 
              if (dij<dk[0]) { /* point closer than largest of current candidates -- add to heap */
                dk[0] = dij; /* add distance to heap */
                ik[0] = ind[j]; /* and corresponding index to heap index */
                if (*k>1) update_heap(dk,ik,*k); /* update heap so it still obeys heap ordering */  
              } /* end of point addition */
            } /* done the one or two points in this box */
          } /* finished with this small box */
        } /* finished with possible candiate box */
      } /* end of else branch */
    } /* todo list end */
    /* So now the dk, ik contain the distances and indices of the k nearest neighbours */
    for (j=0;j<*k;j++) { /* copy to output matrices */
      dist[i + j * *n] = dk[j];
      ni[i + j * *n] = ik[j];
    }     
  } /* end of points loop (i) */
 
  free(dk);
  free(ik);
  free(x);
  *n = pcount;
}

void k_nn(double *X,double *dist,double *a,int *ni,int *n,int *d,int *k,int *get_a) {
/* NOTE: no tie handling... impractical without!

   X is an n by d matrix. Each row is the location of a point
   in some Euclidean d-space.
   Find k nearest neighbours in X of all points in X. 
   ni and dist are both n by k. each row of ni contains the neighbour list.
   Each row of dist is contains the corresponding distances. 
   if get_a is non zero, then volumes of kd boxes are associated with each point
   and returned in a.
   
   Some R test code...
   cd ~simon/mgcv-related/sparse-smooth
   R CMD SHLIB kd-tree.c  
   R 
   
   dyn.load("kd-tree.so")
   set.seed(2)
   n <- 100;d <- 2;k <- 5
   X <- matrix(runif(n*d),n,d)
   dist <- matrix(0,n,k)
   system.time(oo <- .C("k_nn",X=as.double(X),dist=as.double(dist),a=as.double(1:n),ni=as.integer(dist),
                   n=as.integer(n),d=as.integer(d),k=as.integer(k),get.a=as.integer(1)))
   oo$n/n^2 ## efficiency

   dist1 <- dist <- matrix(oo$dist,n,k)
   ni1 <- ni <- matrix(oo$ni+1,n,k)
   ## checking code...
   for (i in 1:n) {
     Xi <- t(t(X)-X[i,])
     di <- rowSums(Xi^2)^.5
     oi <- order(di)
     ni1[i,] <- (1:n)[oi[2:(k+1)]]
     dist1[i,] <- di[ni1[i,]]
     oi <- order(dist[i,])
     dist[i,] <- dist[i,oi]
     ni[i,] <- ni[i,oi]
   }
   range(ni-ni1)
   range(dist-dist1)

*/
  kdtree_type kd; 
  kd_tree(X,n,d,&kd); /* set up the tree */ 
  if (get_a) p_area(a,X,kd,*n,*d);
  k_nn_work(kd,X,dist,ni,n,d,k);
  free_kdtree(kd);
}


void kba_nn(double *X,double *dist,double *a,int *ni,int *n,int *d,int *k,
            int *get_a,double *cut_off) {
/* Obtains a roughly balanced set of 2d + k nearish neighbours. Idea is to take nearest neighbour 
   from kd box immediately above or below each point's box, in each dimension, and to add 
   k further points from the nearest neigbours not already included.  
   For each point:
   1. get 2d+k nearest neighbours.
   2. get 2d balanced neighbours.
   3. find k nearest neighbours in set 1 that are not in set 2.
   Step 3 can go through nearest looking for self and largest.  
*/
  int ii,i,j,nn,d2k,bi,bj,max_i,q,n1,n2;
  double dx,*x,max_dist,d1,d2,maxnd;
  kdtree_type kd; 
  kd_tree(X,n,d,&kd); /* set up the tree */ 
  if (get_a) p_area(a,X,kd,*n,*d);
  d2k = 2 * *d + *k;
  nn = *n; /* follwoing modifies n!!*/
  k_nn_work(kd,X,dist,ni,&nn,d,&d2k); /* get 2d+k nearest neighbours */
  x = (double *)calloc((size_t) *d,sizeof(double));
  for (i=0;i<*n;i++) { /* work through points */
    for (j=0;j<*d;j++) x[j] = X[i + j * *n];
    bi = which_box(&kd,i);
    for (j=0;j<*d;j++) { /* get the balanced neighbours */
     /* lower neighbour... */ 
     if (kd.box[bi].lo[j]!=-kd.huge) { /* then there is a neigbour in this direction */
        dx = (x[j] - kd.box[bi].lo[j])*1.001;
        x[j] -= dx;
        bj = xbox(&kd,x); /* box below bi on axis j*/
        x[j] += dx;
        /* now find point closest to point i in box bj */
      //  Rprintf("bj = %d  p0 = %d  ",bj,kd.box[bj].p0);
        n1 = kd.ind[kd.box[bj].p0];
     //   Rprintf("calling ijdist(%d,%d)...\n",i,n1);
        d1 = ijdist(i,n1,X,*n,*d);
        if (kd.box[bj].p1>kd.box[bj].p0) { 
          n2 = kd.ind[kd.box[bj].p1];
          d2 = ijdist(i,n2,X,*n,*d);
          if (d2<d1) { d1=d2;n1=n2;}
        }
        /* now put n1 into neighbour list in place of furthest neighbour not 
           itself an already computed balanced box point (later computed points 
            are no problem) */
       
        max_dist=0.0;
        max_i=0;
        maxnd=0.0;
        for (q=0;q < d2k;q++) {
          ii = i + *n * q;
          if (dist[ii]>maxnd) maxnd = dist[ii];
          if (ni[ii] == n1) { /* point is already in neighbour set */
            ni[ii] == -(n1+1);    /* signal to ignore for replacement */
            max_i = -1; /* signal that no replacement needed */
            break; 
          }
          if (ni[ii]>0&&dist[ii]>max_dist) { 
            max_dist = dist[ii];max_i=ii;
          }
        }
        if (max_i>=0 && d1 < *cut_off * maxnd) { /* replace furthest replacable item with n1 */
          ni[max_i] = -(n1+1); /* signal not to replace later */
          dist[max_i] = d1;
        }
      }
      /* upper neighbour ... */
      if (kd.box[bi].hi[j]!=kd.huge) { /* then there is a neigbour in this direction */
        dx = (kd.box[bi].hi[j] - x[j])*1.001;
        x[j] += dx;
        bj = xbox(&kd,x); /* box above bi on axis j*/
        x[j] -= dx;
        /* now get nearest point to i from box bj */
        n1 = kd.ind[kd.box[bj].p0];
        d1 = ijdist(i,n1,X,*n,*d);
        if (kd.box[bj].p1>kd.box[bj].p0) { 
          n2 = kd.ind[kd.box[bj].p1];
          d2 = ijdist(i,n2,X,*n,*d);
          if (d2<d1) { d1=d2;n1=n2;}
        }
        /* now put n1 into neighbour list in place of furthest neighbour not 
           itself an already computed balanced box point (later computed points 
           are no problem) */
      
        max_dist=0.0;
        max_i=0;
        maxnd=0.0;
        for (q=0;q < d2k;q++) {
          ii = i + *n * q;
          if (dist[ii]>maxnd) maxnd = dist[ii];
          if (ni[ii] == n1) { /* point is already in neighbour set */
            ni[ii] == -(n1+1);    /* signal to ignore for replacement */
            max_i = -1; /* signal that no replacement needed */
            break; 
          }
          if (ni[ii]>0&&dist[ii]>max_dist) { 
            max_dist = dist[ii];max_i=ii;
          }
        }
        if (max_i>=0 && d1 < *cut_off * maxnd) { /* replace furthest replacable item with n1 */
          ni[max_i] = -(n1+1); /* signal not to replace later */
          dist[max_i] = d1;
        }
      }
    } /* collected balanced neighbours */
    /* finally reset the negative indices to positive */    
    for (q=0;q < d2k;q++) {
       ii = i + *n * q;
       if (ni[ii]<0) ni[ii] = -ni[ii] - 1; 
    }
  }
  free(x); free_kdtree(kd);
}


void sparse_penalty(double *X,int *n,int *d,double *D,int *K,int *k,int *m,int *a_weight) {
/* X is n by d, and each row of X contains the location of a point. 
   There are no repeat points in X.

   D is n by m*k, where m is the number of derivatives in the penalty
     and k is the number of neighbours (*including self*) required to 
     approximate the derivatives. D[i + m*n*l + n*j] is the coefficient 
     associated with the neighbour referenced by K[i + n*j] and the lth
     derivative.   

   Set up is general to allow for future extension of this routine, but currently 
   only the d==2, m=3, k=6 TPS like case is dealt with here. 

*/
  int i,j,*ni;
  double *M, /* matrix mapping derivatives to function values */
    *dist,*area,cut_off=5;

  M = (double *)calloc((size_t) *k * *k,sizeof(double));

  (*k)--;

  ni = (int *)calloc((size_t) *n * *k,sizeof(int)); /* neighbour index list */
  dist = (double *)calloc((size_t) *n * *k,sizeof(double)); /* corresponding distances */
  area = (double *)calloc((size_t) *n,sizeof(double)); /* area associated with each point */

  /* get a balanced version of the nearest neighbours */
  kba_nn(X,dist,area,ni,n,d,k,a_weight,&cut_off);
  for (i=0;i<*n;i++) { /* work through all points */ 
    M[0] = 1.0;for (j=1;j<6;j++) M[j*6] = 0.0;
    for (j=1;i<6;i++) {
      M[j] = 1.0;
      ii = ni[i + (j-1) * *n]; /* neighnour index */
      x = X[ii] - X[i];     
      z = X[ii + *n] - X[i + *n];
      M[j + 6] = x;
      M[j + 12] = z;
      M[j + 18] = x*x/2;
      M[j + 24] = z*z/2;
      M[j + 30] = x*z/2;
    }
    /* Let g = [f,f_x,f_z,f_xx,f_zz,f_xz], then f -> Mg as neighbours
       approach point i. Now invert M, to estimate g using g = M^{-1}f */
    
    /* call mgcv_svd_full to pseudoinvert M */
  }
  

}



