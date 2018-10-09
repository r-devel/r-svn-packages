/* (c) Simon N. Wood (2015) Released under GPL2 */

/* Routines to work with discretized covariate models 
   Data structures:
   * X contains sub matrices making up each term 
   * nx is number of sub-matrices.
   * jth matrix is m[j] by p[j]
   * k contains the index vectors for each matrix in Xb.
   * ts[i] is starting matrix of ith term    
   * dt[i] is number of matrices making up ith term.
   * n_terms is the number of terms. 
   * Q contains reparameterization matrices for each term
   * qs[i] is starting address within Q for ith repara matrix
   R CMD SHLIB discrete2.c NOTE: *Must* put CFLAGS = -fopenmp -O2 in .R/Makevars 

   Here is an interesting bit of code for converting an index kk
   running over the upper triangle of an nt by nt matrix to rows 
   and columns (first row corresponds to kk = 0, 1, 2 ,...)
   i=kk;r=0;while (i >= *nt-r) { i -= *nt - r; r++;} c = r + i; 
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <Rconfig.h>
#include "mgcv.h"
#ifdef OPENMP_ON
#include <omp.h>
#endif

void Cdgemv(char *trans, int *m, int *n, double *alpha, double *a, int *lda,
	    double *x, int *incx, double *beta, double *y, int *incy) {
/* miserably vanilla implementation of dgemv to plug in in place of BLAS call for
   debugging purposes (written to examine openblas thread safety problem) */
  int i,j,q;
  double *yp,*ap,*xp;
  if (*trans == 'T') q = *n; else q =*m; /* length of y */
  if (*alpha==0) { /* matrix part does not contribute */
    for (i=0;i<q;i++,y += *incy) *y *= *beta;
    return;
  } else *beta /= *alpha;
  if (*trans == 'N') { /* y <- alpha a x + beta y */
    for (yp=y,i=0;i<*m;i++,a++,yp += *incy) *yp = *beta * *yp + *a * *x;
    x += *incx;
    for (j=1;j<*n;j++,x += *incx) 
      for (ap=a + *lda * j,yp=y,i=0;i<*m;i++,ap++,yp += *incy) *yp += *ap * *x;
  } else { /* y <- alpha a' x + beta y */
    for (yp=y,i=0;i<*n;i++,yp++) {
      *yp *= *beta;
      for (ap=a + *lda * i,xp=x,j=0;j<*m;j++,ap++,xp += *incx) *yp += *ap * *xp;
    }  
  }
  for (i=0;i<q;i++,y += *incy) *y *= *alpha;
} /* Cdgemv */ 


/* basic extraction operations */ 

void singleXj(double *Xj,double *X, int *m, int *k, int *n,int *j) {
/* Extract a column j of matrix stored in compact form in X, k into Xj. 
   X has m rows. k is of length n. ith row of result is Xj = X[k(i),j]
   (an n vector). This function is O(n). Thread safe.
*/
  double *pe; 
  X += *m * *j; /* shift to start of jth column */ 
  for (pe = Xj + *n;Xj < pe;Xj++,k++) *Xj = X[*k]; 
} /* singleXj */

void tensorXj(double *Xj, double *X, int *m, int *p,int *dt, 
              int *k, int *n, int *j, int *kstart,int *koff) {
/* Extract a column j of tensor product term matrix stored in compact 
   form in X, k into Xj. There are dt sub matrices in Xj. The ith is 
   m[i] by p[i]. There are dt index n - vectors stacked end on end in 
   k. The ith component (starting at 0) has index vector at column 
   kstart[i] + *koff of k. This function is O(n*dt)

   This routine performs pure extraction only if Xj is a vector of 1s on 
   entry. Otherwise the jth column is multiplied element wise by the 
   contents of Xj on entry. Thread safe.
*/ 
  int q=1,l,i,jp,*kp;
  double *p0,*p1,*M;
  p1 = Xj + *n; /* pointer for end of Xj */
  for (i = 0;i < *dt;i++) q *= p[i];
  jp = *j; 
  for (i = 0;i < *dt; i++) {
    q /= p[i]; /* update q */
    l = jp/q; /* column of current marginal */
    jp = jp%q;
    M = X + m[i] * l; /* M now points to start of col l of ith marginal model matrix */
    kp = k + (kstart[i] + *koff) * (ptrdiff_t) *n;
    for (p0=Xj;p0<p1;p0++,kp++) *p0 *= M[*kp];
    //k += *kjump * (ptrdiff_t) *n; /*k cols relating to different marginals may be separated by kjump cols*/
    X += m[i] * p[i]; /* move to the next marginal matrix */
  }
} /* tensorXj */

void singleXty(double *Xy,double *temp,double *y,double *X, int *m,int *p, int *k, int *n, int *add) {
/* forms X'y for a matrix stored in packed form in X, k, with ith row 
   X[k[i],]. X as supplied is m by p, while k is an n vector.
   Xy and temp are respectively p and m vectors and do not need to be cleared on entry.
   Thread safe.
*/
  double *p0,*p1,done=1.0,dzero=0.0; 
  char trans = 'T',ntrans = 'N';
  int one=1;
  for (p0=temp,p1 = p0 + *m;p0<p1;p0++) *p0 = 0.0;
  for (p1=y + *n;y<p1;y++,k++) temp[*k] += *y;
  if (*add) dzero = 1.0;
  F77_CALL(dgemv)(&trans, m, p,&done,X,m,temp,&one,&dzero,Xy,&one);
  // Cdgemv(&trans, m, p,&done,X,m,temp,&one,&dzero,Xy,&one);
  /* dgemm call equivalent to dgemv call above */ 
  //F77_CALL(dgemm)(&trans,&ntrans,p,&one,m,&done,X,m,temp,m,&dzero,Xy,m); 
} /*singleXty*/

void tensorXty(double *Xy,double *work,double *work1, double *y,double *X, 
               int *m, int *p,int *dt,int *k, int *n, int *add, int *kstart,int *koff) {
/* forms X'y for a matrix stored as a compact tensor product in X, k.
   There are dt maginal matrices packed in X, the ith being m[i] by p[i].
   y and work are n vectors. work does not need to be cleared on entry.
   work1 is an m vector where m is the dimension of the final marginal.
   Note: constraint not dealt with here. Thread safe.
   k[,kstart[i] + *koff] is the index for the ith term (starting at i=0)
*/
  int pb=1,i,j,pd;
  double *p1,*yn,*p0,*M; 
  yn = y + *n;
  M = X;
  for (i=0;i<*dt-1;i++) { 
    pb *= p[i]; 
    M += m[i] * p[i]; /* shift M to final marginal */ 
  }
  pd = p[*dt-1];
  for (i=0;i<pb;i++) {
    /* extract the ith column of M_0 * M_1 * ... * M_{d-2} * y, where 
       '*' is the row tensor product here */
    for (p0=y,p1=work;p0<yn;p0++,p1++) *p1 = *p0; /* copy y to work */ 
    j = *dt - 1; 
    tensorXj(work,X,m,p,&j,k,n,&i,kstart,koff);
    /* now form M'work */
    singleXty(Xy+i*pd,work1,work,M,m + *dt-1,&pd, k + (ptrdiff_t) *n * (kstart[j] + *koff) ,n,add);
  }
} /* tensorXty */

void singleXb(double *f,double *work,double *X,double *beta,int *k,int *m, int *p,int *n,
              int *kstart, int *kstop) {
/* Forms X beta, where beta is a p - vector, and X is stored in compact form in 
   X, k, with ith row sum_j X[k[i,j],]. X is stored as an m by p matrix. k is the matrix
   containing the indices. indices for this term start at column kstart and end 
   at column kstop-1 (of k). often there is only one j. 
   work is an m vector. Thread safe.
*/
  char trans='N';
  double done=1.0,dzero=0.0,*p1,*fp;
  int one=1,j;
  F77_CALL(dgemv)(&trans, m, p,&done,X,m,beta,&one,&dzero,work,&one);
  //Cdgemv(&trans, m, p,&done,X,m,beta,&one,&dzero,work,&one);
  /* dgemm call equivalent to dgemv call above */ 
  //F77_CALL(dgemm)(&trans,&trans,m,&one,p,&done,X,m,beta,p,&dzero,work,m); 
  p1 = f + *n;
  fp = f;
  k += *kstart * (ptrdiff_t) *n;
  for (;fp < p1;fp++,k++) *fp = work[*k];
  for (j=1;j < *kstop - *kstart;j++) for (fp=f;fp < p1;fp++,k++) *fp += work[*k];
} /* singleXb */

void tensorXb(double *f,double *X, double *C,double *work, double *beta,
              int *m, int *p,int *dt,int *k, int *n,double *v,int *qc,
              int *kstart, int *kstop) {
/* for X* beta where X* is a tensor product term with nt marginal model matrices
   stored in X. ith such is m[i] by p[i]. 
   work is an n vector. C is an m[d] by pb working matrix.  
   v is a vector such that if Q=I-vv' and Z is Q with the first column dropped then 
   Z is a null space basis for the identifiability constraint. If *qc <= 0 then
   no constraint is applied. Thread safe.
   k[,kstart[i]...kstop[i]-1] are the indices for the ith component of this term 
   (starting at i=0).
*/ 
  char trans='N';
  int pb=1,md,*kp,*kd,pd,i,j,q;
  double *M,done=1.0,dzero=0.0,*p0,*p1,*p2,*p3,*pf,*pc,x;
  M = X;
  for (i=0;i<*dt-1;i++) {
    pb *= p[i];
    M += m[i] * (ptrdiff_t) p[i]; /* final marginal */ 
  }
  md = m[*dt - 1];
  pd = p[*dt -1];
  //kd = k + (*dt-1) * (ptrdiff_t) *n; /* index vector for final term */
  kd = k + kstart[*dt-1] * (ptrdiff_t) *n;
  /* form work = M B, where vec(B) = beta */
  if (*qc<=0) { /* no constraint supplied */
    F77_CALL(dgemm)(&trans,&trans,&md,&pb, &pd, &done,
		    M,&md,beta,&pd,&dzero,C,&md); 
  } else { /* there is a constraint matrix */
    /* first map supplied beta to unconstrained parameterization */ 
    j = pb * pd; /* total number of coeffs - length of unconstrained beta */
    
    *work = 0.0;x=0.0;
    for (p0=work+1,p1=p0+j-1,p2=beta,p3=v+1;p0<p1;p0++,p2++,p3++) { 
      *p0 = *p2; 
      x += *p0 * *p3; /* v'beta where beta padded with extra zero at start */ 
    }
    for (p0=work,p1=p0+j,p2=v;p0<p1;p0++,p2++) *p0 -= *p2 * x; /* (I-vv')(0,beta')' */

    /*F77_CALL(dgemv)(&trans, &j, qc,&done,Q,&j,beta,&one,&dzero,work,&one); old when Q full matrix */
    F77_CALL(dgemm)(&trans,&trans,&md,&pb, &pd, &done,
		    M,&md,work,&pd,&dzero,C,&md);
  }
  p1 = work + *n;
  for (pf=f,p0=f + *n;pf<p0;pf++) *pf = 0.0;
  for (q = 0;q < *kstop - *kstart;q++) /* loop over possible multiple indices */
  for (i=0;i<pb;i++) {
    /* extract ith col of truncated tensor product (final marginal omitted) */
    for (p0=work;p0<p1;p0++) *p0 = 1.0; /* set work to 1 */ 
    j = *dt - 1; 
    tensorXj(work,X,m,p,&j,k,n,&i,kstart,&q);
    pc = C + i * md; /* ith col of C */
    for (kp=kd + q  *(ptrdiff_t) *n,pf=f,p0=work;p0<p1;p0++,pf++,kp++) {
      *pf += pc[*kp] * *p0;
    }
  }
} /* tensorXb */

void Xbd(double *f,double *beta,double *X,int *k,int *ks, int *m,int *p, int *n, 
	 int *nx, int *ts, int *dt, int *nt,double *v,int *qc,int *bc) {
/* Forms f = X beta for X stored in the packed form described in function XWX
   bc is number of cols of beta and f... 
   LIMITED thread safety. Allocates/frees memory using R_chk_calloc/R_chk_free (usually),
   protected within critical sections. Not safe to use with other routines allocating memory using these routines
   within a parallel section. Safe to use as the only allocator of memory within a parallel section. 
*/
  ptrdiff_t *off,*voff;
  int i,j,q,*pt,*tps,dC=0,c1,first;
  double *f0,*pf,*p0,*p1,*p2,*C=NULL,*work,maxp=0,maxrow=0;
  /* obtain various indices */
#pragma omp critical (xbdcalloc)
  { pt = (int *) CALLOC((size_t)*nt,sizeof(int)); /* the term dimensions */
    off = (ptrdiff_t *) CALLOC((size_t)*nx+1,sizeof(ptrdiff_t)); /* offsets for X submatrix starts */
    voff = (ptrdiff_t *) CALLOC((size_t)*nt+1,sizeof(ptrdiff_t)); /* offsets for v subvector starts */
    tps = (int *) CALLOC((size_t)*nt+1,sizeof(int)); /* the term starts in param vector or XWy */
  }
  for (q=i=0;i< *nt; i++) { /* work through the terms */
    for (j=0;j<dt[i];j++,q++) { /* work through components of each term */
      off[q+1] = off[q] + p[q] * (ptrdiff_t) m[q]; /* submatrix start offsets */
      if (maxrow<m[q]) maxrow=m[q];
      if (j>0 && j==dt[i]-1) {
        c1 = pt[i] * (ptrdiff_t) m[q]; 
        if (c1>dC) dC = c1; /* dimension of working matrix C */
      }
      if (j==0) pt[i] = p[q]; else pt[i] *= p[q]; /* term dimension */
    } 
    if (qc[i]>0) voff[i+1] = voff[i] + pt[i]; else voff[i+1] = voff[i]; /* start of ith v matrix */
    if (maxp < pt[i]) maxp = pt[i];
    if (qc[i]<=0) tps[i+1] = tps[i] + pt[i]; /* where ith terms starts in param vector */ 
    else tps[i+1] = tps[i] + pt[i] - 1; /* there is a tensor constraint to apply - reducing param count*/
  }
  /* now form the product term by term... */ 
  i = *n; if (i<maxp) i=maxp; if (i<maxrow) i=maxrow;
#pragma omp critical (xbdcalloc)
  { pf=f0 = (double *)CALLOC((size_t)*n,sizeof(double));
    work = (double *)CALLOC((size_t)i,sizeof(double));
    if (dC) C = (double *)CALLOC((size_t)dC,sizeof(double));
  }
  for (j=0;j < *bc;j++) { /* loop over columns of beta */
    first = 1;
    for (i=0;i < *nt;i++) { /* work through terms */  
        if (first) f0 = f; /* result written straight to f for i==0 */
        if (dt[i]==1) singleXb(f0,work,X+off[ts[i]],beta+tps[i],k /*+ (ptrdiff_t) *n *ts[i]*/,
			       m+ts[i], p+ts[i],n,ks + ts[i], ks + ts[i] + *nx); 
        else tensorXb(f0,X+off[ts[i]],C,work, beta+tps[i],m+ts[i], p+ts[i],dt+i,k 
		      /*+ (ptrdiff_t) *n * ts[i]*/,n,v+voff[i], qc+i,ks + ts[i], ks + ts[i] + *nx);
        if (!first) {
          for (p0=f,p1=f + *n,p2=f0;p0<p1;p0++,p2++) *p0 += *p2; /* f <- f + f0 */     
        } else { 
          f0=pf; /* restore f0 */
          first=0;
        }
    } /* term loop */
    f += *n;beta += tps[*nt]; /* move on to next column */
  } /* col beta loop */
#pragma omp critical (xbdcalloc)
  { if (dC) FREE(C);
    FREE(work);FREE(f0);
    FREE(pt);FREE(off);FREE(voff);FREE(tps);
  }
} /* Xb */

void diagXVXt(double *diag,double *V,double *X,int *k,int *ks,int *m,int *p, int *n, 
	      int *nx, int *ts, int *dt, int *nt,double *v,int *qc,int *pv,int *nthreads) {
/* Forms diag(XVX') where X is stored in the compact form described in XWXd.
   V is a pv by pv matrix. 
   Parallelization is by splitting the columns of V into nthreads subsets.
   Currently inefficient. Could be speeded up by a factor of 2, by supplying a
   square root of V in place of V.
*/
  double *xv,*dc,*p0,*p1,*p2,*p3,*ei,*xi;
  ptrdiff_t bsj,bs,bsf,i,j,kk;
  int one=1;
  #ifndef OPENMP_ON
  *nthreads = 1;
  #endif
  if (*nthreads<1) *nthreads = 1;
  if (*nthreads > *pv) *nthreads = *pv;
  xv = (double *) CALLOC((size_t) *nthreads * *n,sizeof(double)); /* storage for cols of XV */
  xi = (double *) CALLOC((size_t) *nthreads * *n,sizeof(double)); /* storage for cols of X */
  ei = (double *) CALLOC((size_t) *nthreads * *pv,sizeof(double)); /* storage for identity matrix cols */
  dc = (double *) CALLOC((size_t) *nthreads * *n,sizeof(double)); /* storage for components of diag */
  if (*nthreads>1) {
    bs = *pv / *nthreads;
    while (bs * *nthreads < *pv) bs++;
    while (bs * *nthreads - bs >= *pv) (*nthreads)--;
    bsf = *pv - (bs * *nthreads - bs); 
  } else {
    bsf = bs = *pv;
  }
  #ifdef OPENMP_ON
  #pragma omp parallel for private(j,bsj,i,kk,p0,p1,p2,p3) num_threads(*nthreads)
  #endif  
  for (j=0;j < *nthreads;j++) {
    if (j == *nthreads - 1) bsj = bsf; else bsj = bs;
    for (i=0;i<bsj;i++) { /* work through this block's columns */
      kk = j * bs + i;
      ei[j * *pv + kk] = 1;if (i>0) ei[j * *pv + kk - 1] = 0;
      /* Note thread safety of XBd means this must be only memory allocator in this section*/
      Xbd(xv + j * *n,V + kk * *pv,X,k,ks,m,p,n,nx,ts,dt,nt,v,qc,&one); /* XV[:,kk] */
      Xbd(xi + j * *n,ei + j * *pv,X,k,ks,m,p,n,nx,ts,dt,nt,v,qc,&one); /* X[:,kk] inefficient, but deals with constraint*/
      p0 = xi + j * *n;p1=xv + j * *n;p2 = dc + j * *n;p3 = p2 + *n;
      for (;p2<p3;p0++,p1++,p2++) *p2 += *p0 * *p1; /* elementwise product of XV[:,kk] X[:,kk] */
    } 
  } /* parallel loop end */
  /* sum the contributions from the different threads into diag... */
  for (p0=diag,p1=p0+ *n,p2=dc;p0<p1;p0++,p2++) *p0 = *p2;
  for (i=1;i< *nthreads;i++) for (p0=diag,p1=p0+ *n;p0<p1;p0++,p2++) *p0 += *p2;
  FREE(xv);FREE(dc);FREE(xi);FREE(ei);
} /* diagXVXt */


void XWyd(double *XWy,double *y,double *X,double *w,int *k,int *ks, int *m,int *p, int *n, 
         int *nx, int *ts, int *dt, int *nt,double *v,int *qc,
         int *ar_stop,int *ar_row,double *ar_weights) {
/* NOT thread safe */
  double *Wy,*p0,*p1,*p2,*p3,*Xy0,*work,*work1,x;
  ptrdiff_t i,j,*off,*voff;
  int *tps,maxm=0,maxp=0,one=1,zero=0,*pt,add,q;
  if (*ar_stop>=0) { /* model has AR component, requiring sqrt(weights) */
    for (p0 = w,p1 = w + *n;p0<p1;p0++) *p0 = sqrt(*p0);
  }
  /* obtain various indices */
  pt = (int *) CALLOC((size_t)*nt,sizeof(int)); /* the term dimensions */
  off = (ptrdiff_t *) CALLOC((size_t)*nx+1,sizeof(ptrdiff_t)); /* offsets for X submatrix starts */
  voff = (ptrdiff_t *) CALLOC((size_t)*nt+1,sizeof(ptrdiff_t)); /* offsets for v subvector starts */
  tps = (int *) CALLOC((size_t)*nt+1,sizeof(int)); /* the term starts in param vector or XWy */
  for (q=i=0;i< *nt; i++) { /* work through the terms */
    for (j=0;j<dt[i];j++,q++) { /* work through components of each term */
      off[q+1] = off[q] + p[q]*m[q]; /* submatrix start offsets */
      if (j==0) pt[i] = p[q]; else pt[i] *= p[q]; /* term dimension */
      if (maxm<m[q]) maxm=m[q];
    } 
    if (qc[i]>0) voff[i+1] = voff[i] + pt[i]; else voff[i+1] = voff[i];  /* start of ith Q matrix */
    if (maxp < pt[i]) maxp=pt[i];
    if (qc[i]<=0) tps[i+1] = tps[i] + pt[i]; /* where ith terms starts in param vector */ 
    else tps[i+1] = tps[i] + pt[i] - 1; /* there is a tensor constraint to apply - reducing param count*/
  }
  Xy0 =  (double *) CALLOC((size_t)maxp,sizeof(double));
  work =  (double *) CALLOC((size_t)*n,sizeof(double));
  work1 = (double *) CALLOC((size_t)maxm,sizeof(double));
  /* apply W to y */
  Wy = (double *) CALLOC((size_t)*n,sizeof(double)); /* Wy */
  for (p0=Wy,p1=Wy + *n,p2=w;p0<p1;p0++,y++,p2++) *p0 = *y * *p2; 
  if (*ar_stop>=0) { /* AR components present (weights are sqrt, therefore) */
    rwMatrix(ar_stop,ar_row,ar_weights,Wy,n,&one,&zero,work);
    rwMatrix(ar_stop,ar_row,ar_weights,Wy,n,&one,&one,work); /* transpose of transform applied */
    for (p0=w,p1=w + *n,p2=Wy;p0<p1;p0++,p2++) *p2 *= *p0; /* sqrt weights again */
  }
  /* now loop through terms applying the components of X'...*/
  for (i=0;i<*nt;i++) { /* term loop */ 
   add=0; 
   if (dt[i]>1) { /* it's a tensor */
      //tensorXty(Xy0,work,work1,Wy,X+off[ts[i]],m+ts[i],p+ts[i],dt+i,k+ts[i] * (ptrdiff_t) *n,n);
      for (q=0;q<ks[ts[i] + *nx]-ks[ts[i]];q++) {  /* loop through index columns */
        tensorXty(Xy0,work,work1,Wy,X+off[ts[i]],m+ts[i],p+ts[i],dt+i,k,n,&add,ks+ts[i],&q);
        add=1;
      }
      if (qc[i]>0) { /* there is a constraint to apply Z'Xy0: form Q'Xy0 and discard first row... */
        /* Q' = I - vv' */
        for (x=0.0,p0=Xy0,p1=p0 + pt[i],p2=v+voff[i];p0<p1;p0++,p2++) x += *p0 * *p2; /* x = v'Xy0 */
        p0=XWy + tps[i];p1 = p0 + pt[i]-1;p2 = v+voff[i] + 1;p3=Xy0+1;
        for (;p0<p1;p0++,p2++,p3++) *p0 = *p3 - x * *p2; /* (I-vv')Xy0 less first element */
      } else { /* straight copy */
        for (p0=Xy0,p1=p0+pt[i],p2=XWy+tps[i];p0<p1;p0++,p2++) *p2 = *p0;
      }
    } else { /* it's a singleton */
      //singleXty(XWy+tps[i],work1,Wy,X+off[ts[i]], m+ts[i],p+ts[i], k+ts[i] * (ptrdiff_t) *n,n);
      for (q=ks[ts[i]];q<ks[ts[i] + *nx];q++) { /* loop through index columns */  
        singleXty(XWy+tps[i],work1,Wy,X+off[ts[i]], m+ts[i],p+ts[i], k + q * (ptrdiff_t) *n,n,&add);
        add=1;
      }
    }
  } /* term loop */
  FREE(Wy); FREE(Xy0); FREE(work); FREE(work1); 
  FREE(pt); FREE(off); FREE(voff); FREE(tps);
} /* XWy */

void XWXd(double *XWX,double *X,double *w,int *k,int *ks, int *m,int *p, int *n, int *nx, int *ts, 
          int *dt, int *nt,double *v,int *qc,int *nthreads,int *ar_stop,int *ar_row,double *ar_weights) {
/* This version uses open MP by column, within sub-block.

   An alternative in which one thread did each sub-block scaled poorly, largely 
   because the work is very uneven between sub-blocks. 

   Forms Xt'WXt when Xt is divided into blocks of columns, each stored in compact form
   using arguments X and k. 
   * 'X' contains 'nx' blocks, the ith is an m[i] by p[i] matrix containing the unique rows of 
     the ith marginal model matrix. 
   * There are 'nt' model terms. Each term is made up of one or more maginal model matrices.
   * The jth term starts at block ts[j] of X, and has dt[j] marginal matrices. The terms model
     matrix is the row tensor product of its full (n row) marginals.
   * The index vectors converting the unique row matrices to full marginal matrices are in 
     'k', an n-row matrix of integers. Conceptually if Xj and kj represent the jth unique 
      row matrix and index vector then the ith row of the corresponding full marginal matrix 
      is Xj[kj[i],], but things are more complicated when each full term matrix is actually 
      the sum of several matrices (summation convention).
   * To handle the summation convention, each marginal matrix can have several index vectors. 
     'ks' is an nx by 2 matrix giving the columns of k corresponding to the ith marginal 
     model matrix. Specifically columns ks[i,1]:(ks[i,2]-1) of k are the index vectors for the ith 
     marginal. All marginals corresponding to one term must have the same number of index columns.
     The full model matrix for the jth term is constucted by summing over q the full 
     model matrices corresponding to the qth index vectors for each of its marginals.    
   * For example the exression for the full model matrix of the jth term is...
  
     X^full_j = sum_q prod_i X_{ts[j]+i}[k[,ks[i]+q],]  

     - q runs from 0 to ks[i,2] - ks[i,1] - 1; i runs from 0 to dt[j] - 1.
         
   Tensor product terms may have constraint matrices Z, which post multiply the tensor product 
   (typically imposing approximate sum-to-zero constraints). Actually Z is Q with the first column 
   dropped where Q =  I - vv'. qc[i]==0 for singleton terms.  

   AR models are handled via the 3 ar_* arrays. ar_stop[0] < 0 signals no AR. 
  
*/  
  int r,c,i,j,q,*pt,*pd,a,b,*tps,ptot,maxp=0,maxm=0,pa,pb,kk,dk,*start,one=1,zero=0,add; 
  ptrdiff_t *off,*voff;
  double *p0,*p1,*p2, *Xi, *Xj, *temp,*tempn,*xwx,*xwx0,
    *XiB,*XjB,*tempB,*tempnB,*x0,*x1,x;
  #ifndef OPENMP_ON
  *nthreads = 1;
  #endif
  if (*nthreads<1) *nthreads = 1;
  start = (int *) CALLOC((size_t)(*nthreads+1),sizeof(int));
  XiB = (double *) CALLOC((size_t)(*n * *nthreads),sizeof(double)); /* column of X */ 
  XjB = (double *) CALLOC((size_t)(*n * *nthreads),sizeof(double)); /* column of X */
  tempnB = (double *) CALLOC((size_t)(*n * *nthreads),sizeof(double));
  pt = (int *) CALLOC((size_t)*nt,sizeof(int)); /* the term dimensions */
  pd = (int *) CALLOC((size_t)*nt,sizeof(int)); /* storage for last marginal size */
  off = (ptrdiff_t *) CALLOC((size_t)*nx+1,sizeof(ptrdiff_t)); /* offsets for X submatrix starts */
  voff = (ptrdiff_t *) CALLOC((size_t)*nt+1,sizeof(ptrdiff_t)); /* offsets for v subvectors starts */
  tps = (int *) CALLOC((size_t)*nt+1,sizeof(int)); /* the term starts */
  if (*ar_stop>=0) { /* model has AR component, requiring sqrt(weights) */
    for (p0 = w,p1 = w + *n;p0<p1;p0++) *p0 = sqrt(*p0);
  }
  for (q=i=0;i< *nt; i++) { /* work through the terms */
    for (j=0;j<dt[i];j++,q++) {
      if (j==dt[i]-1) pd[i] = p[q]; /* the relevant dimension for deciding product ordering */
      off[q+1] = off[q] + p[q] * (ptrdiff_t) m[q]; /* submatrix start offsets */
      if (j==0) pt[i] = p[q]; else pt[i] *= p[q]; /* term dimension */
      if (maxm<m[q]) maxm=m[q];
    } 
    if (qc[i]>0) voff[i+1] = voff[i] + pt[i]; else voff[i+1] = voff[i]; /* start of ith v vector */
    if (maxp<pt[i]) maxp=pt[i];
    if (qc[i]<=0) tps[i+1] = tps[i] + pt[i]; /* where ith terms starts in param vector */ 
    else tps[i+1] = tps[i] + pt[i] - 1; /* there is a tensor constraint to apply - reducing param count*/
  }
  tempB = (double *) CALLOC((size_t) maxm * *nthreads,sizeof(double));
  xwx = (double *) CALLOC((size_t) maxp * maxp,sizeof(double)); /* working cross product storage */
  xwx0 = (double *) CALLOC((size_t) maxp * maxp,sizeof(double)); /* working cross product storage */
  ptot = tps[*nt]; /* total number of parameters */
  for (r=0;r < *nt;r++) for (c=r;c< *nt;c++) { /* the block loop */
    if (pd[r]>pd[c]) { /* Form Xr'WXc */
      a=r;b=c;
    } else { /* Form Xc'WXr */
      a=c;b=r; 
    }
    /* split cols between threads... */  
    dk = pt[b] / *nthreads; //rk = pt[b] % *nthreads;
    if (dk * *nthreads < pt[b]) dk++;
    start[0]=0; 
    for (i=0;i<*nthreads;i++) { 
      start[i+1] = start[i] + dk;
      if (start[i+1]>pt[b]) start[i+1]=pt[b];
    }
    #ifdef OPENMP_ON
    #pragma omp parallel private(Xi,Xj,i,q,add,temp,tempn,p0,p1,p2) num_threads(*nthreads)
    #endif 
    { /* begin parallel section */
      #ifdef OPENMP_ON
      #pragma omp for
      #endif
      for (kk=0;kk<*nthreads;kk++) { 
        /* allocate thread specific storage... */
        temp = tempB + kk * (ptrdiff_t) maxm;
        Xi = XiB + kk * (ptrdiff_t) *n; Xj = XjB + kk * (ptrdiff_t) *n;
        tempn = tempnB + kk * (ptrdiff_t) *n;
        for (i=start[kk];i<start[kk+1];i++) { /* loop over columns of Xb */ 
          /* extract Xb[:,i]... */
          if (ks[ts[b]]==ks[ts[b] + *nx]-1) { /* single index column */
            if (dt[b]>1) { /* tensor */
              for (p0=Xi,p1=p0+*n;p0<p1;p0++) *p0 = 1.0; /* set Xi to 1 */
	      tensorXj(Xi,X+off[ts[b]], m + ts[b] , p + ts[b], dt+b, 
	      		    k, n,&i,ks + ts[b],&zero);
            } else { /* singleton */
              singleXj(Xi,X+off[ts[b]], m + ts[b], k + (ptrdiff_t)*n * ks[ts[b]] /* ts[b]*/, n,&i);
            }
          } else { /* need to sum over multiple cols of k */
            for (q = 0;q<ks[ts[b] + *nx]-ks[ts[b]];q++) { /* k column loop */
              if (dt[b]>1) { /* tensor */
                for (p0=Xj,p1=p0+*n;p0<p1;p0++) *p0 = 1.0; /* set Xj to 1 */
                tensorXj(Xj,X+off[ts[b]], m + ts[b] , p + ts[b], dt+b, 
			    k, n,&i,ks+ts[b],&q);
              } else { /* singleton */
                singleXj(Xj,X+off[ts[b]], m + ts[b], k + (ptrdiff_t)*n * (q + ks[ts[b]]) /* ts[b]*/, n,&i);
              }
              if (q == 0) for (p0=Xj,p1=p0+*n,p2=Xi;p0<p1;p0++,p2++) *p2 = *p0;
              else for (p0=Xj,p1=p0+*n,p2=Xi;p0<p1;p0++,p2++) *p2 += *p0;
            } /* k column loop */
          }
          /* apply W to Xi - for AR models this is sqrt(W) */
          for (p0=w,p1=w + *n,p2=Xi;p0<p1;p0++,p2++) *p2 *= *p0; 
          if (*ar_stop>=0) { /* AR components present (weights are sqrt, therefore) */
            rwMatrix(ar_stop,ar_row,ar_weights,Xi,n,&one,&zero,tempn);
            rwMatrix(ar_stop,ar_row,ar_weights,Xi,n,&one,&one,tempn); /* transpose of transform applied */
            for (p0=w,p1=w + *n,p2=Xi;p0<p1;p0++,p2++) *p2 *= *p0; /* sqrt weights again */
          }
          /* now form Xa'WXb[:,i]... */ 
          add=0; 
          for (q=0;q<ks[ts[a] + *nx]-ks[ts[a]];q++) {
            if (dt[a]>1) { /* tensor */
              tensorXty(xwx + i * pt[a],tempn,temp,Xi,X+off[ts[a]],m+ts[a],p+ts[a],
	      			    dt+a,k, n,&add,ks+ts[a],&q);
            } else { /* singleton */
              singleXty(xwx + i * pt[a],temp,Xi,X+off[ts[a]],m+ts[a],p+ts[a],k + (ptrdiff_t)*n * (q + ks[ts[a]]),n,&add);
            }
            add = 1; /* for q>0 accumulate result */
          }  
        } /* loop over columns of Xb */
      }
      /* so now xwx contains pt[a] by pt[b] matrix Xa'WXb */
    } /* end parallel section */

    /* if Xb is tensor, may need to apply constraint */
    if (dt[a]>1&&qc[a]>0) { /* first term is a tensor with a constraint */
      x0=x1=xwx; /* pointers to columns of xwx */ 
      /* col by col form (I-vv')xwx, dropping first row... */
      for (j=0;j<pt[b];j++) {
        for (x=0.0,p0=x0,p1=x0+pt[a],p2=v+voff[a];p0<p1;p0++,p2++) x += *p0 * *p2;
        for (p2=v+voff[a]+1,p1=x0+pt[a],x0++;x0<p1;x1++,x0++,p2++) *x1 = *x0 - *p2 *x;
      }
      pa = pt[a] -1;
    } else pa = pt[a];
    if (dt[b]>1&&qc[b]>0) { /* second term is a tensor with a constraint */
      /* copy xwx to xwx0 */
      for (p0=xwx,p1=p0 + pt[b] * (ptrdiff_t) pa,p2=xwx0;p0<p1;p0++,p2++) *p2 = *p0;
      /* row by row form xwx(I-vv') dropping first col... */
      for (j=0;j<pa;j++) {
        for (x=0.0,p0=xwx0+j,p1=v+voff[b],p2=p1+pt[b];p1<p2;p0 += pa,p1++) x += *p0 * *p1; 
        x1 = xwx + j;
        for (p0=xwx0+j+pa,p1=v+voff[b],p2=p1+pt[b],p1++;p1<p2;x1+=pa,p1++,p0+=pa) *x1 = *p0 - *p1 *x; 
      }
      pb=pt[b]-1;
    } else pb=pt[b];
    /* copy result into overall XWX*/  
  
    if (pd[r]>pd[c]) { /* xwx = Xr'WXc */
      for (i=0;i<pa;i++) for (j=0;j<pb;j++) 
      XWX[i+tps[r]+(j+tps[c]) * (ptrdiff_t)ptot] = XWX[j+tps[c]+(i+tps[r]) * (ptrdiff_t)ptot] 
                                                 = xwx[i + pa * (ptrdiff_t)j];
    } else { /* Form xwx = Xc'WXr */
      for (i=0;i<pa;i++) for (j=0;j<pb;j++) 
      XWX[i+tps[c]+(j+tps[r])*(ptrdiff_t)ptot] = XWX[j+tps[r]+(i+tps[c])*(ptrdiff_t)ptot] 
                                               = xwx[i + pa * (ptrdiff_t)j];
    }
  } /* end of block loop */
  FREE(start);
  FREE(XiB); FREE(XjB); FREE(tempnB); FREE(pt); FREE(pd); FREE(off);
  FREE(voff); FREE(tps); FREE(tempB); FREE(xwx); FREE(xwx0);
} /* XWXd */



void XWXij0(double *XWX,int i,int j, double *X,int *k, int *ks, int *m, int *p,int nx,int n,int *ts, int *dt,
	   int nt, double *w,double *ws, int tri,ptrdiff_t *off,double *work) {
/* Forms product X_i'WX_j where X_i and X_j are stored in compact form. 
   * W = diag(w) if tri!=0 and is tridiagonal otherwise, with super in ws and sub in ws + n-1;     
   * off[i] is offset to start of ith matrix (out of nx)  
   If Xk is a tensor product term, let dXk denote row tensor product of all but it's final 
   marginal. 
   
   Workspace: Let mi and mj index the final marginal of term i and j, then the dim of work needed is
              2*n + 3*m[mi] + I(n>m[mi]*m[mj])*m[mi]*m[mj] + max(m[mi]*p[mj],m[mj]*p[mi])

   ORIGINAL

*/
  int si,sj,ri,rj,jm,im,r,c,kk,ddt,ddtj,koff,*K,*Ki,*Kj,nxwx,pim,pjm,
    ii,jj,rfac,q,t,s,ddti,tensi,tensj,acc_w,alpha;
  ptrdiff_t mim,mjm; /* avoid integer overflow in large pointer calculations */ 
  double x,*wl,*dXi,*dXj,*pdXj,*Xt,*Xi,*Xj,*Cq,*Dq,done=1.0,dzero=0.0,
         *C,*D,*W,*wb,*wbs,*wbl;
  char trans = 'T',ntrans = 'N';
  si = ks[ts[i]+nx]-ks[ts[i]]; /* number of terms in summation convention for i */
  if (tri) wl = ws + n - 1; /* sub-diagonal */
  /* compute number of columns in dXi/ number of rows of blocks in product */
  for (ri=1,kk=ts[i];kk<ts[i]+dt[i]-1;kk++) ri *= p[kk];
  im = ts[i]+dt[i]-1; /* the index of the final marginal for term i */
  mim = (ptrdiff_t) m[im];
  nxwx = p[im]*ri; /* rows of XWX = cols of Xi */
  /* Allocate work space for dXi(n), dXj(n) and initialze pdXj*/
  dXi = work;work += n;pdXj = dXj = work;work += n;
  if (dt[i]==1&&dt[j]==1&&m[ts[i]]==n&&m[ts[j]]==n) { /* both sub matrices are dense */
    //Rprintf("de ");
    jm = ts[j];
    mjm = (ptrdiff_t) m[jm];
    pim = p[im];pjm = p[jm];
    for (r=0;r<pim;r++) {
      Xi = X + off[im] + mim * r;
      if (i==j) for (c=r;c<pjm;c++) { /* symmetric case */
	Xj = X + off[jm] + mjm * c;  
	if (tri) {
	  x = Xi[0]*(Xj[0]*w[0]+Xj[1]*ws[0]);
	  for (kk=1;kk<n-1;kk++) x += Xi[kk]*(wl[kk-1]*Xj[kk-1] + w[kk]*Xj[kk] + ws[kk]*Xj[kk+1]);
          x += Xi[n-1]*(wl[n-2]*Xj[n-2]+w[n-1]*Xj[n-1]);				   
	} else for (x=0.0,kk=0;kk<n;kk++) x += Xi[kk]*Xj[kk]*w[kk];
	XWX[c + nxwx * r] = XWX[r + nxwx * c] = x;
      }	else for (c=0;c<pjm;c++) {
        Xj = X + off[jm] + mjm * c;  
	if (tri) {
	  x = Xi[0]*(Xj[0]*w[0]+Xj[1]*ws[0]);
	  for (kk=1;kk<n-1;kk++) x += Xi[kk]*(wl[kk-1]*Xj[kk-1] + w[kk]*Xj[kk] + ws[kk]*Xj[kk+1]);
          x += Xi[n-1]*(wl[n-2]*Xj[n-2]+w[n-1]*Xj[n-1]);				   
	} else for (x=0.0,kk=0;kk<n;kk++) x += Xi[kk]*Xj[kk]*w[kk];
	XWX[r + nxwx * c] = x;
      }	  
    }
  } else if (!tri && i==j && si==1) { /* simplest setup - just accumulate diagonal */
    //Rprintf("wb ");
    /* Allocate space for wb(m[im]), wbs and wbl*/
    wb = work; work += mim;// wbs = work; work += m[i]; wbl = work; work += m[i];
    
    if (dt[i]>1) { /* tensor */
      ddt = dt[i]-1; /* number of marginals, exluding final */
      koff = 0; /* only one index vector per marginal, so no offset */
    }  
    for (r=0;r<ri;r++) { /* loop over block rows of product */
      if (dt[i]>1) { /* extract col r of dXi */
	for (kk=0;kk<n;kk++) dXi[kk] = 1.0;
	tensorXj(dXi, X + off[ts[i]], m + ts[i], p + ts[i],&ddt, 
		 k, &n, &r, ks + ts[i],&koff);
      }
      for (c=r;c<ri;c++) { /* block column loop - note symmetry */
	if (dt[i]>1) { /* extract col c of dXi */
          if (r!=c) {
            dXj=pdXj;for (kk=0;kk<n;kk++) dXj[kk] = 1.0;
	    tensorXj(dXj, X + off[ts[i]], m + ts[i], p + ts[i],&ddt, 
		     k, &n, &c, ks + ts[i],&koff);
	  } else dXj = dXi;
	}
	/* Clear work space to zero... */
	for (ii=0;ii<mim;ii++) wb[ii]=0.0;
        //if (tri) for (ii=0;ii<m[i]-1;ii++) wbs[ii] = wbl[ii] = 0.0; 
	K = k + ks[im] * n; /* index for final margin */
	/* Accumulate the weights ... */
	if (dt[i]>1) {
	  for (kk=0;kk<n;kk++) wb[K[kk]] += dXi[kk]*dXj[kk]*w[kk];
	} else { /* singleton */
	  for (kk=0;kk<n;kk++) wb[K[kk]] += w[kk]; 
	}
	/* Now form the Xi'WXi... */
	Xt = X + off[im]; /* final marginal model matrix */
	pim = p[im]; /* cols of Xt */
	if (r!=c)  for (jj=0;jj<pim;jj++) for (ii=0;ii<pim;ii++) { /* block need not be symmetric */
	   Xi = Xt + mim * ii; /* iith col of Xt */
           Xj = Xt + mim * jj; /* jjth col of Xt */
	   for (x=0.0,kk=0;kk<mim;kk++) x += wb[kk]*Xj[kk]*Xi[kk];
	   XWX[c*pim+jj+(r*pim +ii)*nxwx] = XWX[r*pim+ii+(c*pim +jj)*nxwx] = x; 
	} else for (ii=0;ii<pim;ii++) for (jj=ii;jj<pim;jj++) { /* diagonal and symmetric */ 
          Xi = Xt + mim * ii; /* iith col of Xt */
          Xj = Xt + mim * jj; /* jjth col of Xt */
	  for (x=0.0,kk=0;kk<mim;kk++) x += wb[kk]*Xj[kk]*Xi[kk];
	  XWX[r*pim+ ii + (c*pim+jj)*nxwx] = XWX[r*pim + jj + (c*pim + ii)*nxwx] =
	  XWX[(r*pim+ ii)*nxwx + c*pim+jj] = XWX[(r*pim + jj)*nxwx + c*pim + ii] = x;
	}					  
      } /* block col loop */ 
    }
  } else { /* general case */
    sj = ks[ts[j]+nx]-ks[ts[j]]; /* number of terms in summation convention for j */
    for (rj=1,kk=ts[j];kk<ts[j]+dt[j]-1;kk++) rj *= p[kk];
    jm = ts[j]+dt[j]-1; /* the index of the final marginal for term j */
    ddti = dt[i]-1;ddtj = dt[j]-1; /* number of marginals excluding final */
    if (ddti) tensi = 1; else tensi = 0; /* is term i a tensor? */
    if (ddtj) tensj = 1; else tensj = 0;
    mjm = (ptrdiff_t) m[jm];
    //Rprintf("n = %d mim = %ld mjm = %ld ",n,mim,mjm);
    if (n>mjm*mim) acc_w = 1; else acc_w = 0; /* accumulate \bar W or \bar W X_j / \bar W'X_i)? */
    if (acc_w) {
      if (p[im]*mim*mjm + p[im]*p[jm]*mjm > mim*mjm*p[jm] + p[im]*p[jm]*mim) rfac=0; else rfac=1;
      /* Allocate storage for W (mim*mjm) */
      W = work; work += mim*mjm;
    } else {
      /* now establish whether to form left product, D, or right product C */
      if (tensi) ii = 2; else ii = 1;if (tensj) ii++;
      if (tri) alpha = ii*3 + 3; else alpha = ii + 1; /* ~ ops per iteration of accumulation loop */ 
      if (alpha*si*sj*n*p[im]+mjm*p[im]*p[jm]<alpha*si*sj*n*p[jm]+mim*p[im]*p[jm]) rfac = 0; else rfac=1;
      if (mim == n) rfac = 0; else if (mjm == n) rfac = 1; /* make absolutely sure we do not form n by p*m product */
    } /* end of accumulation storage allocation */
    if (rfac) {
      /* Allocate storge for C mim by p[jm] */
      C = work; work += mim * p[jm];
    } else {
      /* Allocate storage for D mjm by p[im] */
      D = work; work += mjm * p[im]; 
    }	
    for (r=0;r<ri;r++) { /* loop over block rows of product */
      for (c=0;c<rj;c++) { /* block column loop*/
	if (acc_w) for (kk=0;kk<mim*mjm;kk++) W[kk] = 0.0; /* clear W */
	else if (rfac) for (kk=0;kk<mim*p[jm];kk++) C[kk] = 0.0; /* clear C */
	else for (kk=0;kk<mjm*p[im];kk++) D[kk] = 0.0; /* clear D */
        for (s=0;s<si;s++) for (t=0;t<sj;t++) { /* summation convention loop */
	  if (tensi) { /* extract col r of dXi according to sth set of index vectors */
	    for (kk=0;kk<n;kk++) dXi[kk] = 1.0;
	    tensorXj(dXi, X + off[ts[i]], m + ts[i], p + ts[i],&ddti, 
		 k, &n, &r, ks + ts[i],&s);
          }
	  if (tensj) { /* extract col c of dXj according to tth set of index vectors */
	    for (kk=0;kk<n;kk++) dXj[kk] = 1.0;
            tensorXj(dXj, X + off[ts[j]], m + ts[j], p + ts[j],&ddtj, 
		 k, &n, &c, ks + ts[j],&t);
	  }
	  Ki = k + (ks[im]+s) * n; /* index for final margin of i */
          Kj = k + (ks[jm]+t) * n; /* index for final margin of j */
	  if (acc_w) { /* weight accumulation */
	    //Rprintf("W ");
            if (tensi&&tensj) { /* i and j are tensors */
	      for (kk=0;kk<n;kk++) W[Ki[kk]+mim*Kj[kk]] += w[kk]*dXi[kk]*dXj[kk];
	      if (tri) for (kk=0;kk<n-1;kk++) {
		  W[Ki[kk]+mim*Kj[kk+1]] += ws[kk]*dXi[kk]*dXj[kk+1];
		  W[Ki[kk+1]+mim*Kj[kk]] += wl[kk]*dXi[kk+1]*dXj[kk];
	      }	  
	    } else if (tensi) { /* only i is tensor */
              for (kk=0;kk<n;kk++) W[Ki[kk]+mim*Kj[kk]] += w[kk]*dXi[kk];
	      if (tri) for (kk=0;kk<n-1;kk++) {
		  W[Ki[kk]+mim*Kj[kk+1]] += ws[kk]*dXi[kk];
		  W[Ki[kk+1]+mim*Kj[kk]] += wl[kk]*dXi[kk+1];
	      }	 
	    } else if (tensj) { /* only j is tensor */
	      for (kk=0;kk<n;kk++) W[Ki[kk]+mim*Kj[kk]] += w[kk]*dXj[kk];
	      if (tri) for (kk=0;kk<n-1;kk++) {
		  W[Ki[kk]+mim*Kj[kk+1]] += ws[kk]*dXj[kk+1];
		  W[Ki[kk+1]+mim*Kj[kk]] += wl[kk]*dXj[kk];
	      }	  
	    } else { /* i and j are singletons */
              for (kk=0;kk<n;kk++) W[Ki[kk]+mim*Kj[kk]] += w[kk];
	      if (tri) for (kk=0;kk<n-1;kk++) {
		  W[Ki[kk]+mim*Kj[kk+1]] += ws[kk];
		  W[Ki[kk+1]+mim*Kj[kk]] += wl[kk];
	      }	
	    }	    
	  } else if (rfac) { /* direct accumulation of right factor C (mim by p[jm]) */
	    //Rprintf("C ");
            for (q=0;q<p[jm];q++) {
	      Cq = C + q * mim;Xj = X + off[jm] + q * mjm; /* qth cols of C and X */ 
	      if (tensi&&tensj) { /* X_i and X_j are tensors */
	        for (kk=0;kk<n;kk++) Cq[Ki[kk]] += w[kk]*dXi[kk]*dXj[kk]*Xj[Kj[kk]];
	        if (tri) for (kk=0;kk<n-1;kk++) {
		  Cq[Ki[kk]] += ws[kk]*dXi[kk]*dXj[kk+1]*Xj[Kj[kk+1]];
		  Cq[Ki[kk+1]] += wl[kk]*dXi[kk+1]*dXj[kk]*Xj[Kj[kk]];
	        } 
	      } else if (tensi) { /* only X_i is tensor */
		for (kk=0;kk<n;kk++) Cq[Ki[kk]] += w[kk]*dXi[kk]*Xj[Kj[kk]];
	        if (tri) for (kk=0;kk<n-1;kk++) {
		  Cq[Ki[kk]] += ws[kk]*dXi[kk]*Xj[Kj[kk+1]];
		  Cq[Ki[kk+1]] += wl[kk]*dXi[kk+1]*Xj[Kj[kk]];
	        }
	      } else if (tensj) { /* only X_j is tensor */
	        for (kk=0;kk<n;kk++) Cq[Ki[kk]] += w[kk]*dXj[kk]*Xj[Kj[kk]];
	        if (tri) for (kk=0;kk<n-1;kk++) {
		  Cq[Ki[kk]] += ws[kk]*dXj[kk+1]*Xj[Kj[kk+1]];
		  Cq[Ki[kk+1]] += wl[kk]*dXj[kk]*Xj[Kj[kk]];
	        }
	      } else { /* both singletons */
		for (kk=0;kk<n;kk++) Cq[Ki[kk]] += w[kk]*Xj[Kj[kk]];
	        if (tri) for (kk=0;kk<n-1;kk++) {
		  Cq[Ki[kk]] += ws[kk]*Xj[Kj[kk+1]];
		  Cq[Ki[kk+1]] += wl[kk]*Xj[Kj[kk]];
	        }
	      }
	    } /* q loop */
	  } else { /* direct accumulation of left factor D (mjm by p[im]) = \bar W' X_i */
	    //Rprintf("D ");
            for (q=0;q<p[im];q++) {
	      Dq = D + q * mjm;Xi = X + off[im] + q * mim; /* qth cols of C and Xi */ 
	      if (tensi&&tensj) { /* X_i and X_j are tensors */
	        for (kk=0;kk<n;kk++) Dq[Kj[kk]] += w[kk]*dXi[kk]*dXj[kk]*Xi[Ki[kk]];
	        if (tri) for (kk=0;kk<n-1;kk++) {
		  Dq[Kj[kk]] += wl[kk]*dXi[kk+1]*dXj[kk]*Xi[Ki[kk+1]];
		  Dq[Kj[kk+1]] += ws[kk]*dXi[kk]*dXj[kk+1]*Xi[Ki[kk]];
	        } 
	      } else if (tensi) { /* only X_i is tensor */
	        for (kk=0;kk<n;kk++) Dq[Kj[kk]] += w[kk]*dXi[kk]*Xi[Ki[kk]];
	        if (tri) for (kk=0;kk<n-1;kk++) {
		  Dq[Kj[kk]] += wl[kk]*dXi[kk+1]*Xi[Ki[kk+1]];
		  Dq[Kj[kk+1]] += ws[kk]*dXi[kk]*Xi[Ki[kk]];
	        }
	      } else if (tensj) { /* only X_j is tensor */
	        for (kk=0;kk<n;kk++) Dq[Kj[kk]] += w[kk]*dXj[kk]*Xi[Ki[kk]];
	        if (tri) for (kk=0;kk<n-1;kk++) {
		  Dq[Kj[kk]] += wl[kk]*dXj[kk]*Xi[Ki[kk+1]];
		  Dq[Kj[kk+1]] += ws[kk]*dXj[kk+1]*Xi[Ki[kk]];
	        } 
	      } else { /* both singletons */
	        for (kk=0;kk<n;kk++) Dq[Kj[kk]] += w[kk]*Xi[Ki[kk]];
	        if (tri) for (kk=0;kk<n-1;kk++) {
		  Dq[Kj[kk]] += wl[kk]*Xi[Ki[kk+1]];
		  Dq[Kj[kk+1]] += ws[kk]*Xi[Ki[kk]];
	        }
	      }
	    } /* q loop */
	  }  
	} /* end of summation convention loop */
	if (acc_w) { /* form X_im' \bar W X_j from \bar W */
	  if (rfac) { /* form C = \bar W X_j first */
	  /* dgemm(char *transa,char *transb,int *m,int *n,int *k,double *alpha,double *A,
                   int *lda, double *B, int *ldb, double *beta,double *C,int *ldc) 
             transa/b = 'T' or 'N' for A/B transposed or not. C = alpha op(A) op(B) + beta C,
             where op() is transpose or not. C is m by n. k is cols of op(A). ldx is rows of X
             in calling routine (to allow use of sub-matrices) */   
            F77_CALL(dgemm)(&ntrans,&ntrans,m+im,p+jm,m+jm,&done,W,m+im,X+off[jm],m+jm,&dzero,C,m+im); 
	  } else { /* form D = \bar W' X_i first */
            F77_CALL(dgemm)(&trans,&ntrans,m+jm,p+im,m+im,&done,W,m+im,X+off[im],m+im,&dzero,D,m+jm);
	  }  
	}
	if (rfac) { /* Xi'C direct to the r,c p[im] by p[jm] block of XWX */
          F77_CALL(dgemm)(&trans,&ntrans,p+im,p+jm,m+im,&done,X+off[im],m+im,C,m+im,&dzero,XWX+r*p[im]+c*p[jm]*nxwx,&nxwx);
	} else { /* D'Xj (same block of XWX) */
          F77_CALL(dgemm)(&trans,&ntrans,p+im,p+jm,m+jm,&done,D,m+jm,X+off[jm],m+jm,&dzero,XWX+r*p[im]+c*p[jm]*nxwx,&nxwx);
	}  
      } /* block col loop */
    } /* block loop */  
  } /* general case */
} /* XWXij0 */


void XWXijw(double *XWX,int i,int j, double *X,int *k, int *ks, int *m, int *p,int nx,int n,int *ts, int *dt,
	   int nt, double *w,double *ws, int tri,ptrdiff_t *off,double *work) {
/* Forms product X_i'WX_j where X_i and X_j are stored in compact form. 
   * W = diag(w) if tri!=0 and is tridiagonal otherwise, with super in ws and sub in ws + n-1;     
   * off[i] is offset to start of ith matrix (out of nx)  
   If Xk is a tensor product term, let dXk denote row tensor product of all but it's final 
   marginal. 
   
   Workspace: Let mi and mj index the final marginal of term i and j, then the dim of work needed is
              2*n + 3*m[mi] + I(n>m[mi]*m[mj])*m[mi]*m[mj] + max(m[mi]*p[mj],m[mj]*p[mi])

   EXPERIMENTAL - works

*/
  int si,sj,ri,rj,jm,im,r,c,kk,ddt,ddtj,koff,*K,*Ki,*Kj,pim,pjm,
    ii,jj,rfac,q,t,s,ddti,tensi,tensj,acc_w,alpha,*Kik,*Kjk,nxwx;
  ptrdiff_t mim,mjm; /* avoid integer overflow in large pointer calculations */ 
  double x,*wl,*dXi,*dXj,*pdXj,*Xt,*Xi,*Xj,*Cq,*Dq,done=1.0,dzero=0.0,
    *C,*D,*W,*wb,*wbs,*wbl,*p0,*p1,*p2,*p3;
  char trans = 'T',ntrans = 'N';
  si = ks[ts[i]+nx]-ks[ts[i]]; /* number of terms in summation convention for i */
  if (tri) wl = ws + n - 1; /* sub-diagonal */
  /* compute number of columns in dXi/ number of rows of blocks in product */
  for (ri=1,kk=ts[i];kk<ts[i]+dt[i]-1;kk++) ri *= p[kk];
  im = ts[i]+dt[i]-1; /* the index of the final marginal for term i */
  mim = (ptrdiff_t) m[im];
  nxwx = p[im] * ri; /* rows of XWX = cols of Xi */
  /* Allocate work space for dXi(n), dXj(n) and initialze pdXj*/
  dXi = work;work += n;pdXj = dXj = work;work += n;
  if (dt[i]==1&&dt[j]==1&&m[ts[i]]==n&&m[ts[j]]==n) { /* both sub matrices are dense */
    //Rprintf("de ");
    jm = ts[j];
    mjm = (ptrdiff_t) m[jm];
    pim = p[im];pjm = p[jm];
    for (r=0;r<pim;r++) {
      Xi = X + off[im] + mim * r;
      if (i==j) for (c=r;c<pjm;c++) { /* symmetric case */
	Xj = X + off[jm] + mjm * c;  
	if (tri) {
	  x = Xi[0]*(Xj[0]*w[0]+Xj[1]*ws[0]);
	  for (kk=1;kk<n-1;kk++) x += Xi[kk]*(wl[kk-1]*Xj[kk-1] + w[kk]*Xj[kk] + ws[kk]*Xj[kk+1]);
          x += Xi[n-1]*(wl[n-2]*Xj[n-2]+w[n-1]*Xj[n-1]);				   
	} else for (x=0.0,p2=Xi,p3=Xj,p0=w,p1=w+n;p0<p1;p0++,p2++,p3++) x += *p0 * *p2 * *p3;
	  //for (x=0.0,kk=0;kk<n;kk++) x += Xi[kk]*Xj[kk]*w[kk];
	XWX[c + (ptrdiff_t) nxwx * r] = XWX[r + (ptrdiff_t) nxwx * c] = x;
      }	else for (c=0;c<pjm;c++) {
        Xj = X + off[jm] + mjm * c;  
	if (tri) {
	  x = Xi[0]*(Xj[0]*w[0]+Xj[1]*ws[0]);
	  for (kk=1;kk<n-1;kk++) x += Xi[kk]*(wl[kk-1]*Xj[kk-1] + w[kk]*Xj[kk] + ws[kk]*Xj[kk+1]);
          x += Xi[n-1]*(wl[n-2]*Xj[n-2]+w[n-1]*Xj[n-1]);				   
	} else for (x=0.0,p2=Xi,p3=Xj,p0=w,p1=w+n;p0<p1;p0++,p2++,p3++) x += *p0 * *p2 * *p3;
	  //for (x=0.0,kk=0;kk<n;kk++) x += Xi[kk]*Xj[kk]*w[kk];
	XWX[r + (ptrdiff_t)nxwx * c] = x;
      }	  
    }
  } else if (!tri && i==j && si==1) { /* simplest setup - just accumulate diagonal */
    //Rprintf("wb ");
    /* Allocate space for wb(m[im]), wbs and wbl*/
    wb = work; work += mim;// wbs = work; work += m[i]; wbl = work; work += m[i];
    
    if (dt[i]>1) { /* tensor */
      ddt = dt[i]-1; /* number of marginals, exluding final */
      koff = 0; /* only one index vector per marginal, so no offset */
    }  
    for (r=0;r<ri;r++) { /* loop over block rows of product */
      if (dt[i]>1) { /* extract col r of dXi */
	for (kk=0;kk<n;kk++) dXi[kk] = 1.0;
	tensorXj(dXi, X + off[ts[i]], m + ts[i], p + ts[i],&ddt, 
		 k, &n, &r, ks + ts[i],&koff);
      }
      for (c=r;c<ri;c++) { /* block column loop - note symmetry */
	if (dt[i]>1) { /* extract col c of dXi */
          if (r!=c) {
            dXj=pdXj;for (kk=0;kk<n;kk++) dXj[kk] = 1.0;
	    tensorXj(dXj, X + off[ts[i]], m + ts[i], p + ts[i],&ddt, 
		     k, &n, &c, ks + ts[i],&koff);
	  } else dXj = dXi;
	}
	/* Clear work space to zero... */
	for (ii=0;ii<mim;ii++) wb[ii]=0.0;
        //if (tri) for (ii=0;ii<m[i]-1;ii++) wbs[ii] = wbl[ii] = 0.0; 
	K = k + ks[im] * n; /* index for final margin */
	/* Accumulate the weights ... */
	if (dt[i]>1) {
	  for (kk=0;kk<n;kk++) wb[K[kk]] += dXi[kk]*dXj[kk]*w[kk];
	  //for (p0=w,p1=w+n,p2=dXi,p3=dXj;p0<p1;p0++,p2++,p3++,K++) wb[*K] += *p0 * *p2 * *p3; 
	} else { /* singleton */
	  for (kk=0;kk<n;kk++) wb[K[kk]] += w[kk];
          //for (p0=w,p1=w+n;p0<p1;p0++,K++) wb[*K] += *p0;
	}
	/* Now form the Xi'WXi... */
	Xt = X + off[im]; /* final marginal model matrix */
	pim = p[im]; /* cols of Xt */
	if (r!=c)  for (jj=0;jj<pim;jj++) for (ii=0;ii<pim;ii++) { /* block need not be symmetric */
	   Xi = Xt + mim * ii; /* iith col of Xt */
           Xj = Xt + mim * jj; /* jjth col of Xt */
	   for (x=0.0,kk=0;kk<mim;kk++) x += wb[kk]*Xj[kk]*Xi[kk];
	   XWX[c*pim+jj+(r*pim +ii)* (ptrdiff_t) nxwx] = XWX[r*pim+ii+(c*pim +jj)*(ptrdiff_t)nxwx] = x; 
	} else for (ii=0;ii<pim;ii++) for (jj=ii;jj<pim;jj++) { /* diagonal and symmetric */ 
          Xi = Xt + mim * ii; /* iith col of Xt */
          Xj = Xt + mim * jj; /* jjth col of Xt */
	  for (x=0.0,kk=0;kk<mim;kk++) x += wb[kk]*Xj[kk]*Xi[kk];
	  XWX[r*pim+ ii + (c*pim+jj)*(ptrdiff_t)nxwx] = XWX[r*pim + jj + (c*pim + ii)*(ptrdiff_t)nxwx] =
	  XWX[(r*pim+ ii)*(ptrdiff_t)nxwx + c*pim+jj] = XWX[(r*pim + jj)*(ptrdiff_t)nxwx + c*pim + ii] = x;
	}					  
      } /* block col loop */ 
    }
  } else { /* general case */
    sj = ks[ts[j]+nx]-ks[ts[j]]; /* number of terms in summation convention for j */
    for (rj=1,kk=ts[j];kk<ts[j]+dt[j]-1;kk++) rj *= p[kk];
    jm = ts[j]+dt[j]-1; /* the index of the final marginal for term j */
    ddti = dt[i]-1;ddtj = dt[j]-1; /* number of marginals excluding final */
    if (ddti) tensi = 1; else tensi = 0; /* is term i a tensor? */
    if (ddtj) tensj = 1; else tensj = 0;
    mjm = (ptrdiff_t) m[jm];
    //Rprintf("n = %d mim = %ld mjm = %ld ",n,mim,mjm);
    if (n>mjm*mim) acc_w = 1; else acc_w = 0; /* accumulate \bar W or \bar W X_j / \bar W'X_i)? */
    if (acc_w) {
      if (p[im]*mim*mjm + p[im]*p[jm]*mjm > mim*mjm*p[jm] + p[im]*p[jm]*mim) rfac=0; else rfac=1;
      /* Allocate storage for W (mim*mjm) */
      W = work; work += mim*mjm;
    } else {
      /* now establish whether to form left product, D, or right product C */
      if (tensi) ii = 2; else ii = 1;if (tensj) ii++;
      if (tri) alpha = ii*3 + 3; else alpha = ii + 1; /* ~ ops per iteration of accumulation loop */ 
      if (alpha*si*sj*n*p[im]+mjm*p[im]*p[jm]<alpha*si*sj*n*p[jm]+mim*p[im]*p[jm]) rfac = 0; else rfac=1;
      if (mim == n) rfac = 0; else if (mjm == n) rfac = 1; /* make absolutely sure we do not form n by p*m product */
    } /* end of accumulation storage allocation */
    if (rfac) {
      /* Allocate storge for C mim by p[jm] */
      C = work; work += mim * p[jm];
    } else {
      /* Allocate storage for D mjm by p[im] */
      D = work; work += mjm * p[im]; 
    }	
    for (r=0;r<ri;r++) { /* loop over block rows of product */
      for (c=0;c<rj;c++) { /* block column loop*/
	if (acc_w) for (kk=0;kk<mim*mjm;kk++) W[kk] = 0.0; /* clear W */
	else if (rfac) for (kk=0;kk<mim*p[jm];kk++) C[kk] = 0.0; /* clear C */
	else for (kk=0;kk<mjm*p[im];kk++) D[kk] = 0.0; /* clear D */
        for (s=0;s<si;s++) for (t=0;t<sj;t++) { /* summation convention loop */
	  if (tensi) { /* extract col r of dXi according to sth set of index vectors */
	    for (kk=0;kk<n;kk++) dXi[kk] = 1.0;
	    tensorXj(dXi, X + off[ts[i]], m + ts[i], p + ts[i],&ddti, 
		 k, &n, &r, ks + ts[i],&s);
          }
	  if (tensj) { /* extract col c of dXj according to tth set of index vectors */
	    for (kk=0;kk<n;kk++) dXj[kk] = 1.0;
            tensorXj(dXj, X + off[ts[j]], m + ts[j], p + ts[j],&ddtj, 
		 k, &n, &c, ks + ts[j],&t);
	  }
	  Ki = k + (ks[im]+s) * n; /* index for final margin of i */
          Kj = k + (ks[jm]+t) * n; /* index for final margin of j */
	  if (acc_w) { /* weight accumulation */
	    //Rprintf("W ");
            if (tensi&&tensj) { /* i and j are tensors */
	      //for (kk=0;kk<n;kk++) W[Ki[kk]+mim*Kj[kk]] += w[kk]*dXi[kk]*dXj[kk];
	      for (Kik=Ki,Kjk=Kj,p0=w,p1=w+n,p2=dXi,p3=dXj;p0<p1;p0++,p2++,p3++,Kik++,Kjk++) W[*Kik + mim * *Kjk] += *p0 * *p2 * *p3;
	      if (tri) for (kk=0;kk<n-1;kk++) {
		  W[Ki[kk]+mim*Kj[kk+1]] += ws[kk]*dXi[kk]*dXj[kk+1];
		  W[Ki[kk+1]+mim*Kj[kk]] += wl[kk]*dXi[kk+1]*dXj[kk];
	      }	  
	    } else if (tensi) { /* only i is tensor */
              //for (kk=0;kk<n;kk++) W[Ki[kk]+mim*Kj[kk]] += w[kk]*dXi[kk];
	      for (Kik=Ki,Kjk=Kj,p0=w,p1=w+n,p2=dXi;p0<p1;p0++,p2++,Kik++,Kjk++) W[*Kik + mim * *Kjk] += *p0 * *p2;
	      if (tri) for (kk=0;kk<n-1;kk++) {
		  W[Ki[kk]+mim*Kj[kk+1]] += ws[kk]*dXi[kk];
		  W[Ki[kk+1]+mim*Kj[kk]] += wl[kk]*dXi[kk+1];
	      }	 
	    } else if (tensj) { /* only j is tensor */
	      //for (kk=0;kk<n;kk++) W[Ki[kk]+mim*Kj[kk]] += w[kk]*dXj[kk];
	      for (Kik=Ki,Kjk=Kj,p0=w,p1=w+n,p3=dXj;p0<p1;p0++,p3++,Kik++,Kjk++) W[*Kik + mim * *Kjk] += *p0  * *p3;
	      if (tri) for (kk=0;kk<n-1;kk++) {
		  W[Ki[kk]+mim*Kj[kk+1]] += ws[kk]*dXj[kk+1];
		  W[Ki[kk+1]+mim*Kj[kk]] += wl[kk]*dXj[kk];
	      }	  
	    } else { /* i and j are singletons */
              //for (kk=0;kk<n;kk++) W[Ki[kk]+mim*Kj[kk]] += w[kk];
	      for (Kik=Ki,Kjk=Kj,p0=w,p1=w+n;p0<p1;p0++,Kik++,Kjk++) W[*Kik + mim * *Kjk] += *p0;
	      if (tri) for (kk=0;kk<n-1;kk++) {
		  W[Ki[kk]+mim*Kj[kk+1]] += ws[kk];
		  W[Ki[kk+1]+mim*Kj[kk]] += wl[kk];
	      }	
	    }	    
	  } else if (rfac) { /* direct accumulation of right factor C (mim by p[jm]) */
	    //Rprintf("C ");
            for (q=0;q<p[jm];q++) {
	      Cq = C + q * mim;Xj = X + off[jm] + q * mjm; /* qth cols of C and X */ 
	      if (tensi&&tensj) { /* X_i and X_j are tensors */
	        // for (kk=0;kk<n;kk++) Cq[Ki[kk]] += w[kk]*dXi[kk]*dXj[kk]*Xj[Kj[kk]];
		for (Kik=Ki,Kjk=Kj,p0=w,p1=w+n,p2=dXi,p3=dXj;p0<p1;p0++,Kik++,Kjk++,p2++,p3++) Cq[*Kik] += *p0 * *p2 * *p3 * Xj[*Kjk];
	        if (tri) for (kk=0;kk<n-1;kk++) {
		  Cq[Ki[kk]] += ws[kk]*dXi[kk]*dXj[kk+1]*Xj[Kj[kk+1]];
		  Cq[Ki[kk+1]] += wl[kk]*dXi[kk+1]*dXj[kk]*Xj[Kj[kk]];
	        } 
	      } else if (tensi) { /* only X_i is tensor */
		//for (kk=0;kk<n;kk++) Cq[Ki[kk]] += w[kk]*dXi[kk]*Xj[Kj[kk]];
		for (Kik=Ki,Kjk=Kj,p0=w,p1=w+n,p2=dXi;p0<p1;p0++,Kik++,Kjk++,p2++) Cq[*Kik] += *p0 * *p2  * Xj[*Kjk];
	        if (tri) for (kk=0;kk<n-1;kk++) {
		  Cq[Ki[kk]] += ws[kk]*dXi[kk]*Xj[Kj[kk+1]];
		  Cq[Ki[kk+1]] += wl[kk]*dXi[kk+1]*Xj[Kj[kk]];
	        }
	      } else if (tensj) { /* only X_j is tensor */
	        //for (kk=0;kk<n;kk++) Cq[Ki[kk]] += w[kk]*dXj[kk]*Xj[Kj[kk]];
		for (Kik=Ki,Kjk=Kj,p0=w,p1=w+n,p3=dXj;p0<p1;p0++,Kik++,Kjk++,p3++) Cq[*Kik] += *p0 * *p3 * Xj[*Kjk];
	        if (tri) for (kk=0;kk<n-1;kk++) {
		  Cq[Ki[kk]] += ws[kk]*dXj[kk+1]*Xj[Kj[kk+1]];
		  Cq[Ki[kk+1]] += wl[kk]*dXj[kk]*Xj[Kj[kk]];
	        }
	      } else { /* both singletons */
		//for (kk=0;kk<n;kk++) Cq[Ki[kk]] += w[kk]*Xj[Kj[kk]];
		for (Kik=Ki,Kjk=Kj,p0=w,p1=w+n;p0<p1;p0++,Kik++,Kjk++) Cq[*Kik] += *p0 * Xj[*Kjk];
	        if (tri) for (kk=0;kk<n-1;kk++) {
		  Cq[Ki[kk]] += ws[kk]*Xj[Kj[kk+1]];
		  Cq[Ki[kk+1]] += wl[kk]*Xj[Kj[kk]];
	        }
	      }
	    } /* q loop */
	  } else { /* direct accumulation of left factor D (mjm by p[im]) = \bar W' X_i */
	    //Rprintf("D ");
            for (q=0;q<p[im];q++) {
	      Dq = D + q * mjm;Xi = X + off[im] + q * mim; /* qth cols of C and Xi */ 
	      if (tensi&&tensj) { /* X_i and X_j are tensors */
	        // for (kk=0;kk<n;kk++) Dq[Kj[kk]] += w[kk]*dXi[kk]*dXj[kk]*Xi[Ki[kk]];
		for (Kik=Ki,Kjk=Kj,p0=w,p1=w+n,p2=dXi,p3=dXj;p0<p1;p0++,Kik++,Kjk++,p2++,p3++) Dq[*Kjk] += *p0 * *p2 * *p3 * Xi[*Kik];
	        if (tri) for (kk=0;kk<n-1;kk++) {
		  Dq[Kj[kk]] += wl[kk]*dXi[kk+1]*dXj[kk]*Xi[Ki[kk+1]];
		  Dq[Kj[kk+1]] += ws[kk]*dXi[kk]*dXj[kk+1]*Xi[Ki[kk]];
	        } 
	      } else if (tensi) { /* only X_i is tensor */
	        //for (kk=0;kk<n;kk++) Dq[Kj[kk]] += w[kk]*dXi[kk]*Xi[Ki[kk]];
		for (Kik=Ki,Kjk=Kj,p0=w,p1=w+n,p2=dXi;p0<p1;p0++,Kik++,Kjk++,p2++) Dq[*Kjk] += *p0 * *p2  * Xi[*Kik];
	        if (tri) for (kk=0;kk<n-1;kk++) {
		  Dq[Kj[kk]] += wl[kk]*dXi[kk+1]*Xi[Ki[kk+1]];
		  Dq[Kj[kk+1]] += ws[kk]*dXi[kk]*Xi[Ki[kk]];
	        }
	      } else if (tensj) { /* only X_j is tensor */
	        //for (kk=0;kk<n;kk++) Dq[Kj[kk]] += w[kk]*dXj[kk]*Xi[Ki[kk]];
		for (Kik=Ki,Kjk=Kj,p0=w,p1=w+n,p3=dXj;p0<p1;p0++,Kik++,Kjk++,p3++) Dq[*Kjk] += *p0 * *p3 * Xi[*Kik];
	        if (tri) for (kk=0;kk<n-1;kk++) {
		  Dq[Kj[kk]] += wl[kk]*dXj[kk]*Xi[Ki[kk+1]];
		  Dq[Kj[kk+1]] += ws[kk]*dXj[kk+1]*Xi[Ki[kk]];
	        } 
	      } else { /* both singletons */
	        //for (kk=0;kk<n;kk++) Dq[Kj[kk]] += w[kk]*Xi[Ki[kk]];
		for (Kik=Ki,Kjk=Kj,p0=w,p1=w+n;p0<p1;p0++,Kik++,Kjk++) Dq[*Kjk] += *p0 * Xi[*Kik];
	        if (tri) for (kk=0;kk<n-1;kk++) {
		  Dq[Kj[kk]] += wl[kk]*Xi[Ki[kk+1]];
		  Dq[Kj[kk+1]] += ws[kk]*Xi[Ki[kk]];
	        }
	      }
	    } /* q loop */
	  }  
	} /* end of summation convention loop */
	if (acc_w) { /* form X_im' \bar W X_j from \bar W */
	  if (rfac) { /* form C = \bar W X_j first */
	  /* dgemm(char *transa,char *transb,int *m,int *n,int *k,double *alpha,double *A,
                   int *lda, double *B, int *ldb, double *beta,double *C,int *ldc) 
             transa/b = 'T' or 'N' for A/B transposed or not. C = alpha op(A) op(B) + beta C,
             where op() is transpose or not. C is m by n. k is cols of op(A). ldx is rows of X
             in calling routine (to allow use of sub-matrices) */   
            F77_CALL(dgemm)(&ntrans,&ntrans,m+im,p+jm,m+jm,&done,W,m+im,X+off[jm],m+jm,&dzero,C,m+im); 
	  } else { /* form D = \bar W' X_i first */
            F77_CALL(dgemm)(&trans,&ntrans,m+jm,p+im,m+im,&done,W,m+im,X+off[im],m+im,&dzero,D,m+jm);
	  }  
	}
	if (rfac) { /* Xi'C direct to the r,c p[im] by p[jm] block of XWX */
          F77_CALL(dgemm)(&trans,&ntrans,p+im,p+jm,m+im,&done,X+off[im],m+im,C,m+im,&dzero,XWX+r*p[im]+c*p[jm]*(ptrdiff_t)nxwx,&nxwx);
	} else { /* D'Xj (same block of XWX) */
          F77_CALL(dgemm)(&trans,&ntrans,p+im,p+jm,m+jm,&done,D,m+jm,X+off[jm],m+jm,&dzero,XWX+r*p[im]+c*p[jm]*(ptrdiff_t)nxwx,&nxwx);
	}  
      } /* block col loop */
    } /* block loop */  
  } /* general case */
} /* XWXijw */  

void XWXd0w(double *XWX,double *X,double *w,int *k,int *ks, int *m,int *p, int *n, int *nx, int *ts, 
          int *dt, int *nt,double *v,int *qc,int *nthreads,int *ar_stop,int *ar_row,double *ar_weights) {
/* essentially a driver routine for XWXij implementing block oriented cross products
   works if XWXijw above, used
*/   
  int *pt, *pd,i,q,j,maxm=0,maxp=0,maxmp=0,tri,r,c,pa,pb,*tps,ptot,*B,*C,*R,N,kk,tid=0;
  ptrdiff_t *off,*voff,mmp;
  double *work,*ws,*xwx,*xwx0,*x0,*x1,*p0,*p1,*p2,x,*cost;
  #ifndef OPENMP_ON
  *nthreads = 1;
  #endif
  if (*nthreads<1) *nthreads = 1;
  pt = (int *) CALLOC((size_t)*nt,sizeof(int)); /* the term dimensions */
  pd = (int *) CALLOC((size_t)*nt,sizeof(int)); /* storage for last marginal size */
  off = (ptrdiff_t *) CALLOC((size_t)*nx+1,sizeof(ptrdiff_t)); /* offsets for X submatrix starts */
  voff = (ptrdiff_t *) CALLOC((size_t)*nt+1,sizeof(ptrdiff_t)); /* offsets for v subvectors starts */
  tps = (int *) CALLOC((size_t)*nt+1,sizeof(int)); /* the term starts */
  for (q=i=0;i< *nt; i++) { /* work through the terms */
    for (j=0;j<dt[i];j++,q++) {
      if (j==dt[i]-1) pd[i] = p[q]; /* the relevant dimension for deciding product ordering */
      off[q+1] = off[q] + p[q] * (ptrdiff_t) m[q]; /* submatrix start offsets */
      if (j==0) pt[i] = p[q]; else pt[i] *= p[q]; /* term dimension */
      if (maxm<m[q] && m[q] < *n) maxm=m[q]; /* most rows of a non-dense sub-matrix */
      if (maxmp<p[q]) maxmp=p[q];
    } 
    if (qc[i]>0) voff[i+1] = voff[i] + pt[i]; else voff[i+1] = voff[i]; /* start of ith v vector */
    if (maxp<pt[i]) maxp=pt[i];
    if (qc[i]<=0) tps[i+1] = tps[i] + pt[i]; /* where ith terms starts in param vector */ 
    else tps[i+1] = tps[i] + pt[i] - 1; /* there is a tensor constraint to apply - reducing param count*/
  }
  q = 3 * *n + maxm + maxm * maxmp; /* note that we never allocate a W accumulation matrix with more than n elements */
  work = (double *)CALLOC((size_t)q * *nthreads,sizeof(double));
  mmp = maxp;mmp = mmp*mmp;
  xwx = (double *)CALLOC((size_t)mmp * *nthreads,sizeof(double)); /* r,c block storage */
  xwx0 = (double *) CALLOC((size_t) mmp * *nthreads,sizeof(double)); /* working cross product storage */
  ptot = tps[*nt]; /* total number of parameters */
  if (*ar_stop>=0) { /* model has AR component*/
    for (p0 = w,p1 = w + *n;p0<p1;p0++) *p0 = sqrt(*p0); /* sqrt weights */
    /* ar_weights[0,2,4,6,...,2*n-2] is the ld of the square root of the tri-diagonal AR weight matrix
       ar_weights[1,3,5,...,2*n-1] is the sub-diagonal. */ 
    ws = (double *)CALLOC((size_t) 2 * *n -2,sizeof(double)); /* super and sub diagonals */
    //for (i=0;i<*n-1;i++) ws[i+*n-1] = ws[i] = ar_weights[1+2*i]*ar_weights[2*i]*w[i]*w[i+1];
    //w[0] *= ar_weights[0]*ar_weights[0]*w[0];
    //for (i=1;i<*n;i++) w[i] *= (ar_weights[2*i-1]*ar_weights[2*i-1]+ar_weights[2*i]*ar_weights[2*i])*w[i]; 
    for (i=0;i<*n-1;i++) ws[i+*n-1] = ws[i] = ar_weights[2*i+1] * ar_weights[2*i+2] * w[i+1] * w[i];
    for (i=0;i<*n-1;i++) w[i] *= (ar_weights[2*i+1]*ar_weights[2*i+1]+ar_weights[2*i]*ar_weights[2*i])*w[i]; 
    i = *n-1;w[i] *= w[i]*ar_weights[2*i]*ar_weights[2*i];
    /* so now w contains the leading diagonal of the tri-diagonal weight matrix,
       ws is the super diagonal and ws + n - 1 is the sub-diagonal. */
    tri = 1;
  } else tri = 0;
  N = ((*nt + 1) * *nt)/2; B = (int *) CALLOC((size_t)N,sizeof(int));
  C = (int *) CALLOC((size_t)N,sizeof(int)); R = (int *) CALLOC((size_t)N,sizeof(int));
  cost = (double *) CALLOC((size_t)N,sizeof(double)); /* cost per block */
  for (kk=r=0;r < *nt;r++) for (c=r;c< *nt;c++,kk++) {
    B[kk]=kk;R[kk]=r;C[kk]=c;
    /* compute rough cost per block figures */
    if (m[r]*m[c] < *n) { /* weight accumulation */
      cost[kk] = m[r]*m[c]*(double) pt[r]*pt[c]/pd[c];
      x = m[r]*m[c]*(double) pt[c]*pt[r]/pd[r]; if (x<cost[kk]) cost[kk] = x;
    } else { /* direct accumulation */
      cost[kk] = *n * pt[r]*pt[c]/pd[c];;
      x = *n * pt[c]*pt[r]/pd[r]; if (x < cost[kk]) cost[kk] = x;
    }
  }
  revsort(cost,B,N); /* R reverse sort on cost, to re-order B - see R.h*/
  //Rprintf("nt = %d N = %d\n",*nt,N);
  #ifdef OPENMP_ON
  #pragma omp parallel for private(kk,r,c,tid,x0,x1,j,x,p0,p1,p2,pa,pb,i) num_threads(*nthreads) schedule(dynamic)
  #endif
  for (kk=0;kk<N;kk++) { /* the block loop */
    r = R[B[kk]];c=C[B[kk]]; /* set up allows blocks to be computed in any order by re-arranging B */
    #ifdef OPENMP_ON
    tid = omp_get_thread_num();
    #endif
    //Rprintf("r=%d c=%d tid = %d  ",r,c,tid); // dodgy not ts
    XWXijw(xwx + mmp * tid,r,c,X,k,ks,m,p,*nx,*n,ts, dt,*nt,w,ws,tri,off,work + tid*q); /* compute r,c block */
    /* if Xr is tensor, may need to apply constraint */
    if (dt[r]>1&&qc[r]>0) { /* first term is a tensor with a constraint */
      x0=x1=xwx; /* pointers to columns of xwx */ 
      /* col by col form (I-vv')xwx, dropping first row... */
      for (j=0;j<pt[c];j++) {
        for (x=0.0,p0=x0,p1=x0+pt[r],p2=v+voff[r];p0<p1;p0++,p2++) x += *p0 * *p2;
        for (p2=v+voff[r]+1,p1=x0+pt[r],x0++;x0<p1;x1++,x0++,p2++) *x1 = *x0 - *p2 * x;
      }
      pa = pt[r] -1;
    } else pa = pt[r];
    if (dt[c]>1&&qc[c]>0) { /* Xc term is a tensor with a constraint */
      /* copy xwx to xwx0 */
      for (p0=xwx+mmp*tid,p1=p0 + pt[c] * (ptrdiff_t) pa,p2=xwx0+mmp*tid;p0<p1;p0++,p2++) *p2 = *p0;
      /* row by row form xwx(I-vv') dropping first col... */
      for (j=0;j<pa;j++) {
        for (x=0.0,p0=xwx0+mmp*tid+j,p1=v+voff[c],p2=p1+pt[c];p1<p2;p0 += pa,p1++) x += *p0 * *p1; 
        x1 = xwx+mmp*tid + j;
        for (p0=xwx0+mmp*tid+j+pa,p1=v+voff[c],p2=p1+pt[c],p1++;p1<p2;x1+=pa,p1++,p0+=pa) *x1 = *p0 - *p1 *x; 
      }
      pb=pt[c]-1;
    } else pb=pt[c];
    /* copy result into overall XWX*/
    p0 = xwx + mmp*tid;
    for (i=0;i<pa;i++) for (j=0;j<pb;j++) 
      XWX[i+tps[r]+(j+tps[c]) * (ptrdiff_t)ptot] = XWX[j+tps[c]+(i+tps[r]) * (ptrdiff_t)ptot] 
                                                 = p0[i + pa * (ptrdiff_t)j];
  }    
  FREE(pt);FREE(pd);FREE(off);FREE(voff);FREE(tps);FREE(work); if (tri) FREE(ws);FREE(xwx);FREE(xwx0);
  FREE(B);FREE(R);FREE(C);FREE(cost);
} /* XWXd0w */ 


void XWXij(double *XWX,int i,int j,int r,int c, double *X,int *k, int *ks, int *m, int *p,int nx,int n,int *ts, int *dt,
	   int nt, double *w,double *ws, int tri,ptrdiff_t *off,double *work,int nxwx) {
/* Forms product X_i'WX_j where X_i and X_j are stored in compact form. 
   * W = diag(w) if tri!=0 and is tridiagonal otherwise, with super in ws and sub in ws + n-1;     
   * off[i] is offset to start of ith matrix (out of nx)  
   If Xk is a tensor product term, let dXk denote row tensor product of all but it's final 
   marginal. 
   
   Workspace: Let mi and mj index the final marginal of term i and j, then the dim of work needed is
              2*n + 3*m[mi] + I(n>m[mi]*m[mj])*m[mi]*m[mj] + max(m[mi]*p[mj],m[mj]*p[mi])

   EXPERIMENTAL

*/
  int si,sj,ri,rj,jm,im,kk,ddt,ddtj,koff,*K,*Ki,*Kj,pim,pjm,
    ii,jj,rfac,q,t,s,ddti,tensi,tensj,acc_w,alpha,*Kik,*Kjk;//,nxwx;
  ptrdiff_t mim,mjm; /* avoid integer overflow in large pointer calculations */ 
  double x,*wl,*dXi,*dXj,*pdXj,*Xt,*Xi,*Xj,*Cq,*Dq,done=1.0,dzero=0.0,
    *C,*D,*W,*wb,*wbs,*wbl,*p0,*p1,*p2,*p3;
  char trans = 'T',ntrans = 'N';
  si = ks[ts[i]+nx]-ks[ts[i]]; /* number of terms in summation convention for i */
  if (tri) wl = ws + n - 1; /* sub-diagonal */
  /* compute number of columns in dXi/ number of rows of blocks in product */
  for (ri=1,kk=ts[i];kk<ts[i]+dt[i]-1;kk++) ri *= p[kk];
  im = ts[i]+dt[i]-1; /* the index of the final marginal for term i */
  mim = (ptrdiff_t) m[im];
  //nxwx = p[im] * ri; /* rows of XWX = cols of Xi */
  /* Allocate work space for dXi(n), dXj(n) and initialze pdXj*/
  dXi = work;work += n;pdXj = dXj = work;work += n;
  if (dt[i]==1&&dt[j]==1&&m[ts[i]]==n&&m[ts[j]]==n) { /* both sub matrices are dense - NOTE: what if there was a summation convention term in this case? */
    //Rprintf("de ");
    jm = ts[j];
    mjm = (ptrdiff_t) m[jm];
    pim = p[im];pjm = p[jm];
    for (ii=0;ii<pim;ii++) {
      Xi = X + off[im] + mim * ii;
      if (i==j) for (jj=ii;jj<pjm;jj++) { /* symmetric case */
	Xj = X + off[jm] + mjm * jj;  
	if (tri) {
	  x = Xi[0]*(Xj[0]*w[0]+Xj[1]*ws[0]);
	  for (kk=1;kk<n-1;kk++) x += Xi[kk]*(wl[kk-1]*Xj[kk-1] + w[kk]*Xj[kk] + ws[kk]*Xj[kk+1]);
          x += Xi[n-1]*(wl[n-2]*Xj[n-2]+w[n-1]*Xj[n-1]);				   
	} else for (x=0.0,p2=Xi,p3=Xj,p0=w,p1=w+n;p0<p1;p0++,p2++,p3++) x += *p0 * *p2 * *p3;
	  //for (x=0.0,kk=0;kk<n;kk++) x += Xi[kk]*Xj[kk]*w[kk];
	XWX[jj + (ptrdiff_t) nxwx * ii] = XWX[ii + (ptrdiff_t) nxwx * jj] = x;
      }	else for (jj=0;jj<pjm;jj++) {
        Xj = X + off[jm] + mjm * jj;  
	if (tri) {
	  x = Xi[0]*(Xj[0]*w[0]+Xj[1]*ws[0]);
	  for (kk=1;kk<n-1;kk++) x += Xi[kk]*(wl[kk-1]*Xj[kk-1] + w[kk]*Xj[kk] + ws[kk]*Xj[kk+1]);
          x += Xi[n-1]*(wl[n-2]*Xj[n-2]+w[n-1]*Xj[n-1]);				   
	} else for (x=0.0,p2=Xi,p3=Xj,p0=w,p1=w+n;p0<p1;p0++,p2++,p3++) x += *p0 * *p2 * *p3;
	  //for (x=0.0,kk=0;kk<n;kk++) x += Xi[kk]*Xj[kk]*w[kk];
	XWX[ii + (ptrdiff_t)nxwx * jj] = x;
      }	  
    }
  } else if (!tri && i==j && si==1) { /* simplest setup - just accumulate diagonal */
    //Rprintf("wb ");
    /* Allocate space for wb(m[im]), wbs and wbl*/
    wb = work; work += mim;// wbs = work; work += m[i]; wbl = work; work += m[i];
    
    if (dt[i]>1) { /* tensor */
      ddt = dt[i]-1; /* number of marginals, exluding final */
      koff = 0; /* only one index vector per marginal, so no offset */
    }  
    //for (r=0;r<ri;r++) { /* loop over block rows of product */
      if (dt[i]>1) { /* extract col r of dXi */
	for (kk=0;kk<n;kk++) dXi[kk] = 1.0;
	tensorXj(dXi, X + off[ts[i]], m + ts[i], p + ts[i],&ddt, 
		 k, &n, &r, ks + ts[i],&koff);
      }
      //for (c=r;c<ri;c++) { /* block column loop - note symmetry */
	if (dt[i]>1) { /* extract col c of dXi */
          if (r!=c) {
            dXj=pdXj;for (kk=0;kk<n;kk++) dXj[kk] = 1.0;
	    tensorXj(dXj, X + off[ts[i]], m + ts[i], p + ts[i],&ddt, 
		     k, &n, &c, ks + ts[i],&koff);
	  } else dXj = dXi;
	}
	/* Clear work space to zero... */
	for (ii=0;ii<mim;ii++) wb[ii]=0.0;
        //if (tri) for (ii=0;ii<m[i]-1;ii++) wbs[ii] = wbl[ii] = 0.0; 
	K = k + ks[im] * n; /* index for final margin */
	/* Accumulate the weights ... */
	if (dt[i]>1) {
	  for (kk=0;kk<n;kk++) wb[K[kk]] += dXi[kk]*dXj[kk]*w[kk];
	  //for (p0=w,p1=w+n,p2=dXi,p3=dXj;p0<p1;p0++,p2++,p3++,K++) wb[*K] += *p0 * *p2 * *p3; 
	} else { /* singleton */
	  for (kk=0;kk<n;kk++) wb[K[kk]] += w[kk];
          //for (p0=w,p1=w+n;p0<p1;p0++,K++) wb[*K] += *p0;
	}
	/* Now form the Xi'WXi... */
	Xt = X + off[im]; /* final marginal model matrix */
	pim = p[im]; /* cols of Xt */
	if (r!=c)  for (jj=0;jj<pim;jj++) for (ii=0;ii<pim;ii++) { /* block need not be symmetric */
	   Xi = Xt + mim * ii; /* iith col of Xt */
           Xj = Xt + mim * jj; /* jjth col of Xt */
	   for (x=0.0,kk=0;kk<mim;kk++) x += wb[kk]*Xj[kk]*Xi[kk];
	   XWX[c*pim+jj+(r*pim +ii)* (ptrdiff_t) nxwx] = XWX[r*pim+ii+(c*pim +jj)*(ptrdiff_t)nxwx] = x; 
	} else for (ii=0;ii<pim;ii++) for (jj=ii;jj<pim;jj++) { /* diagonal and symmetric */ 
          Xi = Xt + mim * ii; /* iith col of Xt */
          Xj = Xt + mim * jj; /* jjth col of Xt */
	  for (x=0.0,kk=0;kk<mim;kk++) x += wb[kk]*Xj[kk]*Xi[kk];
	  XWX[r*pim+ ii + (c*pim+jj)*(ptrdiff_t)nxwx] = XWX[r*pim + jj + (c*pim + ii)*(ptrdiff_t)nxwx] =
	  XWX[(r*pim+ ii)*(ptrdiff_t)nxwx + c*pim+jj] = XWX[(r*pim + jj)*(ptrdiff_t)nxwx + c*pim + ii] = x;
	}					  
	//} /* block col loop */ 
	//}
  } else { /* general case */
    sj = ks[ts[j]+nx]-ks[ts[j]]; /* number of terms in summation convention for j */
    for (rj=1,kk=ts[j];kk<ts[j]+dt[j]-1;kk++) rj *= p[kk];
    jm = ts[j]+dt[j]-1; /* the index of the final marginal for term j */
    ddti = dt[i]-1;ddtj = dt[j]-1; /* number of marginals excluding final */
    if (ddti) tensi = 1; else tensi = 0; /* is term i a tensor? */
    if (ddtj) tensj = 1; else tensj = 0;
    mjm = (ptrdiff_t) m[jm];
    //Rprintf("n = %d mim = %ld mjm = %ld ",n,mim,mjm);
    if (n>mjm*mim) acc_w = 1; else acc_w = 0; /* accumulate \bar W or \bar W X_j / \bar W'X_i)? */
    if (acc_w) {
      if (p[im]*mim*mjm + p[im]*p[jm]*mjm > mim*mjm*p[jm] + p[im]*p[jm]*mim) rfac=0; else rfac=1;
      /* Allocate storage for W (mim*mjm) */
      W = work; work += mim*mjm;
    } else {
      /* now establish whether to form left product, D, or right product C */
      if (tensi) ii = 2; else ii = 1;if (tensj) ii++;
      if (tri) alpha = ii*3 + 3; else alpha = ii + 1; /* ~ ops per iteration of accumulation loop */ 
      if (alpha*si*sj*n*p[im]+mjm*p[im]*p[jm]<alpha*si*sj*n*p[jm]+mim*p[im]*p[jm]) rfac = 0; else rfac=1;
      if (mim == n) rfac = 0; else if (mjm == n) rfac = 1; /* make absolutely sure we do not form n by p*m product */
    } /* end of accumulation storage allocation */
    if (rfac) {
      /* Allocate storge for C mim by p[jm] */
      C = work; work += mim * p[jm];
    } else {
      /* Allocate storage for D mjm by p[im] */
      D = work; work += mjm * p[im]; 
    }	
    //for (r=0;r<ri;r++) { /* loop over block rows of product */
    //  for (c=0;c<rj;c++) { /* block column loop*/
	if (acc_w) for (kk=0;kk<mim*mjm;kk++) W[kk] = 0.0; /* clear W */
	else if (rfac) for (kk=0;kk<mim*p[jm];kk++) C[kk] = 0.0; /* clear C */
	else for (kk=0;kk<mjm*p[im];kk++) D[kk] = 0.0; /* clear D */
        for (s=0;s<si;s++) for (t=0;t<sj;t++) { /* summation convention loop */
	  if (tensi) { /* extract col r of dXi according to sth set of index vectors */
	    for (kk=0;kk<n;kk++) dXi[kk] = 1.0;
	    tensorXj(dXi, X + off[ts[i]], m + ts[i], p + ts[i],&ddti, 
		 k, &n, &r, ks + ts[i],&s);
          }
	  if (tensj) { /* extract col c of dXj according to tth set of index vectors */
	    for (kk=0;kk<n;kk++) dXj[kk] = 1.0;
            tensorXj(dXj, X + off[ts[j]], m + ts[j], p + ts[j],&ddtj, 
		 k, &n, &c, ks + ts[j],&t);
	  }
	  Ki = k + (ks[im]+s) * n; /* index for final margin of i */
          Kj = k + (ks[jm]+t) * n; /* index for final margin of j */
	  if (acc_w) { /* weight accumulation */
	    //Rprintf("W ");
            if (tensi&&tensj) { /* i and j are tensors */
	      //for (kk=0;kk<n;kk++) W[Ki[kk]+mim*Kj[kk]] += w[kk]*dXi[kk]*dXj[kk];
	      for (Kik=Ki,Kjk=Kj,p0=w,p1=w+n,p2=dXi,p3=dXj;p0<p1;p0++,p2++,p3++,Kik++,Kjk++) W[*Kik + mim * *Kjk] += *p0 * *p2 * *p3;
	      if (tri) for (kk=0;kk<n-1;kk++) {
		  W[Ki[kk]+mim*Kj[kk+1]] += ws[kk]*dXi[kk]*dXj[kk+1];
		  W[Ki[kk+1]+mim*Kj[kk]] += wl[kk]*dXi[kk+1]*dXj[kk];
	      }	  
	    } else if (tensi) { /* only i is tensor */
              //for (kk=0;kk<n;kk++) W[Ki[kk]+mim*Kj[kk]] += w[kk]*dXi[kk];
	      for (Kik=Ki,Kjk=Kj,p0=w,p1=w+n,p2=dXi;p0<p1;p0++,p2++,Kik++,Kjk++) W[*Kik + mim * *Kjk] += *p0 * *p2;
	      if (tri) for (kk=0;kk<n-1;kk++) {
		  W[Ki[kk]+mim*Kj[kk+1]] += ws[kk]*dXi[kk];
		  W[Ki[kk+1]+mim*Kj[kk]] += wl[kk]*dXi[kk+1];
	      }	 
	    } else if (tensj) { /* only j is tensor */
	      //for (kk=0;kk<n;kk++) W[Ki[kk]+mim*Kj[kk]] += w[kk]*dXj[kk];
	      for (Kik=Ki,Kjk=Kj,p0=w,p1=w+n,p3=dXj;p0<p1;p0++,p3++,Kik++,Kjk++) W[*Kik + mim * *Kjk] += *p0  * *p3;
	      if (tri) for (kk=0;kk<n-1;kk++) {
		  W[Ki[kk]+mim*Kj[kk+1]] += ws[kk]*dXj[kk+1];
		  W[Ki[kk+1]+mim*Kj[kk]] += wl[kk]*dXj[kk];
	      }	  
	    } else { /* i and j are singletons */
              //for (kk=0;kk<n;kk++) W[Ki[kk]+mim*Kj[kk]] += w[kk];
	      for (Kik=Ki,Kjk=Kj,p0=w,p1=w+n;p0<p1;p0++,Kik++,Kjk++) W[*Kik + mim * *Kjk] += *p0;
	      if (tri) for (kk=0;kk<n-1;kk++) {
		  W[Ki[kk]+mim*Kj[kk+1]] += ws[kk];
		  W[Ki[kk+1]+mim*Kj[kk]] += wl[kk];
	      }	
	    }	    
	  } else if (rfac) { /* direct accumulation of right factor C (mim by p[jm]) */
	    //Rprintf("C ");
            for (q=0;q<p[jm];q++) {
	      Cq = C + q * mim;Xj = X + off[jm] + q * mjm; /* qth cols of C and X */ 
	      if (tensi&&tensj) { /* X_i and X_j are tensors */
	        // for (kk=0;kk<n;kk++) Cq[Ki[kk]] += w[kk]*dXi[kk]*dXj[kk]*Xj[Kj[kk]];
		for (Kik=Ki,Kjk=Kj,p0=w,p1=w+n,p2=dXi,p3=dXj;p0<p1;p0++,Kik++,Kjk++,p2++,p3++) Cq[*Kik] += *p0 * *p2 * *p3 * Xj[*Kjk];
	        if (tri) for (kk=0;kk<n-1;kk++) {
		  Cq[Ki[kk]] += ws[kk]*dXi[kk]*dXj[kk+1]*Xj[Kj[kk+1]];
		  Cq[Ki[kk+1]] += wl[kk]*dXi[kk+1]*dXj[kk]*Xj[Kj[kk]];
	        } 
	      } else if (tensi) { /* only X_i is tensor */
		//for (kk=0;kk<n;kk++) Cq[Ki[kk]] += w[kk]*dXi[kk]*Xj[Kj[kk]];
		for (Kik=Ki,Kjk=Kj,p0=w,p1=w+n,p2=dXi;p0<p1;p0++,Kik++,Kjk++,p2++) Cq[*Kik] += *p0 * *p2  * Xj[*Kjk];
	        if (tri) for (kk=0;kk<n-1;kk++) {
		  Cq[Ki[kk]] += ws[kk]*dXi[kk]*Xj[Kj[kk+1]];
		  Cq[Ki[kk+1]] += wl[kk]*dXi[kk+1]*Xj[Kj[kk]];
	        }
	      } else if (tensj) { /* only X_j is tensor */
	        //for (kk=0;kk<n;kk++) Cq[Ki[kk]] += w[kk]*dXj[kk]*Xj[Kj[kk]];
		for (Kik=Ki,Kjk=Kj,p0=w,p1=w+n,p3=dXj;p0<p1;p0++,Kik++,Kjk++,p3++) Cq[*Kik] += *p0 * *p3 * Xj[*Kjk];
	        if (tri) for (kk=0;kk<n-1;kk++) {
		  Cq[Ki[kk]] += ws[kk]*dXj[kk+1]*Xj[Kj[kk+1]];
		  Cq[Ki[kk+1]] += wl[kk]*dXj[kk]*Xj[Kj[kk]];
	        }
	      } else { /* both singletons */
		//for (kk=0;kk<n;kk++) Cq[Ki[kk]] += w[kk]*Xj[Kj[kk]];
		for (Kik=Ki,Kjk=Kj,p0=w,p1=w+n;p0<p1;p0++,Kik++,Kjk++) Cq[*Kik] += *p0 * Xj[*Kjk];
	        if (tri) for (kk=0;kk<n-1;kk++) {
		  Cq[Ki[kk]] += ws[kk]*Xj[Kj[kk+1]];
		  Cq[Ki[kk+1]] += wl[kk]*Xj[Kj[kk]];
	        }
	      }
	    } /* q loop */
	  } else { /* direct accumulation of left factor D (mjm by p[im]) = \bar W' X_i */
	    //Rprintf("D ");
            for (q=0;q<p[im];q++) {
	      Dq = D + q * mjm;Xi = X + off[im] + q * mim; /* qth cols of C and Xi */ 
	      if (tensi&&tensj) { /* X_i and X_j are tensors */
	        // for (kk=0;kk<n;kk++) Dq[Kj[kk]] += w[kk]*dXi[kk]*dXj[kk]*Xi[Ki[kk]];
		for (Kik=Ki,Kjk=Kj,p0=w,p1=w+n,p2=dXi,p3=dXj;p0<p1;p0++,Kik++,Kjk++,p2++,p3++) Dq[*Kjk] += *p0 * *p2 * *p3 * Xi[*Kik];
	        if (tri) for (kk=0;kk<n-1;kk++) {
		  Dq[Kj[kk]] += wl[kk]*dXi[kk+1]*dXj[kk]*Xi[Ki[kk+1]];
		  Dq[Kj[kk+1]] += ws[kk]*dXi[kk]*dXj[kk+1]*Xi[Ki[kk]];
	        } 
	      } else if (tensi) { /* only X_i is tensor */
	        //for (kk=0;kk<n;kk++) Dq[Kj[kk]] += w[kk]*dXi[kk]*Xi[Ki[kk]];
		for (Kik=Ki,Kjk=Kj,p0=w,p1=w+n,p2=dXi;p0<p1;p0++,Kik++,Kjk++,p2++) Dq[*Kjk] += *p0 * *p2  * Xi[*Kik];
	        if (tri) for (kk=0;kk<n-1;kk++) {
		  Dq[Kj[kk]] += wl[kk]*dXi[kk+1]*Xi[Ki[kk+1]];
		  Dq[Kj[kk+1]] += ws[kk]*dXi[kk]*Xi[Ki[kk]];
	        }
	      } else if (tensj) { /* only X_j is tensor */
	        //for (kk=0;kk<n;kk++) Dq[Kj[kk]] += w[kk]*dXj[kk]*Xi[Ki[kk]];
		for (Kik=Ki,Kjk=Kj,p0=w,p1=w+n,p3=dXj;p0<p1;p0++,Kik++,Kjk++,p3++) Dq[*Kjk] += *p0 * *p3 * Xi[*Kik];
	        if (tri) for (kk=0;kk<n-1;kk++) {
		  Dq[Kj[kk]] += wl[kk]*dXj[kk]*Xi[Ki[kk+1]];
		  Dq[Kj[kk+1]] += ws[kk]*dXj[kk+1]*Xi[Ki[kk]];
	        } 
	      } else { /* both singletons */
	        //for (kk=0;kk<n;kk++) Dq[Kj[kk]] += w[kk]*Xi[Ki[kk]];
		for (Kik=Ki,Kjk=Kj,p0=w,p1=w+n;p0<p1;p0++,Kik++,Kjk++) Dq[*Kjk] += *p0 * Xi[*Kik];
	        if (tri) for (kk=0;kk<n-1;kk++) {
		  Dq[Kj[kk]] += wl[kk]*Xi[Ki[kk+1]];
		  Dq[Kj[kk+1]] += ws[kk]*Xi[Ki[kk]];
	        }
	      }
	    } /* q loop */
	  }  
	} /* end of summation convention loop */
	if (acc_w) { /* form X_im' \bar W X_j from \bar W */
	  if (rfac) { /* form C = \bar W X_j first */
	  /* dgemm(char *transa,char *transb,int *m,int *n,int *k,double *alpha,double *A,
                   int *lda, double *B, int *ldb, double *beta,double *C,int *ldc) 
             transa/b = 'T' or 'N' for A/B transposed or not. C = alpha op(A) op(B) + beta C,
             where op() is transpose or not. C is m by n. k is cols of op(A). ldx is rows of X
             in calling routine (to allow use of sub-matrices) */   
            F77_CALL(dgemm)(&ntrans,&ntrans,m+im,p+jm,m+jm,&done,W,m+im,X+off[jm],m+jm,&dzero,C,m+im); 
	  } else { /* form D = \bar W' X_i first */
            F77_CALL(dgemm)(&trans,&ntrans,m+jm,p+im,m+im,&done,W,m+im,X+off[im],m+im,&dzero,D,m+jm);
	  }  
	}
	if (rfac) { /* Xi'C direct to the r,c p[im] by p[jm] block of XWX */
          F77_CALL(dgemm)(&trans,&ntrans,p+im,p+jm,m+im,&done,X+off[im],m+im,C,m+im,&dzero,XWX+r*p[im]+c*p[jm]*(ptrdiff_t)nxwx,&nxwx);
	} else { /* D'Xj (same block of XWX) */
          F77_CALL(dgemm)(&trans,&ntrans,p+im,p+jm,m+jm,&done,D,m+jm,X+off[jm],m+jm,&dzero,XWX+r*p[im]+c*p[jm]*(ptrdiff_t)nxwx,&nxwx);
	}  
	//} /* block col loop */
	//} /* block loop */  
  } /* general case */
} /* XWXij */  

void XWXd0(double *XWX,double *X,double *w,int *k,int *ks, int *m,int *p, int *n, int *nx, int *ts, 
          int *dt, int *nt,double *v,int *qc,int *nthreads,int *ar_stop,int *ar_row,double *ar_weights) {
/* essentially a driver routine for XWXij implementing block oriented cross products
   
   This version has the looping over sub-blocks, associated with tensor product terms, located in this routine
   to get better load balancing.
   
   Requires XWX to be over-sized on entry - namely n.params + n.terms by n.params + n.terms instead of
   n.params by n.params.
*/   
  int *pt, *pd,i,q,j,si,maxm=0,maxp=0,maxmp=0,tri,r,c,rb,cb,rt,ct,pa,pb,*tps,*tpsu,ptot,*b,*B,*C,*R,*sb,N,kk,kb,tid=0,nxwx=0;
  ptrdiff_t *off,*voff,mmp;
  double *work,*ws,*Cost,*cost,
    //*xwx,*xwx0,
    *x0,*x1,*p0,*p1,*p2,x;
  #ifndef OPENMP_ON
  *nthreads = 1;
  #endif
  if (*nthreads<1) *nthreads = 1;
  pt = (int *) CALLOC((size_t)*nt,sizeof(int)); /* the term dimensions */
  pd = (int *) CALLOC((size_t)*nt,sizeof(int)); /* storage for last marginal size */
  off = (ptrdiff_t *) CALLOC((size_t)*nx+1,sizeof(ptrdiff_t)); /* offsets for X submatrix starts */
  voff = (ptrdiff_t *) CALLOC((size_t)*nt+1,sizeof(ptrdiff_t)); /* offsets for v subvectors starts */
  tps = (int *) CALLOC((size_t)*nt+1,sizeof(int)); /* the term starts */
  tpsu = (int *) CALLOC((size_t)*nt+1,sizeof(int)); /* the unconstrained term starts */
  for (q=i=0;i< *nt; i++) { /* work through the terms */
    for (j=0;j<dt[i];j++,q++) {
      if (j==dt[i]-1) pd[i] = p[q]; /* the relevant dimension for deciding product ordering */
      off[q+1] = off[q] + p[q] * (ptrdiff_t) m[q]; /* submatrix start offsets */
      if (j==0) pt[i] = p[q]; else pt[i] *= p[q]; /* term dimension */
      if (maxm<m[q] && m[q] < *n) maxm=m[q]; /* most rows of a non-dense sub-matrix */
      if (maxmp<p[q]) maxmp=p[q];
    } 
    if (qc[i]>0) voff[i+1] = voff[i] + pt[i]; else voff[i+1] = voff[i]; /* start of ith v vector */
    if (maxp<pt[i]) maxp=pt[i];
    if (qc[i]<=0) tps[i+1] = tps[i] + pt[i]; /* where ith terms starts in param vector */ 
    else tps[i+1] = tps[i] + pt[i] - 1; /* there is a tensor constraint to apply - reducing param count*/
    tpsu[i+1] = tpsu[i] + pt[i]; /* where ith term starts in unconstrained param vector */ 
  }
  q = 3 * *n + maxm + maxm * maxmp; /* note that we never allocate a W accumulation matrix with more than n elements */
  work = (double *)CALLOC((size_t)q * *nthreads,sizeof(double));
  mmp = maxp;mmp = mmp*mmp;
  //xwx = (double *)CALLOC((size_t)mmp * *nthreads,sizeof(double)); /* r,c block storage */
  //xwx0 = (double *) CALLOC((size_t) mmp * *nthreads,sizeof(double)); /* working cross product storage */
  ptot = tps[*nt]; /* total number of parameters */
  nxwx = tpsu[*nt];
  if (*ar_stop>=0) { /* model has AR component*/
    for (p0 = w,p1 = w + *n;p0<p1;p0++) *p0 = sqrt(*p0); /* sqrt weights */
    /* ar_weights[0,2,4,6,...,2*n-2] is the ld of the square root of the tri-diagonal AR weight matrix
       ar_weights[1,3,5,...,2*n-1] is the sub-diagonal. */ 
    ws = (double *)CALLOC((size_t) 2 * *n -2,sizeof(double)); /* super and sub diagonals */
    //for (i=0;i<*n-1;i++) ws[i+*n-1] = ws[i] = ar_weights[1+2*i]*ar_weights[2*i]*w[i]*w[i+1];
    //w[0] *= ar_weights[0]*ar_weights[0]*w[0];
    //for (i=1;i<*n;i++) w[i] *= (ar_weights[2*i-1]*ar_weights[2*i-1]+ar_weights[2*i]*ar_weights[2*i])*w[i]; 
    for (i=0;i<*n-1;i++) ws[i+*n-1] = ws[i] = ar_weights[2*i+1] * ar_weights[2*i+2] * w[i+1] * w[i];
    for (i=0;i<*n-1;i++) w[i] *= (ar_weights[2*i+1]*ar_weights[2*i+1]+ar_weights[2*i]*ar_weights[2*i])*w[i]; 
    i = *n-1;w[i] *= w[i]*ar_weights[2*i]*ar_weights[2*i];
    /* so now w contains the leading diagonal of the tri-diagonal weight matrix,
       ws is the super diagonal and ws + n - 1 is the sub-diagonal. */
    tri = 1;
  } else tri = 0;
  N = ((*nt + 1) * *nt)/2; B = (int *) CALLOC((size_t)N,sizeof(int));
  C = (int *) CALLOC((size_t)N,sizeof(int)); R = (int *) CALLOC((size_t)N,sizeof(int));
  sb = (int *) CALLOC((size_t)N+1,sizeof(int)); /* at which sub-block does block start */
  Cost = (double *)CALLOC((size_t) N,sizeof(double));
  sb[0] = 0; /* start of first sub-block of block 0 */
  for (kk=r=0;r < *nt;r++) for (c=r;c< *nt;c++,kk++) {
      //B[kk]=kk;
    R[kk]=r;C[kk]=c;
    si = ks[ts[r] + *nx] - ks[ts[r]]; /* number of terms in summation convention for i */
    if (r==c && si==1 && !tri) {
      i = pt[r]/pd[r];
      i = i*(i+1)/2;
    } else {
      i = (pt[r]/pd[r])*(pt[c]/pd[c]);
    }  
    sb[kk+1] = sb[kk] + i; /* next block starts at this sub-block */
    /* compute rough cost per sub block figures */
    if (m[r]*m[c] < *n) { /* weight accumulation */
       Cost[kk] = m[r]*m[c]*(double) pd[c];
       x = m[r]*m[c]*(double) pd[r]; if (x<Cost[kk]) Cost[kk] = x;
    } else { /* direct accumulation */
      Cost[kk] = *n * pd[c];;
      x = *n * pd[r]; if (x < Cost[kk]) Cost[kk] = x;
    }
  }
  b = (int *) CALLOC((size_t)sb[N],sizeof(int));
  B = (int *) CALLOC((size_t)sb[N],sizeof(int));
  cost = (double *)CALLOC((size_t)sb[N],sizeof(double));
  for (kb=0,i=0;i<sb[N];i++) {
    b[i]=i;while (i>=sb[kb+1]) kb++; /* kb is main block */
    rb = R[kb];cb=C[kb];
    cost[i] = Cost[kb];
    B[i] = kb; 
  }  
  revsort(cost,b,sb[N]); /* R reverse sort on cost, to re-order b - see R.h*/
  //Rprintf("nt = %d N = %d\n",*nt,N);
  //for (i=0;i<=N;i++) Rprintf("%d  ",sb[i]);Rprintf("sb\n");
  
  //for (i=0;i<= *nt;i++) Rprintf("%d  ",tpsu[i]);Rprintf("tpsu\n");
  /* In what follows rb and cb are the whole term row column indices. r and c are the sub blocks within 
     the cross-product between two terms. The sub blocks arise when we have tensor product terms. The cleaner 
     design in which the sub-blocks are dealt with in XWXij does not load balance so well, hence this design.*/ 
  kb=0;
  #ifdef OPENMP_ON
#pragma omp parallel for private(j,kb,kk,r,c,rb,cb,rt,ct,tid,i) num_threads(*nthreads) schedule(dynamic)
  #endif
  for (j=0;j<sb[N];j++) { /* the block loop */
    kk = b[j];kb=B[kk];
    //while (kk>=sb[kb+1]) kb++; /* kb is main block */
    rb = R[kb];cb=C[kb]; /* set up allows blocks to be computed in any order by re-arranging B */
    i = kk - sb[kb]; /* sub-block index */
    rt = pt[rb]/pd[rb]; ct = pt[cb]/pd[cb]; /* total rows and cols of sub-blocks */
    /* compute sub-row and column */
    if (sb[kb+1]-sb[kb]<rt*ct) { /* symmetric upper half only needed */ 
      r=0; while (i >= rt - r) { i -= rt - r;r++;}
      c = i + r;
    } else {
      r = i / ct;
      c = i % ct;
    }
    //Rprintf("kb = %d rb = %d cb = %d r = %d c = %d\n",kb,rb,cb,r,c);
    #ifdef OPENMP_ON
    tid = omp_get_thread_num(); /* needed for providing thread specific work space to XWXij */
    #endif
    //Rprintf("r=%d c=%d tid = %d  ",r,c,tid); // dodgy not ts
    XWXij(XWX+tpsu[rb] +  (ptrdiff_t) nxwx * tpsu[cb],rb,cb,r,c,X,k,ks,m,p,*nx,*n,ts, dt,*nt,w,ws,
	  tri,off,work + tid * (ptrdiff_t) q,nxwx); /* compute r,c block */
    /* NOTE: above will write directly to oversized XWX, then have constraints applied post-hoc. */ 
  } /* block loop */

  /* now XWX contains the unconstrained X'WX, but the constraints have to be applied to blocks involving tensor products */
  //Rprintf("\ncleaning up!\n");
  for (r=0;r < *nt;r++) for (c=r;c< *nt;c++) {
    /* if Xr is tensor, may need to apply constraint */
    if (dt[r]>1&&qc[r]>0) { /* first term is a tensor with a constraint */
      /* col by col form (I-vv')xwx, dropping first row... */
      for (j=0;j<pt[c];j++) { /* loop over columns */
        x0 = XWX + tpsu[r] + (tpsu[c]+j) * (ptrdiff_t) nxwx; /* jth col of raw block */
	x1 = XWX + tps[r] + (tpsu[c]+j) * (ptrdiff_t) nxwx;   /* jth col of constrained block */
	for (x=0.0,p0=x0,p1=x0+pt[r],p2=v+voff[r];p0<p1;p0++,p2++) x+= *p0 * *p2; 
	for (p2=v+voff[r]+1,p1=x0+pt[r],x0++;x0<p1;x1++,x0++,p2++) *x1 = *x0 - *p2 * x;
      }
      pa = pt[r]-1;
      //x0=x1=xwx; /* pointers to columns of xwx */ 
      /* col by col form (I-vv')xwx, dropping first row... */
      //for (j=0;j<pt[c];j++) {
      //  for (x=0.0,p0=x0,p1=x0+pt[r],p2=v+voff[r];p0<p1;p0++,p2++) x += *p0 * *p2;
      //  for (p2=v+voff[r]+1,p1=x0+pt[r],x0++;x0<p1;x1++,x0++,p2++) *x1 = *x0 - *p2 * x;
      //}
    } else {
      pa = pt[r];
      if (tpsu[r]!=tps[r]) { /* still need to shift rows upwards */
        for (j=0;j<pt[c];j++) { /* loop over columns */
	  x0 = XWX + tpsu[r] + (tpsu[c]+j) * (ptrdiff_t) nxwx; /* jth col of raw block */
	  x1 = XWX + tps[r] + (tpsu[c]+j) * (ptrdiff_t) nxwx;   /* jth col of shifted block */
	  for (p1 = x0 + pt[r];x0<p1;x0++,x1++) *x1 = *x0;
	}  
      }	
    }  
    if (dt[c]>1&&qc[c]>0) { /* Xc term is a tensor with a constraint */
      /* row by row form xwx(I-vv') dropping first col... */
      for (j=0;j<pa;j++) { /* work down rows */
        x0 = XWX + tps[r] + j + tpsu[c] * (ptrdiff_t) nxwx; /* jth col of raw block */
	x1 = XWX + tps[r] + j + tps[c] * (ptrdiff_t) nxwx;   /* jth col of constrained block */
	for (x=0.0,p0=x0,p1=v+voff[c],p2=p1+pt[c];p1<p2;p0 += nxwx,p1++) x += *p0 * *p1;
	for (p1=v+voff[c]+1,x0+= nxwx;p1<p2;x1 += nxwx,x0+= nxwx,p1++) *x1 = *x0 - *p1 *x;
      }
    } else if (tpsu[c]!=tps[c]) { /* still need to shift cols leftwards */
      for (j=0;j<pa;j++) { /* work down rows */
	x0 = XWX + tps[r] + j + tpsu[c] * (ptrdiff_t) nxwx; /* jth col of raw block */
	x1 = XWX + tps[r] + j + tps[c] * (ptrdiff_t) nxwx;   /* jth col of shifted block */
	for (p1 = x0 + nxwx * pt[c];x0<p1;x0 += nxwx,x1 += nxwx) *x1 = *x0;
      }
    }  
  }

  if (ptot<nxwx) row_squash(XWX,ptot,nxwx,ptot); /* drop the now redundant trailing rows */
  up2lo(XWX,ptot); /* copy upper triangle to lower */
  
  FREE(pt);FREE(pd);FREE(off);FREE(voff);FREE(tps);FREE(tpsu);FREE(work); if (tri) FREE(ws);//FREE(xwx);FREE(xwx0);
  FREE(B);FREE(R);FREE(C);FREE(sb);FREE(Cost);FREE(cost);FREE(b);
} /* XWXd0 */ 


