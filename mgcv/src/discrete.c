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
#ifdef SUPPORT_OPENMP
#include <omp.h>
#endif



/* basic extraction operations */ 

void singleXj(double *Xj,double *X, int *m, int *k, int *n,int *j) {
/* Extract a column j of matrix stored in compact form in X, k into Xj. 
   X has m rows. k is of length n. ith row of result is Xj = X[k(i),j]
   (an n vector). This function is O(n).
*/
  double *pe; 
  X += *m * *j; /* shift to start of jth column */ 
  for (pe = Xj + *n;Xj < pe;Xj++,k++) *Xj = X[*k]; 
} /* singleXj */

void tensorXj(double *Xj, double *X, int *m, int *p,int *dt, 
              int *k, int *n, int *j) {
/* Extract a column j of tensor product term matrix stored in compact 
   form in X, k into Xj. There are dt sub matrices in Xj. The ith is 
   m[i] by p[i]. There are dt index n - vectors stacked end on end in 
   k. This function is O(n*dt)

   This routine performs pure extraction only if Xj is a vector of 1s on 
   entry. Otherwise the jth column is multiplied element wise by the 
   contents of Xj on entry. 
*/ 
  int q=1,l,i,jp;
  double *p0,*p1,*M;
  p1 = Xj + *n; /* pointer for end of Xj */
  for (i = 0;i < *dt;i++) q *= p[i];
  jp = *j; 
  for (i = 0;i < *dt; i++) {
    q /= p[i]; /* update q */
    l = jp/q; /* column of current marginal */
    jp = jp%q;
    M = X + m[i] * l; /* M now points to start of col l of ith marginal model matrix */
    for (p0=Xj;p0<p1;p0++,k++) *p0 *= M[*k];
    X += m[i] * p[i]; /* move to the next marginal matrix */
  }
} /* tensorXj */

void singleXty(double *Xy,double *temp,double *y,double *X, int *m,int *p, int *k, int *n) {
/* forms X'y for a matrix stored in packed form in X, k, with ith row 
   X[k[i],]. X as supplied is m by p, while k is an n vector.
   Xy and temp are respectively p and m vectors and do not need to be cleared on entry.
*/
  double *p0,*p1,done=1.0,dzero=0.0; 
  char trans = 'T';
  int one=1;
  for (p0=temp,p1 = p0 + *m;p0<p1;p0++) *p0 = 0.0;
  for (p1=y + *n;y<p1;y++,k++) temp[*k] += *y;
  F77_CALL(dgemv)(&trans, m, p,&done,X,m,temp,&one,&dzero,Xy,&one);
} /*singleXty*/

void tensorXty(double *Xy,double *work,double *work1, double *y,double *X, 
               int *m, int *p,int *dt,int *k, int *n) {
/* forms X'y for a matrix stored as a compact tensor product in X, k.
   There are dt maginal matrices packed in X, the ith being m[i] by p[i].
   y and work are n vectors. work does not need to be cleared on entry.
   work1 is an m vector where m is the dimension of the final marginal.
   Note: constraint not dealt with here.
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
     * is the row tensor product here */
    for (p0=y,p1=work;p0<yn;p0++,p1++) *p1 = *p0; /* copy y to work */ 
    j = *dt - 1; 
    tensorXj(work,X,m,p,&j,k,n,&i);
    /* now form M'work */
    singleXty(Xy+i*pd,work1,work,M,m + *dt-1,&pd, k + *n * (*dt-1),n);
  }
} /* tensorXty */

void singleXb(double *f,double *work,double *X,double *beta,int *k,int *m, int *p,int *n) {
/* Forms X beta, where peta is a p - vector, and X is stored in compact form in 
   X, k, with ith row X[k[i],]. X is stored as an m by p matrix. k is an n vector.
   work is an m vector.
*/
  char trans='N';
  double done=1.0,dzero=0.0,*p1;
  int one=1;
  F77_CALL(dgemv)(&trans, m, p,&done,X,m,beta,&one,&dzero,work,&one);
  p1 = f + *n;
  for (;f < p1;f++,k++) *f = work[*k];
} /* singleXb */

void tensorXb(double *f,double *X, double *C,double *work, double *beta,
              int *m, int *p,int *dt,int *k, int *n,double *Q,int *qc) {
/* for X* beta where X* is a tensor product term with nt marginal model matrices
   stored in X. ith such is m[i] by p[i]. 
   work is an n vector. C is m[d] by pb.  
   Q is a (prod_i p[i]) by qc constraint matrix
*/ 
  char trans='N';
  int pb=1,md,*kp,*kd,pd,i,j,one=1;
  double *M,done=1.0,dzero=0.0,*p0,*p1,*pf,*pc;
  M = X;
  for (i=0;i<*dt-1;i++) {
    pb *= p[i];
    M += m[i]*p[i]; /* final marginal */ 
  }
  md = m[*dt - 1];
  pd = p[*dt -1];
  kd = k + (*dt-1) * *n;
  /* form work = M B, where vec(B) = beta */
  if (*qc<=0) { /* no constraint supplied */
    F77_CALL(dgemm)(&trans,&trans,&md,&pb, &pd, &done,
		    M,&md,beta,&pd,&dzero,C,&md); 
  } else { /* their is a constraint matrix */
    /* first map supplied beta to unconstrained parameterization */ 
    j = pb*pd;
    F77_CALL(dgemv)(&trans, &j, qc,&done,Q,&j,beta,&one,&dzero,work,&one);
    F77_CALL(dgemm)(&trans,&trans,&md,&pb, &pd, &done,
		    M,&md,work,&pd,&dzero,C,&md);
  }
  p1 = work + *n;
  for (pf=f,p0=f + *n;pf<p0;pf++) *pf = 0.0;
  for (i=0;i<pb;i++) {
    /* extract ith col of truncated tensor product (final marginal ommited) */
    for (p0=work;p0<p1;p0++) *p0 = 1.0; /* set work to 1 */ 
    j = *dt - 1; 
    tensorXj(work,X,m,p,&j,k,n,&i);
    pc = C + i * md; /* ith col of C */
    for (kp=kd,pf=f,p0=work;p0<p1;p0++,pf++,kp++) {
      *pf += pc[*kp] * *p0;
    }
  }
} /* tensorXb */


void XWy(double *XWy,double *y,double *X,double *w,int *k, int *m,int *p, int *n, 
         int *nx, int *ts, int *dt, int *nt,double *Q,int *qc) {
  double *Wy,*p0,*p1,*p2,done=1.0,dzero=0.0,*Xy0,*work,*work1;
  int q,i,j,*pt,*off,*Qoff,*tps,maxm,maxp,one=1;
  char trans='T';
  /* obtain various indices */
  pt = (int *) R_chk_calloc((size_t)*nt,sizeof(int)); /* the term dimensions */
  off = (int *) R_chk_calloc((size_t)*nx+1,sizeof(int)); /* offsets for X submatrix starts */
  Qoff = (int *) R_chk_calloc((size_t)*nt+1,sizeof(int)); /* offsets for Q submatrix starts */
  tps = (int *) R_chk_calloc((size_t)*nt+1,sizeof(int)); /* the term starts */
  for (q=i=0;i< *nt; i++) { /* work through the terms */
    for (j=0;j<dt[i];j++,q++) { /* work through components of each term */
      off[q+1] = off[q] + p[q]*m[q]; /* submatrix start offsets */
      if (j==0) pt[i] = p[q]; else pt[i] *= p[q]; /* term dimension */
      if (maxm<m[q]) maxm=m[q];
    } 
    Qoff[i+1] = Qoff[i] + qc[i] * pt[i]; /* start of ith Q matrix */
    if (maxp<pt[i]) maxp=pt[i];
    if (qc[i]<=0) tps[i+1] = tps[i] + pt[i]; /* where ith terms starts in param vector */ 
    else tps[i+1] = tps[i] + qc[i]; /* there is a tensor constraint to apply - reducing param count*/
  }
  Xy0 =  (double *) R_chk_calloc((size_t)maxp,sizeof(double));
  work =  (double *) R_chk_calloc((size_t)*n,sizeof(double));
  work1 = (double *) R_chk_calloc((size_t)maxm,sizeof(double));
  /* apply W to y - for AR models this needs to be expanded */
  Wy = (double *) R_chk_calloc((size_t)*n,sizeof(double)); /* Wy */
  for (p0=Wy,p1=Wy + *n;p0<p1;p0++,y++,w++) *p0 = *y * *w;
  /* now loop through terms applying the components of X'...*/
  for (i=0;i<*nt;i++) {
    if (dt[i]>1) { /* it's a tensor */
      tensorXty(Xy0,work,work1,Wy,X+off[ts[i]],m+ts[i],p+ts[i],dt+i,k+ts[i] * *n,n);
      if (qc[i]>0) { /* there is a constraint to apply */
        F77_CALL(dgemv)(&trans, pt+i, qc+i,&done,Q+Qoff[i],pt+i,Xy0,&one,&dzero,XWy + tps[i],&one);
      } else { /* straight copy */
        for (p0=Xy0,p1=p0+pt[i],p2=XWy+tps[i];p0<p1;p0++,p2++) *p2 = *p0;
      }
    } else { /* it's a singleton */
      //singleXty(double *Xy,double *temp,double *y,double *X, int *m,int *p, int *k, int *n)
      singleXty(XWy+tps[i],work1,Wy,X+off[ts[i]], m+ts[i],p+ts[i], k+ts[i] * *n,n);
    }
  }
  R_chk_free(Wy); R_chk_free(Xy0); R_chk_free(work); R_chk_free(work1); 
  R_chk_free(pt); R_chk_free(off); R_chk_free(Qoff); R_chk_free(tps);
} /* XWy */

void XWX(double *XWX,double *X,double *w,int *k, int *m,int *p, int *n, int *nx, int *ts, int *dt, int *nt,
         double *Q,int *qc,int *nthreads) {
/* This version uses open MP by column, within sub-block.

   An alternative in which one thread did each sub-block scaled poorly, largely 
   because the work is very uneven between sub-blocks. 

   Forms Xt'WXt when Xt is divided into blocks of columns, each stored in compact form
   using arguments X and k. 'X' contains 'nx' blocks, the ith is an m[i] by p[i] matrix. 
   k contains nx index vectors. If kj and Xj represent the jth sub-matrix and index vector 
   then the ith row of the corresponding full sub-matrix is Xj[kj[i],]. The blocks of Xt 
   need not be directly given by terms like Xj[kj[i],], some are actually given by row 
   tensors of several such matrices. Conceptually Xt has nt < nx blocks (terms) and the ith 
   block starts at sub-matriux  ts[i] of X, and is made up of the row tensor produce of 
   dt[i] consecutive sub-matrices. 
    
   Tensor product terms may have constraint matrices Q, which post multiply the tensor product 
   (typically imposing approximate sum-to-zero constraints). qc contains the number of columns 
   Q for each of the nt terms. It should be zero for a singleton term or a tensor product with 
   no constraints. Otherwise qc contains the number of columns of Q - typically one less than 
   the total number of coeffs for the tensor prodoct. 
  
   NOTE: ideally tensor product constraints would be implemented using single householder 
   rotations.  
*/  
  int r,c,i,j,q,*pt,*pd,*off,a,b,*tps,ptot,maxp=0,maxm=0,*Qoff,pa,pb,kk,dk,rk,*start; 
  double done=1.0,dzero=0.0,*p0,*p1,*p2, *Xi,*temp,*tempn,*xwx,*xwx0,
         *XiB,*tempB,*tempnB;
  char trans='T',not_trans='N';
  #ifndef SUPPORT_OPENMP
  *nthreads = 1;
  #endif
  if (*nthreads<1) *nthreads = 1;
  start = (int *) R_chk_calloc((size_t)(*nthreads+1),sizeof(int));
  XiB = (double *) R_chk_calloc((size_t)(*n * *nthreads),sizeof(double)); /* column of X */
  tempnB = (double *) R_chk_calloc((size_t)(*n * *nthreads),sizeof(double));
  pt = (int *) R_chk_calloc((size_t)*nt,sizeof(int)); /* the term dimensions */
  pd = (int *) R_chk_calloc((size_t)*nt,sizeof(int)); /* storage for last marginal size */
  off = (int *) R_chk_calloc((size_t)*nx+1,sizeof(int)); /* offsets for X submatrix starts */
  Qoff = (int *) R_chk_calloc((size_t)*nt+1,sizeof(int)); /* offsets for Q submatrix starts */
  tps = (int *) R_chk_calloc((size_t)*nt+1,sizeof(int)); /* the term starts */
  for (q=i=0;i< *nt; i++) { /* work through the terms */
    for (j=0;j<dt[i];j++,q++) {
      if (j==dt[i]-1) pd[i] = p[q]; /* the relevant dimension for deciding product ordering */
      off[q+1] = off[q] + p[q]*m[q]; /* submatrix start offsets */
      if (j==0) pt[i] = p[q]; else pt[i] *= p[q]; /* term dimension */
      if (maxm<m[q]) maxm=m[q];
    } 
    Qoff[i+1] = Qoff[i] + qc[i] * pt[i]; /* start of ith Q matrix */
    if (maxp<pt[i]) maxp=pt[i];
    if (qc[i]<=0) tps[i+1] = tps[i] + pt[i]; /* where ith terms starts in param vector */ 
    else tps[i+1] = tps[i] + qc[i]; /* there is a tensor constraint to apply - reducing param count*/
  }
  tempB = (double *) R_chk_calloc((size_t)(maxm * *nthreads),sizeof(double));
  xwx = (double *) R_chk_calloc((size_t)(maxp*maxp),sizeof(double)); /* working cross product storage */
  xwx0 = (double *) R_chk_calloc((size_t)(maxp*maxp),sizeof(double)); /* working cross product storage */
  ptot = tps[*nt]; /* total number of parameters */
  for (r=0;r < *nt;r++) for (c=r;c< *nt;c++) { /* the block loop */
    if (pd[r]>pd[c]) { /* Form Xr'WXc */
      a=r;b=c;
    } else { /* Form Xc'WXr */
      a=c;b=r; 
    }
    /* split cols between threads... */  
    dk = pt[b] / *nthreads; rk = pt[b] % *nthreads;
    if (dk * *nthreads < pt[b]) dk++;start[0]=0; 
    for (i=0;i<*nthreads;i++) { 
      start[i+1] = start[i] + dk;
      if (start[i+1]>pt[b]) start[i+1]=pt[b];
    }
    #ifdef SUPPORT_OPENMP
    #pragma omp parallel private(Xi,i,temp,tempn,p0,p1,p2) num_threads(*nthreads)
    #endif 
    { /* begin parallel section */
      #ifdef SUPPORT_OPENMP
      #pragma omp for
      #endif
      for (kk=0;kk<*nthreads;kk++) { 
        /* allocate thread specific storage... */
        temp = tempB + kk * maxm;
        Xi = XiB + kk * *n;
        tempn = tempnB + kk * *n;
        for (i=start[kk];i<start[kk+1];i++) { /* loop over columns of Xb */ 
          /* extract Xb[:,i]... */
          if (dt[b]>1) { /* tensor */
            for (p0=Xi,p1=p0+*n;p0<p1;p0++) *p0 = 1.0; /* set Xi to 1 */
            tensorXj(Xi,X+off[ts[b]], m + ts[b] , p + ts[b], dt+b, 
                     k + *n * ts[b], n,&i);
          } else { /* singleton */
            singleXj(Xi,X+off[ts[b]], m + ts[b], k + *n * ts[b], n,&i);
          }
          /* apply W to Xi - for AR models this needs to be expanded */
          for (p0=w,p1=w + *n,p2=Xi;p0<p1;p0++,p2++) *p2 *= *p0;
          /* now form Xa'WXb[:,i]... */  
          if (dt[a]>1) { /* tensor */
            tensorXty(xwx + i * pt[a],tempn,temp,Xi,X+off[ts[a]],m+ts[a],p+ts[a],
                      dt+a,k+ *n * ts[a], n);
          } else { /* singleton */
            singleXty(xwx + i * pt[a],temp,Xi,X+off[ts[a]],m+ts[a],p+ts[a],k + *n * ts[a],n);
          }  
        } /* loop over columns of Xb */
      }
      /* so now xwx contains pt[a] by pt[b] matrix Xa'WXb */
    } /* end parallel section */

    /* if Xb is tensor, may need to apply constraint */
    if (dt[a]>1&&qc[a]>0) { /* first term is a tensor with a constraint */
      /* copy xwx to xwx0 */
      for (p0=xwx,p1=p0 + pt[b]*pt[a],p2=xwx0;p0<p1;p0++,p2++) *p2 = *p0;
      /* apply Qa' */
      F77_CALL(dgemm)(&trans,&not_trans,qc+a,pt+b,pt+a, &done,
		    Q+Qoff[a],pt+a,xwx0,pt+a,&dzero,xwx,qc+a);
      pa = qc[a]; /* rows of xwx */
    } else pa = pt[a];
    if (dt[b]>1&&qc[b]>0) { /* second term is a tensor with a constraint */
      /* copy xwx to xwx0 */
      for (p0=xwx,p1=p0 + pt[b]*pa,p2=xwx0;p0<p1;p0++,p2++) *p2 = *p0;
      /* apply Qa from right */
      F77_CALL(dgemm)(&not_trans,&not_trans,&pa,qc+b,pt+b, &done,
		      xwx0,&pa, Q+Qoff[b],pt+b,&dzero,xwx,&pa);
      pb = qc[b]; /* rows of xwx */
    } else pb=pt[b];
    /* copy result into overall XWX*/  
  
    if (pd[r]>pd[c]) { /* xwx = Xr'WXc */
      for (i=0;i<pa;i++) for (j=0;j<pb;j++) 
      XWX[i+tps[r]+(j+tps[c])*ptot] = XWX[j+tps[c]+(i+tps[r])*ptot] = xwx[i + pa * j];
    } else { /* Form xwx = Xc'WXr */
      for (i=0;i<pa;i++) for (j=0;j<pb;j++) 
      XWX[i+tps[c]+(j+tps[r])*ptot] = XWX[j+tps[r]+(i+tps[c])*ptot] = xwx[i + pa * j];
    }
  } /* end of block loop */
  R_chk_free(start);
  R_chk_free(XiB); R_chk_free(tempnB); R_chk_free(pt); R_chk_free(pd); R_chk_free(off);
  R_chk_free(Qoff); R_chk_free(tps); R_chk_free(tempB); R_chk_free(xwx); R_chk_free(xwx0);
  //R_chk_free(); R_chk_free(); R_chk_free(); R_chk_free(); R_chk_free();
} /* XWX */

