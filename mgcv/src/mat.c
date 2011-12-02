/* Convenient C wrappers for calling LAPACK and LINPACK 

   See R-x/include/R_ext/Lapack.h...
*/
#include "mgcv.h"
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
/*#include <dmalloc.h>*/

void mgcv_mmult0(double *A,double *B,double *C,int *bt,int *ct,int *r,int *c,int *n)
/* This code doesn't rely on the BLAS...
 
   Forms r by c product of B and C, transposing each according to bt and ct.
   n is the common dimension of the two matrices, which are stored in R 
   default column order form. Algorithm is inner loop optimized in each case 
   (i.e. inner loops only update pointers by 1, rather than jumping).
   
  
   r<-1000;c<-1000;n<-1000
   A<-matrix(0,r,c);B<-matrix(rnorm(r*n),n,r);C<-matrix(rnorm(n*c),c,n)
   system.time(
   A<-matrix(.C("mgcv_mmult",as.double(A),as.double(B),as.double(C),as.integer(1),as.integer(1),
          as.integer(r),as.integer(c),as.integer(n))[[1]],r,c))
   range(as.numeric(t(B)%*%t(C)-A))
   A<-matrix(0,r,c);B<-matrix(rnorm(r*n),n,r);C<-matrix(rnorm(n*c),n,c)
   system.time(
   A<-matrix(.C("mgcv_mmult",as.double(A),as.double(B),as.double(C),as.integer(1),as.integer(0),
          as.integer(r),as.integer(c),as.integer(n))[[1]],r,c))
   range(as.numeric(t(B)%*%C-A))
   A<-matrix(0,r,c);B<-matrix(rnorm(r*n),r,n);C<-matrix(rnorm(n*c),c,n)
   system.time(
   A<-matrix(.C("mgcv_mmult",as.double(A),as.double(B),as.double(C),as.integer(0),as.integer(1),
          as.integer(r),as.integer(c),as.integer(n))[[1]],r,c))
   range(as.numeric(B%*%t(C)-A))
   A<-matrix(0,r,c);B<-matrix(rnorm(r*n),r,n);C<-matrix(rnorm(n*c),n,c)
   system.time(
   A<-matrix(.C("mgcv_mmult",as.double(A),as.double(B),as.double(C),as.integer(0),as.integer(0),
          as.integer(r),as.integer(c),as.integer(n))[[1]],r,c))
   range(as.numeric(B%*%C-A))  
*/

{ double xx,*bp,*cp,*cp1,*cp2,*cp3,*ap,*ap1;
  int br,cr,i,j;
  if (*bt)
  { if (*ct) /* A=B'C' */
    { /* this one is really awkward: have to use first row of C' as working storage
         for current A row; move most slowly through B */
      for (i=0;i<*r;i++) /* update A *row-wise*, one full sweep of C per i */    
      { cp1 = C + *c;
        /* back up row 1 of C' (in A), and initialize current row of A (stored in row 1 of C') */
        xx = *B;
        for (ap=A,cp=C;cp<cp1;cp++,ap+= *r) { *ap = *cp; *cp *= xx;}  
        B++;cp2=cp1;
        for (j=1;j< *n;j++,B++) 
	    for (xx= *B,cp=C;cp<cp1;cp++,cp2++) *cp += xx * *cp2;      
        /* row i of A is now in row 1 of C', need to swap back */
        for (ap=A,cp=C;cp<cp1;cp++,ap+= *r) { xx = *ap; *ap = *cp; *cp = xx;}
        A++;
      } 
    } else /* A=B'C - easiest case: move most slowly through A*/
    { br= *n;cr= *n;cp2 = C + *c * cr;
      for (ap=A,cp1=C;cp1< cp2;cp1+=cr) for (bp=B,i=0;i< *r;i++,ap++)  
      { for (xx=0.0,cp=cp1,cp3=cp1+ *n;cp< cp3;cp++,bp++) xx += *cp * *bp; /* B[k+br*i]*C[k+cr*j];*/
        *ap=xx;
      }
    }
  } else
  { if (*ct) /* A=BC' - update A column-wise, move most slowly through C (but in big steps)*/
    { cp = C;
      for (j=0;j< *c;j++) /* go through columns of A, one full B sweep per j */
      { ap1 = A + *r;C=cp;xx = *C;bp=B;
        for (ap=A;ap<ap1;ap++,bp++) *ap = xx * *bp ;
        C += *c;
        for (i=1;i< *n;i++,C+= *c) 
	for (xx=*C,ap=A;ap<ap1;ap++,bp++) *ap += xx * *bp;
        A=ap1;cp++;
      }
    } else /* A=BC - update A column-wise, moving most slowly through C */
    { for (j=0;j< *c;j++) /* go through columns of A, one full B sweep per j */
      { ap1 = A + *r;xx = *C;bp=B;
        for (ap=A;ap<ap1;ap++,bp++) *ap = xx * *bp ;
        C++;
        for (i=1;i< *n;i++,C++) 
	for (xx=*C,ap=A;ap<ap1;ap++,bp++) *ap += xx * *bp;
        A=ap1;
      }
    }
  }
} /* end mgcv_mmult0*/

void mgcv_mmult(double *A,double *B,double *C,int *bt,int *ct,int *r,int *c,int *n) {
  /* BLAS version A is c, B is a, C is b, bt is transa ct is transb 
     r is m, c is n, n is k,
  */
  char transa='N',transb='N';
  int lda,ldb,ldc;
  double alpha=1.0,beta=0.0;
  if (*bt) { /* so B is n by r */
    transa = 'T';  
    lda = *n;
  } else  lda = *r; /* B is r by n */

  if (*ct) { /* C is c by n */ 
    transb = 'T';
    ldb = *c;
  } else ldb = *n; /* C is n by c */
  
  ldc = *r;
  F77_NAME(dgemm)(&transa,&transb,r,c,n, &alpha,
		B, &lda,C, &ldb,&beta, A, &ldc);
} /* end mgcv_mmult */

void getXtX0(double *XtX,double *X,int *r,int *c)
/* form X'X (nearly) as efficiently as possible - BLAS free*/
{ double *p0,*p1,*p2,*p3,*p4,x;
  int i,j;
  for (p0=X,i=0;i<*c;i++,p0 += *r) 
  for (p1=X,j=0;j<=i;j++,p1 += *r) {
    for (x=0.0,p2=p0,p3=p1,p4=p0 + *r;p2<p4;p2++,p3++) x += *p2 * *p3;    
    XtX[i + j * *c] = XtX[j + i * *c] = x;
  }
}

void getXtX(double *XtX,double *X,int *r,int *c)
/* form X'X (nearly) as efficiently as possible - uses BLAS*/
{ double alpha=1.0,beta=0.0;
  int i,j;
  char uplo = 'L',trans='T';
  F77_NAME(dsyrk)(&uplo,&trans,c, r, & alpha,X,r,&beta,XtX,c);
  /* fill in upper triangle from lower */
  for (i=0;i<*c;i++) 
  for (j=0;j<i;j++)  XtX[j + i * *c] = XtX[i + j * *c];
 
}


void getXtWX0(double *XtWX, double *X,double *w,int *r,int *c,double *work)
/* Original BLAS free version...
   forms X'WX as efficiently as possible, where W = diag(w)
   and X is an r by c matrix stored column wise. 
   work should be an r-vector (longer is no problem).
*/ 
{ int i,j;
  double *p,*p1,*p2,*pX0,*pX1,xx;
  pX0=X;
  for (i=0;i< *c;i++) { 
    p2 = work + *r;
    for (p=w,p1=work;p1<p2;p++,p1++,pX0++) *p1 = *pX0 * *p; 
    for (pX1=X,j=0;j<=i;j++) {
      for (xx=0.0,p=work;p<p2;p++,pX1++) xx += *p * *pX1;
      XtWX[i * *c + j] = XtWX[j * *c + i] = xx;
    }
  }
}

void getXtWX(double *XtWX, double *X,double *w,int *r,int *c,double *work)
/* BLAS version...
   forms X'WX as efficiently as possible, where W = diag(w)
   and X is an r by c matrix stored column wise. 
   work should be an r-vector (longer is no problem).
*/ 
{ int i,j,one=1;
  char trans='T';
  double *p,*p1,*p2,*pX0,xx,*w2,alpha=1.0,beta=0.0;
  pX0=X;
  w2 = XtWX; /* use first column as work space until end */
  for (i=0;i< *c;i++) { 
    p2 = work + *r;
    /* form work=WX[,i]... */
    for (p=w,p1=work;p1<p2;p++,p1++,pX0++) *p1 = *pX0 * *p;  
    /* Now form X[,1:i]'work ... */
    j = i+1; /* number of columns of X to use */
    F77_NAME(dgemv)(&trans, r, &j,&alpha,X, r,work,&one,&beta,w2, &one);
    if (i==0) xx = w2[0]; /* save the 0,0 element of XtWX (since its in use as workspace) */
    else for (j=0;j<=i;j++) XtWX[i * *c + j] = w2[j];
  }
  /* now fill in the other triangle */
  *XtWX = xx;
  for (i=0;i< *c;i++) for (j=0;j<i;j++) XtWX[j * *c + i] =  XtWX[i * *c + j];
}


void getXtMX(double *XtMX,double *X,double *M,int *r,int *c,double *work)
/* BLAS free version (currently no BLAS version as this is only used on 
   small matrices).
   forms X'MX as efficiently as possible, where M is a symmetric matrix
   and X is an r by c matrix. X and M are stored column wise. 
   work should be an r-vector (longer is no problem).
*/

{ int i,j;
  double *p,*p1,*p2,*pX0,*pX1,xx,*pM;
  pX0=X;
  for (i=0;i< *c;i++) { 
    /* first form MX[:,i] */
    p2 = work + *r;pM=M;
    for (p1=work;p1<p2;pM++,p1++) *p1 = *pX0 * *pM;pX0++;
    for (j=1;j< *r;j++,pX0++) 
    for (p1=work;p1<p2;pM++,p1++) *p1 += *pX0 * *pM;
    /* now form ith row and column of X'MX */
    for (pX1=X,j=0;j<=i;j++) {
      for (xx=0.0,p=work;p<p2;p++,pX1++) xx += *p * *pX1;
      XtMX[i * *c + j] = XtMX[j * *c + i] = xx;
    }
  }
}



void mgcv_chol(double *a,int *pivot,int *n,int *rank)
/* a stored in column order, this routine finds the pivoted choleski decomposition of matrix a 
   library(mgcv)
   X<-matrix(rnorm(16),4,4);D<-X%*%t(X)
   rank<-0
   er<-.C("mgcv_chol",as.double(D),as.integer(rep(0,4)),as.integer(4),as.integer(rank))
   rD<-matrix(er[[1]],4,4);piv<-er[[2]]
   chol(D,pivot=TRUE);rD;piv
   n<-length(piv);ind<-1:n;
   ind[piv]<-1:n
   rD<-rD[,ind]
   L<-mroot(D)
   D;t(rD)%*%rD;L%*%t(L)
*/
{ double *work,*p1,*p2,*p;
  int piv=1;
  work=(double *)calloc((size_t) *n,sizeof(double));
  F77_NAME(dchdc)(a,n,n,work,pivot,&piv,rank);
  /* zero stuff below the leading diagonal */
  for (p2=a+ *n,p1=a+1;p2<a+ *n * *n;p1+= *n+1,p2+= *n) for (p=p1;p<p2;p++) *p=0.0;
  free(work);
}

void mroot(double *A,int *rank,int *n)
/* finds the minimum rank or supplied rank square root of n by n matrix  A by pivoted choleski 
   decomposition returned in A is B such that B'B = A (this is different from R routine mroot)

   R testing code....
   library(mgcv)
   X<-matrix(rnorm(12),4,3);D<-X%*%t(X)
   rank<-0
   er<-.C("mroot",as.double(D),as.integer(rank),as.integer(4))
   rD<-matrix(er[[1]],er[[2]],4)
   D;t(rD)%*%rD
   
*/
{ int *pivot,erank,i,j;
  double *pi,*pj,*p0,*p1,*B;
  pivot=(int *)calloc((size_t)*n,sizeof(int));
  mgcv_chol(A,pivot,n,&erank);
  if (*rank<=0) *rank = erank;
  /* unscramble the returned choleski factor */ 
  B=calloc((size_t)(*n * *n),sizeof(double));
  /* copy A to B and zero A */
  for (p0=A,p1=B,i=0;i< *n;i++,p0+= *n,p1+= *n)
  for (pi=p0,pj=p1;pi<=p0+i;pi++,pj++) { *pj = *pi; *pi=0.0;}
  /* now copy B into A undoing pivoting */
  for (p0=B,i=0;i< *n;p0+= *n,i++) 
  for (j=pivot[i],pj=p0,pi=A + (j-1)* *n;pj<=p0+i;pi++,pj++) *pi = *pj;
  /* remove trailing rows from choleski factor */
  for (pi=A,p0=A,i=0;i< *n;i++,p0+= *n)
  for (pj=p0;pj<p0+ *rank;pj++,pi++) *pi = *pj; 
  free(pivot);free(B);
}


void mgcv_svd(double *x,double *u,double *d,int *r,int *c)
/* call LA_PACK svd routine to form x=UDV'. Assumes that V' not wanted, but that full square U is required.
*/
{ const char jobu='A',jobvt='N';
  int lda,ldu,ldvt=1,lwork;
  int info;
  double *vt=NULL,work1,*work;
  ldu=lda= *r;
  lwork=-1;
  /* workspace query */
  F77_NAME(dgesvd)(&jobu,&jobvt, r, c, x, &lda, d, u, &ldu, vt,&ldvt,
  		   &work1, &lwork, &info);
  lwork=(int)floor(work1);
  if (work1-lwork>0.5) lwork++;
  work=(double *)calloc((size_t)lwork,sizeof(double));
  /* actual call */
  F77_NAME(dgesvd)(&jobu,&jobvt, r, c, x, &lda, d, u, &ldu, vt,&ldvt,
  		   work, &lwork, &info);
  free(work);
}

void mgcv_svd_full(double *x,double *vt,double *d,int *r,int *c)
/* call LA_PACK svd routine to form x=UDV'. U returned in x. V' returned in vt.
   assumed r >= c. U is r by c. D is length c. V is c by c.

# Here is R test code.....
library(mgcv)
n<-4;q<-3
X<-matrix(rnorm(n*q),n,q)
um<-.C("mgcv_svd_full",as.double(X),double(q*q),double(q),as.integer(n),as.integer(q),
        PACKAGE="mgcv")
er<-La.svd(X)
matrix(um[[1]],n,q);er$u
um[[3]];er$d
matrix(um[[2]],q,q);er$v
*/
{ const char jobu='O',jobvt='A';
  int lda,ldu,ldvt,lwork;
  int info;
  double work1,*work,*u=NULL;
  ldu=lda= *r;ldvt = *c;
  lwork=-1;
  /* workspace query */
  F77_NAME(dgesvd)(&jobu,&jobvt, r, c, x, &lda, d, u, &ldu, vt,&ldvt,
  		   &work1, &lwork, &info);
  lwork=(int)floor(work1);
  if (work1-lwork>0.5) lwork++;
  work=(double *)calloc((size_t)lwork,sizeof(double));
  /* actual call */
  F77_NAME(dgesvd)(&jobu,&jobvt, r, c, x, &lda, d, u, &ldu, vt,&ldvt,
  		   work, &lwork, &info);
  free(work);
}


void mgcv_td_qy(double *S,double *tau,int *m,int *n, double *B,int *left,int *transpose)
/* Multiplies m by n matrix B by orthogonal matrix returned from mgcv_tri_diag and stored in 
   S, tau. B is overwritten with result. 
    
   Note that this is a bit inefficient if really only a few rotations matter!

   Calls LAPACK routine dormtr 
*/
{ char trans='N',side='R',uplo='U';
  int nq,lwork=-1,info;
  double *work,work1;
  if (*left) { side = 'L';nq = *m;} else nq = *n;
  if (*transpose) trans = 'T';
  /* workspace query ... */
  F77_NAME(dormtr)(&side,&uplo,&trans,m,n,S,&nq,tau,B,m,&work1,&lwork,&info);
  lwork=(int)floor(work1);if (work1-lwork>0.5) lwork++;
  work=(double *)calloc((size_t)lwork,sizeof(double));
  /* actual call ... */
  F77_NAME(dormtr)(&side,&uplo,&trans,m,n,S,&nq,tau,B,m,work,&lwork,&info);
  free(work);
}

void mgcv_tri_diag(double *S,int *n,double *tau)
/* Reduces symmetric n by n matrix S to tridiagonal form T 
   by similarity transformation. Q'SQ = T, where Q is an
   orthogonal matrix. Only the upper triangle of S is actually
   used.

   On exit the diagonal and superdiagonal of T are written in 
   the corresponding position in S. The elements above the first 
   superdiagonal, along with tau, store the householder reflections 
   making up Q. 

   Note that this is not optimally efficient if actually only a 
   few householder rotations are needed because S does not have 
   full rank.

   The routine calls dsytrd from LAPACK.
     
*/
{ int lwork=-1,info;
  double *work,work1,*e,*d;
  char uplo='U';
  d = (double *)calloc((size_t)*n,sizeof(double));
  e = (double *)calloc((size_t)*n-1,sizeof(double));
  /* work space query... */
  F77_NAME(dsytrd)(&uplo,n,S,n,d,e,tau,&work1,&lwork,&info);
  lwork=(int)floor(work1);if (work1-lwork>0.5) lwork++;
  work=(double *)calloc((size_t)lwork,sizeof(double));
  /* Actual call... */
  F77_NAME(dsytrd)(&uplo,n,S,n,d,e,tau,work,&lwork,&info);
  free(work);free(d);free(e);
}


void mgcv_backsolve0(double *R,int *r,int *c,double *B,double *C, int *bc) 
/* BLAS free version
   Finds C = R^{-1} B where R is the c by c matrix stored in the upper triangle 
   of r by c argument R. B is c by bc. (Possibility of non square argument
   R facilitates use with output from mgcv_qr). This is just a standard back 
   substitution loop.
*/  
{ int i,j,k;
  double x,*pR,*pC;
  for (j=0;j<*bc;j++) { /* work across columns of B & C */
    for (i = *c-1;i>=0;i--) { /* work up each column of B & C */
      x = 0.0;
      /* for (k=i+1;k<*c;k++) x += R[i + *r * k] * C[k + j * *c]; ...following replaces...*/
      pR = R + i + (i+1) * *r;pC = C + j * *c + i + 1;
      for (k=i+1;k<*c;k++,pR+= *r,pC++) x += *pR * *pC;      
      C[i + j * *c] = (B[i + j * *c] - x)/R[i + *r * i];
    }
  }
}

void mgcv_backsolve(double *R,int *r,int *c,double *B,double *C, int *bc) 
/* BLAS version
   Finds C = R^{-1} B where R is the c by c matrix stored in the upper triangle 
   of r by c argument R. B is c by bc. (Possibility of non square argument
   R facilitates use with output from mgcv_qr). This is just a standard back 
   substitution loop.
*/  
{ double *pR,*pC,alpha=1.0;
  char side='L',uplo='U',transa='N',diag='N';
  for (pC=C,pR=pC+ *bc * *c;pC<pR;pC++,B++) *pC = *B; /* copy B to C */
  F77_NAME(dtrsm)(&side,&uplo,&transa, &diag,c, bc, &alpha,R, r,C,c);
}


void mgcv_forwardsolve0(double *R,int *r,int *c,double *B,double *C, int *bc) 
/* BLAS free version
   Finds C = R^{-T} B where R is the c by c matrix stored in the upper triangle 
   of r by c argument R. B is c by bc. (Possibility of non square argument
   R facilitates use with output from mgcv_qr). This is just a standard forward 
   substitution loop.
*/  
{ int i,j,k;
  double x;
  for (j=0;j<*bc;j++) { /* work across columns of B & C */
    for (i = 0;i< *c;i++) { /* work down each column of B & C */
      x=0.0;
      for (k=0;k<i;k++) x+= C[k + j * *c] * R[k + i * *r]; 
      C[i + j * *c] = (B[i + j * *c] - x) / R[i + i * *r]; 
    }
  }
}

void mgcv_forwardsolve(double *R,int *r,int *c,double *B,double *C, int *bc) 
/* BLAS version
   Finds C = R^{-T} B where R is the c by c matrix stored in the upper triangle 
   of r by c argument R. B is c by bc. (Possibility of non square argument
   R facilitates use with output from mgcv_qr). This is just a standard forward 
   substitution loop.
*/  
{ double *pR,*pC,alpha=1.0;
  char side='L',uplo='U',transa='T',diag='N';
  for (pC=C,pR=pC+ *bc * *c;pC<pR;pC++,B++) *pC = *B; /* copy B to C */
  F77_NAME(dtrsm)(&side,&uplo,&transa, &diag,c, bc, &alpha,R, r,C,c);
}


void mgcv_qr(double *x, int *r, int *c,int *pivot,double *tau)
/* call LA_PACK to get pivoted QR decomposition of x
   tau is an array of length min(r,c)
   pivot is array of length c, zeroed on entry, pivoting order on return.
   On exist upper triangle of x is R. Below upper triangle plus tau 
   represent reflectors making up Q.
   pivoting is always performed (not just when matrix is rank deficient), so
   leading diagonal of R is in descending order of magnitude.
   library(mgcv)
   r<-4;c<-3
   X<-matrix(rnorm(r*c),r,c)
   pivot<-rep(1,c);tau<-rep(0,c)
   um<-.C("mgcv_qr",as.double(X),as.integer(r),as.integer(c),as.integer(pivot),as.double(tau))
   qr.R(qr(X));matrix(um[[1]],r,c)[1:c,1:c]
*/
{ int info,lwork=-1,*ip;
  double work1,*work;
  /* workspace query */
  F77_NAME(dgeqp3)(r,c,x,r,pivot,tau,&work1,&lwork,&info);
  lwork=(int)floor(work1);if (work1-lwork>0.5) lwork++;
  work=(double *)calloc((size_t)lwork,sizeof(double));
   /* actual call */
  F77_NAME(dgeqp3)(r,c,x,r,pivot,tau,work,&lwork,&info); 
  free(work);
  /*if (*r<*c) lwork= *r; else lwork= *c;*/ 
  for (ip=pivot;ip < pivot + *c;ip++) (*ip)--;
  /* ... for 'tis C in which we work and not the 'cursed Fortran... */
  
}


void mgcv_qrqy(double *b,double *a,double *tau,int *r,int *c,int *k,int *left,int *tp)
/* applies k reflectors of Q of a QR decomposition to r by c matrix b.
   Apply Q from left if left!=0, right otherwise.
   Transpose Q only if tp!=0.
   Information about Q has been returned from mgcv_qr, and is stored in tau and 
   below the leading diagonal of a.
   library(mgcv)
   r<-4;c<-3
   X<-matrix(rnorm(r*c),r,c)
   qrx<-qr(X)
   pivot<-rep(1,c);tau<-rep(0,c)
   um<-.C("mgcv_qr",a=as.double(X),as.integer(r),as.integer(c),as.integer(pivot),tau=as.double(tau))
   y<-1:4;left<-1;tp<-0;cy<-1
   er<-.C("mgcv_qrqy",as.double(y),as.double(um$a),as.double(um$tau),as.integer(r),as.integer(cy),as.integer(c),
        as.integer(left),as.integer(tp),PACKAGE="mgcv")
   er[[1]];qr.qy(qrx,y)
   
*/

{ char side='L',trans='N';
  int lda,lwork=-1,info;
  double *work,work1;
 
  if (! *left) { side='R';lda = *c;} else lda= *r;
  if ( *tp) trans='T'; 
  /* workspace query */
  F77_NAME(dormqr)(&side,&trans,r,c,k,a,&lda,tau,b,r,&work1,&lwork,&info);
  lwork=(int)floor(work1);if (work1-lwork>0.5) lwork++;
  work=(double *)calloc((size_t)lwork,sizeof(double));
  /* actual call */
  F77_NAME(dormqr)(&side,&trans,r,c,k,a,&lda,tau,b,r,work,&lwork,&info); 
  free(work);
   
}


void update_qr(double *Q,double *R,int *n, int *q,double *lam, int *k)
/* Let X=QR where X is n by q and R is q by q upper triangular and Q is n by q.
   A single element extra row, x, is to be appended to X and Q and R updated 
   accordingly. x is zero except for kth element lam.

   Let Q* be the full orthogonal matrix of which Q is the upper left portion, then 
   
   [X] = [Q* 0][R]
   [x]   [0  1][0]
               [x]

   The rhs of the above can be bought into standard QR form by application of givens rotations 
   from the left to the augmented R matrix and the corresponding inverse rotations from the right 
   to the augmented Q* matrix. The rotations from the right applied to the augmented Q* have the 
   effect of rotating columns of Q into the final column of the augmented matrix and vice-versa.

   Since the columns between Q and the final column are not involved, only Q and R need to be updated 
   here, the rest of Q* being irrelevant. This routine does not augment the Q by an extra row, it is 
   assumed that the calling function only requires the update of the input rows.

   All matrices are assumed to be packed in column order. Some very minor further optimizations could 
   be added (e.g. using fact that working and most of x are zero at first iteration), but it's really 
   unlikely to yield a worthwhile saving.

   Some R code for testing the routine:

library(mgcv)
n<-4;q<-3
X<-matrix(rnorm(n*q),n,q)
#X[,q-1]<-X[,q]
qrx<-qr(X,tol=0)
Q<-qr.Q(qrx);R<-qr.R(qrx);lam<-1
um<-.C("update_qr",as.double(Q),as.double(R),as.integer(n),as.integer(q),
        as.double(lam),as.integer(q-1),PACKAGE="mgcv")
R1<-matrix(um[[2]],q,q);Q1<-matrix(um[[1]],n,q)
Xa<-matrix(0,n+1,q)
Xa[1:n,]<-X;Xa[n+1,q]<-lam
qr.R(qr(Xa,tol=0))   

*/

{ double *x,*work,c,s,r,x0,x1,m,*xip,*xjp,*riip,*rijp,*Qp,*wp;
  x=(double *)calloc((size_t)*q,sizeof(double)); 
  work=(double *)calloc((size_t)*n,sizeof(double)); /* working extra column of Q */
  x[*k] = *lam;
  /* conceptually i runs from k to q in the following loop */
  for (Qp=Q+ *k * *n,riip=R+ *k * *q + *k,xip=x+ *k ;xip< x+ *q;xip++,riip+= *q+1)
  { /* rotate x[i] into R[i,i], using over/underflow proof rotator */
    x0= *xip; /* x[i] */
    x1= *riip; /* R[1 * *q + i] */
    m = fabs(x0);s=fabs(x1); if (s>m) m=s;
    x1/=m;x0/=m;
    r=sqrt(x0*x0+x1*x1);
    c=x1/r;s=x0/r;
    *riip=m*r;/* *xip=0.0; but never neaded*/    
    /* conceptually j runs from i+1 to q in the following loop */
    for (rijp=riip + *q,xjp=xip+1;xjp<x+ *q;xjp++,rijp+= *q) /* apply rotation to x[j], R[i,j] */
    { x1= *rijp;    /* R[j* *q+i] */
      /* *xjp == x[j] */
      *rijp = c*x1 - s* *xjp;
      *xjp = s*x1 + c* *xjp;
    }
    /* apply transpose of rotation from right to Q */
    /* conceptually j runs from 0 to n in following loop */
    for (wp=work;wp<work+ *n;wp++,Qp++)
    { x1 = *Qp;   /*Q[i* *n+j]*/
      /* *wp == work[j] */
      *Qp = c*x1-s* *wp;
      *wp = s*x1+c* *wp;
    }
  }
  free(x);free(work);
}

void mgcv_symeig(double *A,double *ev,int *n,int *use_dsyevd,int *get_vectors,
                 int *descending)

/* gets eigenvalues and (optionally) vectors of symmetric matrix A. 
   
   Two alternative underlying LAPACK routines can be used, 
   either dsyevd (slower, robust) or dsyevr (faster, seems less robust). 
   Vectors returned in columns of A, values in ev (ascending).
   
   ******************************************************
   *** Eigenvalues are returned  in *ascending* order ***
   *** unless descending is set to be non-zero        ***
   ******************************************************

   Testing R code....
   library(mgcv)
   n<-4;A<-matrix(rnorm(n*n),n,n);A<-A%*%t(A);d<-array(0,n)
   er<-eigen(A)
   um<-.C("mgcv_symeig",as.double(A),as.double(d),as.integer(n),
           as.integer(1),as.integer(1),PACKAGE="mgcv")
   er$vectors;matrix(um[[1]],n,n)
   er$values;um[[2]]
*/  

{ char jobz='V',uplo='U',range='A'; 
  double work1,*work,dum1=0,abstol=0.0,*Z,*dum2,x,*p;
  int lwork = -1,liwork = -1,iwork1,info,*iwork,dumi=0,n_eval=0,*isupZ,i;
  if (*get_vectors) jobz='V'; else jobz='N';
  if (*use_dsyevd)
  { F77_NAME(dsyevd)(&jobz,&uplo,n,A,n,ev,&work1,&lwork,&iwork1,&liwork,&info);
    lwork=(int)floor(work1);if (work1-lwork>0.5) lwork++;
    work=(double *)calloc((size_t)lwork,sizeof(double));
    liwork = iwork1;iwork= (int *)calloc((size_t)liwork,sizeof(int));
    F77_NAME(dsyevd)(&jobz,&uplo,n,A,n,ev,work,&lwork,iwork,&liwork,&info);
    free(work);free(iwork);
  } else
  { Z=(double *)calloc((size_t)(*n * *n),sizeof(double)); /* eigen-vector storage */
    isupZ=(int *)calloc((size_t)(2 * *n),sizeof(int)); /* eigen-vector support */
    F77_NAME(dsyevr)(&jobz,&range,&uplo,
		   n,A,n,&dum1,&dum1,&dumi,&dumi,
		   &abstol,&n_eval,ev, 
    		     Z,n,isupZ, &work1,&lwork,&iwork1,&liwork,&info);
    lwork=(int)floor(work1);if (work1-lwork>0.5) lwork++;
    work=(double *)calloc((size_t)lwork,sizeof(double));
    liwork = iwork1;iwork= (int *)calloc((size_t)liwork,sizeof(int));
    F77_NAME(dsyevr)(&jobz,&range,&uplo,
		   n,A,n,&dum1,&dum1,&dumi,&dumi,
		   &abstol,&n_eval,ev, 
    		     Z,n,isupZ, work,&lwork,iwork,&liwork,&info);
    free(work);free(iwork);
    
    if (*descending) for (i=0;i<*n/2;i++) { /* reverse the eigenvalues */
      x = ev[i]; ev[i] = ev[*n-i-1];ev[*n-i-1] = x;
    }

    if (*get_vectors) {  /* copy vectors back into A */
      if (*descending) { /* need to reverse order */
        dum2 = Z + *n * (*n-1);
        for (work=dum2;work>=Z;work -= *n)
	  for (p=work;p<work + *n;p++,A++) *A = *p; 
      } else { /* simple copy */
        dum2 = Z + *n * *n;
        for (work=Z;work<dum2;work++,A++) *A = *work;
      }
    }
    free(Z);free(isupZ);
  }

 }

void mgcv_trisymeig(double *d,double *g,double *v,int *n,int getvec,int descending) 
/* Find eigen-values and vectors of n by n symmetric tridiagonal matrix 
   with leading diagonal d and sub/super diagonals g. 
   eigenvalues returned in d, and eigenvectors in columns of v, if
   getvec!=0. If *descending!=0 then eigenvalues returned in descending order,
   otherwise ascending. eigen-vector order corresponds.  

   Routine is divide and conquer followed by inverse iteration. 

   dstevd could be used instead, with just a name change.
   dstevx may be faster, but needs argument changes.
*/ 
		    
{ char compz;
  double *work,work1,x,*dum1,*dum2;
  int ldz=0,info,lwork=-1,liwork=-1,*iwork,iwork1,i,j;

  if (getvec) { compz='I';ldz = *n;} else { compz='N';ldz=0;}

  /* workspace query first .... */
  F77_NAME(dstedc)(&compz,n,
		   d, g, /* lead and su-diag */
		   v, /* eigenvectors on exit */  
                   &ldz, /* dimension of v */
		   &work1, &lwork,
		   &iwork1, &liwork, &info);

   lwork=(int)floor(work1);if (work1-lwork>0.5) lwork++;
   work=(double *)calloc((size_t)lwork,sizeof(double));
   liwork = iwork1;
   iwork= (int *)calloc((size_t)liwork,sizeof(int));

   /* and the actual call... */
   F77_NAME(dstedc)(&compz,n,
		   d, g, /* lead and su-diag */
		   v, /* eigenvectors on exit */  
                   &ldz, /* dimension of v */
		   work, &lwork,
		   iwork, &liwork, &info);

   if (descending) { /* need to reverse eigenvalues/vectors */
     for (i=0;i<*n/2;i++) { /* reverse the eigenvalues */
       x = d[i]; d[i] = d[*n-i-1];d[*n-i-1] = x;
       dum1 = v + *n * i;dum2 = v + *n * (*n-i-1); /* pointers to heads of cols to exchange */
       for (j=0;j<*n;j++,dum1++,dum2++) { /* work down columns */
         x = *dum1;*dum1 = *dum2;*dum2 = x;
       }
     }
   }

   free(work);free(iwork);
   *n=info; /* zero is success */
}



void Rlanczos(double *A,double *U,double *D,int *n, int *m, int *lm,double *tol) {
/* Prototype faster lanczos_spd for calling from R.
   A is n by n symmetric matrix. Let k = m + max(0,lm).
   U is n by k and D is a k-vector.
   m is the number of upper eigenvalues required and lm the number of lower.
   If lm<0 then the m largest magnitude eigenvalues (and their eigenvectors)
   are returned 

   Matrices are stored in R (and LAPACK) format (1 column after another).

   ISSUE: 1. eps_stop tolerance is set *very* tight.
          2. Currently all eigenvectors of Tj are found, although only the next unconverged one
             is really needed. Might be better to be more selective using dstein from LAPACK. 
          3. Basing whole thing on dstevx might be faster
          4. Is random start vector really best? convergence seems very slow. Might be better to 
             use e.g. a sine wave, and simply change its frequency if it seems to be in null space.
             Demmel (1997) suggests using a random vector, to avoid any chance of orthogonality with
             an eigenvector!
        
*/
  int biggest=0,f_check,i,k,kk,ok,l,j,vlength=0,neg_conv,pos_conv,ni,pi,neg_closed,pos_closed,converged,incx=1;
  double **q,*v=NULL,bt,xx,yy,*a,*b,*d,*g,*z,*err,*p0,*p1,*Ap,*zp,*qp,normTj,eps_stop=DOUBLE_EPS,max_err,alpha=1.0,beta=0.0;
  unsigned long jran=1,ia=106,ic=1283,im=6075; /* simple RNG constants */
  const char uplo='U';
  eps_stop = *tol; 

  if (*lm<0) { biggest=1;*lm=0;} /* get m largest magnitude eigen-values */
  f_check = (*m + *lm)/2; /* how often to get eigen_decomp */
  if (f_check<10) f_check =10;
  kk = (int) floor(*n/10); if (kk<1) kk=1;  
  if (kk<f_check) f_check = kk;

  q=(double **)calloc((size_t)(*n+1),sizeof(double *));

  /* "randomly" initialize first q vector */
  q[0]=(double *)calloc((size_t)*n,sizeof(double));
  b=q[0];bt=0.0;
  for (i=0;i < *n;i++)   /* somewhat randomized start vector!  */
  { jran=(jran*ia+ic) % im;   /* quick and dirty generator to avoid too regular a starting vector */
    xx=(double) jran / (double) im -0.5; 
    b[i]=xx;xx=-xx;bt+=b[i]*b[i];
  } 
  bt=sqrt(bt); 
  for (i=0;i < *n;i++) b[i]/=bt;

  /* initialise vectors a and b - a is leading diagonal of T, b is sub/super diagonal */
  a=(double *)calloc((size_t) *n,sizeof(double));
  b=(double *)calloc((size_t) *n,sizeof(double));
  g=(double *)calloc((size_t) *n,sizeof(double));
  d=(double *)calloc((size_t) *n,sizeof(double)); 
  z=(double *)calloc((size_t) *n,sizeof(double));
  err=(double *)calloc((size_t) *n,sizeof(double));
  for (i=0;i< *n;i++) err[i]=1e300;
  /* The main loop. Will break out on convergence. */
  for (j=0;j< *n;j++) 
  { /* form z=Aq[j]=A'q[j], the O(n^2) step ...  */
    /*for (Ap=A,zp=z,p0=zp+*n;zp<p0;zp++) 
      for (*zp=0.0,qp=q[j],p1=qp+*n;qp<p1;qp++,Ap++) *zp += *Ap * *qp;*/
    /*  BLAS versions y := alpha*A*x + beta*y, */
    F77_NAME(dsymv)(&uplo,n,&alpha,
		A,n,
		q[j],&incx,
		&beta,z,&incx);
    /* Now form a[j] = q[j]'z.... */
    for (xx=0.0,qp=q[j],p0=qp+*n,zp=z;qp<p0;qp++,zp++) xx += *qp * *zp;
    a[j] = xx;
 
    /* Update z..... */
    if (!j)
    { /* z <- z - a[j]*q[j] */
      for (zp=z,p0=zp+*n,qp=q[j];zp<p0;qp++,zp++) *zp -= xx * *qp;
    } else
    { /* z <- z - a[j]*q[j] - b[j-1]*q[j-1] */
      yy=b[j-1];      
      for (zp=z,p0=zp + *n,qp=q[j],p1=q[j-1];zp<p0;zp++,qp++,p1++) *zp -= xx * *qp + yy * *p1;    
  
      /* Now stabilize by full re-orthogonalization.... */
      


      for (i=0;i<=j;i++) 
      { /* form xx= z'q[i] */
        /*for (xx=0.0,qp=q[i],p0=qp + *n,zp=z;qp<p0;zp++,qp++) xx += *zp * *qp;*/
        xx = -F77_NAME(ddot)(n,z,&incx,q[i],&incx); /* BLAS version */
        /* z <- z - xx*q[i] */
        /*for (qp=q[i],zp=z;qp<p0;qp++,zp++) *zp -= xx * *qp;*/
        F77_NAME(daxpy)(n,&xx,q[i],&incx,z,&incx); /* BLAS version */
      } 
      
      /* exact repeat... */
      for (i=0;i<=j;i++) 
      { /* form xx= z'q[i] */
        /* for (xx=0.0,qp=q[i],p0=qp + *n,zp=z;qp<p0;zp++,qp++) xx += *zp * *qp; */
        xx = -F77_NAME(ddot)(n,z,&incx,q[i],&incx); /* BLAS version */
        /* z <- z - xx*q[i] */
        /* for (qp=q[i],zp=z;qp<p0;qp++,zp++) *zp -= xx * *qp; */
        F77_NAME(daxpy)(n,&xx,q[i],&incx,z,&incx); /* BLAS version */
      } 
      /* ... stabilized!! */
    } /* z update complete */
    

    /* calculate b[j]=||z||.... */
    for (xx=0.0,zp=z,p0=zp+*n;zp<p0;zp++) xx += *zp * *zp;b[j]=sqrt(xx); 
  
    /*if (b[j]==0.0&&j< *n-1) ErrorMessage(_("Lanczos failed"),1);*/ /* Actually this isn't really a failure => rank(A)<l+lm */
    /* get q[j+1]      */
    if (j < *n-1)
    { q[j+1]=(double *)calloc((size_t) *n,sizeof(double));
      for (xx=b[j],qp=q[j+1],p0=qp + *n,zp=z;qp<p0;qp++,zp++) *qp = *zp/xx; /* forming q[j+1]=z/b[j] */
    }

    /* Now get the spectral decomposition of T_j.  */

    if (((j>= *m + *lm)&&(j%f_check==0))||(j == *n-1))   /* no  point doing this too early or too often */
    { for (i=0;i<j+1;i++) d[i]=a[i]; /* copy leading diagonal of T_j */
      for (i=0;i<j;i++) g[i]=b[i]; /* copy sub/super diagonal of T_j */   
      /* set up storage for eigen vectors */
      if (vlength) free(v); /* free up first */
      vlength=j+1; 
      v = (double *)calloc((size_t)vlength*vlength,sizeof(double));
 
      /* obtain eigen values/vectors of T_j in O(j^2) flops */
    
      kk = j + 1;
      mgcv_trisymeig(d,g,v,&kk,1,1);
      /* ... eigenvectors stored one after another in v, d[i] are eigenvalues */

      /* Evaluate ||Tj|| .... */
      normTj=fabs(d[0]);if (fabs(d[j])>normTj) normTj=fabs(d[j]);

      for (k=0;k<j+1;k++) /* calculate error in each eigenvalue d[i] */
      { err[k]=b[j]*v[k * vlength + j]; /* bound on kth e.v. is b[j]* (jth element of kth eigenvector) */
        err[k]=fabs(err[k]);
      }
      /* and check for termination ..... */
      if (j >= *m + *lm)
      { max_err=normTj*eps_stop;
        if (biggest) { /* getting m largest magnitude eigen values */
          /* Finished only when the smallest element of the positive converged set and the 
             smallest element of the negative converged set are both not in the largest magnitude 
             set, or these sets are finished, and the largest magnitude set is of size *m, or when the total number 
             converged equals the matrix dimension */
	  pos_closed=neg_closed=0;
          neg_conv=0;i=j; while (i>=0&&err[i]<max_err&&d[i]<0.0) { neg_conv++;i--;}
          if (i>=0&&err[i]<max_err) neg_closed=1; /* all negatives found */
	  pos_conv=0;i=0; while (i<=j&&err[i]<max_err&&d[i]>=0.0) { i++;pos_conv++;}
          if (i<=j&&err[i]<max_err) pos_closed=1; /* all positives found */
 
          if (neg_conv+pos_conv >= *m) { /* some chance of having finished */
            pi=0;ni=0; /* counters for how many of neg and pos converged to include in largest set */
            while (pi+ni < *m) {
              if ((d[pi] > -d[j-ni]&&pi<pos_conv)||ni>=neg_conv) pi++; else ni++;
            }
            /* now pi and ni are the number of terms to include in the largest m 
               from the +ve and -ve converged sets */
            converged=1;
            if (!pos_closed&&pi==pos_conv) converged=0; /* don't know that there is not a larger value to come */            
            if (!neg_closed&&ni==neg_conv) converged=0; /* ditto */
          } else converged = 0;
 
          if (converged) {
            *m = pi;
            *lm = ni;
            j++;break;
          }
        } else
        { ok=1;
          for (i=0;i < *m;i++) if (err[i]>max_err) ok=0;
          for (i=j;i > j - *lm;i--) if (err[i]>max_err) ok=0;
          if (ok) 
          { j++;break;}
        }
      }
    }
  }
  /* At this stage, complete construction of the eigen vectors etc. */
  
  /* Do final polishing of Ritz vectors and load va and V..... */
  /*  for (k=0;k < *m;k++) // create any necessary new Ritz vectors 
  { va->V[k]=d[k];
    for (i=0;i<n;i++) 
    { V->M[i][k]=0.0; for (l=0;l<j;l++) V->M[i][k]+=q[l][i]*v[k][l];}
  }*/

  /* assumption that U is zero on entry! */

  for (k=0;k < *m;k++) /* create any necessary new Ritz vectors */
  { D[k]=d[k];
    for (l=0;l<j;l++)
    for (xx=v[l + k * vlength],p0=U + k * *n,p1 = p0 + *n,qp=q[l];p0<p1;p0++,qp++) *p0 += *qp * xx;
  }

  for (k= *m;k < *lm + *m;k++) /* create any necessary new Ritz vectors */
  { kk=j-(*lm + *m - k); /* index for d and v */
    D[k]=d[kk];
    for (l=0;l<j;l++)
    for (xx=v[l + kk * vlength],p0=U + k * *n,p1 = p0 + *n,qp=q[l];p0<p1;p0++,qp++) *p0 += *qp * xx;
  }
 
  /* clean up..... */
  free(a);
  free(b);
  free(g);
  free(d);
  free(z);
  free(err);
  if (vlength) free(v);
  for (i=0;i< *n+1;i++) if (q[i]) free(q[i]);free(q);  
  *n = j; /* number of iterations taken */
} /* end of Rlanczos */

