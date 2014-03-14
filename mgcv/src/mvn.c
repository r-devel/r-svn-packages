/* (c) Simon N Wood. 2014. Released under GPL2. 
  
  likelihood and derivative evaluation for multivariate Gaussian 
  additive models.

  R CMD SHLIB mvn.c to compile.

*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>
#include <R_ext/BLAS.h>
#include "mgcv.h"

void mvn_ll(double *y,double *X,double *beta,int *n,int *lpi, /* note zero indexing */
            int *m,double *ll,double *lb,int *nt) {
/* inputs:
    * 'y' is an m by n matrix, each column of which is a m-dimensional observation of a 
      multivariate normal r.v.
    * 'X' is a sequence of model matrices. The first (0th) model matrix runs from columns 0 to lpi[0]-1,
      the jth from cols lpi[j-1] to lpi[j]-1. lpi indexing starts from 0!!
    * 'beta' is a parameter vector corresponding to X. The m*(m+1)/2 elements starting at lpi[m] are the 
      parameters of the Choleki factor of the precision matrix.
    * nt is number of threads to use.     
  
    outputs:
    * 'll' is the evaluated log likelihood.
    * 'lb' is the grad vector 
*/
  double *R,*theta,ldetR,*Xl,*bl,oned=1.0,zerod=0.0,*p,*p1,*p2,*p3,xx,zz,
    *v1,*v2,*mu,*Rymu;
  int i,j,k,l,pl,one=1,bt,ct;
  const char not_trans='N';
  /* Create the Choleski factor of the precision matrix */
  R = (double *)R_chk_calloc((size_t)*m * *m,sizeof(double));
  theta = beta + lpi[*m]; /* parameters of R */
  ldetR = 0.0; /* log|R| */
  for (k=0,i=0;i<*m;i++) { 
    R[i + *m + i] = exp(theta[k]);ldetR += theta[k];k++; 
    for (j=i+1;j<*m;j++) { R[i + *m + j] = theta[k];k++;}
  }  
  /* obtain y - mu */
  mu  = (double *)R_chk_calloc((size_t)*n,sizeof(double));
  for (l=0;l<*m;l++) { /* loop through components */
    if (l==0) { Xl = X;pl = lpi[0];bl=beta;} /* Xl is lth model matrix with pl columns, coef vec bl */ 
    else { Xl = X + *n * lpi[l-1];pl = lpi[l]-lpi[l-1];bl = beta + lpi[l-1];}   
    F77_CALL(dgemv)(&not_trans,n,&pl,&oned,Xl,n, bl, &one,&zerod, mu, &one); /* BLAS call for mu = Xl bl */
    /* now subtract mu from relevant component of y */
    for (p=mu,p1= mu + *n,p2=y+k;p<p1;p++,p2 += *m) *p2 -= *p;
  }
  R_chk_free(mu);
  /* so y now contains y-mu */
  
  /* R(y-mu) is required repeatedly... */
  Rymu =  (double *)R_chk_calloc((size_t)*n * *m,sizeof(double));
  bt=0;ct=0;mgcv_pmmult(Rymu,R,y,&bt,&ct,m,n,m,nt);  
  /* compute the log likelihood */
  for (*ll=0.0,p=Rymu,p1=p + *n * *m;p<p1;p++) *ll += *p * *p;
  *ll = - *ll/2 + ldetR * *n;  
  
  /* now the grad vector */
  
  p = lb;
  /* first the derivatives w.r.t. the coeffs of the linear predictors */
  for (l=0;l<*m;l++) { /* work through dimensions */
    if (l==0) { Xl = X;pl = lpi[0];bl=beta;} /* Xl is lth model matrix with pl columns, coef vec bl */ 
    else { Xl = X + *n * lpi[l-1];pl = lpi[l]-lpi[l-1];bl = beta + lpi[l-1];} 
    for (i=0;i<pl;i++) { /* work through coefs for this dimension */
      *p = 0.0; /* clear lb */
      for (j=0;j<*n;j++) {
        xx = Xl[j + *n * i];
        for (p1=R + l * *m,p2 = p1 + l,p3 = Rymu + *m *j;p1<p2;p1++,p3++) *p += xx * *p1 * *p3; 
      }
      p++;
    }
  }
  /* now the derivatives w.r.t. the parameters, theta, of R */ 
  v1 =  (double *)R_chk_calloc((size_t)*m,sizeof(double));
  v2 =  (double *)R_chk_calloc((size_t)*m,sizeof(double));
  for (k=0,i=0;i<*m;i++) { /* i is row */ 
    /* get tr(R^{-1}R R_theta^k) */
    v1[i] = xx = exp(theta[k]); /* the non-zero element of R_theta^i at i,i */;
    mgcv_backsolve(R,m,m,v1,v2,&one);
    v1[i] = 0.0;
    *p += *n * v2[i]; /* v2[i] contains tr(R^{-1}R R_theta^k) */
    k++; /* increment the theta index */ 
    /* quadratic form involves only ith dimension */
    for (zz=0.0,l=0,p1 = Rymu+i,p2=y+i;l<*n;l++,p1 += *m,p2 += *m) zz += *p1 * *p2 * xx;
    *p += -zz/2;p++;
    for (j=i+1;j<*m;j++) { /* j is col */
      v1[i] = 1.0;/* the non-zero element of R_theta^k at i,j */
      mgcv_backsolve(R,m,m,v1,v2,&one);
      v1[i] = 0.0;
      *p += *n * v2[j]; /* v2[j] contains tr(R^{-1}R R_theta^k) */
      k++; /* increment the theta index */ 
      for (zz=0.0,l=0,p1 = Rymu+i,p2=y+i;l<*n;l++,p1 += *m,p2 += *m) zz += *p1 * *p2;
      *p += -zz/2;p++;
    }
  }  
  R_chk_free(v1);R_chk_free(v2);
  R_chk_free(R);R_chk_free(Rymu);
} /* mvn_ll */


