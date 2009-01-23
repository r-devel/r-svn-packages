/* Copyright (C) 2007,2008,2009 Simon N. Wood  simon.wood@r-project.org

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
(www.gnu.org/copyleft/gpl.html)

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
USA. */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define ANSI
/*#define DEBUG*/
#include "matrix.h"
#include "mgcv.h"

void getXtWX(double *XtWX, double *X,double *w,int *r,int *c,double *work)
/* forms X'WX as efficiently as possible, where W = diag(w)
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

void getXtMX(double *XtMX,double *X,double *M,int *r,int *c,double *work)
/* forms X'MX as efficiently as possible, where M is a symmetric matrix
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

double trBtAB(double *A,double *B,int *n,int*m) 
/* form tr(B'AB) where A is n by n and B is n by m, m < n,
   basic point is that this is sum_ijk A_ik B_ij B_kj
 */
{ double tr=0.0,x,*p,*p1,*p2;
  int j,k;
  for (j=0;j<*m;j++)
    for (k=0;k<*n;k++) {
      p = A + *n * k;p2 = p + *n;
      p1 = B + *n * j;
      x = B[k + j * *n];
      for (;p<p2;p++,p1++) tr += *p * *p1 * x;   
    }
  return(tr);
}


void multSk(double *y,double *x,int *xcol,int k,double *rS,int *rSncol,int *q,double *work)
/* function to form y = Sk x, where a square root of Sk 
   is packed somewhere inside rS. x must be q by xcol. The 
   kth square root is q by rSncol[k]. The square roots are packed 
   one after another columnwise (R default).
   
   work and y must be the same dimension as x.
*/
{ int i,off,nc,bt,ct;
  double *rSk;
  off=0; /* the start of the kth square root */
  for (i=0;i<k;i++) off += *q * rSncol[i];
  rSk = rS + off; /* pointer to the kth square root */
  nc = rSncol[k];
  bt=1;ct=0;
  mgcv_mmult(work,rSk,x,&bt,&ct,&nc,xcol,q);
  bt=0;
  mgcv_mmult(y,rSk,work,&bt,&ct,q,xcol,&nc);
}

double diagABt(double *d,double *A,double *B,int *r,int *c)
/* obtain diag(AB') as efficiently as possible, and return tr(AB') A and B are
   r by c stored column-wise.
*/
{ int j;
  double tr,*pa,*pb,*p1,*pd;
  for (pa=A,pb=B,p1=pa + *r,pd=d;pa<p1;pa++,pb++,pd++) *pd = *pa * *pb;
  for (j=1;j < *c;j++)
  for (p1=d + *r,pd=d;pd<p1;pa++,pb++,pd++) *pd += *pa * *pb;
  /* d now contains diag(AB') */
  for (tr=0.0,pd=d,p1=d + *r;pd < p1;pd++) tr += *pd;
  return(tr);
}


double trAB(double *A,double *B,int *n, int *m)
/* Get tr(AB) where A is n by m and B is m by n 
*/ 
{ double *p,*pa,*pb,tr=0.0;
  int i,j;
  for (pa=A,pb=B,j=0;j<*m;j++,pb++)
    for (p=pb,i=0;i<*n;i++,p += *m,pa++) tr+= *p * *pa;
  return(tr);
}


void get_bSb(double *bSb,double *bSb1, double *bSb2,double *sp,double *E,
             double *rS,int *rSncol,int *Encol, int *q,int *M,double *beta, 
             double *b1, double *b2,int *deriv)
/* Routine to obtain beta'Sbeta and its derivatives w.r.t. the log smoothing 
   parameters, this is part of REML calculation... 
   
   S= EE'.

   b1 and b2 contain first and second derivatives of q-vector beta w.r.t. 
   \pho_k. They are packed as follows....

   * b1 will contain dbeta/d\rho_0, dbeta/d\rho_1 etc. So, for example, dbeta_i/d\rho_j
     (indices starting at zero) is located in b1[q*j+i].
   
   * b2 will contain d^2beta/d\rho_0d\rho_0, d^2beta/d\rho_1d\rho_0,... but rows will not be
     stored if they duplicate an existing row (e.g. d^2beta/d\rho_0d\rho_1 would not be 
     stored as it already exists and can be accessed by interchanging the sp indices).
     So to get d^2beta_k/d\rho_id\rho_j: 
     i)   if i<j interchange the indices
     ii)  off = (j*m-(j+1)*j/2+i)*q (m is number of sp's) 
     iii) v2[off+k] is the required derivative.       

*/
{ double *Sb,*Skb,*work,*work1,*p1,*p0,*p2,xx;
  int i,j,bt,ct,one=1,m,k,rSoff,mk,km; 
  
  work = (double *)calloc((size_t)*q,sizeof(double)); 
  Sb = (double *)calloc((size_t)*q,sizeof(double));
  bt=1;ct=0;mgcv_mmult(work,E,beta,&bt,&ct,Encol,&one,q);
  bt=0;ct=0;mgcv_mmult(Sb,E,work,&bt,&ct,q,&one,Encol); /* S \hat \beta */

  for (*bSb=0.0,i=0;i<*q;i++) *bSb += beta[i] * Sb[i]; /* \hat \beta' S \hat \beta */

  if (*deriv <=0) {free(work);free(Sb);return;}

  work1 = (double *)calloc((size_t)*q,sizeof(double));
  Skb = (double *)calloc((size_t)*M * *q,sizeof(double));
 
  for (p1=Skb,rSoff=0,i=0;i<*M;i++) { /* first part of first derivatives */
     /* form S_k \beta * sp[i]... */
     bt=1;ct=0;mgcv_mmult(work,rS + rSoff ,beta,&bt,&ct,rSncol+i,&one,q);
     for (j=0;j<rSncol[i];j++) work[j] *= sp[i]; 
     bt=0;ct=0;mgcv_mmult(p1,rS + rSoff ,work,&bt,&ct,q,&one,rSncol+i);
     rSoff += *q * rSncol[i];

     /* now the first part of the first derivative */
     for (xx=0.0,j=0;j<*q;j++,p1++) xx += beta[j] * *p1;
     bSb1[i] = xx; 
  }


  if (*deriv>1)  for (m=0;m < *M;m++) { /* Hessian */
     bt=1;ct=0;mgcv_mmult(work1,E,b1+m * *q,&bt,&ct,Encol,&one,q);
     bt=0;ct=0;mgcv_mmult(work,E,work1,&bt,&ct,q,&one,Encol);  /* S dbeta/drho_m */

    for (k=m;k < *M;k++) {
      km=k * *M + m;mk=m * *M + k;  /* second derivatives needed */
      /* d2beta'/drho_k drho_m S beta */
      for (xx=0.0,p0=Sb,p1=Sb + *q;p0<p1;p0++,b2++) xx += *b2 * *p0;
      bSb2[km] = 2*xx; 
       
      /* dbeta'/drho_k S d2beta/drho_m */
      for (xx=0.0,p0=b1+k* *q,p1=p0 + *q,p2=work;p0<p1;p0++,p2++) xx += *p2 * *p0;
      bSb2[km] += 2*xx;

      /* dbeta'/drho_k S_m beta sp[m] */
      for (xx=0.0,p0=Skb + k* *q,p1=p0 + *q,p2= b1+m* *q;p0<p1;p0++,p2++) xx += *p2 * *p0;
      bSb2[km] += 2*xx;
 
      /* dbeta'/drho_m S_k beta sp[k] */
      for (xx=0.0,p0=Skb + m* *q,p1=p0 + *q,p2= b1+k* *q;p0<p1;p0++,p2++) xx += *p2 * *p0;
      bSb2[km] += 2*xx;

      if (k==m) bSb2[km] += bSb1[k]; else bSb2[mk] = bSb2[km];
    }
  } /* done Hessian */

  /* Now finish off the first derivatives */
  bt=1;ct=0;mgcv_mmult(work,b1,Sb,&bt,&ct,M,&one,q);
  for (i=0;i<*M;i++) bSb1[i] += 2*work[i];
  
  free(Sb);free(work);free(Skb);free(work1);

}

double frobenius_norm(double *X,int *r, int *c)
/* The Frobenius norm of r by c matrix X. Interestingly, this gives an 
   upper bound on the two norm (largest singular value). 
*/
{ double fnorm=0.0,*p1;
  int n;
  n = *r * *c;
  for (p1=X+n;X<p1;X++) fnorm += *X * *X;
  return(sqrt(fnorm));
}

double qr_ldet_inv(double *X,int *r,double *Xi,int *get_inv) 
/* Obtains the log|X| and the inverse of X (r by r), by pivoted QR decomposition. 
   The inverse is returned (unpivoted) in Xi. 
   The function returns log|X| as its value.
*/
{ double *tau,ldet,*p,*Qt;
  int *pivot,i,TRUE=1,j;
  /* Allocated working storage ...*/
  pivot = (int *)calloc((size_t)*r,sizeof(int));
  tau = (double *)calloc((size_t)*r,sizeof(double));
  
  mgcv_qr(X,r,r,pivot,tau); /* get QR=X itself */

  /* evaluate log|X| = sum_i log(|R_ii|) ...*/
  for (ldet=0.0,p=X,i=0;i<*r;i++,p += *r+1) ldet += log(fabs(*p));
  
  if (*get_inv) {
  /* Now get the inverse of X. X^{-1} = R^{-1}Q' */
    Qt = (double *)calloc((size_t)*r * *r,sizeof(double));
    for (p=Qt,i=0;i<*r;i++,p += *r+1) *p = 1.0;
    mgcv_qrqy(Qt,X,tau,r,r,r,&TRUE,&TRUE); /* Extracting the orthogonal factor Q' */

    mgcv_backsolve(X,r,r,Qt,Xi,r); /* Now Xi contains the row pivoted inverse of X */

    /* Finally unpivot Xi. 
       pivot[i] is the unpivoted row that pivoted row i should end up in  
    */

    for (p=Xi,j=0;j<*r;j++,p += *r) { /* work column by colum using tau as working storage */

      for (i=0;i<*r;i++) tau[pivot[i]] = p[i]; /* ith row of pivoted -> pivot[i] row of unpivoted */
      for (i=0;i<*r;i++) p[i] = tau[i];        /* store unpivoted column in Xi */

    }
    free(Qt);
  } /* end if (*get_inv) */
  free(pivot);free(tau);
  return(ldet);
}

void get_detS2(double *sp,double *sqrtS, int *rSncol, int *q,int *M, int * deriv, 
               double *det, double *det1, double *det2, double *d_tol,
               double *r_tol,int *fixed_penalty)
/* Routine to evaluate log|S| and its derivatives wrt log(sp), in a stable manner, using 
   an orthogonal transformation strategy based on QR decomposition.

   Inputs are:
   `sp' the array of smoothing parameters.
   `sqrtS' the `M' square root penalty matrices. The ith is `q' by `rSncol[i]'. They are 
        packed one after the other. 
   `deriv' is the order of derivatives required. 0,1 or 2.
   `d_tol' is the tolerance to use for grouping dominant terms. 
   `r_tol' (<< d_tol) is the tolerance used for rank determination.
   `fixed_penalty' non-zero indicates that there is a fixed component of 
          total penalty matrix S, the square root of which is in the final 
          q * rSncol[M+1] elements of sqrtS.                 

   Outputs are:
   `det' the log determinant.
   `det1' M-array of derivatives of log det wrt log sp. 
   `det2' M by M Hessian of log det wrt log sp.   
   
*/
{ double *R,*work,*tau,*rS1,*rS2, *S,*Si,*Sb,*B,*Sg,*p,*p1,*p2,*p3,*p4,*frob,max_frob,x,*spf,Rcond;
  int *pivot,iter,i,j,k,bt,ct,rSoff,K,Q,Qr,*gamma,*gamma1,*alpha,r,max_col,Mf,tot_col=0,left,tp;

  if (*fixed_penalty) { 
    Mf = *M + 1;  /* total number of components, including fixed one */
    spf = (double *)calloc((size_t)Mf,sizeof(double));
    for (i=0;i<*M;i++) spf[i]=sp[i];spf[*M]=1.0; /* includes sp for fixed term */
  } 
  else {spf=sp;Mf = *M;} /* total number of components, including fixed one */

  /* Create working copies of sqrtS, which can be modified:
     rS1 is repeatedly orthogonally transformed, while rS2 is row pivoted. 
  */
  if (*deriv) { /* only need to modify if derivatives needed */
    for (j=i=0;i<Mf;i++) j += rSncol[i];tot_col=j;
    j *= *q;
    rS1 = (double *)calloc((size_t) j,sizeof(double));
    rS2 = (double *)calloc((size_t) j,sizeof(double));
    for (p=rS1,p3=rS2,p1=rS1+j,p2=sqrtS;p<p1;p++,p2++,p3++) *p3 = *p = *p2;
  }
  /* Explicitly form the Si (stored in a single block), so S_i is stored
     in Si + i * q * q (starting i from 0). As iteration progresses,
     blocks are shrunk -- always q by Q */
  p = Si = (double *)calloc((size_t)*q * *q * Mf,sizeof(double));
  max_col = *q; /* need enough storage just in case square roots are over-sized */
  for (rSoff=i=0;i<Mf;p+= *q * *q,rSoff+=rSncol[i],i++) {
    bt=0;ct=1;mgcv_mmult(p,sqrtS+rSoff * *q,sqrtS+rSoff * *q,&bt,&ct,q,q,rSncol+i);
    if (rSncol[i]>max_col) max_col=rSncol[i];
  }

 
  /* Initialize the sub-dominant set gamma and the counters */
  K = 0;Q = *q;
  frob =  (double *)calloc((size_t)Mf,sizeof(double)); 
  gamma = (int *)calloc((size_t)Mf,sizeof(int));  /* terms remaining to deal with */
  gamma1 = (int *)calloc((size_t)Mf,sizeof(int)); /* new gamma */
  alpha = (int *)calloc((size_t)Mf,sizeof(int));  /* dominant terms */
  for (i=0;i<Mf;i++) gamma[i] = 1; /* no terms dealt with initially */
  
  /* Other storage... */

  S = (double *) calloc((size_t) Q * Q,sizeof(double)); /* Transformed S (total) */
  Sb = (double *) calloc((size_t) Q * Q,sizeof(double)); /* summation storage */
  pivot = (int *)calloc((size_t) Q,sizeof(int)); /* pivot storage */
  tau = (double *) calloc((size_t) Q,sizeof(double)); /* working storage */  
  work = (double *)calloc((size_t)(4 * Q),sizeof(double));

  Sg = (double *) calloc((size_t) Q * Q,sizeof(double)); /* summation storage */
  B = (double *) calloc((size_t) Q * max_col,sizeof(double)); /* Intermediate storage */
  R = (double *) calloc((size_t) Q * Q,sizeof(double)); /* storage for unpivoted QR factor */
  /* Start the main orthogonal transform loop */
  iter =0;
  while(1) {
    iter ++;

  /* Find the Frobenius norms of the Si in set gamma */
    max_frob=0.0;
    for (p=Si,i=0;i<Mf;i++,p += *q * Q) 
      if (gamma[i]) { /* don't bother if already dealt with */ 
        frob[i] = frobenius_norm(p,q,&Q);
        if (frob[i] *spf[i] >max_frob) max_frob=frob[i]  * spf[i];
    }
  /* Find sets alpha and gamma' */
    for (i=0;i<Mf;i++) {
      if (gamma[i]) { /* term is still to be dealt with */
        if (frob[i]  * spf[i] > max_frob * *d_tol) { 
          alpha[i] = 1;gamma1[i] = 0; /* deal with it now */
        } else {
          alpha[i] = 0;gamma1[i] = 1; /* put it off */ 
        }
      } else { /* wasn't in gamma, so not in alpha or gamma1 */
        alpha[i] = gamma1[i] = 0;
      }
    }

  /* Form the scaled sum of the Si in alpha and get its rank by pivoted QR
     and condition estimation...
  */
    for (p=Sb,p1=p + *q * Q;p<p1;p++) *p=0.0; /* clear Sb */
    for (p=Si,i=0;i<Mf;i++,p += *q * Q) if (alpha[i]) { 
      x = frob[i];
      for (p1=p,p2=Sb,p3=p + *q * Q;p1<p3;p1++,p2++) *p2 += *p1 / x;
    } 
    for (i=0;i<*q;i++) pivot[i]=0; 
    mgcv_qr(Sb, &Q, q ,pivot,tau); /* obtain pivoted QR decomposition of Sb */
    /* Now obtain the rank, r, of Sb (see Golub and van Loan, 1996, p.129 & p.260)... */ 
    r = Q;
    R_cond(Sb,&Q,&r,work,&Rcond);
    while (*r_tol * Rcond > 1) { r--;R_cond(Sb,&Q,&r,work,&Rcond);}
    Qr = Q-r;  

    /* ...  r is the rank of Sb, or any other positively weighted sum over alpha */

    /*  printf("\n iter = %d,  rank = %d,   Q = %d",iter,r,Q);
    printf("\n gamma = ");for (i=0;i<Mf;i++) printf(" %d",gamma[i]);
    printf("\n alpha = ");for (i=0;i<Mf;i++) printf(" %d",alpha[i]);
    printf("\n gamma1 = ");for (i=0;i<Mf;i++) printf(" %d",gamma1[i]);*/
   

  /* If Q==r then terminate (form S first if it's the first iteration) */
    
    if (Q==r) { 
      if (iter==1 ) { /* form S */
        for (p=Si,i=0;i<Mf;i++,p += Q*Q) { 
          x = spf[i];
          for (p1=p,p2=S,p3=p+Q*Q;p1<p3;p1++,p2++) *p2 += *p1 * x;
        }
        break; 
      } else break; /* just use current S */ 
    } /* end if (Q==r) */

  /* Form the dominant term and QR-decompose it */
    for (p=Sb,p1=p + *q * Q;p<p1;p++) *p = 0.0; /* clear Sb */
    for (p=Si,i=0;i<Mf;i++,p += *q * Q) if (alpha[i]) { /* summing S[[i]]*sp[i] over i in alpha */
      x = spf[i];
      for (p1=p,p2=Sb,p3=p+ *q * Q;p1<p3;p1++,p2++) *p2 += *p1 * x;
    }
    for (i=0;i<*q;i++) pivot[i]=0; 
    mgcv_qr(Sb, &Q, q ,pivot,tau); /* obtain pivoted QR decomposition of Sb */
    
  /* unpivot R, which means that no further pivoting is needed */
    for (p=R,p1=R + *q * r;p<p1;p++) *p=0.0; /* clear R */
    for (i=0;i<r;i++) for (j=i;j<*q;j++) R[i + pivot[j] * r] = Sb[i + j * Q]; 

   /* DEBUG ONLY... */
    /*  printf("\npivot = ");for (j=0;j<*q;j++) printf("%d ",pivot[j]);
    printf("Current R...\n");
    for (i=0;i<r;i++) { for (j=0;j<*q;j++) printf("%7.2g  ",Sb[i + Q *j]); printf("\n");} */

  /* Form the sum over the elements in gamma1, Sg */

    for (p=Sg,p1=p + *q * Q;p<p1;p++) *p=0.0; /* clear Sg */
    for (p=Si,i=0;i<Mf;i++,p += *q * Q) if (gamma1[i]) { /* summing S[[i]]*sp[i] over i in gamma1 */
      x = spf[i];
      for (p1=p,p2=Sg,p3=p+ *q * Q;p1<p3;p1++,p2++) *p2 += *p1 * x;
    } 

  /* Form S' the orthogonal transform of S */
  
  /* Form Q'Sg... */ 
    left=1;tp=1; 
    mgcv_qrqy(Sg,Sb,tau,&Q,q,&Q,&left,&tp);

  /* copy transformed Sg into remainder of transformed S */
    for (i=0;i<Q;i++) for (j=0;j<*q;j++) S[i+K + j * *q] = Sg[i + j * Q];
  /* and add R in the appropriate place ... */ 
    for (i=0;i<r;i++) for (j=0;j<*q;j++) S[i+K + j * *q] += R[i + j * r];       


  /* transform remaining S_i in gamma1 */ 
    for (p1=p=Si,i=0;i<Mf;i++,p += *q * Q,p1 += *q *Qr) if (gamma1[i]) {
      left=1;tp=1;
      mgcv_qrqy(p,Sb,tau,&Q,q,&Q,&left,&tp); /* Q'Si */
      p2=p+r;p3=p1;
      for (j=0;j<*q;j++,p2+=r) for (k=0;k<Qr;k++,p3++,p2++) *p3 = *p2; /* copy to correct place */ 
    } 
 
  /* Transform the square roots of Si */
    if (*deriv) { /* transformed rS1 only needed for derivatives */
      /* copy last Q rows of rS1 into rS2 */
      for (i=0;i<Q;i++) for (j=0;j<tot_col;j++) rS2[i+Q*j] = rS1[K + i + *q * j];
      /* pre-multiply rS2 by Q */ 
      left=1;tp=1;
      mgcv_qrqy(rS2,Sb,tau,&Q,&tot_col,&Q,&left,&tp); /* Q'rS2 */
      
      /* copy rS2 into last Q rows of rS1 */
      for (i=0;i<Q;i++) for (j=0;j<tot_col;j++) rS1[K + i + *q * j] = rS2[i+Q*j];
    
      /* zero the last Qr rows of the rS1 in alpha */
      for (p=rS1,k=0;k<Mf;p +=rSncol[k] * *q, k++) if (alpha[k]) {
        for (i=K+r;i<*q;i++) for (j=0;j<rSncol[k];j++) p[i + j * *q] = 0.0;
      }
    }


    /* DEBUG ONLY... */
    /* printf("Current S...\n");
       for (i=0;i<*q;i++) { for (j=0;j<*q;j++) printf("%7.2g  ",S[i + *q *j]); printf("\n");}*/
 
  /* Update K, Q and gamma */   
    K = K + r; Q = Qr;
    for (i=0;i<Mf;i++) gamma[i] = gamma1[i];
  } /* end of Orthogonal Transform Loop */

  /* transpose S */
  for (i=0;i<*q;i++) for (j=0;j<*q;j++) R[i + *q * j] = S[j + *q * i]; 

  /* DEBUG ONLY... */
  /* printf("Final S...\n");
     for (i=0;i<*q;i++) { for (j=0;j<*q;j++) printf("%7.2g  ",S[i + *q *j]); printf("\n");}*/

  /* Now get the determinant and inverse of the transformed S (stored in B) */
  *det = qr_ldet_inv(R,q,B,deriv); /* R=S' here */
  /* finally, the derivatives, based on transformed S inverse and transformed square roots */  
  
  if (*deriv) { /* get the first derivatives */
    /* first accumulate S^{-T} sqrtS into Si */
    bt=0;ct=0;mgcv_mmult(Si,B,sqrtS,&bt,&ct,q,&tot_col,q);
    /* Now get the required derivatives */
    for (p=Si,p1=rS1,i=0;i<*M;p += *q *rSncol[i],p1+= *q *rSncol[i],i++) {
      for (x=0.0,p2=p1,p3=p,p4=p1+ *q*rSncol[i];p2<p4;p2++,p3++) x += *p2 * *p3; 
      det1[i] = x*sp[i]; /* tr(S^{-1}S_i) */
    }
  }
  
  if (*deriv==2) { /* get second derivatives, as well */
    for (p=Si,p1=rS2,p2 = p1 + *q * tot_col;p1<p2;p1++,p++) *p1 = *p; /* copy S^{-1} sqrtS into rS2 */   

    /* loop through creating S^{-1} S_i and storing in Si...*/
    for (p1=Si,p=rS2,p2=rS1,i=0;i<*M;p2+= *q * rSncol[i], p += *q *rSncol[i],i++,p1 += *q * *q) {
      bt=0;ct=1;mgcv_mmult(p1,p,p2,&bt,&ct,q,q,rSncol+i);
    }

    /* DEBUG ONLY...
       for (i=0;i<*M;i++) { for (x=0.0,j=0;j<*q;j++) x += Si[i* *q * *q + j + j* *q];det1[i]=x*sp[i];}*/

    for (i=0;i<*M;i++) for (j=i;j<*M;j++) 
      det2[i + *M * j] = det2[j + *M * i] = -sp[i]*sp[j]*trAB(Si + *q * *q *i,Si + *q * *q *j,q,q);
    for (i=0;i<*M;i++) det2[i + *M * i] += det1[i];
  }
  free(R);
  free(work);
  free(frob);
  free(gamma);
  free(gamma1);
  free(alpha);
  free(S);
  free(Sb);
  free(Sg);
  if (*deriv) { free(rS1);free(rS2);}
  if (*fixed_penalty) {free(spf);}
  free(Si);
  free(B);
  free(pivot);free(tau);
} /* end of get_detS3 */




void get_detS2a(double *sp,double *sqrtS, int *rSncol, int *q,int *M, int * deriv, 
               double *det, double *det1, double *det2, double *d_tol,
               double *r_tol,int *fixed_penalty)
/* Routine to evaluate log|S| and its derivatives wrt log(sp), in a stable manner, using 
   a similarity transform strategy.

   Inputs are:
   `sp' the array of smoothing parameters.
   `sqrtS' the `M' square root penalty matrices. The ith is `q' by `rSncol[i]'. They are 
        packed one after the other. 
   `deriv' is the order of derivatives required. 0,1 or 2.
   `d_tol' is the tolerance to use for grouping dominant terms. 
   `r_tol' (<< d_tol) is the tolerance used for rank determination.
   `fixed_penalty' non-zero indicates that there is a fixed component of 
          total penalty matrix S, the square root of which is in the final 
          q * rSncol[M+1] elements of sqrtS.                 

   Outputs are:
   `det' the log determinant.
   `det1' M-array of derivatives of log det wrt log sp. 
   `det2' M by M Hessian of log det wrt log sp.   
   
*/
{ double *rS, *Un, *U, *S,*Si,*Sb,*B,*C,*Sg,*p,*p1,*p2,*p3,*frob,*ev,max_frob,x,*spf;
  int iter,i,j,k,bt,ct,rSoff,K,Q,Qr,*gamma,*gamma1,*alpha,TRUE=1,FALSE=0,r,max_col,Mf;

  if (*fixed_penalty) { 
    Mf = *M + 1;  /* total number of components, including fixed one */
    spf = (double *)calloc((size_t)Mf,sizeof(double));
    for (i=0;i<*M;i++) spf[i]=sp[i];spf[*M]=1.0; /* includes sp for fixed term */
  } 
  else {spf=sp;Mf = *M;} /* total number of components, including fixed one */

  /* Create a working copy of sqrtS, which can be modified  */
  if (*deriv) { /* only need to modify if derivatives needed */
    for (j=i=0;i<Mf;i++) j += rSncol[i];j *= *q;
    rS = (double *)calloc((size_t) j,sizeof(double));
    for (p=rS,p1=rS+j,p2=sqrtS;p<p1;p++,p2++) *p = *p2;
  }
  /* Explicitly form the Si (stored in a single block), so S_i is stored
     in Si + i * q * q (starting i from 0). As iteration progresses,
     blocks are shrunk -- always Q by Q */
  p = Si = (double *)calloc((size_t)*q * *q * Mf,sizeof(double));
  max_col = *q; /* need enough storage just in case square roots are over-sized */
  for (rSoff=i=0;i<Mf;p+= *q * *q,rSoff+=rSncol[i],i++) {
    bt=0;ct=1;mgcv_mmult(p,sqrtS+rSoff * *q,sqrtS+rSoff * *q,&bt,&ct,q,q,rSncol+i);
    if (rSncol[i]>max_col) max_col=rSncol[i];
  }

 
  /* Initialize the sub-dominant set gamma and the counters */
  K = 0;Q = *q;
  frob =  (double *)calloc((size_t)Mf,sizeof(double)); 
  gamma = (int *)calloc((size_t)Mf,sizeof(int));  /* terms remaining to deal with */
  gamma1 = (int *)calloc((size_t)Mf,sizeof(int)); /* new gamma */
  alpha = (int *)calloc((size_t)Mf,sizeof(int));  /* dominant terms */
  for (i=0;i<Mf;i++) gamma[i] = 1; /* no terms dealt with initially */
  
  /* Other storage... */

  S = (double *) calloc((size_t) Q * Q,sizeof(double)); /* Transformed S (total) */
  U=Sb = (double *) calloc((size_t) Q * Q,sizeof(double)); /* summation storage */
  Sg = (double *) calloc((size_t) Q * Q,sizeof(double)); /* summation storage */
  ev = (double *) calloc((size_t) Q,sizeof(double));     /* eigenvalue storage */
  B = (double *) calloc((size_t) Q * max_col,sizeof(double)); /* Intermediate storage */
  C = (double *) calloc((size_t) Q * max_col,sizeof(double)); /* Intermediate storage */

  /* Start the main similarity transform loop */
  iter =0;
  while(1) {
    iter ++;

  /* Find the Frobenius norms of the Si in set gamma */
    max_frob=0.0;
    for (p=Si,i=0;i<Mf;i++,p += Q * Q) 
      if (gamma[i]) { /* don't bother if already dealt with */ 
        frob[i] = frobenius_norm(p,&Q,&Q);
        if (frob[i] *spf[i] >max_frob) max_frob=frob[i]  * spf[i];
    }
  /* Find sets alpha and gamma' */
    for (i=0;i<Mf;i++) {
      if (gamma[i]) { /* term is still to be dealt with */
        if (frob[i]  * spf[i] > max_frob * *d_tol) { 
          alpha[i] = 1;gamma1[i] = 0; /* deal with it now */
        } else {
          alpha[i] = 0;gamma1[i] = 1; /* put it off */ 
        }
      } else { /* wasn't in gamma, so not in alpha or gamma1 */
        alpha[i] = gamma1[i] = 0;
      }
    }

  /* Form the scaled sum of the Si in alpha and eigen-decompose it to get its rank */
    for (p=Sb,p1=p+Q*Q;p<p1;p++) *p=0.0; /* clear Sb */
    for (p=Si,i=0;i<Mf;i++,p += Q*Q) if (alpha[i]) { 
      x = frob[i];
      for (p1=p,p2=Sb,p3=p+Q*Q;p1<p3;p1++,p2++) *p2 += *p1 / x;
    } 
    mgcv_symeig(Sb,ev,&Q,&FALSE,&FALSE,&FALSE); /* get eigenvalues (ascending) of scaled sum over alpha */
    
    

    r=1;
    while(r<Q&&(ev[Q-r-1]>ev[Q-1] * *r_tol)) r++; 
    /* ...  r is the rank of Sb, or any other positively weighted sum over alpha */

    /* printf("\n iter = %d,  rank = %d,   Q = %d",iter,r,Q);
    printf("\n gamma = ");for (i=0;i<Mf;i++) printf(" %d",gamma[i]);
    printf("\n alpha = ");for (i=0;i<Mf;i++) printf(" %d",alpha[i]);
    printf("\n gamma1 = ");for (i=0;i<Mf;i++) printf(" %d",gamma1[i]);
    printf("\n ev = ");for (i=0;i<Q;i++) printf("  %g",ev[i]);*/

  /* If Q==r then terminate (form S first if it's the first iteration) */
    
    if (Q==r) { 
      if (iter==1 ) { /* form S */
        for (p=Si,i=0;i<Mf;i++,p += Q*Q) { 
          x = spf[i];
          for (p1=p,p2=S,p3=p+Q*Q;p1<p3;p1++,p2++) *p2 += *p1 * x;
        }
        break; 
      } else break; /* just use current S */ 
    } /* end if (Q==r) */

  /* Form the dominant term and eigen-decompose it */
    for (p=Sb,p1=p+Q*Q;p<p1;p++) *p = 0.0; /* clear Sb */
    for (p=Si,i=0;i<Mf;i++,p += Q*Q) if (alpha[i]) { /* summing S[[i]]*sp[i] over i in alpha */
      x = spf[i];
      for (p1=p,p2=Sb,p3=p+Q*Q;p1<p3;p1++,p2++) *p2 += *p1 * x;
    } 
    mgcv_symeig(Sb,ev,&Q,&FALSE,&TRUE,&TRUE); /* get eigen decomposition of dominant term (ev descending) */
    
  /* .... U points to Sb, which now contains eigen-vectors*/

  /* Form the sum over the elements in gamma1, Sg */

    for (p=Sg,p1=p+Q*Q;p<p1;p++) *p=0.0; /* clear Sg */
    for (p=Si,i=0;i<Mf;i++,p += Q*Q) if (gamma1[i]) { /* summing S[[i]]*sp[i] over i in gamma1 */
      x = spf[i];
      for (p1=p,p2=Sg,p3=p+Q*Q;p1<p3;p1++,p2++) *p2 += *p1 * x;
    } 

  /* Form S' the similarity transformed S */
    if (K>0) { /* deal with upper right component B */
      /* first copy out B into C */ 
      for (j=0;j<Q;j++) for (i=0;i<K;i++) C[i + K * j] = S[i + *q * (j+K)];
      /* Now form BU (store in B)*/
      bt=0;ct=0;mgcv_mmult(B,C,U,&bt,&ct,&K,&Q,&Q);
      /* Replace B into S */
      for (j=0;j<Q;j++) for (i=0;i<K;i++) S[i + *q * (j+K)]= S[j + K + *q * i] = B[i + K * j];
    }

    /* Now deal with the lower right component, C */
    /* U'SgU  */
    bt=0;ct=0;mgcv_mmult(B,Sg,U,&bt,&ct,&Q,&Q,&Q); /* SgU is in B */
    bt=1;ct=0;mgcv_mmult(C,U,B,&bt,&ct,&Q,&Q,&Q);  /* U'SgU is in C */ 
    for (i=0;i<r;i++) C[i+i * Q] += ev[i];  /* Adding in the (truly) non zero eigen-values */
   
    /* Now copy C back into right part of S' */
    for (j=0;j<Q;j++) for (i=0;i<Q;i++) S[i + K + *q * (j+K)] = C[i + Q * j];
    
  /* Transform the square roots of Si in alpha and gamma1 (Can leave fixed term alone - not needed)*/
    if (*deriv) { /* transformed rS_i only needed for derivatives */
      for (p=rS,k=0;k<*M;p += rSncol[k] * *q,k++) if (alpha[k]) {  /* p points to the square root of S_i */    
        /* extract the part of rS_i to be modified */
        for (i=0;i<Q;i++) for (j=0;j<rSncol[k];j++) C[i + Q * j] = p[i + K + *q * j];
        bt=1;ct=0;mgcv_mmult(B,U,C,&bt,&ct,&r,rSncol+k,&Q); 
        for (i=0;i<r;i++) for (j=0;j<rSncol[k];j++) p[i + K + *q * j] = B[i + r * j];
        for (i=K+r;i<K+Q;i++) for (j=0;j<rSncol[k];j++) p[i + *q * j] = 0.0;
      } else if (gamma1[k]) { 
        for (i=0;i<Q;i++) for (j=0;j<rSncol[k];j++) C[i + Q * j] = p[i + K + *q * j];
        bt=1;ct=0;mgcv_mmult(B,U,C,&bt,&ct,&Q,rSncol+k,&Q);
        for (i=0;i<Q;i++) for (j=0;j<rSncol[k];j++) p[i + K + *q * j] = B[i + Q * j];
      }
    }

  /* Transform the Si in gamma' */
    Qr = Q - r;Un = U + r * Q;
    for (p1=p=Si,i=0;i<Mf;i++,p += Q*Q,p1 +=Qr*Qr) if (gamma1[i]) { /* p points to old Si, and p1 to new */
      bt=1;ct=0;mgcv_mmult(B,Un,p,&bt,&ct,&Qr,&Q,&Q);
      bt=0;ct=0;mgcv_mmult(p1,B,Un,&bt,&ct,&Qr,&Qr,&Q); 
    }
  /* Update K, Q and gamma */   
    K = K + r; Q = Qr;
    for (i=0;i<Mf;i++) gamma[i] = gamma1[i];
  } /* end of Similarity Transfrom Loop */

  /* Now get the determinant and inverse of the transformed S (stored in B) */
  *det = qr_ldet_inv(S,q,B,deriv);
  /* finally, the dervivatives, based on transformed S inverse and transformed square roots */  
  
  if (*deriv) { /* get the first derivatives */
    for (p=rS,i=0;i<*M;p += *q *rSncol[i],i++) det1[i] = trBtAB(B,p,q,rSncol+i)*sp[i]; /* tr(S^{-1}S_i) */
  }
  
  if (*deriv==2) { /* get second derivatives, as well */
    for (p1=Si,p=rS,i=0;i<*M;p += *q *rSncol[i],i++,p1 += *q * *q) { /* loop through creating S^{-1} S_i and storing in Si*/
      bt=0;ct=0;mgcv_mmult(C,B,p,&bt,&ct,q,rSncol+i,q);
      bt=0;ct=1;mgcv_mmult(p1,C,p,&bt,&ct,q,q,rSncol+i);
    }
    for (i=0;i<*M;i++) for (j=i;j<*M;j++) 
      det2[i + *M * j] = det2[j + *M * i] = -sp[i]*sp[j]*trAB(Si + *q * *q *i,Si + *q * *q *j,q,q);
    for (i=0;i<*M;i++) det2[i + *M * i] += det1[i];
  }
  free(frob);
  free(gamma);
  free(gamma1);
  free(alpha);
  free(S);
  free(Sb);
  free(Sg);
  if (*deriv) { free(rS);}
  if (*fixed_penalty) {free(spf);}
  free(Si);
  free(ev);
  free(B);
  free(C);
}




void get_ddetXWXpS(double *det1,double *det2,double *P,double *K,double *sp,
             double *rS,int *rSncol,double *Tk,double *Tkm,int *n,int *q,int *r,int *M,int *deriv)

/* obtains derivatives of |X'WX + S| wrt the log smoothing parameters, as required for REML. 
   The determinant itself has to be obtained during intial decompositions: see gdi().

   * P is q by r
   * K is n by r
  
   * this routine assumes that sp contains smoothing parameters, rather than log smoothing parameters.
 
   * Note that P and K are as in Wood (2008) JRSSB 70, 495-518.

*/

{ double *diagKKt,xx,*KtTK,*PtrSm,*PtSP,*trPtSP,*work,*pdKK,*p1;
  int m,k,bt,ct,j,one=1,km,mk,rSoff,deriv2,max_col;
  if (*deriv==2) deriv2=1; else deriv2=0;
  /* obtain diag(KK') */ 
  if (*deriv) {
    diagKKt = (double *)calloc((size_t)*n,sizeof(double));
    xx = diagABt(diagKKt,K,K,n,r); 
  } else { /* nothing to do */
      return;
  }
  /* set up work space */
  work =  (double *)calloc((size_t)*n,sizeof(double));

  /* now loop through the smoothing parameters to create K'TkK */
  if (deriv2) {
    KtTK = (double *)calloc((size_t)(*r * *r * *M),sizeof(double));
    for (k=0;k < *M;k++) {
      j = k * *r * *r;
      getXtWX(KtTK+ j,K,Tk + k * *n,n,r,work);
    }
  } else { KtTK=(double *)NULL;} /* keep compiler happy */
  
  /* start first derivative */ 
  bt=1;ct=0;mgcv_mmult(det1,Tk,diagKKt,&bt,&ct,M,&one,n); /* tr(TkKK') */ 

  /* Finish first derivative and create create P'SmP if second derivs needed */
  
  max_col = *q;
  for (j=0;j<*M;j++) if (max_col<rSncol[j]) max_col=rSncol[j]; /* under ML can have q < max(rSncol) */

  PtrSm = (double *)calloc((size_t)(*r * max_col ),sizeof(double)); /* storage for P' rSm */
  trPtSP = (double *)calloc((size_t) *M,sizeof(double));

  if (deriv2) {
    PtSP = (double *)calloc((size_t)(*M * *r * *r ),sizeof(double));
  } else { PtSP = (double *) NULL;}

  for (rSoff=0,m=0;m < *M;m++) { /* loop through penalty matrices */
     bt=1;ct=0;mgcv_mmult(PtrSm,P,rS+rSoff * *q,&bt,&ct,r,rSncol+m,q);
     rSoff += rSncol[m];
     trPtSP[m] = sp[m] * diagABt(work,PtrSm,PtrSm,r,rSncol+m); /* sp[m]*tr(P'S_mP) */ 
     det1[m] += trPtSP[m]; /* completed first derivative */
     if (deriv2) { /* get P'S_mP */
       bt=0;ct=1;mgcv_mmult(PtSP+ m * *r * *r,PtrSm,PtrSm,&bt,&ct,r,r,rSncol+m);
     }
  }

  /* Now accumulate the second derivatives */

  if (deriv2) for (m=0;m < *M;m++) for (k=m;k < *M;k++){
     km=k * *M + m;mk=m * *M + k;
     /* tr(Tkm KK') */
     for (xx=0.0,pdKK=diagKKt,p1=pdKK + *n;pdKK<p1;pdKK++,Tkm++) xx += *Tkm * *pdKK;
     det2[km] = xx;

     /* - tr(KTkKK'TmK) */
     det2[km] -= diagABt(work,KtTK + k * *r * *r,KtTK+ m * *r * *r,r,r);

     /* sp[k]*tr(P'S_kP) */
     if (k==m) det2[km] += trPtSP[m];

     /* -sp[m]*tr(K'T_kKP'S_mP) */
     det2[km] -= sp[m]*diagABt(work,KtTK + k * *r * *r,PtSP + m * *r * *r,r,r);
     
     /* -sp[k]*tr(K'T_mKP'S_kP) */
     det2[km] -= sp[k]*diagABt(work,KtTK + m * *r * *r,PtSP + k * *r * *r,r,r);
 
     /* -sp[m]*sp[k]*tr(P'S_kPP'S_mP) */
     det2[km] -= sp[m]*sp[k]*diagABt(work,PtSP + k * *r * *r,PtSP + m * *r * *r,r,r);

     det2[mk] = det2[km];     
  }


 
  /* free up some memory */
  if (deriv2) {free(PtSP);free(KtTK);}
  free(diagKKt);free(work);
  free(PtrSm);free(trPtSP);

}


void get_trA2(double *trA,double *trA1,double *trA2,double *P,double *K,double *sp,
	      double *rS,int *rSncol,double *Tk,double *Tkm,double *w,int *n,int *q,int *r,int *M,int *deriv)

/* obtains trA and its first two derivatives wrt the log smoothing parameters 
   * P is q by r
   * K is n by r
   * U1 is q by r
   * this routine assumes that sp contains smoothing parameters, rather than log smoothing parameters.

   * If deriv is 0 then only tr(A) is obtained here.
   * This version uses only K and P, and is for the case where expressions involve weights which
     are reciprocal varainces, not the squares of weights which are reciprocal standard deviations.


*/

{ double *diagKKt,*diagKKtKKt,xx,*KtTK,*KtTKKtK,*KKtK,*KtK,*work,*pTk,*pTm,*pdKKt,*pdKKtKKt,*p0,*p1,*p2,*p3,*pd,
    *PtrSm,*PtSP,*KPtrSm,*diagKPtSPKt,*diagKPtSPKtKKt,*PtSPKtK, *KtKPtrSm, *KKtKPtrSm,*Ip,*IpK;
  int i,m,k,bt,ct,j,one=1,km,mk,rSoff,deriv2,neg_w=0;
  if (*deriv==2) deriv2=1; else deriv2=0;
  /* Get the sign array for negative w_i */
  Ip = (double *)calloc((size_t)*n,sizeof(double));
  for (p0=w,p1=p0+ *n,p2=Ip;p0<p1;p0++,p2++) if (*p0 < 0) {*p2 = -1.0;neg_w=1;} else *p2 = 1.0;

  /* obtain tr(A) and diag(A) = diag(KK'Ip) */ 
  diagKKt = (double *)calloc((size_t)*n,sizeof(double));
  *trA = diagABt(diagKKt,K,K,n,r); 
  if (neg_w) { /* correct trA */
    for (*trA=0.0,p0=diagKKt,p1=p0 + *n,p2=Ip;p0<p1;p0++,p2++) *trA += *p2 * *p0;
  }
  if (!*deriv) {
    free(Ip);free(diagKKt);
    return;
  }

  /* set up work space */
  work =  (double *)calloc((size_t)*n,sizeof(double));
  /* Get K'IpK and KK'IpK  */
  KtK = (double *)calloc((size_t)*r * *r,sizeof(double));
  if (neg_w) { 
    IpK = (double *)calloc((size_t) *r * *n,sizeof(double));
    for (p0=IpK,p3=K,i=0;i<*r;i++) 
      for (p1=Ip,p2=p1 + *n;p1<p2;p1++,p0++,p3++) *p0 = *p1 * *p3; 
  } else IpK = K;
  bt=1;ct=0;mgcv_mmult(KtK,K,IpK,&bt,&ct,r,r,n);  
  KKtK = (double *)calloc((size_t)*n * *r,sizeof(double));
  bt=0;ct=0;mgcv_mmult(KKtK,K,KtK,&bt,&ct,n,r,r);  

  /* obtain diag(KK'KK') */
  diagKKtKKt = (double *)calloc((size_t)*n,sizeof(double));
  xx = diagABt(diagKKtKKt,KKtK,K,n,r);
 
  /* now loop through the smoothing parameters to create K'TkK and K'TkKK'K */
  if (deriv2) {
    KtTK = (double *)calloc((size_t)(*r * *r * *M),sizeof(double));
    KtTKKtK = (double *)calloc((size_t)(*r * *r * *M),sizeof(double));
    for (k=0;k < *M;k++) {
      j = k * *r * *r;
      getXtWX(KtTK+ j,K,Tk + k * *n,n,r,work);
      bt=ct=0;mgcv_mmult(KtTKKtK + k * *r * *r ,KtTK + j,KtK,&bt,&ct,r,r,r);
    }
  } else { KtTK=KtTKKtK=(double *)NULL;}
  
  /* evaluate first and last terms in first derivative of tr(F) */
  bt=1;ct=0;mgcv_mmult(trA1,Tk,diagKKt,&bt,&ct,M,&one,n); /* tr(KK'Tk) */ 
  bt=1;ct=0;mgcv_mmult(work,Tk,diagKKtKKt,&bt,&ct,M,&one,n); /* tr(KK'TkKK') */
  for (i=0;i<*M;i++) trA1[i] +=  - work[i];

  /* now evaluate terms in Hessian of tr(F) which depend on what's available so far */
  if (deriv2) for (m=0;m < *M;m++) for (k=m;k < *M;k++){
     km=k * *M + m;mk=m * *M + k;

     /* tr(KK'Tkm  - KK'TkKK') */
     for (xx=0.0,pdKKt=diagKKt,pdKKtKKt=diagKKtKKt,p1=pdKKt + *n;pdKKt<p1;pdKKt++,pdKKtKKt++,Tkm++) 
          xx += *Tkm * (*pdKKt - *pdKKtKKt);
     trA2[km] = xx;

     /* -2 tr(K'TkKK'TmK)*/
     trA2[km] -= 2*diagABt(work,KtTK + k * *r * *r,KtTK+ m * *r * *r,r,r);

     /* 2 tr(K'TkKK'TmKK'K) -- needs correction*/
     xx = 2*diagABt(work,KtTK+k * *r * *r,KtTKKtK+m * *r * *r,q,r);
    
     trA2[km] += xx;

     trA2[mk] = trA2[km];     
  }

  /* free up some memory */
  if (deriv2) {free(KtTKKtK);free(KtTK);} 

  free(diagKKtKKt);free(diagKKt);

  /* create KP'rSm, KK'KP'rSm and P'SmP */
  PtrSm = (double *)calloc((size_t)(*r * *q ),sizeof(double)); /* transient storage for P' rSm */
  KPtrSm = (double *)calloc((size_t)(*n * *q),sizeof(double)); /* transient storage for K P' rSm */
  diagKPtSPKt = (double *)calloc((size_t)(*n * *M),sizeof(double));
  if (deriv2) {
    PtSP = (double *)calloc((size_t)(*M * *r * *r ),sizeof(double));
    PtSPKtK = (double *)calloc((size_t)(*M * *r * *r ),sizeof(double));
    KtKPtrSm = (double *)calloc((size_t)(*r * *q),sizeof(double));/* transient storage for K'K P'rSm */ 
    KKtKPtrSm = (double *)calloc((size_t)(*n * *q),sizeof(double));/* transient storage for K'K P'rSm */ 
    diagKPtSPKtKKt = (double *)calloc((size_t)(*n * *M),sizeof(double));
  } else { PtSP=KtKPtrSm=diagKPtSPKtKKt=(double *)NULL; }
  for (rSoff=0,m=0;m < *M;m++) {
    bt=1;ct=0;mgcv_mmult(PtrSm,P,rS+rSoff * *q,&bt,&ct,r,rSncol+m,q);
    bt=0;ct=0;mgcv_mmult(KPtrSm,K,PtrSm,&bt,&ct,n,rSncol+m,r); 
    if (deriv2) {
      bt=0;ct=0;mgcv_mmult(KtKPtrSm,KtK,PtrSm,&bt,&ct,r,rSncol+m,r); 
      bt=0;ct=1;mgcv_mmult(PtSP+ m * *r * *r,PtrSm,PtrSm,&bt,&ct,r,r,rSncol+m);
    
      bt=0;ct=0;mgcv_mmult(KKtKPtrSm,KKtK,PtrSm,&bt,&ct,n,rSncol+m,r);      
      bt=0;ct=1;mgcv_mmult(PtSPKtK + m * *r * *r,PtrSm,KtKPtrSm,&bt,&ct,r,r,rSncol+m); 
      xx = diagABt(diagKPtSPKtKKt+ m * *n,KPtrSm,KKtKPtrSm,n,rSncol+m);
    }
    rSoff += rSncol[m];
    xx = sp[m] * diagABt(diagKPtSPKt+ m * *n,KPtrSm,KPtrSm,n,rSncol+m);
       if (neg_w) { /* have to correct xx for negative w_i */
      for (xx=0.0,p0=diagKPtSPKt+m * *n,p1=p0 + *n,p2=Ip;p0<p1;p0++,p2++) xx += *p0 * *p2;
      xx *= sp[m];
    }
    trA1[m] -= xx; /* finishing trA1 */
    if (deriv2) trA2[m * *M + m] -=xx; /* the extra diagonal term of trA2 */
  }
  if (!deriv2) { /* trA1 finished, so return */
    free(PtrSm);free(KPtrSm);free(diagKPtSPKt);
    free(work);free(KtK);free(KKtK);
    return;
  }
  /* now use these terms to finish off the Hessian of tr(F) */ 
   for (m=0;m < *M;m++) for (k=m;k < *M;k++){
     km=k * *M + m;mk=m * *M + k;

     /* 2 sp[m] tr(KK'TkKP'SmPK') */
     pTk = Tk + k * *n;
     for (xx=0.0,pd = diagKPtSPKtKKt + m * *n,p1=pd + *n;pd < p1;pd++,pTk++) xx += *pd * *pTk;
     trA2[km] += 2*sp[m] *xx;

     /* 2 sp[k] tr(KK'TmKP'SkPK') */
     pTm = Tk + m * *n;
     for (xx=0.0,pd = diagKPtSPKtKKt + k * *n,p1=pd + *n;pd < p1;pd++,pTm++) xx += *pd * *pTm;
     trA2[km] += 2*sp[k] *xx;
     
     /* - sp[m] tr(TkKP'SmPK') */
     pTk = Tk + k * *n;
     for (xx=0.0,pd = diagKPtSPKt + m * *n,p1=pd + *n;pd < p1;pd++,pTk++) xx += *pd * *pTk;
     trA2[km] -= sp[m] *xx;
     
     /* - sp[k] tr(TmKP'SkPK') */
     pTm = Tk + m * *n;
     for (xx=0.0,pd = diagKPtSPKt + k * *n,p1=pd + *n;pd < p1;pd++,pTm++) xx += *pd * *pTm;
     trA2[km] -= sp[k] *xx;

     /* 2 sp[m] sp[k] tr(KP'SkPP'SmPK') */
     trA2[km] += 2 * sp[k]*sp[m]*diagABt(work,PtSP + m * *r * *r,PtSPKtK + k * *r * *r,q,r);
      
     trA2[mk] =trA2[km];
   } 
   /* clear up */
   free(PtrSm);free(KPtrSm);free(PtSP);free(KtKPtrSm);free(diagKPtSPKt);
   free(diagKPtSPKtKKt);free(work);free(KtK);free(KKtK);free(PtSPKtK);free(KKtKPtrSm);
   free(Ip);if (neg_w) free(IpK);  
}



void get_trA(double *trA,double *trA1,double *trA2,double *U1,double *KU1t,double *P,double *K,double *sp,
             double *rS,int *rSncol,double *Tk,double *Tkm,int *n,int *q,int *r,int *M,int *deriv)

/* obtains trA and its first two derivatives wrt the log smoothing parameters 
   * P is q by r
   * K is n by r
   * U1 is q by r
   * this routine assumes that sp contains smoothing parameters, rather than log smoothing parameters.

   * If deriv is 0 then only tr(A) is obtained here.
   * Note that P and K are as in Wood (2008) JRSSB 70, 495-518, but the expressions used here
     are slightly different. They are as efficient, but are from before I noticed that U1 
     could be eliminated. 

   Used by Wood 2008. Now superceded by get_trA2.
*/

{ double *diagA,*diagAA,xx,*KtTK,*U1KtTK,*work,*pTk,*pTm,*pdA,*pdAA,*p1,*pd,
         *PtrSm,*U1PtrSm,*U1PtSP,*KPtrSm,*KU1tU1PtrSm,*diagBtSB,*diagBtSBA;
  int i,m,k,bt,ct,j,one=1,km,mk,rSoff,deriv2;
  if (*deriv==2) deriv2=1; else deriv2=0;
  /* obtain tr(A) and diag(A) */ 
  if (*deriv) {
    diagA = (double *)calloc((size_t)*n,sizeof(double));
    *trA = diagABt(diagA,K,K,n,r);
  } else { /* then only tr(A) is required so return now*/
      for (xx=0.0,i=0,j=i+ *q * *r;i<j;i++,U1++) xx+= *U1 * *U1;
      *trA = xx;
      return;
  }
  /* set up work space */
  work =  (double *)calloc((size_t)*n,sizeof(double));
  /* obtain diag(AA) */
  diagAA = (double *)calloc((size_t)*n,sizeof(double));
  xx = diagABt(diagAA,KU1t,KU1t,n,q);
  /* now loop through the smoothing parameters to create K'TkK and U1K'TkK */
  if (deriv2) {
    KtTK = (double *)calloc((size_t)(*r * *r * *M),sizeof(double));
    U1KtTK = (double *)calloc((size_t)(*q * *r * *M),sizeof(double));
    for (k=0;k < *M;k++) {
      j = k * *r * *r;
      getXtWX(KtTK+ j,K,Tk + k * *n,n,r,work);
      bt=ct=0;mgcv_mmult(U1KtTK+ k * *q * *r ,U1,KtTK + j,&bt,&ct,q,r,r);
    }
  } else { KtTK=U1KtTK=(double *)NULL;}
  
  /* evaluate first term in first derivative of tr(A) */
  bt=1;ct=0;mgcv_mmult(trA1,Tk,diagA,&bt,&ct,M,&one,n); /* tr(TkA) */ 
  bt=1;ct=0;mgcv_mmult(work,Tk,diagAA,&bt,&ct,M,&one,n); /* tr(ATkA) */
  for (i=0;i<*M;i++) trA1[i] = 2*(trA1[i] - work[i]);
  
  /* now evaluate terms in Hessian of tr(A) which depend on what's available so far */
  if (deriv2) for (m=0;m < *M;m++) for (k=m;k < *M;k++){
     km=k * *M + m;mk=m * *M + k;
     /* 2tr(Tkm A - ATkmA) */
     for (xx=0.0,pdA=diagA,pdAA=diagAA,p1=pdA + *n;pdA<p1;pdA++,pdAA++,Tkm++) xx += *Tkm * (*pdA - *pdAA);
     trA2[km] = 2*xx;

     /* 4tr(TkTmA - ATmTkA) */
     pTk = Tk + k * *n;pTm = Tk + m * *n;
     for (xx=0.0,pdA=diagA,pdAA=diagAA,p1=pdA + *n;pdA<p1;pdA++,pdAA++,pTk++,pTm++) 
     xx += *pTk * *pTm  * (*pdA - *pdAA);
     trA2[km] += 4*xx; 

     /* -4 tr(TkATmA + TmATkA) */
     trA2[km] -= 8*diagABt(work,KtTK + k * *r * *r,KtTK+ m * *r * *r,r,r);

     /* 8 tr(ATkATmA) */
     trA2[km] += 8*diagABt(work,U1KtTK+k * *q * *r,U1KtTK+m * *q * *r,q,r);
 
     trA2[mk] = trA2[km];     
  }

  /* free up some memory */
  if (deriv2) {free(U1KtTK);free(KtTK);}

  free(diagAA);free(diagA);

  /* create KP'rSm, KU1tU1P'rSm and U1P'SmP */
  PtrSm = (double *)calloc((size_t)(*r * *q ),sizeof(double)); /* transient storage for P' rSm */
  KPtrSm = (double *)calloc((size_t)(*n * *q),sizeof(double)); /* transient storage for K P' rSm */
  diagBtSB = (double *)calloc((size_t)(*n * *M),sizeof(double));
  if (deriv2) {
    U1PtrSm = (double *)calloc((size_t)(*q * *q ),sizeof(double)); /* transient storage for U1 P' rSm */
    U1PtSP = (double *)calloc((size_t)(*M * *q * *r ),sizeof(double));
    KU1tU1PtrSm = (double *)calloc((size_t)(*n * *q),sizeof(double));/* transient storage for K U1'U1 P'rSm */ 
    diagBtSBA = (double *)calloc((size_t)(*n * *M),sizeof(double));
  } else {U1PtrSm=U1PtSP=KU1tU1PtrSm=diagBtSBA=(double *)NULL; }
  for (rSoff=0,m=0;m < *M;m++) {
    bt=1;ct=0;mgcv_mmult(PtrSm,P,rS+rSoff * *q,&bt,&ct,r,rSncol+m,q);
    bt=0;ct=0;mgcv_mmult(KPtrSm,K,PtrSm,&bt,&ct,n,rSncol+m,r); 
    if (deriv2) {
      bt=0;ct=0;mgcv_mmult(U1PtrSm,U1,PtrSm,&bt,&ct,q,rSncol+m,r); 
      bt=0;ct=1;mgcv_mmult(U1PtSP+ m * *q * *r,U1PtrSm,PtrSm,&bt,&ct,q,r,rSncol+m);
      /* Now do KU1tU1P'rSm, recycling PtrSm as transient storage */
      bt=1;ct=0;mgcv_mmult(PtrSm,U1,U1PtrSm,&bt,&ct,r,rSncol+m,q);
      bt=0;ct=0;mgcv_mmult(KU1tU1PtrSm,K,PtrSm,&bt,&ct,n,rSncol+m,r);      
      xx = diagABt(diagBtSBA+ m * *n,KPtrSm,KU1tU1PtrSm,n,rSncol+m);
    }
    rSoff += rSncol[m];
    xx = sp[m] * diagABt(diagBtSB+ m * *n,KPtrSm,KPtrSm,n,rSncol+m);
    trA1[m] -= xx; /* finishing trA1 */
    if (deriv2) trA2[m * *M + m] -=xx; /* the extra diagonal term of trA2 */
  }
  if (!deriv2) { /* trA1 finished, so return */
    free(PtrSm);free(KPtrSm);free(diagBtSB);
    return;
  }
  /* now use these terms to finish off the Hessian of tr(A) */ 
   for (m=0;m < *M;m++) for (k=m;k < *M;k++){
     km=k * *M + m;mk=m * *M + k;

     /* 4 sp[m] tr(ATkB'SmB) */
     pTk = Tk + k * *n;
     for (xx=0.0,pd = diagBtSBA + m * *n,p1=pd + *n;pd < p1;pd++,pTk++) xx += *pd * *pTk;
     trA2[km] += 4*sp[m] *xx;

     /* 4 sp[k] tr(ATmB'SkB) */
     pTm = Tk + m * *n;
     for (xx=0.0,pd = diagBtSBA + k * *n,p1=pd + *n;pd < p1;pd++,pTm++) xx += *pd * *pTm;
     trA2[km] += 4*sp[k] *xx;
     
     /* -2 sp[m] tr(TkB'SmB) */
     pTk = Tk + k * *n;
     for (xx=0.0,pd = diagBtSB + m * *n,p1=pd + *n;pd < p1;pd++,pTk++) xx += *pd * *pTk;
     trA2[km] -= 2*sp[m] *xx;
     
     /* -2 sp[k] tr(TmB'SkB) */
     pTm = Tk + m * *n;
     for (xx=0.0,pd = diagBtSB + k * *n,p1=pd + *n;pd < p1;pd++,pTm++) xx += *pd * *pTm;
     trA2[km] -= 2*sp[k] *xx;

     /* 2 sp[m] sp[k] tr(B'SmG^{-1}SkB) */
     trA2[km] += 2 * sp[k]*sp[m]*diagABt(work,U1PtSP + m * *q * *r,U1PtSP + k * *q * *r,q,r);
      
     trA2[mk] =trA2[km];
   } 
   /* clear up */
   free(PtrSm);free(U1PtrSm);free(U1PtSP);free(KPtrSm);free(KU1tU1PtrSm);free(diagBtSB);
   free(diagBtSBA);free(work);
   
}


void B1B2zBaseSetup(double *B2z,double *B1z,double *z,double *P,double *K,
           double *KKtz,double *PKtz,double *KPtSPKtz,double *rS,
           int *rSncol,int *n,int *q, int *r,int *M,double *sp,double *work,
           int *deriv)

/* Initializes B1z, B2z and creates
   KKtz, PKtz and KPSPKtz
   work must have dimension of at least 2*n+M*q

   Used by Wood 2008. Now defunct.
*/

{ double *PPtSPKtz,*v1,*v2,*dp,*dp0,*dp1,*pB2z,*pPPtSPKtz,xx;
  int i,k,one=1,m,bt,ct,deriv2;
  /* A. portion out work */
  if (*deriv==2) deriv2=1; else deriv2=0;
  dp=work;
  v1 = dp;dp += *n;  
  v2 = dp;dp += *n;
  pPPtSPKtz=PPtSPKtz = dp; dp += *q * *M;
  /* B. create KKtz and PKtz */
  bt=1;ct=0;mgcv_mmult(v1,K,z,&bt,&ct,r,&one,n);
  bt=0;ct=0;mgcv_mmult(KKtz,K,v1,&bt,&ct,n,&one,r);
  bt=0;ct=0;mgcv_mmult(PKtz,P,v1,&bt,&ct,q,&one,r);
  /* C. loop through sp's creating PP'SkPK'z, KP'SkPK'z and intializing B1z */
  for (k=0;k < *M;k++) {
    multSk(v1,PKtz,&one,k,rS,rSncol,q,v2);
    bt=1;ct=0;mgcv_mmult(v2,P,v1,&bt,&ct,r,&one,q);
    bt=0;ct=0;mgcv_mmult(pPPtSPKtz,P,v2,&bt,&ct,q,&one,r);
    if (deriv2) {
     bt=0;ct=0;mgcv_mmult(KPtSPKtz,K,v2,&bt,&ct,n,&one,r);
     KPtSPKtz += *n; /* move to next slot */
    }
    xx = -sp[k];
    for (i=0;i<*q;i++,B1z++,pPPtSPKtz++) *B1z = xx * *pPPtSPKtz;  
    
  }
  /* D. double loop through sps to set up B2z */
  if (deriv2)
  { pB2z=B2z;
    for (m=0;m < *M;m++)
    for (k=m;k < *M;k++)
    { /* 1. obtain PP'SmPP'SkPK'z */
      pPPtSPKtz = PPtSPKtz + k * *q;
      multSk(v1,pPPtSPKtz,&one,m,rS,rSncol,q,v2);    
      bt=1;ct=0;mgcv_mmult(v2,P,v1,&bt,&ct,r,&one,q);
      bt=0;ct=0;mgcv_mmult(v1,P,v2,&bt,&ct,q,&one,r);
      /* result now in v1, put (sp[m]*sp[k])*v1 into the relevant part of B2z */
      dp1 = v1 + *q;xx=sp[m]*sp[k];
      for (dp=v1,dp0=pB2z;dp<dp1;dp++,dp0++) *dp0 = xx * *dp;

      /* 2. obtain PP'SkPP'SmPK'z (simple k,m interchange of term 1)*/
      pPPtSPKtz = PPtSPKtz + m * *q;
      multSk(v1,pPPtSPKtz,&one,k,rS,rSncol,q,v2);    
      bt=1;ct=0;mgcv_mmult(v2,P,v1,&bt,&ct,r,&one,q);
      bt=0;ct=0;mgcv_mmult(v1,P,v2,&bt,&ct,q,&one,r);
      /* result now in v1, add (sp[m]*sp[k])*v1 to relevant part of B2z */
      dp1 = v1 + *q;xx=sp[m]*sp[k];
      for (dp=v1,dp0=pB2z;dp<dp1;dp++,dp0++) *dp0 += xx * *dp;
    
      /* 3. the final term PP'SkPK' */
      if (m==k) {
	  dp = PPtSPKtz + k * *q;
          dp1 = dp + *q;xx=sp[k];
          for (dp0=pB2z;dp<dp1;dp++,dp0++) *dp0 -= xx * *dp;
      }

      pB2z += *q;
    }
  } /* end of if (*deriv=2) */
}



void B1B2z(double *B2z,double *B1z,double *B2zBase,double *B1zBase,double *z,
           double *Tk,double *Tkm,double *P,double *K,
           double *KKtz,double *PKtz,double *KPtSPKtz,double *rS,
           int *rSncol,int *n,int *q, int *r,int *M,double *sp,double *work,
           int *deriv)
/* Routine to apply first and second derivatives of B to z'.
   Some key dimensions:
   * M is the number of smoothing parameters
   * z is an n-vector
   * Tk is packed as a first derivative structure
   * TKM is packed as a second derivative structure.
   * P is q by r
   * K is n by r
   * KKtz is an n-vector
   * PKtz is a q-vector
   * KPSPKtz[i] is an n-vector
   * the rS are q by rSncol[i] matrices 

   B2zBase and B1zBase are the parts of the derivatives that do not change
   through iteration. 

   * work should contain a block of memory for (2*q+4*n)*M + 2*n doubles.
   * deriv==2 for second derivatives, otherwise first derivs only

   Used by Wood 2008. No defunct.
*/
{ int i,k,bt,ct,one=1,m,deriv2;
  double *dp,*dp0,*dp1,*dp2,*TKKtz,*PKtTKKtz,*KKtTKKtz,*PKtTz,*Tz,*KKtTz,*v1,*v2,
         *pTk,*pTKKtz,*pPKtTKKtz,*pKKtTKKtz,*pPKtTz,*pTz,*pKKtTz,*pTkm,*pB2z,
         xx,*pKKtz,*pKPtSPKtz;
  if (*deriv==2) deriv2=1; else deriv2=0;
  /* A. Farm out workspace */
  dp = work;v1 = dp;dp += *n;v2 = dp; dp += *n;
  TKKtz = dp; dp += *M * *n;
  KKtTKKtz = dp; dp += *M * *n;
  PKtTKKtz = dp; dp += *M * *q;
  PKtTz = dp; dp += *M * *q;
  Tz = dp; dp += *M * *n;
  KKtTz = dp; dp += *M * *n;
  /* B. initialize B2z and B1z to base values */
  dp1 = B1z + *q * *M;
  for (dp=B1z,dp0=B1zBase;dp<dp1;dp++,dp0++) *dp = *dp0;
  if (deriv2) {
    dp1 = B2z + *q * *M * (1 + *M) /2;
    for (dp=B2z,dp0=B2zBase;dp<dp1;dp++,dp0++) *dp = *dp0;
  }
  /* C. Initial loop through smoothing parameters, creating:
        * TKKtz[i], PKtTKKtz[i], KKtTKKtz[i], PKtTz[i],Tz[i],KKtTz[i]
        * B1z
  */
  pTk = Tk;
  pTKKtz=TKKtz;
  pPKtTKKtz = PKtTKKtz;
  pKKtTKKtz=KKtTKKtz;
  pTz = Tz;
  pKKtTz=KKtTz;
  pPKtTz=PKtTz;
  for (k=0;k<*M;k++) { /* loop through smoothing parameters */
    /* form TKKtz[i] */
    dp1 = pTk + *n;  
    for (dp = pTk,dp0=KKtz;dp<dp1;dp++,pTKKtz++,dp0++) *pTKKtz = *dp * *dp0; 
    /* Now form r-vector v1 = KtTKKtz */
    bt=1;ct=0;
    mgcv_mmult(v1,K,pTKKtz - *n,&bt,&ct,r,&one,n);
    /* ... which is the basis for PKtTKKtz */
    bt=0;ct=0;
    mgcv_mmult(pPKtTKKtz,P,v1,&bt,&ct,q,&one,r);
    /* ... and also for KKtTKKtz */
    if (deriv2) mgcv_mmult(pKKtTKKtz,K,v1,&bt,&ct,n,&one,r);
    /* Form Tz... */
    dp1 = pTk + *n;
    for (dp0=z,dp=pTk;dp<dp1;dp++,pTz++,dp0++) *pTz = *dp0 * *dp;
    /* Form r-vector v1 = KtTz */
    bt=1;ct=0;
    mgcv_mmult(v1,K,pTz - *n,&bt,&ct,r,&one,n);
    /* and hence PKtTz */
    bt=0;ct=0;
    mgcv_mmult(pPKtTz,P,v1,&bt,&ct,q,&one,r);
    /* and also KKtTz */
    bt=0;ct=0;
    if (deriv2) mgcv_mmult(pKKtTz,K,v1,&bt,&ct,n,&one,r);
    /* can now update B1z */
    for (i=0;i < *q;i++,B1z++,pPKtTz++,pPKtTKKtz++)
	*B1z += *pPKtTz - 2 * *pPKtTKKtz; 
    /* move pointers to next derivative */ 
    if (deriv2) {
      pKKtTKKtz += *n;
      pKKtTz += *n;
    }
    pTk += *n;
  }
  if (!deriv2) return; /* only first derivatives needed */
  /* D. double loop through smoothing parameters to obtain B2z
  */
  pTkm=Tkm;pB2z=B2z;
  for (m=0;m < *M;m++)
  for (k=m;k < *M;k++)
  { /* 1. obtain PK'TmKK'TkKK'z */
    pKKtTKKtz = KKtTKKtz + k * *n;
    dp = Tk + m * *n; /* pointer to Tm array*/
    dp1 = v1 + *n; /* end of v1 = TmKK'TkKK'z */
    for (dp0=v1;dp0<dp1;dp0++,dp++,pKKtTKKtz++) *dp0 = *dp * *pKKtTKKtz;
    bt=1;ct=0;mgcv_mmult(v2,K,v1,&bt,&ct,r,&one,n);
    bt=0;ct=0;mgcv_mmult(v1,P,v2,&bt,&ct,q,&one,r);
    /* result now in v1, add 4*v1 to relevant part of B2z */
    dp1 = v1 + *q;
    for (dp=v1,dp0=pB2z;dp<dp1;dp++,dp0++) *dp0 += 4 * *dp;
 
    /* 2. obtain PP'SmPK'TkKK'z */ 
    pPKtTKKtz = PKtTKKtz + k * *q;
    multSk(v1,pPKtTKKtz,&one,m,rS,rSncol,q,v2);
    bt=1;ct=0;mgcv_mmult(v2,P,v1,&bt,&ct,r,&one,q);
    bt=0;ct=0;mgcv_mmult(v1,P,v2,&bt,&ct,q,&one,r);
    /* result now in v1 - add to B2z */
    xx = 2 * sp[m];
    dp1 = v1 + *q;
    for (dp=v1,dp0=pB2z;dp<dp1;dp++,dp0++) *dp0 += xx * *dp;

    /* 3. obtain PK'TmTkKK'z */ 
    pTKKtz = TKKtz + k * *n;
    dp0 = Tk + m * *n;dp1 = dp0 + *n;
    for (dp=dp0,dp2=v1;dp < dp1;dp++,pTKKtz++,dp2++) *dp2 = *dp * *pTKKtz;
    bt=1;ct=0;mgcv_mmult(v2,K,v1,&bt,&ct,r,&one,n);
    bt=0;ct=0;mgcv_mmult(v1,P,v2,&bt,&ct,q,&one,r);
    /* result now in v1 - for addition to B2z */
    dp1 = v1 + *q;
    for (dp=v1,dp0=pB2z;dp<dp1;dp++,dp0++) *dp0 -= 4 * *dp;

    /* 4. obtain PK'TkmKK'z */
    pKKtz = KKtz;
    dp1 = v1 + *n;
    for (dp0=pTkm,dp=v1;dp < dp1;dp++,dp0++,pKKtz++) *dp = *dp0 * *pKKtz;
    bt=1;ct=0;mgcv_mmult(v2,K,v1,&bt,&ct,r,&one,n);
    bt=0;ct=0;mgcv_mmult(v1,P,v2,&bt,&ct,q,&one,r);    
    /* add result in v1 to B2z */
    dp1 = v1 + *q;
    for (dp=v1,dp0=pB2z;dp<dp1;dp++,dp0++) *dp0 -= 2 * *dp;

    /* 5. obtain PK'TkKK'TmKK'z (code modified from m,k version of term) */
    pKKtTKKtz = KKtTKKtz + m * *n;
    dp = Tk + k * *n; /* pointer to Tk array*/
    dp1 = v1 + *n; /* end of v1 = TkKK'TmKK'z */
    for (dp0=v1;dp0<dp1;dp0++,dp++,pKKtTKKtz++) *dp0 = *dp * *pKKtTKKtz;
    bt=1;ct=0;mgcv_mmult(v2,K,v1,&bt,&ct,r,&one,n);
    bt=0;ct=0;mgcv_mmult(v1,P,v2,&bt,&ct,q,&one,r);
    /* result now in v1, add 4*v1 to relevant part of B2z */
    dp1 = v1 + *q;
    for (dp=v1,dp0=pB2z;dp<dp1;dp++,dp0++) *dp0 += 4 * *dp;

    /* 6. obtain PK'TkKP'SmPK'z */
    pKPtSPKtz = KPtSPKtz + m * *n;
    dp = Tk + k * *n;dp1 = dp + *n;
    for (dp0=v1;dp<dp1;dp++,dp0++,pKPtSPKtz++) *dp0 = *dp * *pKPtSPKtz;
    bt=1;ct=0;mgcv_mmult(v2,K,v1,&bt,&ct,r,&one,n);
    bt=0;ct=0;mgcv_mmult(v1,P,v2,&bt,&ct,q,&one,r);
    /* result now in v1 - add into B2z */
    xx = 2 * sp[m];
    dp1 = v1 + *q;
    for (dp=v1,dp0=pB2z;dp<dp1;dp++,dp0++) *dp0 += xx * *dp;

    /* 7. obtain PK'TkKK'Tmz  */
    pKKtTz = KKtTz + m * *n;
    dp = Tk + k * *n;dp1 = dp + *n;
    for (dp0=v1;dp < dp1;dp++,dp0++,pKKtTz++) *dp0 = *dp * *pKKtTz;
    bt=1;ct=0;mgcv_mmult(v2,K,v1,&bt,&ct,r,&one,n);
    bt=0;ct=0;mgcv_mmult(v1,P,v2,&bt,&ct,q,&one,r);
    /* result now in v1, add -2*v1 to relevant part of B2z */
    dp1 = v1 + *q;
    for (dp=v1,dp0=pB2z;dp<dp1;dp++,dp0++) *dp0 -= 2 * *dp;

    /* 8. obtain PK'TmKK'Tkz  (code simple modification of term 7) */
    pKKtTz = KKtTz + k * *n;
    dp = Tk + m * *n;dp1 = dp + *n;
    for (dp0=v1;dp < dp1;dp++,dp0++,pKKtTz++) *dp0 = *dp * *pKKtTz;
    bt=1;ct=0;mgcv_mmult(v2,K,v1,&bt,&ct,r,&one,n);
    bt=0;ct=0;mgcv_mmult(v1,P,v2,&bt,&ct,q,&one,r);
    /* result now in v1, add -2*v1 to relevant part of B2z */
    dp1 = v1 + *q;
    for (dp=v1,dp0=pB2z;dp<dp1;dp++,dp0++) *dp0 -= 2 * *dp;

    /* 9. obtain PP'SmPK'Tkz */ 
    pPKtTz = PKtTz + k * *q;
    multSk(v1,pPKtTz,&one,m,rS,rSncol,q,v2);    
    bt=1;ct=0;mgcv_mmult(v2,P,v1,&bt,&ct,r,&one,q);
    bt=0;ct=0;mgcv_mmult(v1,P,v2,&bt,&ct,q,&one,r);
    /* result now in v1, add -sp[m]*v1 to relevant part of B2z */
    dp1 = v1 + *q;xx=sp[m];
    for (dp=v1,dp0=pB2z;dp<dp1;dp++,dp0++) *dp0 -= xx * *dp;

    /* 10. obtain PK'TmTkz */
    pTz = Tz + k * *n;
    dp = Tk + m * *n;dp1 = dp + *n;
    for (dp0 = v1;dp < dp1;dp++,dp0++,pTz++) *dp0 = *dp * *pTz;
    bt=1;ct=0;mgcv_mmult(v2,K,v1,&bt,&ct,r,&one,n);
    bt=0;ct=0;mgcv_mmult(v1,P,v2,&bt,&ct,q,&one,r);
    /* result now in v1, add v1 to relevant part of B2z */
    dp1 = v1 + *q;
    for (dp=v1,dp0=pB2z;dp<dp1;dp++,dp0++) *dp0 += *dp;
        
    /* 11. obtain PK'Tkmz */
    for (dp=pTkm,dp1=pTkm + *n,dp0=v1,dp2=z;dp < dp1;dp++,dp0++,dp2++) *dp0 = *dp * *dp2;
    bt=1;ct=0;mgcv_mmult(v2,K,v1,&bt,&ct,r,&one,n);
    bt=0;ct=0;mgcv_mmult(v1,P,v2,&bt,&ct,q,&one,r);
    /* result now in v1, add v1 to relevant part of B2z */
    dp1 = v1 + *q;
    for (dp=v1,dp0=pB2z;dp<dp1;dp++,dp0++) *dp0 += *dp;

    /* 12. obtain PK'TmKP'SkPK'z */
    pKPtSPKtz = KPtSPKtz + k * *n;
    dp = Tk + m * *n;dp1 = dp + *n;
    for (dp0=v1;dp < dp1;dp++,dp0++,pKPtSPKtz++) *dp0 = *dp * *pKPtSPKtz;
    bt=1;ct=0;mgcv_mmult(v2,K,v1,&bt,&ct,r,&one,n);
    bt=0;ct=0;mgcv_mmult(v1,P,v2,&bt,&ct,q,&one,r);
    /* result now in v1, add 2*sp[k]*v1 to relevant part of B2z */
    dp1 = v1 + *q;xx=2*sp[k];
    for (dp=v1,dp0=pB2z;dp<dp1;dp++,dp0++) *dp0 += xx * *dp;

  

    /* 13. PP'SkPK'TmKK'z */ 
    pPKtTKKtz = PKtTKKtz + m * *q;
    multSk(v1,pPKtTKKtz,&one,k,rS,rSncol,q,v2);    
    bt=1;ct=0;mgcv_mmult(v2,P,v1,&bt,&ct,r,&one,q);
    bt=0;ct=0;mgcv_mmult(v1,P,v2,&bt,&ct,q,&one,r);
    /* result now in v1, add 2*sp[k]*v1 to relevant part of B2z */
    dp1 = v1 + *q;xx=2*sp[k];
    for (dp=v1,dp0=pB2z;dp<dp1;dp++,dp0++) *dp0 += xx * *dp;

 
    /* 14. obtain PP'SkPK'Tmz */
    pPKtTz = PKtTz + m * *q;
    multSk(v1,pPKtTz,&one,k,rS,rSncol,q,v2);    
    bt=1;ct=0;mgcv_mmult(v2,P,v1,&bt,&ct,r,&one,q);
    bt=0;ct=0;mgcv_mmult(v1,P,v2,&bt,&ct,q,&one,r);
    /* result now in v1, add -sp[k]*v1 to relevant part of B2z */
    dp1 = v1 + *q;xx=sp[k];
    for (dp=v1,dp0=pB2z;dp<dp1;dp++,dp0++) *dp0 -= xx * *dp;

  
    /* update complete, move pB2z and pTkm to next case */
    pB2z += *q;
    pTkm += *n;
  } 
}


void getB1z1(double *B1z1,double *z1,double *K,double *P,double *Tk,double *sp,
             double *rS,int *rSncol,int *n,int *r, int *q,int *M,double *work)
/* function to form dB/d \rho_k dz'/\rho_m for k,m=1..M and return the result in
   B1z1 which has dimension q*M^2. The storage order is:
   dB/d \rho_0 dz'/\rho_0, dB/d \rho_0 dz'/\rho_1, ...
   All matrices are stored column-wise, as in R. Conceptual dimensions are:
   * K is n by r
   * P is q by r
   * elements of z1 and Tk are n-vectors.
   * work must be of length 2*(n+q)*M, at least.
   -- Used by Wood 2008. Now defunct.
   
*/
{ double *v1,*v2,*PKtz1,*KKtz1,*dp,*dp0,*dp1,*dp2,*pz1,*pv1,xx;
  int bt,ct,k,j;
  /* A. allocate work*/
  dp=work;
  v1 = dp; dp += *n * *M;
  v2 = dp; dp += *q * *M; /* note limited size!*/
  PKtz1 = dp; dp += *q * *M;
  KKtz1 = dp; dp += *n * *M; 
  
  /* B. obtain PK'z1 and KK'z1  */ 
  bt=1;ct=0;mgcv_mmult(v2,K,z1,&bt,&ct,r,M,n);
  bt=0;ct=0;mgcv_mmult(KKtz1,K,v2,&bt,&ct,n,M,r);
  bt=0;ct=0;mgcv_mmult(PKtz1,P,v2,&bt,&ct,q,M,r);
  
  /* C. loop through the smoothing parameters updating B1z1 */
  for (k=0;k<*M;k++) {
    /* PP'SkPK'z1 */
    multSk(v2,PKtz1,M,k,rS,rSncol,q,v1);
    bt=1;ct=0;mgcv_mmult(v1,P,v2,&bt,&ct,r,M,q);
    bt=0;ct=0;mgcv_mmult(v2,P,v1,&bt,&ct,q,M,r);
    /* v2 now contains result - need to add (-sp[k] times) to appropriate block of B1z1 */
    dp1=B1z1 + *q * *M;xx = -sp[k];
    for (pv1=v2,dp=B1z1;dp < dp1;dp++,pv1++) *dp = xx * *pv1;  
    
    /* PK'Tkz1 */
    dp0 = Tk + k * *n;dp1 = dp0 + *n;
    for (pz1=z1,pv1=v1,j=0;j<*M;j++) 
     for (dp=dp0;dp<dp1;dp++,pz1++,pv1++) *pv1 = *dp * *pz1;  
    bt=1;ct=0;mgcv_mmult(v2,K,v1,&bt,&ct,r,M,n);
    bt=0;ct=0;mgcv_mmult(v1,P,v2,&bt,&ct,q,M,r);
    /* v1 now contains result - need to add to appropriate block of B1z1 */
    dp1=B1z1 + *q * *M;
    for (pv1=v1,dp=B1z1;dp < dp1;dp++,pv1++) *dp += *pv1;

    /* PK'TkKK'z1 */
    dp0 = Tk + k * *n;dp1 = dp0 + *n;
    for (pv1=v1,dp2=KKtz1,j=0;j<*M;j++) 
     for (dp=dp0;dp<dp1;dp++,dp2++,pv1++) *pv1 = *dp * *dp2;
    bt=1;ct=0;mgcv_mmult(v2,K,v1,&bt,&ct,r,M,n);
    bt=0;ct=0;mgcv_mmult(v1,P,v2,&bt,&ct,q,M,r);
    /* v1 now contains result - need to add (-2 times) to appropriate block of B1z1 */
    dp1=B1z1 + *q * *M;
    for (pv1=v1,dp=B1z1;dp < dp1;dp++,pv1++) *dp += -2 * *pv1; 
        
    /* move B1z1 on to next block */
    B1z1 += *q * *M;
  }
}

void rc_prod(double *y,double *z,double *x,int *xcol,int *n)
/* obtains element-wise product of z and each of the  xcol columns of x,
   returning the result in y. z and the columns of x are n-vectors. 
   (equivalent to y = diag(z)%*%x)
*/ 
{ int i;
  double *pz,*pz1;
  pz1 = z + *n;
  for (i=0;i < *xcol;i++) 
      for (pz=z;pz<pz1;pz++,y++,x++) *y = *pz * *x;
}


void Rinv(double *Ri,double *R,int *c,int *r, int *ri)
/* invert c by c upper triangular matrix R, actually stored in upper 
   part of r by c matrix. Result returned in top of  Ri (actually ri by c).
*/

{ int i,j,k,eye;
  double xx,*rc;
  rc=Ri;
  for (i=0;i<*c;i++) {
      for (eye=1,k=i;k>=0;k--) {
	  for (xx=0.0,j=k+1;j < *c;j++) xx += R[k + j * *r] * rc[j];
          rc[k]=(eye-xx)/R[k + k * *r];
          eye=0;
      }
      for (k=i+1;k<*c;k++) rc[k]=0.0;
      rc += *ri;
  }
}


void pearson2(double *P, double *P1, double *P2,
              double *y,double *mu,double *V, double *V1,double *V2,double *g1,double *g2,
              double *p_weights,double *eta1, double *eta2,int n,int M,int deriv, int deriv2)
/* Alternative calculation of the derivatives of the Pearson statistic, which avoids assuming that
   z and w are based on Fisher scoring */
{ double resid,xx,*Pe1,*Pe2,*pp,*p1,*p0,*v2,*Pi1,*Pi2;
  int i,k,m,n_2dCols=0,one=1;
  if (deriv) {
    Pe1 = (double *)calloc((size_t)n,sizeof(double)); /* for dP/deta */
    Pi1 = (double *)calloc((size_t) n * M,sizeof(double)); /* for dPi/drho */
    if (deriv2) { 
      n_2dCols = (M * (1 + M))/2;
      Pe2 = (double *)calloc((size_t)n,sizeof(double)); /* for d2P/deta2 */
      v2 = (double *)calloc((size_t)n,sizeof(double));
      Pi2 = (double *)calloc((size_t)n_2dCols*n,sizeof(double)); /* for d2P_i/drho */
    }
  }
  *P=0.0;
  for (i=0; i < n;i++) {
    resid = y[i]-mu[i];
    xx = resid*p_weights[i]/V[i];
    *P += xx*resid;
    if (deriv) {
      Pe1[i] = - xx* (2 + resid*V1[i])/g1[i];
      if (deriv2) {
        Pe2[i] = - Pe1[i]*g2[i]/g1[i] + 
	  (2*p_weights[i]/V[i]+2*xx*V1[i] - Pe1[i]*V1[i]*g1[i] - xx*resid*(V2[i]-V1[i]*V1[i]))/(g1[i]*g1[i]);
      }
    }
  } /* derivs wrt eta completed */

  if (deriv) { /* transform to derivs wrt rho */
    rc_prod(Pi1,Pe1,eta1,&M,&n); /* Pi1 = dP_i/drho_k done */
      if (deriv2) {  
        rc_prod(Pi2,Pe1,eta2,&n_2dCols,&n);
        for (pp=Pi2,m=0;m < M;m++) for (k=m;k < M;k++) {
	    rc_prod(Pe1,eta1 + n *  m,eta1 + n * k,&one,&n);
            rc_prod(v2,Pe2,Pe1,&one,&n);
            p1=v2 + n;
            for (p0=v2;p0<p1;p0++,pp++) *pp += *p0;        
        } /* Pi2 update completed */
      }
  } /* derivatives of Pi wrt rho completed */

  /* now sum the derivatives over i */

  if (deriv) {
    pp = Pi1;
    for (k=0;k<M;k++) { xx=0.0; for (i=0;i<n;i++,pp++) xx += *pp;P1[k] = xx;}
    if (deriv2) {
        for (pp=Pi2,m=0;m < M;m++) for (k=m;k < M;k++) {
          xx=0.0;
          for (i=0;i<n;i++,pp++) xx += *pp;
          P2[k*M+m] = P2[m*M+k] = xx;
        } 
    }
  } /* end of derivative summation */
  
  /* clear up */
  if (deriv) {
    free(Pe1);free(Pi1);
    if (deriv2) {
      free(Pe2);free(Pi2);free(v2);
    }
  }

} 


void pearson(double *w, double *w1,double *w2,double *z,double *z1, double *z2,
             double *eta,double *eta1,double *eta2,double *P, double *P1, double *P2,
             double *work,int n,int M,int deriv, int deriv2)
/* Function to evaluate the pearson statistic sum_i [w_i(z_i-eta_i)]^2
   and its derivstives wrt the log smoothing parameters. Arrays ending 
   1 or 2 contain 1st or second derivatives of their base name quantity 
   wrt log smoothing parameters. n is length of z, w, eta. M is number 
   of smoothing parameters. Note: Only works with Fisher scoring.
   Is code implementing Wood 2008.
*/
{ double *zeta,*wzeta,*p0,*p1,*p2,*p3,*p4,*zetaSq,*wSqzeta,*wzetaSq,*wSqzetaSq,xx;
  int i,bt,ct,one=1,m,k;
  zeta = work;work +=n;
  wzeta = work;work +=n;
  zetaSq = work;work +=n;
  wSqzeta = work;work +=n;
  wzetaSq = work;work +=n;
  wSqzetaSq = work; work +=n;
  for (p0=zeta,p1=zeta+n,p2=z,p3=eta;p0<p1;p0++,p2++,p3++,zetaSq++) 
  { *p0 = *p2 - *p3; *zetaSq = *p0 * *p0;} /* get z - eta */
   zetaSq -= n;  
  for (*P=0.0,p0=wzeta,p1=wzeta+n,p2=zeta,p3=w;p0<p1;p0++,p2++,p3++,wSqzeta++,wzetaSq++,zetaSq++) 
  { *p0 = *p2 * *p3; /* wzeta= w_i(z_i - eta_i) */
    *P += *p0 * *p0; /*Pearson statistic */ 
    *wSqzeta = *p0 * *p3; /* w_i^2(z_i - eta_i) */
    *wzetaSq = *p3 * *zetaSq; /* w_i (z_i-eta_i)^2 */
  } 
  wSqzeta -= n;zetaSq -=n;wzetaSq -= n;
  if (!deriv) return; /* no derivatives required */
  if (deriv2) {
    for (p0=w,p1=w+n;p0<p1;p0++,wzetaSq++,wSqzetaSq++)  
      *wSqzetaSq = *p0 * *wzetaSq;
    wSqzetaSq -= n;wzetaSq -= n;
  }   
  /* do first derivatives */
  bt=1;ct=0;mgcv_mmult(P1,wzetaSq,w1,&bt,&ct,&one,&M,&n);
  bt=1;ct=0;mgcv_mmult(work,wSqzeta,z1,&bt,&ct,&one,&M,&n);
  for (i=0;i<M;i++) P1[i] += work[i];
  bt=1;ct=0;mgcv_mmult(work,wSqzeta,eta1,&bt,&ct,&one,&M,&n);
  for (i=0;i<M;i++) { P1[i] += -work[i];P1[i]*=2;}
  if (deriv2)  /* get the second derivatives */
  for (m=0;m < M;m++) for (k=m;k < M;k++) {
      xx=0.0;
      for (i=0;i<n;i++,w2++,z2++,eta2++) /* terms involving second derivatives */
	  xx += *w2 * wzetaSq[i] + wSqzeta[i] * (*z2 - *eta2);
      p2=w1+m*n;p3=w1+k*n;
      for (p0=zetaSq,p1=zetaSq+n;p0<p1;p0++,p2++,p3++) xx += *p0 * *p2 * *p3;
      p2=w1+m*n;p3=z1+k*n;p4=eta1+k*n;
      for (p0=wzeta,p1=wzeta+n;p0<p1;p0++,p2++,p3++,p4++) xx+= 2 * *p0 * *p2 * (*p3 - *p4);
      p2=w1+k*n;p3=z1+m*n;p4=eta1+m*n;
      for (p0=wzeta,p1=wzeta+n;p0<p1;p0++,p2++,p3++,p4++) xx+= 2 * *p0 * *p2 * (*p3 - *p4);
      p0=z1+m*n;p1=eta1+m*n;p2=z1+k*n;p3=eta1+k*n;p4=p3+n;
      for (;p3<p4;w++,p0++,p1++,p2++,p3++) xx += *w * *w * (*p0 - *p1) * (*p2 - *p3);
      w-=n;
      P2[m*M+k]=P2[k*M+m]= 2*xx;
  }
}


void ift(double *X,double *P,double *rS,double *beta,double *sp,double *w,
         double *dwdeta,double *b1, double *b2,double *eta1,double *eta2,
         int *n,int *q,int *r, int *M,int *rSncol,int *deriv2)

/* Uses the implicit function theorem to get derivatives of beta wrt rho = log(sp) cheaply
   without iteration, when Newton or Fisher-canonical-link are used... 
   X is n by q
   P is q by r
   there are M smoothing parameters (unlogged) in sp
   beta is a q vector
   b1 is q by M
   b2 is q by n_2dCols 
*/
{ int n_2dCols,i,j,k,one=1,bt,ct;
  double *work,*Skb,*pp,*p0,*p1,*work1;
  work = (double *) calloc((size_t)*n,sizeof(double));
  work1 = (double *) calloc((size_t)*n,sizeof(double));
  Skb = (double *) calloc((size_t)*q,sizeof(double));
  n_2dCols = (*M * (1 + *M))/2;
  for (i=0;i<*M;i++) { /* first derivative loop */
    multSk(Skb,beta,&one,i,rS,rSncol,q,work); /* get S_i \beta */
    for (j=0;j<*q;j++) Skb[j] *= -sp[i]; 
    bt=1;ct=0;mgcv_mmult(work,P,Skb,&bt,&ct,r,&one,q); 
    bt=0;ct=0;mgcv_mmult(b1 + i * *q,P,work,&bt,&ct,q,&one,r); /* so ith col of b1 contains -2sp[i]PP'S_ib */ 
  } /* first derivatives of beta finished */

  bt=0;ct=0;mgcv_mmult(eta1,X,b1,&bt,&ct,n,M,q); /* first deriv of eta */

  if (*deriv2) { /* then second derivatives needed */
    pp = b2;   
    for (i=0;i<*M;i++) for (k=i;k<*M;k++) { 
      p0 = eta1 + *n * i;p1 = eta1 + *n * k;
      /*  for (j=0;j<*n;j++,p0++,p1++) work[j] = - *p0 * *p1 * 2 * w[j] * dwdeta[j];
          .... this version only for use when weights are square root weights,
          following is used instead
      */
      for (j=0;j<*n;j++,p0++,p1++) work[j] = - *p0 * *p1 * dwdeta[j];
      bt=1;ct=0;mgcv_mmult(Skb,X,work,&bt,&ct,q,&one,n); /* X'f */
      multSk(work,b1+k* *q,&one,i,rS,rSncol,q,work1); /* get S_i dbeta/drho_k */
      for (j=0;j<*q;j++) Skb[j] += -sp[i]*work[j];
      multSk(work,b1+i* *q,&one,k,rS,rSncol,q,work1); /* get S_k dbeta/drho_i */
      for (j=0;j<*q;j++) Skb[j] += -sp[k]*work[j];
      bt=1;ct=0;mgcv_mmult(work,P,Skb,&bt,&ct,r,&one,q); 
      bt=0;ct=0;mgcv_mmult(pp,P,work,&bt,&ct,q,&one,r);
      if (i==k) for (j=0;j< *q;j++) pp[j] += b1[i * *q + j];
      pp += *q;
    }

    bt=0;ct=0;mgcv_mmult(eta2,X,b2,&bt,&ct,n,&n_2dCols,q); /* second derivatives of eta */
  }

  free(work);free(Skb);free(work1);
}


double MLpenalty(double *det1,double *det2,double *Tk,double *Tkm,double *U1, double *R,double *Q, int *nind,double *sp,
                 double *rS,int *rSncol,int *q,int *n,int *nn,int *Ms,int *M,int *neg_w,double *rank_tol,int *deriv) {
/* Routine to obtain the version of log|X'WX+S| that applies to ML, rather than REML.
   * U1 is basis for penalty range space
   * Q, R are the QR factors of diag(abs(W))X augmenented by the square root of S
   * nind is the array indexing the locations of the `neg_w' -ve elements of W.
   * q is the number of model coefficients
   * Ms is the penalty null space dimension.
   * nn is the number of rows in Q. 

   Basic task of the routine is to project Hessian of the penalized log likelihood 
   into the range space of the penalty, in order to obtain the correction term that 
   applies for ML. 
*/

  double *RU1,*tau,*work,*Ri,*Qb,Rcond,*K,*P,*IQ,*IQQ,*Vt,
         *d,*p0,*p1,*pd,*p2,*p3,ldetXWXS,ldetI2D;
  int bt,ct,qM,*pivot,rank,ncol,i,j,k,left,tp,*pi,*pi1; 

  qM = *q - *Ms;

  RU1 = (double *)calloc((size_t) *q * qM,sizeof(double));
  bt=0;ct=0;mgcv_mmult(RU1,R,U1,&bt,&ct,q,&qM,q);
 
  /* A pivoted QR decomposition of RU1 is needed next */
  tau=(double *)calloc((size_t)qM,sizeof(double)); /* part of reflector storage */
  pivot=(int *)calloc((size_t)qM,sizeof(int));
  
  mgcv_qr(RU1,q,&qM,pivot,tau); /* RU1 and tau now contain the QR decomposition information */
  /* pivot[i] gives the unpivoted position of the ith pivoted parameter.*/
  
  /* Now obtain rank of new factor R */
  
  work = (double *)calloc((size_t)(4 * qM),sizeof(double));
  rank = qM;
  R_cond(RU1,q,&rank,work,&Rcond);
  while (*rank_tol * Rcond > 1) { rank--;R_cond(RU1,q,&rank,work,&Rcond);}
  free(work);

  /* Ri needed (possibly truncated to rank by rank) */

  Ri =  (double *)calloc((size_t) rank * rank,sizeof(double)); 
  Rinv(Ri,RU1,&rank,q,&rank); /* getting R^{-1} */
  
  /* new Q factor needed explicitly */

  Qb = (double *)calloc((size_t) *q * rank,sizeof(double)); 
  for (i=0;i< rank;i++) Qb[i * *q + i] = 1.0;
  left=1;tp=0;mgcv_qrqy(Qb,RU1,tau,q,&rank,&qM,&left,&tp); /* Q from the QR decomposition */

  free(tau);

  K = (double *)calloc((size_t) *n * rank,sizeof(double));
  P = (double *)calloc((size_t) qM * rank,sizeof(double));

  if (*neg_w) { /* need to deal with -ve weight correction */
    if (*neg_w < *q+1) k = *q+1; else k = *neg_w;
    IQ = (double *)calloc((size_t) k * *q,sizeof(double)); 
    for (i=0;i< *neg_w;i++) { /* Copy the rows of Q corresponding to -ve w_i into IQ */
      p0 = IQ + i;p1 = Q + nind[i];
      for (j=0;j<*q;j++,p0+=k,p1+= *nn) *p0 = *p1;
    }
    /* Note that IQ may be zero padded, for convenience */
    IQQ = (double *)calloc((size_t) k * rank,sizeof(double)); 
    bt=0;ct=0;mgcv_mmult(IQQ,IQ,Qb,&bt,&ct,&k,&rank,q); /* I^-Q_1 \bar Q is k by rank */
    free(IQ);
     
    /* Get the SVD of IQQ */
    Vt = (double *)calloc((size_t) rank * rank,sizeof(double));
    d = (double *)calloc((size_t)  rank,sizeof(double));
    mgcv_svd_full(IQQ,Vt,d,&k,&rank); /* SVD of IQ */
    free(IQQ);
    for (i=0;i<rank;i++) {
      d[i] = 1 - 2*d[i]*d[i];
      if (d[i]<=0) d[i]=0.0; 
      else {
        ldetI2D += log(d[i]); /* log|I-2D^2| */ 
        d[i] = 1/sqrt(d[i]);
      }
    } /* d now contains diagonal of diagonal matrix (I-2D^2)^{-1/2} (possibly pseudoinverse) */
    /* Now form (I-2D^2)^.5 Vt and store in Vt... */
    for (p0=Vt,i=0;i<rank;i++)
    for (p1=d,p2=d+rank;p1<p2;p1++,p0++) *p0 *= *p1;
    
    /* Form K */
    IQ = (double *)calloc((size_t) *n * *q,sizeof(double));
    for (p0=IQ,p1=Q,j=0;j<*q;j++,p1 += *nn) /* copy just Q1 into IQ */
      for (p2 = p1,p3=p1 + *n;p2<p3;p0++,p2++) *p0 = *p2; 
    
    work = (double *)calloc((size_t)*q * rank,sizeof(double));
    bt=0;ct=1;mgcv_mmult(work,Qb,Vt,&bt,&ct,q,&rank,&rank); /* \bar Q V (I - 2D^2)^.5 */

    bt=0;ct=0;mgcv_mmult(K,IQ,work,&bt,&ct,n,&rank,q);
    free(work);
    
    /* Form P */
    bt=0;ct=1;mgcv_mmult(IQ,Ri,Vt,&bt,&ct,&rank,&rank,&rank);
    for (p0=P,p1=IQ,j=0;j<rank;j++,p0+= qM) /* copy R^{-1}V'(I-2D)^{-.5} into first rows of P */
      for (p2=p0,p3 = p2 + rank;p2<p3;p1++,p2++) *p2 = *p1; 
    free(IQ);free(d);free(Vt);   
    
  } else { /* no negative weights, so P and K can be obtained directly */
    ldetI2D = 0.0;
    /* Form K */
    IQ = (double *)calloc((size_t) *n * *q,sizeof(double));
    for (p0=IQ,p1=Q,j=0;j<*q;j++,p1 += *nn) /* copy just Q1 into IQ */
      for (p2 = p1,p3=p1 + *n;p2<p3;p0++,p2++) *p0 = *p2; 
    bt=0;ct=0;mgcv_mmult(K,IQ,Qb,&bt,&ct,n,&rank,q);
    /* Form P */
    for (p0=P,p1=Ri,j=0;j<rank;j++,p0+= qM) /* copy R^{-1} into first rows of P */
    for (p2=p0,p3=p0 + rank;p2<p3;p1++,p2++) *p2 = *p1; 
    free(IQ);
  }

  free(Ri);

  /* Evaluate the required log determinant... */

  for (ldetXWXS=0.0,i=0;i<rank;i++) ldetXWXS += log(fabs(RU1[i + i * *q])); 
  ldetXWXS *= 2;
 
  ldetXWXS += ldetI2D; /* the negative weights correction */
  
  free(RU1);

  /* rS also needs to be transformed and pivoted... */

  ncol = 0;for (i=0;i<*M;i++) ncol += rSncol[i];
  work = (double *)calloc((size_t) qM * ncol,sizeof(double));
  bt=1;ct=0;mgcv_mmult(work,U1,rS,&bt,&ct,&qM,&ncol,q);

  for (p2=rS,pd=work,i=0;i<ncol;i++,pd += qM) /* work across columns of work */
  { for (pi=pivot,pi1=pivot+qM;pi<pi1;pi++,p2++) *p2 = pd[*pi];  /* apply pivot into rS */
  } 
  
  free(work);free(Qb);free(pivot);

  /* Now we have all the ingredients to obtain required derivatives of the log determinant... */
  
  if (*deriv)
    get_ddetXWXpS(det1,det2,P,K,sp,rS,rSncol,Tk,Tkm,n,&qM,&rank,M,deriv);

  free(P);free(K);
  return(ldetXWXS);
}

void gdi(double *X,double *E,double *rS,double *UrS,double *U1,
    double *sp,double *z,double *w,double *mu,double *eta, double *y,
	 double *p_weights,double *g1,double *g2,double *g3,double *g4,double *V0,
	 double *V1,double *V2,double *V3,double *beta,double *D1,double *D2,
    double *P0, double *P1,double *P2,double *trA,
    double *trA1,double *trA2,double *rV,double *rank_tol,double *conv_tol, int *rank_est,
	 int *n,int *q, int *M,int *Mp,int *Encol,int *rSncol,int *deriv,int *use_svd,
	 int *REML,int *fisher,int *fixed_penalty)     
/* Version of gdi, based on derivative ratios and Implicit Function Theorem 
   calculation of the derivatives of beta. Assumption is that Fisher is only used 
   with canonical link, when it is equivalent to Newton anyway.

   This version deals properly with negative weights, which can occur with Newton based 
   PIRLS. In consequence w's in this routine are proportional to reciprocal variances,
   not reciprocal standard deviations.
  
   The function is to be called at convergence of a P-IRLS scheme so that 
   z, w, mu and functions of these can be treated as fixed, and only the 
   derivatives need to be updated.

   All matrices are packed into arrays in column order (i.e. col1, col2,...)
   as in R. 

   All names ending in 1,2 or 3 are derivatives of some sort, with the integer
   indicating the order of differentiation. 

   The arguments of this function point to the following:
   * X is and n by q model matrix.
   * E is a q by Encol square root of the total penalty matrix, so EE'=S
   * rS is a list of square roots of individual penalty matrices, packed
     in one array. The ith such matrix rSi, say, has dimension q by rSncol[i]
     and the ith penalty is [rSi][rSi]'.
   * U1 is an (orthogonal) basis for the penalty range space (q by (q-Mp), where Mp
     is the null space dimension).
   * sp is an M array of smoothing parameters (NOT log smoothing parameters)
   * z and w are n-vectors of the pseudodata and iterative weights
   * p_weights is an n-vector of prior weights (as opposed to the iterative weights in w)
   * mu and y are n-vectors of the fitted values and data.
   * g1,g2,g3,g4 are the n-vectors of the link derivatives 
     Note that g''(mu) g'''(mu) and g''''(mu) are *divided by* g'(mu)
   * V0, V1, V2, V3 are n-vectors of the variance function and first three derivatives,
     Note that V'(mu), V''(mu) & V'''(mu) are divided by V(mu)
   * D1 and D2 are an M-vector and M by M matrix for returning the first 
     and second derivatives of the deviance wrt the log smoothing parameters.
     if *REML is non zero then the derivs will be of the penalized deviance,
     and b'Sb will be returned in conv_tol  
   * trA1 and trA2 are an M-vector and M by M matrix for returning the first 
     and second derivatives of tr(A) wrt the log smoothing parameters.
     If *REML is non zero then the derivatives of the REML penalty are 
     returned instead (with the REML penalty returned in `rank_tol', hack, hack).
   * P0,P1,P2 are for returning the Pearson statistic and its derivatives, or 
     the Pearson scale estimate and derivatives if *REML is non - zero. 
   * rank_est is for returning the estimated rank of the problem.
   * the remaining arguments are the dimensions already refered to except for:
   * deriv, which controls which derivatives are produced:
       deriv==0 for no derivatives: only trA, rV and beta returned
       deriv==1 for first derivatives only
       deriv==2 for gradient and Hessian
     -- on exit contains the number of iteration steps required.   

    * If REML is +ve non-zero, then the REML penalty returned in rank_tol, with it's 
      derivatives in trA1, trA2: it is to be added to the *deviance* to get D_r.
    * If REML is -ve non-zero, then the ML penalty is returned in place of the REML one.
    * non-zero `fisher' indicates that Fisher scoring, rather than full Newton,
      is the basis for iteration. 
    * non-zero `fixed_penalty' inticates that S includes a fixed penalty component,
      the range space projected square root of which is in the final element of `UrS'.
      This information is used by get_detS2().

   The method has 4 main parts:

   1. The initial QR- decomposition and SVD are performed, various quantities which 
      are independent of derivatives are created

   2. IFT used to obtain derivatives of the coefficients wrt the log smoothing 
      parameters. 

   3. Evaluation of the derivatives of the deviance wrt the log smoothing parameters
      (i.e. calculation of D1 and D2)

   4. Evaluation of the derivatives of tr(A) (i.e. trA1 and trA2)

   The method involves first and second derivatives of a number of k-vectors wrt
   log smoothing parameters (\rho), where k is q or n. Consider such a vector, v. 
   
   * v1 will contain dv/d\rho_0, dv/d\rho_1 etc. So, for example, dv_i/d\rho_j
     (indices starting at zero) is located in v1[q*j+i].
   
   * v2 will contain d^2v/d\rho_0d\rho_0, d^2v/d\rho_1d\rho_0,... but rows will not be
     stored if they duplicate an existing row (e.g. d^2v/d\rho_0d\rho_1 would not be 
     stored as it already exists and can be accessed by interchanging the sp indices).
     So to get d^2v_k/d\rho_id\rho_j: 
     i)   if i<j interchange the indices
     ii)  off = (j*m-(j+1)*j/2+i)*q 
     iii) v2[off+k] is the required derivative.       

    

*/
{ double *zz,*WX,*tau,*work,*pd,*p0,*p1,*p2,*p3,*p4,*K=NULL,
         *Ri,*d,*Vt,xx,*b1,*b2,*P,
         *c0,*c1,*c2,*a1,*a2,*eta1,*eta2,
         *PKtz,*v1,*v2,*wi,*w1,*w2,*pw2,*Tk,*Tkm,
         *pb2, *dev_grad,*dev_hess=NULL,Rcond,
    ldetXWXS=0.0,reml_penalty=0.0,bSb=0.0,*R,
    *alpha,*alpha1,*alpha2,*raw,*Q1,*IQ, *U, d_tol;
  int i,j,k,*pivot,ScS,*pi,rank,left,tp,bt,ct,iter=0,m,one=1,n_2dCols,n_b1,n_b2,
    n_eta1,n_eta2,n_work,deriv2,null_space_dim,neg_w=0,*nind,nn,ii,ldetI2D;

  if (*deriv==2) deriv2=1; else deriv2=0;

  d_tol = sqrt(*rank_tol * 100);
  /* first step is to obtain P and K */
  nn= *n + *Encol;
  zz = (double *)calloc((size_t)nn,sizeof(double)); /* storage for z=[sqrt(|W|)z,0] */
  raw = (double *)calloc((size_t) *n,sizeof(double)); /* storage for sqrt(|w|) */
  
  for (i=0;i< *n;i++) 
    if (w[i]<0) { neg_w++;raw[i] = sqrt(-w[i]);} 
    else raw[i] = sqrt(w[i]);

  if (neg_w) {
    nind = (int *)calloc((size_t)neg_w,sizeof(int)); /* index the negative w_i */
    k=0;for (i=0;i< *n;i++) if (w[i]<0) { nind[k]=i;k++;}
  } else { nind = (int *)NULL;}

  for (i=0;i< *n;i++) zz[i] = z[i]*raw[i]; /* form z itself*/

  for (i=0;i<neg_w;i++) { k=nind[i];zz[k] = -zz[k];} 

  WX = (double *) calloc((size_t) ( nn * *q),sizeof(double));
  for (j=0;j<*q;j++) 
  { for (i=0;i<*n;i++) /* form WX */
    { k = i + nn * j;
      WX[k]=raw[i]*X[i + *n *j];
    }
    for (ii=0,i = *n;ii<*Encol;i++,ii++) /* append E' */ 
    { k = i + nn * j;
      WX[k] = E[j + *q * ii];
    }
  } 
  /* get the QR decomposition of WX */
  tau=(double *)calloc((size_t)*q,sizeof(double)); /* part of reflector storage */
 
  pivot=(int *)calloc((size_t)*q,sizeof(int));
  
  mgcv_qr(WX,&nn,q,pivot,tau); /* WX and tau now contain the QR decomposition information */
  /* pivot[i] gives the unpivoted position of the ith pivoted parameter.*/
  
  /* first find the rank of R */
  work = (double *)calloc((size_t)(4 * *q),sizeof(double));
  rank = *q;
  R_cond(WX,&nn,&rank,work,&Rcond);
  while (*rank_tol * Rcond > 1) { rank--;R_cond(WX,&nn,&rank,work,&Rcond);}
  free(work);

  /* printf("q = %d  rank = %d\n",*q,rank);  DEBUG */

  /* Note that in the following the Q is extracted with q rather than rank columns
     this is because ML estimation requires the `full' rather than truncated version.
     Q1 can still be treated as nn by rank, of course.
  */
  Q1 = (double *)calloc((size_t) nn * *q,sizeof(double)); 
  for (i=0;i< *q;i++) Q1[i * nn + i] = 1.0;
  left=1;tp=0;mgcv_qrqy(Q1,WX,tau,&nn,q,q,&left,&tp); /* Q from the QR decomposition */

  Ri =  (double *)calloc((size_t) rank * rank,sizeof(double)); 
  Rinv(Ri,WX,&rank,&nn,&rank); /* getting R^{-1} */
  
  K = (double *)calloc((size_t) *n * rank,sizeof(double));
  P = (double *)calloc((size_t) *q * rank,sizeof(double));

  ldetI2D = 0.0; /* REML determinant correction */

  if (neg_w) { /* then the correction for the negative w_i has to be evaluated */
    if (neg_w < rank+1) k = rank+1; else k = neg_w;
    IQ = (double *)calloc((size_t) k * rank,sizeof(double)); 
    for (i=0;i<neg_w;i++) { /* Copy the rows of Q corresponding to -ve w_i into IQ */
      p0 = IQ + i;p1 = Q1 + nind[i];
      for (j=0;j<rank;j++,p0+=k,p1+= nn) *p0 = *p1;
    }
    /* Note that IQ may be zero padded, for convenience */
    Vt = (double *)calloc((size_t) rank * rank,sizeof(double));
    d = (double *)calloc((size_t)  rank,sizeof(double));
    mgcv_svd_full(IQ,Vt,d,&k,&rank); /* SVD of IQ */
    free(IQ);
    for (i=0;i<rank;i++) {
      d[i] = 1 - 2*d[i]*d[i];
      if (d[i]<=0) d[i]=0.0; 
      else {
        ldetI2D += log(d[i]); /* log|I-2D^2| */ 
        d[i] = 1/sqrt(d[i]);
      }
    } /* d now contains diagonal of diagonal matrix (I-2D^2)^{-1/2} (possibly pseudoinverse) */
    /* Now form (I-2D^2)^.5 Vt and store in Vt... */
    for (p0=Vt,i=0;i<rank;i++)
    for (p1=d,p2=d+rank;p1<p2;p1++,p0++) *p0 *= *p1;

    /* Form K */
    IQ = (double *)calloc((size_t) *n * rank,sizeof(double));
    for (p0=IQ,p1=Q1,j=0;j<rank;j++,p1 += nn) /* copy just Q1 into IQ */
      for (p2 = p1,p3=p1 + *n;p2<p3;p0++,p2++) *p0 = *p2; 
    bt=0;ct=1;mgcv_mmult(K,IQ,Vt,&bt,&ct,n,&rank,&rank);
    /* Form P */
    bt=0;ct=1;mgcv_mmult(IQ,Ri,Vt,&bt,&ct,&rank,&rank,&rank);
    for (p0=P,p1=IQ,j=0;j<rank;j++,p0+= *q) /* copy R^{-1}V'(I-2D)^{-.5} into first rows of P */
      for (p2=p0,p3 = p2 + rank;p2<p3;p1++,p2++) *p2 = *p1; 
    free(IQ);free(d);free(Vt);   
  } else { /* no negative weights so P and K much simpler */
    /* Form K */
    for (p0=K,p1=Q1,j=0;j<rank;j++,p1 += nn) /* copy just Q1 into K */
    for (p2 = p1,p3=p1 + *n;p2<p3;p0++,p2++) *p0 = *p2; 
    /* Form P */
    for (p0=P,p1=Ri,j=0;j<rank;j++,p0+= *q) /* copy R^{-1} into first rows of P */
    for (p2=p0,p3=p0 + rank;p2<p3;p1++,p2++) *p2 = *p1; 
  }
  
  /* At this stage P and K are complete */

   if (*REML>0) {  
      for (ldetXWXS=0.0,i=0;i<rank;i++) ldetXWXS += log(fabs(WX[i + i * nn])); 
      ldetXWXS *= 2;
      ldetXWXS += ldetI2D; /* correction for negative weights */
    }


  /* Apply pivoting to the parameter space - this simply means reordering the rows of Sr and the 
     rS_i, and then unscrambling the parameter vector at the end (along with any covariance matrix)
     pivot[i] gives the unpivoted position of the ith pivoted parameter.
  */

  ScS=0;for (pi=rSncol;pi<rSncol + *M;pi++) ScS+= *pi;  /* total columns of input rS */
  n_work = (4 * *n + 2 * *q) * *M + 2 * *n;
  k = (*M * (1 + *M))/2 * *n;
  if (n_work < k) n_work = k;
  work = (double *)calloc((size_t) n_work,sizeof(double)); /* work space for several routines*/
  p0 = work + *q;
  for (pd=rS,i=0;i<ScS;i++,pd += *q) /* work across columns */
  { for (pi=pivot,p2=work;p2<p0;pi++,p2++) *p2 = pd[*pi];  /* apply pivot into work */
    p3 = pd + *q;
    for (p1=pd,p2=work;p1<p3;p1++,p2++) *p1 = *p2;  /* copy back into rS */
  } /* rS pivoting complete, do E .... */
  p0 = work + *q;
  for (pd=E,i=0;i< *Encol;i++,pd += *q) /* work across columns */
  { for (pi=pivot,p2=work;p2<p0;pi++,p2++) *p2 = pd[*pi];  /* apply pivot into work */
    p3 = pd + *q;
    for (p1=pd,p2=work;p1<p3;p1++,p2++) *p1 = *p2;  /* copy back into E */
  }

  /* Now pivot the penalty range space, U1, iff ML used */

  if (*REML<0) { /* so it's ML */
    p0 = work + *q;j = *q - *Mp;
    for (pd=U1,i=0;i<j;i++,pd += *q) /* work across columns */
    { for (pi=pivot,p2=work;p2<p0;pi++,p2++) *p2 = pd[*pi];  /* apply pivot into work */
      p3 = pd + *q;
      for (p1=pd,p2=work;p1<p3;p1++,p2++) *p1 = *p2;  /* copy back into U1 */
    }
  }

  /* pivot columns of X - can't think of a smart arsed way of doing this*/
  for (i=0;i< *n;i++) {
      for (j=0;j<*q;j++) work[j]=X[pivot[j] * *n + i];
      for (j=0;j<*q;j++) X[j * *n + i] = work[j];
  }
  
  PKtz = (double *)calloc((size_t) *q,sizeof(double)); /* PK'z --- the pivoted coefficients*/
  bt=1;ct=0;mgcv_mmult(work,K,zz,&bt,&ct,&rank,&one,n);
  bt=0;ct=0;mgcv_mmult(PKtz,P,work,&bt,&ct,q,&one,&rank);  


  /************************************************************************************/
  /* free some memory */                    
  /************************************************************************************/
  if (*REML<0) { /* ML is required (rather than REML), and need to save some stuff */
    R = (double *) calloc((size_t) *q * *q ,sizeof(double)); /* save just R (untruncated) */
    for (p0=R,p1=WX,i=0;i<*q;i++,p1+=nn,p0 += *q) 
      for (p2=p1,p3=p0,p4=p0+i;p3<=p4;p3++,p2++) *p3 = *p2;
    
  } else { free(Q1);free(nind); } /* needed later for ML calculation */

  free(raw);free(WX);free(tau);free(Ri);
 
  /************************************************************************************/
  /* The coefficient derivative setup starts here */
  /************************************************************************************/
  /* set up some storage first */
  if (*deriv) {
    n_2dCols = (*M * (1 + *M))/2;
    n_b2 = *q * n_2dCols;
    b2 = (double *)calloc((size_t)n_b2,sizeof(double)); /* 2nd derivs of beta */
   
    n_b1 = *q * *M;
    b1 = (double *)calloc((size_t)n_b1,sizeof(double)); /* 1st derivs of beta */
   
    n_eta1 = *n * *M;
    eta1 = (double *)calloc((size_t)n_eta1,sizeof(double));
    Tk = (double *)calloc((size_t)n_eta1,sizeof(double));
   
    w1 = (double *)calloc((size_t)n_eta1,sizeof(double));

    n_eta2 = *n * n_2dCols;
    eta2 = (double *)calloc((size_t)n_eta2,sizeof(double));
    Tkm = (double *)calloc((size_t)n_eta2,sizeof(double));
  
    w2 = (double *)calloc((size_t)n_eta2,sizeof(double));

 
    v1 = work;v2=work + *n * *M; /* a couple of working vectors */ 
   
    /* Set up constants involved updates (little work => leave readable!)*/
    c0=(double *)calloc((size_t)*n,sizeof(double));
    c1=(double *)calloc((size_t)*n,sizeof(double));  
    c2=(double *)calloc((size_t)*n,sizeof(double));
  
    a1=(double *)calloc((size_t)*n,sizeof(double));  
    a2=(double *)calloc((size_t)*n,sizeof(double));
    alpha=alpha1=alpha2 =(double *)NULL;
    if (*fisher) { /* Fisher scoring updates */
      for (i=0;i< *n;i++) c0[i]=y[i]-mu[i];
      /* c1 = (y-mu)*g''/g' = dz/deta */
      for (i=0;i<*n;i++) c1[i]=g2[i]*c0[i];
      /* d2z/deta2 = c2 = ((y-mu)*(g'''/g'-(g''/g')^2) - g''/g')/g' */
      for (i=0;i<*n;i++) c2[i]=(c0[i]*(g3[i]-g2[i]*g2[i])-g2[i])/g1[i];

      /* set up constants involved in w updates */
      /* dw/deta = - w[i]*(V'/V+2g''/g')/g' */
      for (i=0;i< *n;i++) a1[i] = -  w[i] *(V1[i] + 2*g2[i])/g1[i];
     
      
      /* d2w/deta2 .... */
      for (i=0;i< *n;i++) 
        a2[i] = a1[i]*(a1[i]/w[i]-g2[i]/g1[i]) - w[i]*(V2[i]-V1[i]*V1[i] + 2*g3[i]-2*g2[i]*g2[i])/(g1[i]*g1[i]) ;

    } else { /* full Newton updates */
      
      alpha = (double *) calloc((size_t)*n,sizeof(double));
      alpha1 = (double *) calloc((size_t)*n,sizeof(double));
      alpha2 = (double *) calloc((size_t)*n,sizeof(double));
      for (i=0;i< *n;i++) {
        alpha[i] = 1 + (y[i]-mu[i])*(V1[i]+g2[i]);
        xx = V2[i]-V1[i]*V1[i]+g3[i]-g2[i]*g2[i]; /* temp. storage */
        alpha1[i] = (-(V1[i]+g2[i]) + (y[i]-mu[i])*xx)/alpha[i];
        alpha2[i] = (-2*xx + (y[i]-mu[i])*(V3[i]-3*V1[i]*V2[i]+2*V1[i]*V1[i]*V1[i]+g4[i]-3*g3[i]*g2[i]+2*g2[i]*g2[i]*g2[i]))/alpha[i];
      }
    
      /* end of preliminaries, now setup the multipliers that go forward */
      /* dz/deta... */
      for (i=0;i<*n;i++) c1[i] = 1 - 1/alpha[i] + (y[i]-mu[i])/alpha[i]*(g2[i]-alpha1[i]);
      /* d2z/deta2... */
      for (i=0;i<*n;i++) c2[i] = ((y[i]-mu[i])*(g3[i]-g2[i]*g2[i]-alpha2[i]+alpha1[i]*alpha1[i]) -
                                  (1+(y[i]-mu[i])*alpha1[i])*(g2[i]-alpha1[i]) + alpha1[i])/(alpha[i]*g1[i]);
     
      /* dw/deta ... */
      for (i=0;i<*n;i++) a1[i] = w[i]*(alpha1[i]-V1[i]-2*g2[i])/g1[i];
      /* d2w/deta2... */
      for (i=0;i<*n;i++) a2[i] = a1[i]*(a1[i]/w[i]-g2[i]/g1[i]) - 
                                 w[i]*(alpha1[i]*alpha1[i] - alpha2[i] + V2[i]-V1[i]*V1[i] + 2*g3[i]-2*g2[i]*g2[i])/(g1[i]*g1[i]) ;

      free(alpha);free(alpha1);free(alpha2);
      
    } /* end of full Newton setup */


    /* a useful array for Tk and Tkm */
    wi=(double *)calloc((size_t)*n,sizeof(double)); 
    for (i=0;i< *n;i++) { wi[i]=1/fabs(w[i]);}


    /* get gradient vector and Hessian of deviance wrt coefficients */
    for (i=0;i< *n ;i++) v1[i] = -2*p_weights[i]*(y[i]-mu[i])/(V0[i]*g1[i]);
    dev_grad=(double *)calloc((size_t)*q,sizeof(double));
    bt=1;ct=0;mgcv_mmult(dev_grad,X,v1,&bt,&ct,q,&one,n);
    
    if (deriv2) { /* get hessian of deviance w.r.t. beta */
      for (i=0;i< *n ;i++) v1[i] = 2*w[i];
      dev_hess=(double *)calloc((size_t)(*q * *q),sizeof(double));
      getXtWX(dev_hess,X,v1,n,q,v2);
    } 
  
  } /* end of if (*deriv) */ 
  else { /* keep compilers happy */
    b1=eta1=eta2=c0=c1=c2=(double *)NULL;
    a1=a2=wi=dev_grad=w1=w2=b2=(double *)NULL;
    Tk=Tkm=(double *)NULL;
  }
  /************************************************************************************/
  /* End of the coefficient derivative preparation  */
  /************************************************************************************/


  /************************************************************************************/
  /* Implicit Function Theorem code */
  /************************************************************************************/
  if (*deriv) {
    /* obtain derivatives of beta (b1,b2) and eta (eta1,eta2) using the IFT (a1 = dw/deta) */

    /* Note that PKtz used as pivoted version of beta, but not clear that PKtz really essential if IFT used */

   
    ift(X,P,rS,PKtz,sp,w,a1,b1,b2,eta1,eta2,n,q,&rank,M,rSncol,&deriv2);
  
  
    /* Now use IFT based derivatives to obtain derivatives of W and hence the T_* terms */

    /* get derivatives of w */  
    rc_prod(w1,a1,eta1,M,n); /* w1 = dw/d\rho_k done */
    if (deriv2) {
      rc_prod(w2,a1,eta2,&n_2dCols,n); 
      for (pw2=w2,m=0;m < *M;m++) for (k=m;k < *M;k++) {
        rc_prod(v1,eta1 + *n * m,eta1 + *n * k,&one,n);
        rc_prod(v2,a2,v1,&one,n);
        p1=v2 + *n;
        for (p0=v2;p0<p1;p0++,pw2++) *pw2 += *p0;           
      } /* w2 completed */
    }
    /* get Tk and Tkm */
      
    rc_prod(Tk,wi,w1,M,n); 
    if (deriv2) rc_prod(Tkm,wi,w2,&n_2dCols,n);
    
    /* evaluate gradient and Hessian of deviance */

    bt=1;ct=0;mgcv_mmult(D1,b1,dev_grad,&bt,&ct,M,&one,q); /* gradient of deviance is complete */
      
    if (deriv2) {       
      getXtMX(D2,b1,dev_hess,q,M,v1);
          
      for (pb2=b2,m=0;m < *M;m++) for (k=m;k < *M;k++) { /* double sp loop */
          p1 = dev_grad + *q;  
          for (xx=0.0,p0=dev_grad;p0<p1;p0++,pb2++) xx += *p0 * *pb2;
          D2[k * *M + m] += xx;
          D2[m * *M + k] = D2[k * *M + m];
      } /* Hessian of Deviance is complete !! */
    }

  } /* end of if (*deriv) */

  /* END of IFT */


  /* REML NOTE: \beta'S\beta stuff has to be done here on pivoted versions.
     Store bSb in P0, bSb1 in P1 and bSb2 in P2.
  */
  if (*REML) {
    get_bSb(&bSb,trA1,trA2,sp,E,rS,rSncol,Encol,q,M,PKtz,b1,b2,deriv);
    if (*deriv) for (p2=D2,p1=trA2,i = 0; i< *M;i++) { /* penalized deviance derivs needed */
        D1[i] += trA1[i];
        if (deriv2) for (j=0;j<*M;j++,p1++,p2++) *p2 += *p1;   
    } 
  }
  /* unpivot P into rV and PKtz into beta */

  for (i=0;i< *q;i++) beta[pivot[i]] = PKtz[i];

  for (p1=P,i=0;i < rank; i++) for (j=0;j<*q;j++,p1++) rV[pivot[j] + i * *q] = *p1;
  p0 = rV + *q * rank;p1 = rV + *q * *q;
  for (p2=p0;p2<p1;p2++) *p2 = 0.0; /* padding any trailing columns of rV with zeroes */

  /* Now get the remainder of the REML penalty */
 
  if (*REML) { /* REML or ML */
    /* First log|S|_+ */ 
    U = (double *) calloc((size_t) *q * *q,sizeof(double)); /* matrix for eigen-vectors of S (null space first) */  
    i = *q - *Mp; /* number of columns of UrS */
    /* P0, trA1 and trA2 used here for log det and its derivatives... */
    get_detS2(sp,UrS,rSncol,&i,M,deriv,P0,trA1, trA2,&d_tol,rank_tol,fixed_penalty);
    reml_penalty = - *P0;
    /* DEBUGGING printouts.... */
    /*   printf("\nderiv = %d  |S| = %g  d|S|...\n",*deriv,reml_penalty);
	 for (i=0;i<*M;i++) printf("%g  ",trA1[i]);printf("\n");*/

  } 

  if (*REML>0) { /* It's REML */
    /* Now deal with log|X'WX+S| */   
    reml_penalty += ldetXWXS;
    get_ddetXWXpS(P1,P2,P,K,sp,rS,rSncol,Tk,Tkm,n,q,&rank,M,deriv); /* P1/2 really contain det derivs */
   
    if (*deriv) for (p2=trA2,p1=P2,i = 0; i< *M;i++) { 
      trA1[i] = P1[i] - trA1[i];
      if (deriv2) for (j=0;j<*M;j++,p1++,p2++) *p2 = *p1 - *p2;   
    } 
    free(U); /* not needed */
  } /* So trA1 and trA2 actually contain the derivatives for reml_penalty */

  if (*REML<0) { /* it's ML, and more complicated */
    
    /* get derivs of ML log det in P1 and P2... */

    ldetXWXS = MLpenalty(P1,P2,Tk,Tkm,U1,R,Q1,nind,sp,rS,rSncol,q,n,&nn,Mp,M,&neg_w,rank_tol,deriv);
    
    reml_penalty += ldetXWXS;

    if (*deriv) for (p2=trA2,p1=P2,i = 0; i< *M;i++) { 
      trA1[i] =  P1[i] - trA1[i];
      if (deriv2) for (j=0;j<*M;j++,p1++,p2++) *p2 =  *p1  - *p2;   
    } 
    
    free(R);free(Q1);free(nind);free(U);
  } /* note that rS scrambled from here on... */


  pearson2(P0,P1,P2,y,mu,V0,V1,V2,g1,g2,p_weights,eta1,eta2,*n,*M,*deriv,deriv2);
  
  if (*REML) { /* really want scale estimate and derivatives in P0-P2, so rescale */
    j = *n - null_space_dim;
    *P0 /= j;
    if (*deriv) for (p1 = P1,p2 = P1 + *M;p1<p2;p1++) *p1 /= j;
    if (*deriv>1) for (p1 = P2,p2 = P2 + *M * *M;p1<p2;p1++) *p1 /= j; 
  }

  /* clean up memory, except what's needed to get tr(A) and derivatives 
     
  */ 
  
  free(pivot);free(work);free(PKtz);free(zz);
  
  if (*deriv) {
    free(b1);free(eta1);
    free(eta2);
    free(c0);free(c1);free(c2);
    free(a1);free(a2);free(wi);free(dev_grad);
    free(w1);free(w2);free(b2);

    if (deriv2) { free(dev_hess);}
  }
  
 

  /* Note: the following gets only trA if REML is being used,
           so as not to overwrite the derivatives actually needed  */
  if (*REML) i=0; else i = *deriv;
  get_trA2(trA,trA1,trA2,P,K,sp,rS,rSncol,Tk,Tkm,w,n,q,&rank,M,&i);


  free(P);free(K);
  if (*deriv)
    { free(Tk);free(Tkm);
  }

  if (*REML) {*rank_tol = reml_penalty;*conv_tol = bSb;}

  *deriv = iter; /* the number of iteration steps taken */
} /* end of gdi() */






void gdi2(double *X,double *E,double *rS,
    double *sp,double *z,double *w,double *mu,double *eta, double *y,
	 double *p_weights,double *g1,double *g2,double *g3,double *g4,double *V0,
	 double *V1,double *V2,double *V3,double *beta,double *D1,double *D2,
    double *P0, double *P1,double *P2,double *trA,
    double *trA1,double *trA2,double *rV,double *rank_tol,double *conv_tol, int *rank_est,
	 int *n,int *q, int *M,int *Encol,int *rSncol,int *deriv,int *use_svd,
	 int *REML,int *fisher)     
/* This is the original Wood (2008) code (written 2006). Now superceded by gdi().

   It expects derivatives, 
   rather than derivative ratios in g* and V*, and iterates for derivative of beta.

   Function to iterate for first and second derivatives of the deviance 
   of a GAM fit, and to evaluate the first and second derivatives of
   tr(A). Derivatives are w.r.t. log smoothing parameters.

   The function is to be called at convergence of a P-IRLS scheme so that 
   z, w, mu and functions of these can be treated as fixed, and only the 
   derivatives need to be updated.

   All matrices are packed into arrays in column order (i.e. col1, col2,...)
   as in R. 

   All names ending in 1,2 or 3 are derivatives of some sort, with the integer
   indicating the order of differentiation. 

   The arguments of this function point to the following:
   * X is and n by q model matrix.
   * E is a q by Encol square root of the total penalty matrix, so EE'=S
   * rS is a list of square roots of individual penalty matrices, packed
     in one array. The ith such matrix rSi, say, has dimension q by rSncol[i]
     and the ith penalty is [rSi][rSi]'.
   * sp is an M array of smoothing parameters (NOT log smoothing parameters)
   * z and w are n-vectors of the pseudodata and iterative weights
   * p_weights is an n-vector of prior weights (as opposed to the iterative weights in w)
   * mu and y are n-vectors of the fitted values and data.
   * g1,g2,g3,g4 are the n-vectors of the link derivatives: 
     g'(mu), g''(mu) g'''(mu) and g''''(mu)
   * V0, V1, V2, V3 are n-vectors of the variance function and first two derivatives.
     V0(mu), V'(mu), V''(mu) & V'''(mu) 
   * D1 and D2 are an M-vector and M by M matrix for returning the first 
     and second derivatives of the deviance wrt the log smoothing parameters.
     if *REML is non zero then the derivs will be of the penalized deviance,
     and b'Sb will be returned in conv_tol  
   * trA1 and trA2 are an M-vector and M by M matrix for returning the first 
     and second derivatives of tr(A) wrt the log smoothing parameters.
     If *REML is non zero then the derivatives of the REML penalty are 
     returned instead (with the REML penalty returned in `rank_tol', hack, hack).
   * P0,P1,P2 are for returning the Pearson statistic and its derivatives, or 
     the Pearson scale estimate and derivatives if *REML is non - zero. 
   * rank_est is for returning the estimated rank of the problem.
   * the remaining arguments are the dimensions already refered to except for:
   * deriv, which controls which derivatives are produced:
       deriv==0 for no derivatives: only trA, rV and beta returned
       deriv==1 for first derivatives only
       deriv==2 for gradient and Hessian
     -- on exit contains the number of iteration steps required.   

    * If REML is non-zero, then the REML penalty returned in rank_tol, with it's 
      derivatives in trA1, trA2: it is to be added to the *deviance* to get D_r.
    * non-zero `fisher' indicates that Fisher scoring, rather than full Newton,
      is the basis for iteration. 

   The method has 4 main parts:

   1. The initial QR- decomposition and SVD are performed, various quantities which 
      are independent of derivatives are created

   2. Iteration to obtain the derivatives of the coefficients wrt the log smoothing 
      parameters. 

   3. Evaluation of the derivatives of the deviance wrt the log smoothing parameters
      (i.e. calculation of D1 and D2)

   4. Evaluation of the derivatives of tr(A) (i.e. trA1 and trA2)

   The method involves first and second derivatives of a number of k-vectors wrt
   log smoothing parameters (\rho), where k is q or n. Consider such a vector, v. 
   
   * v1 will contain dv/d\rho_0, dv/d\rho_1 etc. So, for example, dv_i/d\rho_j
     (indices starting at zero) is located in v1[q*j+i].
   
   * v2 will contain d^2v/d\rho_0d\rho_0, d^2v/d\rho_1d\rho_0,... but rows will not be
     stored if they duplicate an existing row (e.g. d^2v/d\rho_0d\rho_1 would not be 
     stored as it already exists and can be accessed by interchanging the sp indices).
     So to get d^2v_k/d\rho_id\rho_j: 
     i)   if i<j interchange the indices
     ii)  off = (j*m-(j+1)*j/2+i)*q 
     iii) v2[off+k] is the required derivative.       

    

*/
{ double *zz,*WX,*tau,*work,*pd,*p0,*p1,*p2,*p3,*K=NULL,*R,*d,*Vt,*V,*U1,*KU1t=NULL,xx,*b1,*b2,*P,
         *c0,*c1,*c2,*a0,*a1,*a2,*B2z,*B2zBase,*B1z,*B1zBase,*eta1,*mu1,*eta2,*KKtz,
         *PKtz,*KPtSPKtz,*v1,*v2,*wi,*wis,*z1,*z2,*zz1,*zz2,*pz2,*w1,*w2,*pw2,*Tk,*Tkm,
         *pb2,*B1z1, *dev_grad,*dev_hess=NULL,diff,mag,*D1_old,*D2_old,Rcond,*tau2,
         ldetXWXS=0.0,reml_penalty=0.0,bSb=0.0,*U,
         *fa,*fa1,*fa2,*fb,*fc,*fc1,*fc2,*fc3,*fd,*fd1,*fd2;
  int i,j,k,*pivot,ScS,*pi,rank,r,left,tp,bt,ct,iter=0,m,one=1,n_2dCols,n_b1,n_b2,
    n_eta1,n_eta2,n_work,ok,deriv2,*pivot2,null_space_dim;

 
  if (*deriv==2) deriv2=1; else deriv2=0;
  zz = (double *)calloc((size_t)*n,sizeof(double)); /* storage for z'=Wz */
  for (i=0;i< *n;i++) zz[i] = z[i]*w[i]; /* form z'=Wz itself*/
  WX = (double *) calloc((size_t) (*n * *q),sizeof(double));
  for (j=0;j<*q;j++) for (i=0;i<*n;i++) /* form WX */
  { k = i + *n * j;
    WX[k]=w[i]*X[k];
  } 
  /* get the QR decomposition of WX */
  tau=(double *)calloc((size_t)*q,sizeof(double)); /* part of reflector storage */
 
  pivot=(int *)calloc((size_t)*q,sizeof(int));
  /* Accuracy can be improved by pivoting on some occasions even though it's not going to be 
     `used' as such here - see Golub and Van Loan (1983) section 6.4. page 169 for reference. */
  mgcv_qr(WX,n,q,pivot,tau); /* WX and tau now contain the QR decomposition information */
  /* Apply pivoting to the parameter space - this simply means reordering the rows of Sr and the 
     rS_i, and then unscrambling the parameter vector at the end (along with any covariance matrix)
     pivot[i] gives the unpivoted position of the ith pivoted parameter.
  */
  ScS=0;for (pi=rSncol;pi<rSncol + *M;pi++) ScS+= *pi;  /* total columns of input rS */
  n_work = (4 * *n + 2 * *q) * *M + 2 * *n;
  k = (*M * (1 + *M))/2 * *n;
  if (n_work < k) n_work = k;
  work = (double *)calloc((size_t) n_work,sizeof(double)); /* work space for several routines*/
  p0 = work + *q;
  for (pd=rS,i=0;i<ScS;i++,pd += *q) /* work across columns */
  { for (pi=pivot,p2=work;p2<p0;pi++,p2++) *p2 = pd[*pi];  /* apply pivot into work */
    p3 = pd + *q;
    for (p1=pd,p2=work;p1<p3;p1++,p2++) *p1 = *p2;  /* copy back into rS */
  } /* rS pivoting complete, do E .... */
  p0 = work + *q;
  for (pd=E,i=0;i< *Encol;i++,pd += *q) /* work across columns */
  { for (pi=pivot,p2=work;p2<p0;pi++,p2++) *p2 = pd[*pi];  /* apply pivot into work */
    p3 = pd + *q;
    for (p1=pd,p2=work;p1<p3;p1++,p2++) *p1 = *p2;  /* copy back into E */
  }
  /* pivot columns of X - can't think of a smart arsed way of doing this*/
  for (i=0;i< *n;i++) {
      for (j=0;j<*q;j++) work[j]=X[pivot[j] * *n + i];
      for (j=0;j<*q;j++) X[j * *n + i] = work[j];
  }
  
  /* Now form the augmented R matrix [R',E']' */
  r = *Encol + *q;
  R=(double *)calloc((size_t)(r * *q),sizeof(double));  
  for (j=0;j< *q;j++) for (i=0;i<=j;i++) R[i+r*j] = WX[i + *n * j];
  for (j=0;j< *q;j++) for (i= *q;i<r;i++) R[i+r*j]=E[j+ (i - *q) * *q ];

  if (*use_svd) {
    /* Get singular value decomposition, and hang the expense */

    d=(double *)calloc((size_t)*q,sizeof(double));
    Vt=(double *)calloc((size_t)(*q * *q),sizeof(double));
    mgcv_svd_full(R,Vt,d,&r,q);  
 
    /* now truncate the svd in order to deal with rank deficiency */
    rank= *q;xx=d[0] * *rank_tol;
    while(d[rank-1]<xx) rank--;
    *rank_est = rank;
    /* REML NOTE: |X'WX+S| is the product of the d's squared */  
    
    if (*REML) { 
        for (ldetXWXS=0.0,i=0;i<rank;i++) ldetXWXS += log(d[i]);    
       ldetXWXS *= 2;
    }

    V = (double *) calloc((size_t)(*q * rank),sizeof(double));
    U1 = (double *) calloc((size_t)(*q * rank),sizeof(double));
    /* produce the truncated V (q by rank): columns dropped so V'V=I but VV'!=I   */
    for (i=0;i< *q;i++) for (j=0;j< rank;j++) V[i+ *q * j]=Vt[j+ *q * i];
    /* produce the truncated U1 (q by rank): rows and columns dropped - no-longer orthogonal */
    for (i=0;i< *q;i++) for (j=0;j< rank;j++) U1[i+ *q * j]=R[i+r*j];
    free(R);free(Vt);
  
    /* At this stage the parameter space is pivoted, and is of dimension `rank' <= *q.
       d=diag(D) and V and U1 are available. Q can be applied via calls to mgcv_qrqy.
       Now obtain P=VD^{-1},K=QU1, KU1' and the other quantities that can be obtained before 
       iteration.
    */
    P=V; /* note: really modifying V here, V can't be used after this point */
    p3=d+rank;
    for (p0=d;p0<p3;p0++) for (i=0;i< *q;i++,P++) *P /= *p0;
    P=V;
    free(d);
  } else { /* use a second pivoted QR */
    tau2=(double *)calloc((size_t)*q,sizeof(double)); /* part of reflector storage */
 
    pivot2=(int *)calloc((size_t)*q,sizeof(int)); /* indexing vector for second pivoting */

    mgcv_qr(R,&r,q,pivot2,tau2); /* R and tau2 now contain the QR decomposition information */

    /* need to get the rank */
    rank = *q;
    R_cond(R,&r,&rank,work,&Rcond);
    while (*rank_tol * Rcond > 1) { rank--;R_cond(R,&r,&rank,work,&Rcond);}
    *rank_est = rank;

    /* REML NOTE: |X'WX+S| is the product of the R[i,i]s squared */  

    if (*REML) { 
      for (ldetXWXS=0.0,i=0;i<rank;i++) ldetXWXS += log(fabs(R[i + i * r])); 
      ldetXWXS *= 2;
    }

    /* Now get P, which is q by rank*/
    V = (double *) calloc((size_t)(*q * rank),sizeof(double));
    Rinv(V,R,&rank,&r,q);
    for (i=rank;i<*q;i++) for (j=0;j<rank;j++) V[i + j * *q]=0.0;
    P=V; /* note: don't re-use V from here on */
    /* finally need U1 */
    Vt = (double *) calloc((size_t)(r * *q),sizeof(double));
    for (p0=Vt,i=0;i<*q;i++,p0 += r+1) *p0 = 1.0; 
    left=1;tp=0;mgcv_qrqy(Vt,R,tau2,&r,q,q,&left,&tp); /* Vt now contains U */
    U1 = (double *) calloc((size_t)(*q * rank),sizeof(double));
    for (i=0;i< *q;i++) for (j=0;j< rank;j++) U1[i+ *q * j]=Vt[i+r*j];
    free(Vt);free(R);free(tau2);
    /* need to unpivot rows of P */
  
    for (j=0;j<rank;j++) {
	for (i=0;i<*q;i++) work[pivot2[i]] = P[i + j * *q];
        for (i=0;i<*q;i++) P[i + j * *q] = work[i]; 
    }
    free(pivot2);   
  }

  PKtz = (double *)calloc((size_t) *q,sizeof(double)); /* PK'z */
  
  if (*deriv) { /* then following O(nq^2) required */
    K = (double *)calloc((size_t) *n * rank,sizeof(double));
    p0=U1;p1=K; /* first q rows of U0 should be U1 */
    for (j=0;j<rank;j++,p1+= *n) { p3=p1 + *q;for (p2=p1;p2<p3;p2++,p0++) *p2 = *p0;} 
    left=1;tp=0;mgcv_qrqy(K,WX,tau,n,&rank,q,&left,&tp); /* QU1 = Q%*%U1, now */

    KU1t = (double *)calloc((size_t) *n * *q,sizeof(double));
    bt=0;ct=1;mgcv_mmult(KU1t,K,U1,&bt,&ct,n,q,&rank);
  } else { /* evaluate coefficients more efficiently */
    /* PKtz is P[U1]'Q'zz */
    left=1;tp=1;mgcv_qrqy(zz,WX,tau,n,&one,q,&left,&tp); /* puts Q'zz in zz */
    bt=1;ct=0;mgcv_mmult(work,U1,zz,&bt,&ct,&rank,&one,q); /* puts [U1]'Q'zz in work */
    /*    bt=1;ct=0;mgcv_mmult(work,K,zz,&bt,&ct,&rank,&one,n);*/
    bt=0;ct=0;mgcv_mmult(PKtz,P,work,&bt,&ct,q,&one,&rank);
  }

  /************************************************************************************/
  /* The coefficient derivative iteration starts here */
  /************************************************************************************/
  /* set up some storage first */
  if (*deriv) {
    n_2dCols = (*M * (1 + *M))/2;
    n_b2 = *q * n_2dCols;
    b2 = (double *)calloc((size_t)n_b2,sizeof(double)); /* 2nd derivs of beta */
    B2zBase = (double *)calloc((size_t)n_b2,sizeof(double)); /* part of 2nd derivs of beta */
    B2z = (double *)calloc((size_t)n_b2,sizeof(double)); /* part of 2nd derivs of beta */
    n_b1 = *q * *M;
    b1 = (double *)calloc((size_t)n_b1,sizeof(double)); /* 1st derivs of beta */
    B1zBase = (double *)calloc((size_t)n_b1,sizeof(double)); /* part of 1st derivs of beta */
    B1z = (double *)calloc((size_t)n_b1,sizeof(double)); /* part of 1st  derivs of beta */
    n_eta1 = *n * *M;
    eta1 = (double *)calloc((size_t)n_eta1,sizeof(double));
    Tk = (double *)calloc((size_t)n_eta1,sizeof(double));
    mu1 = (double *)calloc((size_t)n_eta1,sizeof(double));
    z1 = (double *)calloc((size_t)n_eta1,sizeof(double));
    zz1 = (double *)calloc((size_t)n_eta1,sizeof(double));
    w1 = (double *)calloc((size_t)n_eta1,sizeof(double));

    n_eta2 = *n * n_2dCols;
    eta2 = (double *)calloc((size_t)n_eta2,sizeof(double));
    Tkm = (double *)calloc((size_t)n_eta2,sizeof(double));
    z2 = (double *)calloc((size_t)n_eta2,sizeof(double));
    zz2 = (double *)calloc((size_t)n_eta2,sizeof(double));
    w2 = (double *)calloc((size_t)n_eta2,sizeof(double));

    B1z1 = (double *)calloc((size_t)(*M * *M * *q),sizeof(double));
    KKtz = (double *)calloc((size_t) *n,sizeof(double)); /* KK'z */
  
    KPtSPKtz = (double *)calloc((size_t)(*n * *M),sizeof(double)); /* KP'S_kPK'z */
  
 
    v1 = work;v2=work + *n * *M; /* a couple of working vectors */ 
    /* Now get the iteration independent parts of the derivatives, which are also
       the intial values for the derivatives */
    B1B2zBaseSetup(B2zBase,B1zBase,zz,P,K,KKtz,PKtz,KPtSPKtz,rS,
                 rSncol,n,q,&rank,M,sp,work,deriv);
    /* need to copy B2zBase and B1zBase into b2 and b1 */
    if (deriv2) {
      p1=b2+n_b2;for (p0=b2,p2=B2zBase;p0<p1;p0++,p2++)  *p0 = *p2;  
    }
    p1=b1+n_b1;for (p0=b1,p2=B1zBase;p0<p1;p0++,p2++)  *p0 = *p2;  

    /* Set up constants involved in z (not z'!) updates (little work => leave readable!)*/
    c0=(double *)calloc((size_t)*n,sizeof(double));
    c1=(double *)calloc((size_t)*n,sizeof(double));  
    c2=(double *)calloc((size_t)*n,sizeof(double));
    a0=(double *)calloc((size_t)*n,sizeof(double));
    a1=(double *)calloc((size_t)*n,sizeof(double));  
    a2=(double *)calloc((size_t)*n,sizeof(double));
    fa=fa1=fa2=fb=fc=fc1=fc2=fc3=fd=fd1=fd2=(double *)NULL;
    if (*fisher) { /* Fisher scoring updates */
      for (i=0;i< *n;i++) c0[i]=y[i]-mu[i];
      for (i=0;i<*n;i++) c2[i]=g2[i]/g1[i];
      /* c1 = (y-mu)*g2/g1 */
      for (i=0;i<*n;i++) c1[i]=c2[i]*c0[i];
      /* c2 = (y-mu)*(g3/g1-g2/g1) - g2/g1 */
      for (i=0;i<*n;i++) c2[i]=c0[i]*(g3[i]-g2[i]*g2[i]/g1[i])/g1[i]-c2[i];

#ifdef DEBUG
      printf("\n c0:\n");
	for (i=0;i<*n;i++) printf("  %g",c0[i]);
      printf("\n c1:\n");
	for (i=0;i<*n;i++) printf("  %g",c1[i]);
       printf("\n c2:\n");
	for (i=0;i<*n;i++) printf("  %g",c2[i]);
#endif    

      /* set up constants involved in w updates */
   
      for (i=0;i< *n;i++) a0[i] = - w[i]*w[i]*w[i]*(V1[i]*g1[i]+2*V0[i]*g2[i])/(2*p_weights[i]) ;
      for (i=0;i< *n;i++) a1[i] = 3/w[i];
      for (i=0;i< *n;i++) 
        a2[i] = -w[i]*w[i]*w[i]*(V2[i]*g1[i]+3*V1[i]*g2[i]+2*g3[i]*V0[i])/(g1[i]*2*p_weights[i]);
#ifdef DEBUG
      printf("\n\n\n\n a0:\n");
	 for (i=0;i<*n;i++) printf("  %g",a0[i]);
      printf("\n a1:\n");
	for (i=0;i<*n;i++) printf("  %g",a1[i]);
      printf("\n a2:\n");
	for (i=0;i<*n;i++) printf("  %g",a2[i]);
#endif
    } else { /* full Newton updates */
      fa = (double *) calloc((size_t)*n,sizeof(double));
      fa1 = (double *) calloc((size_t)*n,sizeof(double));
      fa2 = (double *) calloc((size_t)*n,sizeof(double));
      for (i=0;i< *n;i++) {
        fa[i] = 1/g1[i];
        fa1[i] = -g2[i]*fa[i]*fa[i];
        fa2[i] = -2*fa1[i]*g2[i]*fa[i] - g3[i]*fa[i]*fa[i];
      }
      fc = (double *) calloc((size_t)*n,sizeof(double));
      fc1 = (double *) calloc((size_t)*n,sizeof(double));
      fc2 = (double *) calloc((size_t)*n,sizeof(double));
      fc3 = (double *) calloc((size_t)*n,sizeof(double));
      for (i=0;i< *n;i++) {
        fc[i] = V0[i]*g1[i];
        fc1[i] = V1[i]*g1[i] + V0[i]*g2[i];
        fc2[i] = V2[i]*g1[i] + 2*V1[i]*g2[i] + V0[i]*g3[i];
        fc3[i] = V3[i]*g1[i] + 3*(V2[i]*g2[i]+V1[i]*g3[i]) + V0[i]*g4[i];
      }
      fb = (double *) calloc((size_t)*n,sizeof(double));
      for (i=0;i< *n;i++) fb[i] = y[i] - mu[i];
      fd = (double *) calloc((size_t)*n,sizeof(double));
      fd1 = (double *) calloc((size_t)*n,sizeof(double));
      fd2 = (double *) calloc((size_t)*n,sizeof(double));
      for (i=0;i< *n;i++) {
        fd[i] = fa[i]*(1+fb[i]*fc1[i]/fc[i]);
        fd1[i] = fd[i]*(fa1[i]-fa[i]*fc1[i]/fc[i]) + fa[i]*fa[i]*fb[i]*fc2[i]/fc[i];
        fd2[i] = fd1[i]*(fa1[i]-fa[i]*fc1[i]/fc[i]) +
	         fd[i]*fa[i]*(fa2[i] - (fa1[i]*fc1[i]+fa[i]*fc2[i]-fa[i]*fc1[i]*fc1[i]/fc[i])/fc[i]) +
                 2*fa[i]*fa[i]*fa1[i]*fb[i]*fc2[i]/fc[i] +
	         fa[i]*fa[i]*fa[i]*(fb[i]*fc3[i]-fc2[i]-fb[i]*fc2[i]*fc1[i]/fc[i])/fc[i];
      }
      /* end of preliminaries, now setup the multipliers that go forward */
      /* dz/deta... */
      for (i=0;i<*n;i++) c1[i] = 1 - (fa[i]+fb[i]*fd1[i]/fd[i])/fd[i];
      /* d2z/deta2... */
      for (i=0;i<*n;i++) c2[i] = ((2*fa[i]*fd1[i] - fb[i]*fd2[i] + 2*fb[i]*fd1[i]*fd1[i]/fd[i])/fd[i] - fa1[i]*fa[i])/fd[i];
      /* dw/deta... */
      for (i=0;i<*n;i++) a0[i] = 0.5*p_weights[i]*(fd1[i]/sqrt(fd[i]) - sqrt(fd[i])*fc1[i]*fa[i]/fc[i])/sqrt(fc[i]); 
      /* multiplier for dw/deta product in d2w/deta2... */
      for (i=0;i<*n;i++) a1[i] = 1/w[i];
      /* multiplier for remainder term in d2w/deta2... */
      for (i=0;i<*n;i++) a2[i] =  0.5*w[i]*(fd2[i]/fd[i]+fa[i]*(fc1[i]*fc1[i]*fa[i]/fc[1] - fc2[i]*fa[i] - fc1[i]*fa1[i])/fc[i]);
      free(fa);free(fa1);free(fa2);free(fb);free(fc);free(fc1);free(fc2);free(fc3);free(fd);free(fd1);free(fd2);
      
    } /* end of full Newton setup */
    /* some useful arrays for Tk and Tkm */
    wi=(double *)calloc((size_t)*n,sizeof(double));
    wis=(double *)calloc((size_t)*n,sizeof(double));
    for (i=0;i< *n;i++) { wi[i]=1/w[i];wis[i]=wi[i]*wi[i];}

    /* get gradient vector and Hessian of deviance wrt coefficients */
    for (i=0;i< *n ;i++) v1[i] = -2*p_weights[i]*(y[i]-mu[i])/(V0[i]*g1[i]);
    dev_grad=(double *)calloc((size_t)*q,sizeof(double));
    bt=1;ct=0;mgcv_mmult(dev_grad,X,v1,&bt,&ct,q,&one,n);
    
    if (deriv2) { /* get hessian of deviance w.r.t. beta */
      for (i=0;i< *n ;i++) 
      v1[i] = 2*p_weights[i]*
            (1/V0[i] + (y[i]-mu[i])/(V0[i]*V0[i]*g1[i])*(V1[i]*g1[i]+V0[i]*g2[i]))/(g1[i]*g1[i]);
      dev_hess=(double *)calloc((size_t)(*q * *q),sizeof(double));
      getXtWX(dev_hess,X,v1,n,q,v2);
    } 
  
    /* create storage for gradient and Hessian of deviance wrt sp's from previous iteration,
       for convergence testing */
    D1_old =(double *)calloc((size_t)(*M),sizeof(double));
    D2_old =(double *)calloc((size_t)(*M * *M),sizeof(double));

    /* create initial gradient and Hessian of deviance */

    bt=1;ct=0;mgcv_mmult(D1,b1,dev_grad,&bt,&ct,M,&one,q); /* gradient of deviance is complete */
    for (p0=D1,p2=D1_old,p1=D1 + *M;p0<p1;p0++,p2++) *p2 = *p0; /* store D1 for convergence testing */

    if (deriv2) {
      getXtMX(D2,b1,dev_hess,q,M,v1);
          
      for (pb2=b2,m=0;m < *M;m++) for (k=m;k < *M;k++) { /* double sp loop */
         p1 = dev_grad + *q;  
         for (xx=0.0,p0=dev_grad;p0<p1;p0++,pb2++) xx += *p0 * *pb2;
         D2[k * *M + m] += xx;
         D2[m * *M + k] = D2[k * *M + m];
      } /* Hessian of Deviance is complete !! */ 
    
      /* store D2 for convergence testing */
      for (p0=D2,p2=D2_old,p1=D2 + *M * *M;p0<p1;p0++,p2++) *p2 = *p0;
    }
    
  
    /* NOTE: when DEBUG complete, better to store initial D1 and D2 directly in D1_old and D2_old */
    
    for (iter=0;iter<100;iter++) { /* main derivative iteration */ 
      /* get derivatives of eta and mu */
      bt=0;ct=0;mgcv_mmult(eta1,X,b1,&bt,&ct,n,M,q);
      if (deriv2) {
        bt=0;ct=0;mgcv_mmult(eta2,X,b2,&bt,&ct,n,&n_2dCols,q);
     
        p2 = g1 + *n;
        for (p3=mu1,p0=eta1,i=0;i < *M;i++) 
	    for (p1=g1;p1<p2;p1++,p0++,p3++) *p3 = *p0 / *p1;
      }
      /* update the derivatives of z */
      rc_prod(z1,c1,eta1,M,n); /* z1 = dz/d\rho_k done */
      if (deriv2) {  
        rc_prod(z2,c1,eta2,&n_2dCols,n);
        for (pz2=z2,m=0;m < *M;m++) for (k=m;k < *M;k++) {
	    rc_prod(v1,mu1 + *n * m,eta1 + *n * k,&one,n);
            rc_prod(v2,c2,v1,&one,n);
            p1=v2 + *n;
            for (p0=v2;p0<p1;p0++,pz2++) *pz2 += *p0;        
        } /* z2 update completed */
      }
     /* update derivatives of w */  
      rc_prod(w1,a0,eta1,M,n); /* w1 = dw/d\rho_k done */
      if (deriv2) {
        rc_prod(w2,a0,eta2,&n_2dCols,n); 
        for (pw2=w2,m=0;m < *M;m++) for (k=m;k < *M;k++) {
	    rc_prod(v1,eta1 + *n * m,eta1 + *n * k,&one,n);
            rc_prod(v2,a2,v1,&one,n);
            p1=v2 + *n;
            for (p0=v2;p0<p1;p0++,pw2++) *pw2 += *p0;        
            pw2 -= *n;
            rc_prod(v1,w1 + *n * m,w1 + *n * k,&one,n);
            rc_prod(v2,a1,v1,&one,n);
            p1=v2 + *n;
            for (p0=v2;p0<p1;p0++,pw2++) *pw2 += *p0;     
        } /* w2 update completed */
      }
      /* update Tk and Tkm */
      
      rc_prod(Tk,wi,w1,M,n); /* Tk done */
      if (deriv2) {
        rc_prod(Tkm,wi,w2,&n_2dCols,n);
        for (p0=Tkm,m=0;m < *M;m++) for (k=m;k < *M;k++) {
	    rc_prod(v1,w1+k * *n,w1+m * *n,&one,n);
            rc_prod(v2,wis,v1,&one,n);
            p2 = v2 + *n;
            for (p1=v2;p1<p2;p1++,p0++) *p0 -= *p1; 
        } /* Tkm finished */
      } 
      /* update the derivatives of z' (zz) */
      rc_prod(zz1,z,w1,M,n); /* dw_i/d\rho_k z_i */
      rc_prod(v1,w,z1,M,n);  /* dz_i/d\rho_k w_i */
      p2 = v1 + *M * *n;
      for (p0=v1,p1=zz1;p0<p2;p0++,p1++) *p1 += *p0; /*zz1=dz'/d\rho_k finished */
      
      if (deriv2) {
        rc_prod(zz2,z,w2,&n_2dCols,n);
        rc_prod(work,w,z2,&n_2dCols,n); /* NOTE: w2 over-written here! */
        p2 = zz2 + n_2dCols * *n;
        for (p0=zz2,p1=work;p0<p2;p0++,p1++) *p0 += *p1; 
     
        for (pz2=zz2,m=0;m < *M;m++) for (k=m;k < *M;k++) {
           rc_prod(v1,w1+ m * *n ,z1 + k * *n,&one,n);
           rc_prod(v2,w1+ k * *n ,z1 + m * *n,&one,n);
           p1 = v1 + *n;
           for (p0=v1,p2=v2;p0<p1;p0++,p2++,pz2++) *pz2 += *p0 + *p2; 
        } /* zz2 complete */
      }       
      /* update derivatives of \beta */
      
      /* Start with the horrid term B2z, and get B1z at the same time */
     
      B1B2z(B2z,B1z,B2zBase,B1zBase,zz,Tk,Tkm,P,K,KKtz,PKtz,KPtSPKtz,rS,
	rSncol,n,q,&rank,M,sp,work,deriv);
     
      /* now evaluate Bzz1 and Bzz2 and add them to B1z and B2z,
         using w1 and w2 as work-space. */

      bt=1;ct=0;mgcv_mmult(work,K,zz1,&bt,&ct,&rank,M,n);
      bt=0;ct=0;mgcv_mmult(b1,P,work,&bt,&ct,q,M,&rank); /* b1 = B zz1, currently */

      p2 = b1 + *M * *q;
      for (p0=b1,p1=B1z;p0<p2;p0++,p1++) *p0 += *p1; /* b1 complete */

      if (deriv2) {
        bt=1;ct=0;mgcv_mmult(work,K,zz2,&bt,&ct,&rank,&n_2dCols,n);
        bt=0;ct=0;mgcv_mmult(b2,P,work,&bt,&ct,q,&n_2dCols,&rank); /* b2 = B zz2, currently */
        p2 = b2 + n_2dCols * *q;
        for (p0=b2,p1=B2z;p0<p2;p0++,p1++) *p0 += *p1; /* b2 = B2 zz + B zz2, currently */
       
        /* now get the B1 zz1 cross terms by calling getB1z1 */

        getB1z1(B1z1,zz1,K,P,Tk,sp,rS,rSncol,n,&rank,q,M,work); 
       
        /* (dB/d\rho_k dz'/d\rho_m)[i] is in B1Z1[q*m+M*k*q+i] */
      
        pb2 = b2;   
        for (m=0;m < *M;m++) for (k=m;k < *M;k++) { /* double sp loop */
	    p0=B1z1 + *q * m + *M * *q * k;p2 = p0 + *q;
            p1=B1z1 + *q * k + *M * *q * m;
            for (;p0<p2;p0++,p1++,pb2++) *pb2 += *p1 + *p0; 
        } /* b2 complete */
      } 
      /* evaluate gradient and Hessian of deviance (since these are what convergence 
         should be judged on) */
      bt=1;ct=0;mgcv_mmult(D1,b1,dev_grad,&bt,&ct,M,&one,q); /* gradient of deviance is complete */
      
      if (deriv2) {       
        getXtMX(D2,b1,dev_hess,q,M,v1);
          
        for (pb2=b2,m=0;m < *M;m++) for (k=m;k < *M;k++) { /* double sp loop */
          p1 = dev_grad + *q;  
          for (xx=0.0,p0=dev_grad;p0<p1;p0++,pb2++) xx += *p0 * *pb2;
          D2[k * *M + m] += xx;
          D2[m * *M + k] = D2[k * *M + m];
        } /* Hessian of Deviance is complete !! */
      }

      /* now test for convergence */
      ok=1;
      /* test first derivative convergence */ 
      for (diff=mag=0,p0=D1,p2=D1_old,p1=D1 + *M;p0<p1;p0++,p2++) {
         xx = fabs(*p0 - *p2); /* change in derivative */
         if (xx>diff) diff=xx;
         xx = (fabs(*p0) + fabs(*p2))/2; /* size of derivative */
         if (xx>mag) mag=xx; 
      } 
      if (diff > mag * *conv_tol) ok=0;
      /* and do same for second derivatives */
      if (deriv2) {
        for (diff=mag=0,p0=D2,p2=D2_old,p1=D2 + *M * *M;p0<p1;p0++,p2++) {
           xx = fabs(*p0 - *p2); /* change in derivative */
           if (xx>diff) diff=xx;
           xx = (fabs(*p0) + fabs(*p2))/2; /* size of derivative */
           if (xx>mag) mag=xx;
        } 
        if (diff > mag * *conv_tol) ok=0;
      }
      if (ok) break; /* converged */
        
      /* store D1 and D2 for convergence testing */
      for (p0=D1,p2=D1_old,p1=D1 + *M;p0<p1;p0++,p2++) *p2 = *p0;
      if (deriv2) for (p0=D2,p2=D2_old,p1=D2 + *M * *M;p0<p1;p0++,p2++) *p2 = *p0;


    } /* end of main derivative iteration */
  } /* end of if (*deriv) */ 
  else { /* keep compilers happy */
    b1=B1zBase=B1z=eta1=mu1=eta2=B1z1=KKtz=KPtSPKtz=c0=c1=c2=(double *)NULL;
    a0=a1=a2=wi=wis=dev_grad=D1_old=D2_old=z1=z2=zz1=zz2=w1=w2=b2=B2zBase=B2z=(double *)NULL;
    Tk=Tkm=(double *)NULL;
  }
  /************************************************************************************/
  /* End of the coefficient derivative iteration  */
  /************************************************************************************/

  /* REML NOTE: \beta'S\beta stuff has to be done here on pivoted versions.
     Store bSb in P0, bSb1 in P1 and bSb2 in P2.
  */
  if (*REML) {
    get_bSb(&bSb,trA1,trA2,sp,E,rS,rSncol,Encol,q,M,PKtz,b1,b2,deriv);
    if (*deriv) for (p2=D2,p1=trA2,i = 0; i< *M;i++) { /* penalized deviance derivs needed */
        D1[i] += trA1[i];
        if (deriv2) for (j=0;j<*M;j++,p1++,p2++) *p2 += *p1;   
    } 
  }
  /* unpivot P into rV and PKtz into beta */

  for (i=0;i< *q;i++) beta[pivot[i]] = PKtz[i];

  for (p1=P,i=0;i < rank; i++) for (j=0;j<*q;j++,p1++) rV[pivot[j] + i * *q] = *p1;
  p0 = rV + *q * rank;p1 = rV + *q * *q;
  for (p2=p0;p2<p1;p2++) *p2 = 0.0; /* padding any trailing columns of rV with zeroes */

  /* Now get the remainder of the REML penalty */
  if (*REML) {
    /* First deal with log|X'WX+S| */   
    reml_penalty = ldetXWXS;
    get_ddetXWXpS(trA1,trA2,P,K,sp,rS,rSncol,Tk,Tkm,n,q,&rank,M,deriv); /* trA? really contain det derivs */
   
    /* Now log|S|_+ */
    U = (double *) calloc((size_t) *q * *q,sizeof(double));
    free(U);
    reml_penalty -= *P0;
    if (*deriv) for (p2=trA2,p1=P2,i = 0; i< *M;i++) { 
      trA1[i] -= P1[i];
      if (deriv2) for (j=0;j<*M;j++,p1++,p2++) *p2 -= *p1;   
    } 
  } /* So trA1 and trA2 actually contain the derivatives for reml_penalty */

  pearson(w,w1,w2,z,z1,z2,eta,eta1,eta2,P0,P1,P2,work,*n,*M,*deriv,deriv2);

  if (*REML) { /* really want scale estimate and derivatives in P0-P2, so rescale */
    j = *n - null_space_dim;
    *P0 /= j;
    if (*deriv) for (p1 = P1,p2 = P1 + *M;p1<p2;p1++) *p1 /= j;
    if (*deriv>1) for (p1 = P2,p2 = P2 + *M * *M;p1<p2;p1++) *p1 /= j; 
  }

  /* clean up memory, except what's needed to get tr(A) and derivatives 
     Note: Vt and R already freed. P is really V - don't free yet.
  */ 
  
  free(zz);free(WX);free(tau);free(pivot);free(work);free(PKtz);
  
  if (*deriv) {
    free(b1);free(B1zBase);free(B1z);free(eta1);free(mu1);
    free(eta2);free(B1z1);free(KKtz);
    free(KPtSPKtz);free(c0);free(c1);free(c2);free(a0);
    free(a1);free(a2);free(wi);free(wis);free(dev_grad);
    free(D1_old);
    free(D2_old);free(z1);free(z2);free(zz1);free(zz2);
    free(w1);free(w2);free(b2);free(B2zBase);free(B2z);

    if (deriv2) { free(dev_hess);}
  }
  
 

  /* Note: the following gets only trA if REML is being used,
           so as not to overwrite the derivatives actually needed  */
  if (*REML) i=0; else i = *deriv;
  get_trA(trA,trA1,trA2,U1,KU1t,P,K,sp,rS,rSncol,Tk,Tkm,n,q,&rank,M,&i);

  /* clear up the remainder */
  free(U1);free(V);

  if (*deriv)
  { free(Tk);free(Tkm);free(KU1t);free(K);
  }

  if (*REML) {*rank_tol = reml_penalty;*conv_tol = bSb;}

  *deriv = iter; /* the number of iteration steps taken */
} /* end of gdi() */


void R_cond(double *R,int *r,int *c,double *work,double *Rcondition)
/* estimates the condition number of c by c upper triangular matrix 
   R stored in an r by c matrix, stored in column order. 
   work should be of minimum length c * 4. 
   Algorithm is due to Cline, Moler, Stewart and Wilkinson (1979) as reported in
   Golub and Van Loan (1996).
   n <- 100;m <- 10
   X <- matrix(runif(n*m),n,m)
   X[,1] <- X[,m]
   qrx <- qr(X,LAPACK=TRUE)
   R <- qr.R(qrx)
   .C("R_cond",R=as.double(R),r=as.integer(m),c=as.integer(m),
       work=as.double(rep(0,4*m)),Rcond=as.double(1),PACKAGE="mgcv")$Rcond
*/ 
{ double kappa,*pm,*pp,*y,*p,ym,yp,pm_norm,pp_norm,y_inf=0.0,R_inf=0.0;
  int i,j,k;
  /* allocate work */
  pp=work;work+= *c;pm=work;work+= *c;
  y=work;work+= *c;p=work;work+= *c;   
  for (i=0;i<*c;i++) p[i] = 0.0;
  for (k=*c-1;k>=0;k--) {
      yp = (1-p[k])/R[k + *r *k];
      ym = (-1-p[k])/R[k + *r *k];
      for (pp_norm=0.0,i=0;i<k;i++) { pp[i] = p[i] + R[i + *r * k] * yp;pp_norm += fabs(pp[i]);}
      for (pm_norm=0.0,i=0;i<k;i++) { pm[i] = p[i] + R[i + *r * k] * ym;pm_norm += fabs(pm[i]);}
      if (fabs(yp)+pp_norm >= fabs(ym)+pm_norm) {
	  y[k]=yp;
          for (i=0;i<k;i++) p[i] = pp[i];
      } else {
          y[k]=ym;
          for (i=0;i<k;i++) p[i] = pm[i];
      }
      kappa=fabs(y[k]);
      if (kappa>y_inf) y_inf=kappa;
  }
  for (i=0;i<*c;i++) { 
    for (kappa=0.0,j=i;j<*c;j++) kappa += fabs(R[i + *r * j]);  
    if (kappa>R_inf) R_inf = kappa;
  }
  kappa=R_inf*y_inf;
  *Rcondition=kappa;
}


void pls_fit0(double *y,double *X,double *w,double *E,int *n,int *q,int *cE,double *eta,
             double *penalty,double *rank_tol)
/* Fast but stable PLS fitter. Obtains linear predictor, eta, of weighted penalized linear model,
   without evaluating the coefficients, but also returns coefficients in case they are needed. 

   Original version which can not deal with -ve weights since it assumes w_i are the square root 
   weights on entry....

   On return:
   * eta contains the linear predictor
   * penalty is the evaluated penalty
   * the first q elements of y are the coefficients. 
   
   n <- 100;x <- runif(n);w <- rep(1,n)
   X <- model.matrix(~x+I(x^2))
   X[,1] <- X[,3]
   y <- rnorm(n)
   E <- diag(3)
   oo <- .C("pls_fit",y=as.double(y),as.double(X),as.double(w),as.double(E),as.integer(n),
            as.integer(ncol(X)),as.integer(3),eta=as.double(1:n),penalty=as.double(1),
            as.double(.Machine$double.eps*100),PACKAGE="mgcv")
   er <- lm(c(y,rep(0,ncol(E)))~rbind(X,t(E))-1,weights=c(w^2,rep(1,ncol(E))))
   range(fitted(er)[1:n]-oo$eta)
   coef(er);oo$y[1:ncol(X)]
*/

{ int nn,i,j,k,ii,rank,one=1,*pivot,left,tp;
  double *z,*WX,*tau,Rcond,xx,*work;
  nn= *n + *cE;
  z = (double *)calloc((size_t)nn,sizeof(double)); /* storage for z=[Wz,0] */
  for (i=0;i< *n;i++) z[i] = y[i]*w[i]; /* form z itself*/

  WX = (double *) calloc((size_t) ( nn * *q),sizeof(double));
  for (j=0;j<*q;j++) 
  { for (i=0;i<*n;i++) /* form WX */
    { k = i + nn * j;
      WX[k]=w[i]*X[i + *n *j];
    }
    for (ii=0,i = *n;ii<*cE;i++,ii++) /* append E' */ 
    { k = i + nn * j;
      WX[k] = E[j + *q * ii];
    }
  } 
  /* get the QR decomposition of WX */
  tau=(double *)calloc((size_t)*q,sizeof(double)); /* part of reflector storage */
 
  pivot=(int *)calloc((size_t)*q,sizeof(int));
  
  mgcv_qr(WX,&nn,q,pivot,tau); /* WX and tau now contain the QR decomposition information */
  /* pivot[i] gives the unpivoted position of the ith pivoted parameter.*/
  
  /* first find the rank of R */
  work = (double *)calloc((size_t)(4 * *q),sizeof(double));
  rank = *q;
  R_cond(WX,&nn,&rank,work,&Rcond);
  while (*rank_tol * Rcond > 1) { rank--;R_cond(WX,&nn,&rank,work,&Rcond);}
  free(work);

  /* Now get the fitted values X \beta, *without* finding \beta */
  left=1;tp=1;mgcv_qrqy(z,WX,tau,&nn,&one,q,&left,&tp); /* z = Q'z */
  for (i=rank;i<nn;i++) z[i]=0.0;
  for (i=0;i<rank;i++) y[i] = z[i];        /* y = Q'z */ 
  left=1;tp=0;mgcv_qrqy(z,WX,tau,&nn,&one,q,&left,&tp);
  for (i=0;i<*n;i++) eta[i] = z[i]/w[i]; /* the linear predictor */
  for (*penalty=0.0,i=*n;i<nn;i++) *penalty += z[i]*z[i]; /* the penalty term */
  
  /* now find  \hat \beta = R^{-1}Q'z, which are needed if P-IRLS starts to diverge
     in order to be able to evaluate penalty on step reduction */
  
 

  /* now back substitute to find \hat \beta */  
  for (k=rank;k<*q;k++) z[k]=0.0; /* truncated parameters */
  for (k=rank-1;k>=0;k--) {
      for (xx=0.0,j=k+1;j < rank;j++) xx += WX[k + nn * j]*z[j];
      z[k] = (y[k] - xx)/WX[k + nn * k];
  }
  /* unpivot result into y */
  for (i=0;i< *q;i++) y[pivot[i]] = z[i];
 
  /*, should this be required... */
  free(z);free(WX);free(tau);free(pivot);
}


void pls_fit(double *y,double *X,double *w,double *E,int *n,int *q,int *cE,double *eta,
             double *penalty,double *rank_tol)
/* Fast but stable PLS fitter. Obtains linear predictor, eta, of weighted penalized linear model,
   without evaluating the coefficients, but also returns coefficients in case they are needed. 
   
   In this version the w_i are the w_i in \sum_i w_i (y_i - X_i \beta)^2
   rather than being the square root of these. Some w_i may be negative (as
   may occur when using Newton, rather than Fisher updates on IRLS). Note that it is still 
   assumed that any zero weighted data will have been dropped before the call.

   On return:
   
   * if *n is -ve then X'WX+E'E was not +ve definite (which means that the routine should be
     called again with weights based on Fisher scoring).

   otherwise:

   * eta contains the linear predictor
   * penalty is the evaluated penalty
   * the first q elements of y are the coefficients. 
   

   n <- 100;x <- runif(n);w <- runif(100)
   X <- model.matrix(~x+I(x^2))
   X[,1] <- X[,3]
   y <- rnorm(n)
   E <- diag(3)
   oo <- .C("pls_fit",y=as.double(y),as.double(X),as.double(w),as.double(E),as.integer(n),
            as.integer(ncol(X)),as.integer(3),eta=as.double(1:n),penalty=as.double(1),
            as.double(.Machine$double.eps*100),PACKAGE="mgcv")
   er <- lm(c(y,rep(0,ncol(E)))~rbind(X,t(E))-1,weights=c(w,rep(1,ncol(E))))
   range(fitted(er)[1:n]-oo$eta)
   as.numeric(coef(er));oo$y[1:ncol(X)]


   n <- 100;x <- runif(n);w <- runif(100)-.2
   X <- model.matrix(~x+I(x^2))
   y <- rnorm(n)
   E <- diag(3)
   oo <- .C("pls_fit",y=as.double(y),as.double(X),as.double(w),as.double(E),n=as.integer(n),
            as.integer(ncol(X)),as.integer(3),eta=as.double(1:n),penalty=as.double(1),
            as.double(.Machine$double.eps*100),PACKAGE="mgcv")
   beta <- solve((t(X)%*%(w*X)+t(E)%*%E),t(X)%*%(w*y))   
   eta <- X%*%beta

   range(eta-oo$eta)
   as.numeric(beta);oo$y[1:ncol(X)]

   oo$penalty;sum((E%*%beta)^2)
    

*/

{ int nn,i,j,k,ii,rank,one=1,*pivot,left,tp,neg_w=0,*nind,bt,ct;
  double *z,*WX,*tau,Rcond,xx,*work,*Q1,*IQ,*raw,*d,*Vt,*p0,*p1;
  nn= *n + *cE;
  
  z = (double *)calloc((size_t)nn,sizeof(double)); /* storage for z=[sqrt(|W|)z,0] */
  raw = (double *)calloc((size_t) *n,sizeof(double)); /* storage for sqrt(|w|) */
  
  for (i=0;i< *n;i++) 
    if (w[i]<0) { neg_w++;raw[i] = sqrt(-w[i]);} 
    else raw[i] = sqrt(w[i]);

  if (neg_w) {
    nind = (int *)calloc((size_t)neg_w,sizeof(int)); /* index the negative w_i */
    k=0;for (i=0;i< *n;i++) if (w[i]<0) { nind[k]=i;k++;}
  } else { nind = (int *)NULL;}

  for (i=0;i< *n;i++) z[i] = y[i]*raw[i]; /* form z itself*/

  for (i=0;i<neg_w;i++) {k=nind[i];z[k] = -z[k];} 

  WX = (double *) calloc((size_t) ( nn * *q),sizeof(double));
  for (j=0;j<*q;j++) 
  { for (i=0;i<*n;i++) /* form WX */
    { k = i + nn * j;
      WX[k]=raw[i]*X[i + *n *j];
    }
    for (ii=0,i = *n;ii<*cE;i++,ii++) /* append E' */ 
    { k = i + nn * j;
      WX[k] = E[j + *q * ii];
    }
  } 
  /* get the QR decomposition of WX */
  tau=(double *)calloc((size_t)*q,sizeof(double)); /* part of reflector storage */
 
  pivot=(int *)calloc((size_t)*q,sizeof(int));
  
  mgcv_qr(WX,&nn,q,pivot,tau); /* WX and tau now contain the QR decomposition information */
  /* pivot[i] gives the unpivoted position of the ith pivoted parameter.*/
  
  /* first find the rank of R */
  work = (double *)calloc((size_t)(4 * *q),sizeof(double));
  rank = *q;
  R_cond(WX,&nn,&rank,work,&Rcond);
  while (*rank_tol * Rcond > 1) { rank--;R_cond(WX,&nn,&rank,work,&Rcond);}
  free(work);

  if (neg_w) { /* then the correction for the negative w_i has to be evaluated */
    Q1 = (double *)calloc((size_t) nn * rank,sizeof(double)); 
    for (i=0;i<rank;i++) Q1[i * nn + i] = 1.0;
    left=1;tp=0;mgcv_qrqy(Q1,WX,tau,&nn,&rank,q,&left,&tp); /* Q from the QR decomposition */
    if (neg_w < rank+1) k = rank+1; else k = neg_w;
    IQ = (double *)calloc((size_t) k * rank,sizeof(double)); 
    for (i=0;i<neg_w;i++) { /* Copy the rows of Q corresponding to -ve w_i into IQ */
      p0 = IQ + i;p1 = Q1 + nind[i];
      for (j=0;j<rank;j++,p0+=k,p1+= nn) *p0 = *p1;
    }
    free(Q1); 
    /* Note that IQ may be zero padded, for convenience */
    Vt = (double *)calloc((size_t) rank * rank,sizeof(double));
    d = (double *)calloc((size_t)  rank,sizeof(double));
    mgcv_svd_full(IQ,Vt,d,&k,&rank); /* SVD of IQ */
    free(IQ);
    for (i=0;i<rank;i++) {
      d[i] = 1 - 2*d[i]*d[i];
      if (d[i]< - *rank_tol) { /* X'WX not +ve definite, clean up and abort */
        *n = -1;
        free(Vt);free(d);free(pivot);free(tau);free(nind);free(raw);free(z);free(WX);
        return;
      }
      if (d[i]<=0) d[i]=0.0; else d[i] = 1/d[i];
    }
    /* d now contains diagonal of diagonal matrix (I-2D^2)^{-1} (possibly pseudoinverse) */
  } else {Vt = d = (double *)NULL; }
  /* The -ve w_i correction is now complete */

  /* Now get the fitted values X \beta, *without* finding \beta */
  left=1;tp=1;mgcv_qrqy(z,WX,tau,&nn,&one,q,&left,&tp); /* z = Q'z */
  for (i=rank;i<nn;i++) z[i]=0.0;

  if (neg_w) { /* apply the correction factor for negative w_i terms */
     bt=0;ct=0;mgcv_mmult(y,Vt,z,&bt,&ct,&rank,&one,&rank); /* V' Q_1' z */
     for (i=0;i<rank;i++) y[i] *= d[i];
     bt=1;ct=0;mgcv_mmult(z,Vt,y,&bt,&ct,&rank,&one,&rank); /* V (I-2D^2)^{-1} V' Q_1' z */
  }

  for (i=0;i<rank;i++) y[i] = z[i];        /* y = Q'z, or corrected version */ 
  left=1;tp=0;mgcv_qrqy(z,WX,tau,&nn,&one,q,&left,&tp);
  for (i=0;i<*n;i++) eta[i] = z[i]/raw[i]; /* the linear predictor */

  for (*penalty=0.0,i=*n;i<nn;i++) *penalty += z[i]*z[i]; /* the penalty term */
  
  /* now find  \hat \beta = R^{-1}Q'z, which are needed if P-IRLS starts to diverge
     in order to be able to evaluate penalty on step reduction */
  
 

  /* now back substitute to find \hat \beta */  
  for (k=rank;k<*q;k++) z[k]=0.0; /* truncated parameters */
  for (k=rank-1;k>=0;k--) {
      for (xx=0.0,j=k+1;j < rank;j++) xx += WX[k + nn * j]*z[j];
      z[k] = (y[k] - xx)/WX[k + nn * k];
  }
  /* unpivot result into y */
  for (i=0;i< *q;i++) y[pivot[i]] = z[i];
 
  /*, should this be required... */
  free(z);free(WX);free(tau);free(pivot);free(raw);
  if (neg_w) { free(nind);free(d);free(Vt);}
}


