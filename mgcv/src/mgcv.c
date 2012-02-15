/* Source code for mgcv.dll/.so multiple smoothing parameter estimation code,
suitable for interfacing to R 

Copyright (C) 2000-2005 Simon N. Wood  simon.wood@r-project.org

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

    /*#define WINDOWS*/ 

#ifdef WINDOWS
#include <windows.h>   /* For easy self window porting, without neading everything R.h needs */
#include <stdarg.h>
#else
#include <R.h>
#endif
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tprs.h"
#include "mgcv.h"
#include "matrix.h"
#include "qp.h"
#include "general.h"

#define round(a) ((a)-floor(a) <0.5 ? (int)floor(a):(int) floor(a)+1)


#ifdef WINDOWS
/* following routine provides a version of Rprintf that outputs to a file instead 
   of console, to facilitate debugging under windows.....
*/
typedef struct{
double a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20;
} _dodgy_print_storage;

void Rprintf(char *str,...)
{ _dodgy_print_storage ug; 
  FILE *f;
  va_list argptr;
  f=fopen("d:/simon/Rdump.txt","at");
  va_start(argptr,str);
  ug=va_arg(argptr,_dodgy_print_storage);
  fprintf(f,str,ug);
  va_end(argptr);
  fclose(f);
}
#endif

void ErrorMessage(char *msg,int fatal)

{ 
#ifdef WINDOWS
  MessageBox(HWND_DESKTOP,msg,"Info!",MB_ICONEXCLAMATION|MB_OK); 
#else  
if (fatal) error("%s",msg);else warning("%s",msg);
#endif
}

void infobox(char *msg)

{ 
#ifdef WINDOWS
MessageBox(HWND_DESKTOP,msg,"Info!",MB_ICONEXCLAMATION|MB_OK); 
#else
warning("%s",msg);
#endif
}

/* The following are some rather ancient routines used to set up an example
   additive model using regression (cubic) splines, via RGAMsetup(). */
void RUnpackSarray(int m,matrix *S,double *RS)
/* unpacks the R array RS into an array of matrices initialized to the correct dimensions 
   let kk = sum_{i=0}^k S[i].r*S[i].c
   Then the kth matrix starts at element kk of RS and stops at element k(k+1)
   ... let this extracted array be M. S[k].M[i][j]=M[i+S[k].r*j] - in this way we ensure that 
   M can be extracted straight to a matrix in R with 
   A<-matrix(M,S[k].r,S[k].c) 
*/ 
{ int start,i,j,k;
  start=0;
  for (k=0;k<m;k++)
  { for (i=0;i<S[k].r;i++) for (j=0;j<S[k].c;j++) S[k].M[i][j]=RS[start+i+S[k].r*j];
    start += S[k].r*S[k].c;
  }
}

void RPackSarray(int m,matrix *S,double *RS)
/* Packs an array of matrices S[] into an R array RS in the manner described in RUnpackSarray
*/
{ int start,i,j,k;
  start=0;
  for (k=0;k<m;k++)
  { for (i=0;i<S[k].r;i++) for (j=0;j<S[k].c;j++) RS[start+i+S[k].r*j]=S[k].M[i][j];
    start += S[k].r*S[k].c;
  }

}




/********** The following are from spline.c (rather plodding coding) ***********/

/* The next 4 functions are basis functions for the 1st derivative
   representation of a cubic spline */

double b0(x0,x1,x) double x0,x1,x;

/* multiplies function value at x0 */

{ double res,h,xx1;
  if (x<x0) return(1.0);
  if (x>x1) return(0.0);
  h=x1-x0;xx1=x-x1;
  res=2.0*(x-x0+0.5*h)*xx1*xx1/(h*h*h);
  return(res);
}


double b1(x0,x1,x) double x0,x1,x;

/* multiplies function value at x1 */

{ double res,h,xx0;
  if (x<x0) return(0.0);
  if (x>x1) return(1.0);
  h=x1-x0;xx0=x-x0;
  res= -2.0*(x-x1-0.5*h)*xx0*xx0/(h*h*h);
  return(res);
}

double d0(x0,x1,x) double x0,x1,x;

/* multiplies gradient at x0 */

{ double res,h,xx1;
  if (x<x0) res=x-x0;  /* before start of interval - use linear extrapolation */
  else 
  if (x>x1) res=0.0;  /* after end of interval - d0 plays no part */
  else                /* within interval */
  { h=x1-x0;xx1=x-x1;
    res=(x-x0)*xx1*xx1/(h*h);
  }
  return(res);
}

double d1(x0,x1,x) double x0,x1,x;

/* multiplies gradient at x1 */

{ double res,h,xx0;
  if (x<x0) res=0.0; /* d1 plays no part */
  else
  if (x>x1) res=x-x1; /* linear extrapolation after end of interval */
  else                /* in interval */
  { h=x1-x0;xx0=x-x0;
    res=xx0*xx0*(x-x1)/(h*h);
  }
  return(res);
}




matrix getD(h,nak) matrix h;int nak;

/* the matrix mapping the value of the spline to the gradients at the knots.
   nak is true for 'not-a-knot' end conditions at the early end, otherwise
   'natural' end conditions are used. If there are only 2 knots then the spline
   is taken as a straight line if only 1 a constant. */

{ long i,j,n;
  matrix T,D,Res;
  n=h.r+1;
  T=initmat(n,n);D=initmat(n,n);Res=initmat(n,n);
  for (i=0;i<n;i++) for (j=0;j<n;j++)
  { T.M[i][j]=0.0;D.M[i][j]=0.0;}
  if (n==1L)
  { Res.M[0][0]=0.0;
  } else
  if (n==2L)
  { Res.M[0][0]=Res.M[1][0]=-1.0/h.V[0];
    Res.M[0][1]=Res.M[1][1]=1.0/h.V[0];
  } else
  { for (i=0;i<n;i++) T.M[i][i]=2.0;
    for (i=1;i<n-1;i++)
    { T.M[i][i-1]=h.V[i]/(h.V[i]+h.V[i-1]);
      T.M[i][i+1]=1.0-T.M[i][i-1];
      D.M[i][i-1]= -3.0*T.M[i][i-1]/h.V[i-1];
      D.M[i][i+1]=3.0*T.M[i][i+1]/h.V[i];
      D.M[i][i]= -(D.M[i][i+1]+D.M[i][i-1]);
    }
    if (!nak)
    { T.M[0][1]=1.0;D.M[0][0]= -3.0/h.V[0];D.M[0][1]= -D.M[0][0];}
    else
    { T.M[0][1]=2.0*(h.V[0]+h.V[1])/h.V[1];
      D.M[0][0]= -2.0*(3.0*h.V[0]+2.0*h.V[1])/
		(h.V[0]*(h.V[0]+h.V[1]));
      D.M[0][2]=2.0*h.V[0]*h.V[0]/
      (h.V[1]*h.V[1]*(h.V[0]+h.V[1]));
      D.M[0][1]= -D.M[0][0]-D.M[0][2];
    }
    T.M[n-1][n-2]=1.0;D.M[n-1][n-2]= -3.0/h.V[n-2];
    D.M[n-1][n-1]= -D.M[n-1][n-2];
    invert(&T);
    matmult(Res,T,D,0,0);
  }
  freemat(T);freemat(D);
  return(Res);
}

void MonoCon(matrix *A,matrix *b,matrix *x,int control,double lower,double upper ) 

/* gets matrices A and b for constraints of the form Ay>=b ensuring monotonic
   change  of the cubic spline interpolating (x_i,y_i) where h_i=x_{i+1}-x_i 
   
   control indicates type of constraints:
   up=control/4 - 0 for decrease, 1 for increase
   lo=(control-up*4)/2 - 1 for lower bound, 0 no lower bound
   hi=(control-up*4-lo*2) - 1 for upper bound, 0 no upper bound
   control = 4*up+2*lo+hi
*/

{ long i,j,n;
  int up,lo,hi;
  double m;
  matrix h,D;
  h=initmat(x->r-1,1L);
  n=h.r;
  for (i=0;i<n;i++) h.V[i]=x->V[i+1]-x->V[i];
  D=getD(h,0);
  up=control/4;control=control%4;
  lo=control/2;control=control%2;
  hi=control;
  if (up) m= -1.0; else m=1.0; 
  (*A)=initmat(4*n+hi+lo,n+1);
  for (i=0;i<n;i++)
  { for (j=0;j<n+1;j++)
    { if (j==i)
      { A->M[i][j]=(D.M[i][j]+3.0/h.V[i])*m;   /**not certain of d.M update**/
	    A->M[i+n][j]=(D.M[i+1][j]+3.0/h.V[i])*m;
	    A->M[i+2*n][j]=m;
	    A->M[i+3*n][j]= -D.M[i][j]*m;
      } else
      if (j==(i+1))
      { A->M[i][j]=(D.M[i][j]-3.0/h.V[i])*m;
	    A->M[i+n][j]=(D.M[i+1][j]-3.0/h.V[i])*m;
	    A->M[i+2*n][j]= -m;
	    A->M[i+3*n][j]= -D.M[i][j]*m;
      } else
      { A->M[i][j]=D.M[i][j]*m;
	    A->M[i+n][j]=D.M[i+1][j]*m;
	    A->M[i+2*n][j]=0.0;
	    A->M[i+3*n][j]= -D.M[i][j]*m;
      }
    }
  }
  *b = initmat(A->r,1L);
  if (lo)
  { for (j=0;j<n+1;j++) A->M[4*n][j]=0.0;
    if (up) A->M[4*n][0]=1.0; else A->M[4*n][n]=1.0;
    b->V[4*n]=lower;
  }
  if (hi)
  { for (j=0;j<n+1;j++) A->M[4*n][j]=0.0;
    if (up) A->M[4*n+lo][n]=-1.0; else A->M[4*n+lo][0]=-1.0;
    b->V[4*n+lo]=upper;
  }
  freemat(D);
  freemat(h);
}


void tmap(matrix tm,matrix t,double time,int kill)
 /* to release static matrix allocation set kill to 1 otherwise 0 and
	     prepare for a new sequence of knot positions in t*/

/* tm maps values of a function at the t values contained in vector t to
   the value of a spline through those points at 'time' ;tgm does the same
   for the gradient of the spline */

{ static matrix D;static char first=1;
  matrix h;
  double **DM,*dum,*tmV,*DMi,*DMi1,d0v,d1v,x,xx0,xx1,xx02,xx12,h1,h2,h3,b0v,b1v,x0,x1;
  long i,k,tr;
  if (first)
  { first=0;h=initmat(t.r-1,1L);
    for (i=0L;i<t.r-1;i++) h.V[i]=t.V[i+1]-t.V[i];
    D=getD(h,0); /* time trajectories always have natural end conditions */
    freemat(h);
  }
  if (t.r==1L)
  { tm.V[0]=1.0;}
  else
  { DM=D.M;dum=t.V;
    i=0L;tr=t.r-2;dum++;
    while ((time > *dum)&&(i<tr)) {i++;dum++;}
    tr=t.r;tmV=tm.V;
    DMi=DM[i];DMi1=DM[i+1];
    x0=t.V[i];x1=t.V[i+1];
    x=time;
    h1=x1-x0;h2=h1*h1;h3=h2*h1;
    xx0=x-x0;xx1=x-x1;
    xx02=xx0*xx0;xx12=xx1*xx1;
   
    if (x<x0) 
    { d0v=xx0;  
      d1v=0.0;
      b0v=1.0;
      b1v=0.0;
    }
    else if (x>x1) 
    { d0v=0.0;  
      d1v=xx1;
      b0v=0.0;
      b1v=1.0;
    }else                /* within interval */
    { d0v=xx0*xx12/h2;d1v=xx02*xx1/h2;
      b0v=2.0*(xx0+0.5*h1)*xx12/h3;
      b1v= -2.0*(xx1-0.5*h1)*xx02/h3;
    }
    
    for (k=0;k<t.r;k++,tmV++,DMi++,DMi1++)
    { *tmV = *DMi * d0v + *DMi1 * d1v;
    }
    tm.V[i] += b0v;
    tm.V[i+1] += b1v;
   
  }
  if (kill)
  { first=1;
   freemat(D);
  }
}

void tmap2(matrix tm,matrix t,double time,int kill)
 /* to release static matrix allocation set kill to 1 otherwise 0 and
	     prepare for a new sequence of knot positions in t*/

/* tm maps values of a function at the t values contained in vector t to
   the value of a spline through those points at 'time' ;tgm does the same
   for the gradient of the spline */

{ static matrix D;static char first=1;
  matrix h;
  long i,k;
  if (first)
  { first=0;h=initmat(t.r-1,1L);
    for (i=0L;i<t.r-1;i++) h.V[i]=t.V[i+1]-t.V[i];
    D=getD(h,0); /* time trajectories always have natural end conditions */
    freemat(h);
  }
  if (t.r==1L)
  { tm.V[0]=1.0;}
  else
  { i=0L;while((time>t.V[i+1])&&(i<t.r-2)) i++;
    for (k=0;k<t.r;k++)
    tm.V[k]=D.M[i][k]*d0(t.V[i],t.V[i+1],time)+
	    D.M[i+1][k]*d1(t.V[i],t.V[i+1],time);
    tm.V[i]+=b0(t.V[i],t.V[i+1],time);
    tm.V[i+1]+=b1(t.V[i],t.V[i+1],time);
   
  }
  if (kill)
  { first=1;
    freemat(D);
  }
}

void getHBH(HBH,h,nak,rescale) matrix *HBH,h;int nak,rescale;

/* Generates the wiggliness measure matrix for vector h; nak=0 for natural
   end conditions or nak=1 to use the not a knot condition at the lower end;
   set rescale=1 to produce a measure rescaled for the unit interval, set to
   zero otherwise */

{ long n,i,j;
  matrix C,B,BI,H,hn;
  double interval=0.0;
  n=h.r;
  if (rescale)
  { for (i=0;i<h.r;i++) interval+=h.V[i];
    hn=initmat(h.r,1L);
    for (i=0;i<h.r;i++) hn.V[i]=h.V[i]/interval;
  } else hn=h;
  (*HBH)=initmat(n+1,n+1);
  if (!nak)
  { C=initmat(n-1,n+1);
    B=initmat(n-1,n-1);
    H=initmat(n-1,n+1);
    for (i=0;i<n-1;i++)
    { for (j=0;j<n-1;j++)
      { B.M[i][j]=0.0;
	     H.M[i][j]=0.0;
      }
      H.M[i][n-1]=0.0;
      H.M[i][n]=0.0;
    }
    for (i=0;i<n-1;i++)
    { B.M[i][i]=(hn.V[i]+hn.V[i+1])/3.0;
      H.M[i][i]=1.0/hn.V[i];
      H.M[i][i+1]= -1.0/hn.V[i]-1.0/hn.V[i+1];
      H.M[i][i+2]=1.0/hn.V[i+1];
    }
    for (i=0;i<n-2;i++)
    { B.M[i][i+1]=hn.V[i+1]/6.0;
      B.M[i+1][i]=hn.V[i+1]/6.0;
    }
    invert(&B);
    matmult(C,B,H,0,0);
    matmult((*HBH),H,C,1,0);
    freemat(C);freemat(B);freemat(H);
  } else
  { H=initmat(n,n+1);
    BI=initmat(n,n);B=initmat(n,n);
    for (i=0;i<H.r;i++) for (j=0;j<H.c;j++) H.M[i][j]=0.0;
    for (i=1;i<n;i++)
    { H.M[i][i-1]=1.0/hn.V[i-1];H.M[i][i]= -1.0/hn.V[i-1]-1.0/hn.V[i];
      H.M[i][i+1]=1.0/hn.V[i];
    }
    for (i=0;i<n;i++) for (j=0;j<n;j++)
    { BI.M[i][j]=0.0;B.M[i][j]=0.0;}
    for (i=1;i<n;i++)
    { B.M[i][i-1]=hn.V[i-1]/6.0;B.M[i][i]=(hn.V[i-1]+hn.V[i])/3.0;
      if (i<(n-1))
      { B.M[i][i+1]=hn.V[i]/6.0;
	BI.M[i][i+1]=B.M[i][i+1];
      }
      for (j=0;j<2;j++) BI.M[i][j+i-1]=B.M[i][j+i-1];
    }
    B.M[0][0]= -hn.V[1];B.M[0][1]=hn.V[0]+hn.V[1];
    B.M[0][2]= -hn.V[0];
    BI.M[0][0]=hn.V[0]/3.0;BI.M[0][1]=hn.V[0]/6.0;
    C=initmat(n,n);
    invert(&B);
    matmult(C,BI,B,0,0);
    matmult(BI,B,C,1,0);
    freemat(B);freemat(C);
    C=initmat(n,n+1);
    matmult(C,BI,H,0,0);
    matmult((*HBH),H,C,1,0);
    freemat(C);freemat(BI);freemat(H);
  }
  if (rescale) freemat(hn);
}


void getSmooth(S,x,rescale) matrix *S,x;int rescale;

/* gets a natural wiggliness measure for a spline with knots at the elements
   of vector x. Set rescale not zero to pretend that the domain is the unit
   interval. */

{ matrix h;
  long i;
  h=initmat(x.r-1L,1L);
  for (i=0;i<x.r-1;i++) h.V[i]=x.V[i+1]-x.V[i];
  getHBH(S,h,0,rescale);
  freemat(h);
}


void crspline(double *x,int n,int knots,matrix *X,matrix *S, matrix *C, matrix *xp,int control)

/* Sets up a cubic regression spline, by specifying a set of knot locations 
   spread evenly throughout the covariate values, and using cubic hermite
   polynomials to represent the spline.

   This is a very fast way of setting up 1-d regression splines, but the
   basis provided by tprs() will be better. This code gives the basis used
   in mgcv prior to version 0.6.

   The inputs are:
   x - the array of covariate values
   n - the length of x
   knots - the number of knots

   The outputs are:

   X - design matrix 
   S - wiggliness penalty matrix
   C - constraint matrix 
   xp - knot position vector - can also be supplied as input

   control 0 indicates full call as above
   control 1 indicates prediction call: xp *must* be supplied and only X 
             calculated

   so that a penalized regression spline could be fitted by minimising:

       ||Xb-y||^2 + \lambda b'Sb

   subject to Cb=0 if the sum of the values of the spline at the xp values 
   is to be zero - a constraint that is useful in GAM modelling, to ensure
   model identifiability.
*/

{ int j,i,k;
  matrix y,my;
  double dx,xx;
  /* sort x values into order and get list of unique x values .... */
  if (control==0)
  { if (xp->V[0]>=xp->V[1]) /* then knot positions have not been supplied */
    { y=initmat((long)n,1L);
      for (j=0;j<n;j++) y.V[j]=x[j];
      y.r=(long)n;
      sort(y); /* next reduce to list of unique values..... */
      k=0;for (i=0;i<n;i++) if (y.V[k]!=y.V[i]) { k++;y.V[k]=y.V[i];} y.r=(long)k+1;
      dx=(y.r-1)/(knots-1.0);
      /* now place the knots..... */
      xp->V[0]=y.V[0];
      for (i=1;i<knots-1;i++)  /* place knots */
      { xx=dx*i;
        k=(int)floor(xx);
        xx -= k;
        xp->V[i]=(1-xx)*y.V[k]+xx*y.V[k+1];
      } 
      xp->V[knots-1]=y.V[y.r-1];
      freemat(y);
    }
 
    /* create the wiggliness measure matrix...... */
    getSmooth(S,*xp,0);
 
    /* create the constraint matrix ...... */
    *C=initmat(1L,(long)knots);
    for (i=0;i<knots;i++) C->M[0][i]=1.0;
    /* and finally, create the design matrix ...... */
  }
  *X=initmat((long)n,xp->r);
  my=initmat(xp->r,1L);
  for (j=0;j<n;j++)
  { tmap(my,*xp,x[j],0);
    for (i=0;i<my.r;i++) X->M[j][i]=my.V[i];
  }   
  tmap(my,*xp,x[0],1); /* kill matrix allocation in tmap */
 
  freemat(my);
}

void construct_cr(double *x,int *nx,double *k,int *nk,double *X,double *S,double *C,int *control)

/* Routine to be called from R to set up a cubic regression spline basis given
   x - array of x values (un-ordered)
   nx - number of x values
   k - array of/for knot locations. Should be in ascending order - if first two are not then 
       locations are generated automatically
   nk - number of knots
   The routine returns the knot locations in k plus
   X - the model matrix       n by nk
   S - the penalty matrix     nk by nk
   C - the constraint matrix  1 by nk
   
   control is an array of control constant:
   0 indicates full set up as above
   1 indicates prediction call, in which case k *must* be supplied and only X is returned
   
*/

{ matrix Xm,Sm,Cm,xp;
  int i;
  xp=initmat((long)*nk,1L);
  for (i=0;i<xp.r;i++) xp.V[i]=k[i];
  /* following initializes Xm and also Sm and Cm if *control ==0 */
  crspline(x,*nx,*nk,&Xm,&Sm,&Cm,&xp,*control);
  for (i=0;i<xp.r;i++) k[i]=xp.V[i];
  RArrayFromMatrix(X,Xm.r,&Xm);
  freemat(Xm);freemat(xp);
  if (*control==0)
  { RArrayFromMatrix(S,Sm.r,&Sm);
    RArrayFromMatrix(C,Cm.r,&Cm);
    freemat(Sm);freemat(Cm);
  }
} 



void MinimumSeparation(double *gx,double *gy,int *gn,double *dx,double *dy, int *dn,double *dist)
/* For each point gx[i],gy[i] calculates the minimum  Euclidian distance to a point in dx[], dy[].
   These distances are stored in dist. 
*/

{ double sep,xx,yy,*dum,*xdum,*ydum;
  int n,m;
  n = *gn;m = *dn;
  for (dum=dist;dum < dist + n; dum++,gx++,gy++)
  { xx= *gx - *dx;yy = *gy - *dy;*dum = xx*xx + yy*yy; /* first separation */
    for (xdum=dx+1,ydum=dy+1;xdum < dx + m;xdum++,ydum++)
    { xx= *gx - *xdum;yy = *gy - *ydum;sep = xx*xx + yy*yy; /* subsequent separations */
      if (sep < *dum) *dum = sep;
    }
    *dum = sqrt(*dum);
  }
} 



void RuniqueCombs(double *X,int *ind,int *r, int *c)

/* X is a matrix. This routine finds its unique rows and strips out the 
   duplicates. This is useful for finding out the number of unique covariate
   combinations present in a set of data. */

{ matrix B,Xd;
  int i,*ind1;
  B=Rmatrix(X,(long)(*r),(long)(*c));
  Xd=initmat(B.r,B.c+1);
  Xd.c--;mcopy(&B,&Xd);freemat(B);Xd.c++;
  for (i=0;i<Xd.r;i++) Xd.M[i][Xd.c-1]=(double)i;
  ind1=Xd_strip(&Xd);
  for (i=0;i<*r;i++) ind[i] = ind1[i]; /* copy index for return */
  Xd.c--; /* hide index array  */
  RArrayFromMatrix(X,Xd.r,&Xd);  /* NOTE: not sure about rows here!!!! */
  *r = (int)Xd.r; 
  freemat(Xd);free(ind1);
#ifdef MEM_CHECK
  dmalloc_log_unfreed();  dmalloc_verify(NULL);
#endif 
}

void RMonoCon(double *Ad,double *bd,double *xd,int *control,double *lower,double *upper,int *n)
/* obtains coefficient matrices for imposing monotoicity (and optionally bounds) on a 
   cubic regression spline with n knots located as specified in xd. 
   
   control indicates type of constraints:
   up=control/4 - 0 for decrease, 1 for increase
   lo=(control-up*4)/2 - 1 for lower bound, 0 no lower bound
   hi=(control-up*4-lo*2) - 1 for upper bound, 0 no upper bound
   control = 4*up+2*lo+hi
   
   lower and upper are the bounds to impose (ignored if control doesn't
   indicate that they should be used).

   Ad will have 4(n-1)+lo+hi rows and n columns
   bd will have 4(n-1)+lo+hi rows

*/ 

{ int i;
  matrix x,A,b;
  x=initmat((long)*n,1L);
  for (i=0;i<x.r;i++) x.V[i]=xd[i];
  MonoCon(&A,&b,&x,*control,*lower,*upper); 
  RArrayFromMatrix(Ad,A.r,&A);
  RArrayFromMatrix(bd,b.r,&b);
 
  freemat(x);freemat(A);freemat(b);  
#ifdef MEM_CHECK
  dmalloc_log_unfreed();  dmalloc_verify(NULL);
#endif 
}



void RQT(double *A,int *r,int*c)

/* Obtains the QT decomposition of matrix A (stored according to R conventions)
   AQ=[0,T] where T is reverse lower triangular (upper left is zero). r<c and 
   first c-r columns of Q are basis vectors for the null space of A 
   (Q orthogonal). Let this null space basis be Z. It is actually stored as 
   a series of r Householder rotations over the rows of A. Let u_i be the ith
   row of A (A[i,], i>=1) then the last i-1 elements of u_i are zero, while if 
   H_i=(I-u_i u_i') then Q=H_1 H_2 H_3 ...H_r.
   
   The main purpose of this routine *was* to provide a suitable representation 
   of the null space of any equality constraints on the problem addressed by 
   mgcv(). So if the constraints are Cp=0, RQT() was called to get an 
   appropriate null space basis in A. The non-obvious representation usually 
   saves much computing, since there are usually few constraints, resulting in 
   a high dimensional null space - in this case the Householder representation
   is very efficient. 

   However, the current version of mgcv() expects to get the constraint matrix    itself and not the null space. 
*/

{ matrix Q,B;
  B=Rmatrix(A,(long)(*r),(long)(*c));
  Q=initmat(B.r,B.c);
  QT(Q,B,0);
  RArrayFromMatrix(A,(long)(*r),&Q);
  freemat(Q);freemat(B); 
#ifdef MEM_CHECK
  dmalloc_log_unfreed();  dmalloc_verify(NULL);
#endif
}







void RprintM(matrix *A)

{ 
#ifdef WINDOWS
#else
int i,j;
  if (A->c==1L) 
  { for (i=0;i<A->r;i++) 
   Rprintf("%8.3g ",A->V[i]);Rprintf("\n");
  } 
else
 for (i=0;i<A->r;i++)
     { for (j=0;j<A->c;j++) Rprintf("%8.3g ",A->M[i][j]);Rprintf("\n");}
#endif
}

void  RPCLS(double *Xd,double *pd,double *yd, double *wd,double *Aind,double *bd,
            double *Afd,double *Hd,double *Sd,
            int *off,int *dim,double *theta, int *m,int *nar)

/* Interface routine for PCLS the constrained penalized weighted least squares solver.
   nar is an array of dimensions. Let:
   n=nar[0] - number of data
   np=nar[1] - number of parameters
   nai=nar[2] - number of inequality constraints
   naf=nar[3] - number of fixed constraints
   getH=nar[4] - 0 for no hat matrix, 1 to produce one. 
   
   Problem to be solved is:

   minimise      ||W^0.5 (y - Xp)||^2 + p'Bp
   subject to    Ain p >= b  & Af p = "constant"

   where B = \sum_{i=1}^m \theta_i S_i and W=diag(w)

   - in fact S_i are not stored whole - rather the smallest non-zero sub-matrix of each S_i is 
   stored in a densely packed form in S[]: see routines RpackSarray() and RUnpackSarray() for 
   details of the sub-matrix packing. off[i],off[i] is the location within the full S_i to
   insert the sub-matrix actually stored which is of dimension dim[i] by dim[i].

   W = diag(w) 

   on exit p contains the best fit parameter vector. 

*/
{ matrix y,X,p,w,Ain,Af,b,H,*S;
  int n,np,i,*active;
 
  np=nar[1];n=nar[0];
  /* unpack from R into matrices */
  X=Rmatrix(Xd,(long)n,(long)np);
  p=Rmatrix(pd,(long)np,1L);
  y=Rmatrix(yd,(long)n,1L);
  w=Rmatrix(wd,(long)n,1L);
  if (nar[2]>0) Ain=Rmatrix(Aind,(long)nar[2],(long)np); else Ain.r=0L;
  if (nar[3]>0) Af=Rmatrix(Afd,(long)nar[3],(long)np); else Af.r=0L;
  if (nar[2]>0) b=Rmatrix(bd,(long)nar[2],1L);else b.r=0L;
 
  if (*m) S=(matrix *)calloc((size_t) *m,sizeof(matrix));
  else S=&H; /* avoid spurious compiler warning */
  for (i=0;i< *m;i++) S[i]=initmat((long)dim[i],(long)dim[i]);
  RUnpackSarray(*m,S,Sd);
  
  if (nar[4]) H=initmat(y.r,y.r); else H.r=H.c=0L;
  active=(int *)calloc((size_t)(p.r+1),sizeof(int)); /* array for active constraints at best fit active[0] will be  number of them */
  /* call routine that actually does the work */
 
  PCLS(&X,&p,&y,&w,&Ain,&b,&Af,&H,S,off,theta,*m,active);

  /* copy results back into R arrays */ 
  for (i=0;i<p.r;i++) pd[i]=p.V[i];
 
  if (H.r) RArrayFromMatrix(Hd,H.r,&H);
  /* clear up .... */
  free(active);
 
  for (i=0;i< *m;i++) freemat(S[i]);
  if (*m) free(S);
 
  freemat(X);freemat(p);freemat(y);freemat(w);
  if (H.r) freemat(H);
  if (Ain.r) freemat(Ain);
  if (Af.r) freemat(Af);
  if (b.r) freemat(b);
#ifdef MEM_CHECK
  dmalloc_log_unfreed();  dmalloc_verify(NULL);
#endif
}

/*********************************************************************************************/
/* Bug fix and revision record:

1. 20/10/00: Knot placement method in GAMsetup() modified. Previous method had an error, so 
   that when df for a term was close to the number of data, a non-existent covariate value
   (i.e. out of array bound). New code also yields more regular placement, and now deals with 
   repeat values of covariates.
3. 5/1/01: Modified RGAMsetup(), GAMsetup(), gam_map() and RGAMpredict() so that nsdf is now
   total number of non-spline parameters including any constant. Hence R code must now provide 
   column for constant explicitly.
4. 5/1/01: fixed bug in RGAMpredict - standard errors of parametric components of linear predictor
   were wrongly calculated.
5. 30/5/01: GAMsetup re-organised to ease introduction of new tprs basis and multi-dimensional 
   smooths
6. 10/2001: b0,b1,d0,d1 modified so that extrapolation is linear beyond ends of spline as
   it should be for a natural spline.
7. 31/10/01: mgcv.c covariance and edf calculations made more robust - for poorly conditioned 
   cases choleski can fail in calculation of covariance matrix: in these cases use svd instead.
8. 2/11/01: RGAMpredict, RGAMsetup, GAMsetup, and gam_map modified to take array of penalty
   orders, p_order. This allows user explicit control of the order of penalty in spline terms, 
   while still supporting autoselection when p_order[i]==0. 
9. 9/11/01: RGAMpredict modified to allow 5th control option, which returns a matrix mapping 
            params to l.p. vector
10. 9/11/01: New routine RPackSArray and RUnpackSarray so that storage of arrays of penalty 
             matrices is not so wasteful.
11. 12/11/01: UZ and Xu now packed efficiently using above 2 routines.
12. 12/11/01: Routine RPCLS added for solving linearly constrained penalized least squares problems 
              by quadratic programming
13. 13/11/01: Routine RMonoCon added for finding monotonic constraint matrices for cubic regression 
              splines.
14. 5/2/02: GAMsetup and RGAMsetup modified to deal with "by" variables - covariates that multiply a 
            whole smooth term. The centering conditions have not been changed.
15. 6/9/02: Slight modification to gam_map() - terms are not calculated if corresponding by variable 
            is zero. This can save flops in fairly advanced use (e.g. posum package)
16.23/10/02: mgcv modified in order to check that Tr(A) calculations are sensible, and that termwise 
             effective degrees of freedom are calculated correctly. The problem arises with ill-conditioned 
             models when an inversion required for the term-wise effective degrees of freedom can 
             become unstable.
17.23/10/02: Bug in TrA calculation when smoothing parameters supplied. X'X used in place of X'WX - fixed.

18. 24/1/04: RGAMpredict, RGAMsetup, GAMsetup and gam_map deleted, to make way for a more object oriented
             and modular approach to model setup and prediction, based on "smooth objects". Constructor and 
             prediction code added instead.
*/

