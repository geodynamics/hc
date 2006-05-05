
/* 

small matrix inversion routines used by hc_polsol of H & C code


lifted from numerical recipes, but changed the matrix addressing
a[0...2][0...2]

this was moderately well tested.


$Id: hc_matrix.c,v 1.10 2006/03/20 05:32:48 becker Exp $

*/

#ifdef COMPILE_TEST
/* 
   test routine 
*/
#include <stdio.h>
#include <math.h>
#define HC_PREC double
#define FALSE 0
#define HC_EPS_PREC 5e-15
void hc_ludcmp_3x3(HC_PREC [3][3],int *);
void hc_lubksb_3x3(HC_PREC [3][3],int *,HC_PREC *);
void main(void){
  HC_PREC amat[3][3],bvec[3];int i,j,indx[3];
  char fstring[10];
  /* read in A from stdin */
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
    fscanf(stdin,"%lf",&amat[i][j]);
  /* read in x from stdin */
  for(i=0;i<3;i++)
    fscanf(stdin,"%lf",&bvec[i]);
  /* solve A x = b, where b will be modified  */
  hc_ludcmp_3x3(amat,indx);hc_lubksb_3x3(amat,indx,bvec);
  for(i=0;i<3;i++)fprintf(stdout,"%g\n",bvec[i]);
}
#else
/* 

   normal compilation 

*/
#include "hc.h"
#endif

#define NR_TINY 1.0e-20;

void hc_ludcmp_3x3(HC_PREC a[3][3],int *indx)
{
  int i,imax=0,j,k;
  HC_PREC big,dum,sum,temp;
  HC_PREC vv[3];
  
  for (i=0;i < 3;i++) {
    big=0.0;
    for (j=0;j < 3;j++)
      if ((temp = fabs(a[i][j])) > big) 
	big=temp;
    if (fabs(big) < HC_EPS_PREC) {
      fprintf(stderr,"hc_ludcmp_3x3: singular matrix in routine, big: %g\n",
	      big);
      //hc_print_3x3(a,stderr);
      for(j=0;j<3;j++)
	fprintf(stderr,"%g %g %g\n",a[j][0],a[j][1],a[j][2]);
      exit(-1);
    }
    vv[i]=1.0/big;
  }
  for (j=0;j < 3;j++) {
    for (i=0;i < j;i++) {
      sum = a[i][j];
      for (k=0;k < i;k++) 
	sum -= a[i][k] * a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for (i=j;i < 3;i++) {
      sum=a[i][j];
      for (k=0;k < j;k++)
	sum -= a[i][k] * a[k][j];
      a[i][j]=sum;
      if ( (dum = vv[i]*fabs(sum)) >= big) {
	big=dum;
	imax=i;
      }
    }
    if (j != imax) {
      for (k=0;k < 3;k++) {
	dum = a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;
      }
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (fabs(a[j][j]) < HC_EPS_PREC) 
      a[j][j] = NR_TINY;
    if (j != 2) {
      dum=1.0/(a[j][j]);
      for (i=j+1;i < 3;i++) 
	a[i][j] *= dum;
    }
  }
}
#undef NR_TINY
void hc_lubksb_3x3(HC_PREC a[3][3], int *indx, HC_PREC *b)
{
  int i,ii=0,ip,j;
  HC_PREC sum;
  for (i=0;i < 3;i++) {
    ip = indx[i];
    sum = b[ip];
    b[ip]=b[i];
    if (ii)
      for (j=ii-1;j <= i-1;j++) 
	sum -= a[i][j]*b[j];
    else if (fabs(sum)>HC_EPS_PREC) 
      ii = i+1;
    b[i]=sum;
  }
  for (i=2;i>=0;i--) {
    sum=b[i];
    for (j=i+1;j < 3;j++) 
      sum -= a[i][j]*b[j];
    b[i] = sum/a[i][i];
  }
}

