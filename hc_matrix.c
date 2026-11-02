/*

small linear algebra routines for the surface boundary condition
solve in hc_polsol

the system is at most 3 x 3. the matrix is stored as a[3][3], and n
(2 or 3) selects the upper left n x n block

hc_ludcmp_3x3 replaces a by its LU decomposition with partial
(row) pivoting, L has unit diagonal, and indx[] records the row
permutation. hc_lubksb_3x3 then solves A x = b for a right hand side
b, which is overwritten with the solution.

this is a plain implementation of Gaussian elimination with partial
pivoting (e.g. Golub & van Loan, Matrix Computations, Alg. 3.4.1),
it replaced an earlier version that was adapted from Numerical
Recipes

*/
#include "hc.h"

void hc_ludcmp_3x3(HC_PREC a[3][3],int n,int *indx)
{
  int i,j,k,imax;
  HC_PREC big,tmp,fac;

  for(k=0;k < n;k++){
    /*
       find the pivot row for column k
    */
    imax = k;
    big = fabs(a[k][k]);
    for(i=k+1;i < n;i++){
      tmp = fabs(a[i][k]);
      if(tmp > big){
	big = tmp;
	imax = i;
      }
    }
    if(big < HC_EPS_PREC){
      fprintf(stderr,"hc_ludcmp_3x3: singular matrix, column %i, pivot %g\n",
	      k,(double)big);
      for(i=0;i < n;i++){
	for(j=0;j < n;j++){
	  fprintf(stderr,"%g ",(double)a[i][j]);
	}
	fprintf(stderr,"\n");
      }
      exit(-1);
    }
    indx[k] = imax;
    /*
       swap rows k and imax
    */
    if(imax != k){
      for(j=0;j < n;j++){
	tmp = a[k][j];
	a[k][j] = a[imax][j];
	a[imax][j] = tmp;
      }
    }
    /*
       eliminate below the pivot, store the multipliers in the lower
       triangle
    */
    for(i=k+1;i < n;i++){
      fac = a[i][k] / a[k][k];
      a[i][k] = fac;
      for(j=k+1;j < n;j++){
	a[i][j] -= fac * a[k][j];
      }
    }
  }
}

void hc_lubksb_3x3(HC_PREC a[3][3],int n,int *indx,HC_PREC *b)
{
  int i,j;
  HC_PREC tmp,sum;
  /*
     apply the row permutation to b
  */
  for(i=0;i < n;i++){
    if(indx[i] != i){
      tmp = b[i];
      b[i] = b[indx[i]];
      b[indx[i]] = tmp;
    }
  }
  /*
     forward substitution, L y = P b, L has unit diagonal
  */
  for(i=1;i < n;i++){
    sum = b[i];
    for(j=0;j < i;j++){
      sum -= a[i][j] * b[j];
    }
    b[i] = sum;
  }
  /*
     back substitution, U x = y
  */
  for(i=n-1;i >= 0;i--){
    sum = b[i];
    for(j=i+1;j < n;j++){
      sum -= a[i][j] * b[j];
    }
    b[i] = sum / a[i][i];
  }
}
