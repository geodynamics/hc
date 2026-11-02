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

/*

banded LU factorization with partial pivoting and the corresponding
solve, for the multiple shooting poloidal solver (hc_polsol_ms.c)

storage follows the LAPACK general band convention: the matrix A
(n x n) with kl sub- and ku super-diagonals is kept in ab with
leading dimension ldab = 2*kl+ku+1, column major, with

   A(i,j) at ab[(kl+ku+i-j) + j*ldab]

the extra kl rows at the top of ab hold the fill-in from pivoting,
and have to be zero on input (HC_BAND_INDEX gives the storage
location). ipiv[n] returns the row interchanges. return value is 0
on success, j+1 if the pivot in column j is zero.

*/
int hc_band_lu_factor(int n,int kl,int ku,HC_PREC *ab,int *ipiv)
{
  int j,i,c,km,jp,cmax,ldab;
  HC_PREC big,tmp,fac;

  ldab = 2*kl+ku+1;
  for(j=0;j < n;j++){
    km = kl;
    if(km > n-1-j){
      km = n-1-j;
    }
    /*
       pivot search in column j, rows j...j+km
    */
    jp = j;
    big = fabs(ab[HC_BAND_INDEX(j,j,kl,ku,ldab)]);
    for(i=j+1;i <= j+km;i++){
      tmp = fabs(ab[HC_BAND_INDEX(i,j,kl,ku,ldab)]);
      if(tmp > big){
	big = tmp;
	jp = i;
      }
    }
    ipiv[j] = jp;
    if(big == 0.0){
      return j+1;
    }
    /*
       the row interchange affects columns j...j+kl+ku
    */
    cmax = j+kl+ku;
    if(cmax > n-1){
      cmax = n-1;
    }
    if(jp != j){
      for(c=j;c <= cmax;c++){
	tmp = ab[HC_BAND_INDEX(j,c,kl,ku,ldab)];
	ab[HC_BAND_INDEX(j,c,kl,ku,ldab)]  = ab[HC_BAND_INDEX(jp,c,kl,ku,ldab)];
	ab[HC_BAND_INDEX(jp,c,kl,ku,ldab)] = tmp;
      }
    }
    /*
       multipliers
    */
    fac = 1.0/ab[HC_BAND_INDEX(j,j,kl,ku,ldab)];
    for(i=j+1;i <= j+km;i++){
      ab[HC_BAND_INDEX(i,j,kl,ku,ldab)] *= fac;
    }
    /*
       rank one update of the trailing band
    */
    for(c=j+1;c <= cmax;c++){
      tmp = ab[HC_BAND_INDEX(j,c,kl,ku,ldab)];
      if(tmp != 0.0){
	for(i=j+1;i <= j+km;i++){
	  ab[HC_BAND_INDEX(i,c,kl,ku,ldab)] -= ab[HC_BAND_INDEX(i,j,kl,ku,ldab)] * tmp;
	}
      }
    }
  }
  return 0;
}
/*
   solve A x = b given the factorization from hc_band_lu_factor, b is
   overwritten with x
*/
void hc_band_lu_solve(int n,int kl,int ku,HC_PREC *ab,int *ipiv,HC_PREC *b)
{
  int j,i,km,imin,ldab;
  HC_PREC tmp;

  ldab = 2*kl+ku+1;
  /*
     forward: apply interchanges and L
  */
  for(j=0;j < n;j++){
    if(ipiv[j] != j){
      tmp = b[j];
      b[j] = b[ipiv[j]];
      b[ipiv[j]] = tmp;
    }
    km = kl;
    if(km > n-1-j){
      km = n-1-j;
    }
    for(i=j+1;i <= j+km;i++){
      b[i] -= ab[HC_BAND_INDEX(i,j,kl,ku,ldab)] * b[j];
    }
  }
  /*
     backward: U, which has kl+ku super diagonals after pivoting
  */
  for(j=n-1;j >= 0;j--){
    b[j] /= ab[HC_BAND_INDEX(j,j,kl,ku,ldab)];
    imin = j-kl-ku;
    if(imin < 0){
      imin = 0;
    }
    for(i=imin;i < j;i++){
      b[i] -= ab[HC_BAND_INDEX(i,j,kl,ku,ldab)] * b[j];
    }
  }
}
