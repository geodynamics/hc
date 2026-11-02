#include "hc.h"
//
// fourier transform routines as used by rick_sh routines
// based on Rick O'Connell's subroutines, which are modified 
// Numerical Recipes
//
// $Id: rick_fft_c.c,v 1.2 2006/01/22 02:46:15 becker Exp $
//

void rick_cs2ab(SH_RICK_PREC *rdata,int n)
{
  //
  //  Transforms spectral coefficients from cos-sin series to
  //    complex discrete fourier series. Function is real, and
  //    transformed by realft(rdata,n/2,1). Number of data points
  //    is n. Does not recover real component for frequency n/2.
  
  int i;
  SH_RICK_PREC en;
  en = (SH_RICK_PREC)n;
  rdata[0] *= en;
  en/=2.0;
  for(i=2;i<n;i++)
    rdata[i] *= en;
}


void 
rick_ab2cs (rdata, n)
SH_RICK_PREC *rdata;
int n;
{

//  Changes coefficients of complex spectrum of a real function
//   transformed by realft.f to real coefficients of a series
//   of C*cos(m*x)+S*sin(mx). Coefficients are ordered as
//   C(0),S(0),C(1),S(1),C(2),...,C(n/2-1),S(n/2-1). This loses
//   the real part of spectrum for frequency n/2.
//   The number of data points is n, The call to realft is
//    call realft(rdata,n/2,1)
//
  int i;
  SH_RICK_PREC en;
  en = 1.0/(SH_RICK_PREC)n;

  rdata[0] *= en;
  rdata[1] = 0.0;
  en *= 2.0;
  for(i=2;i<n;i++)
    rdata[i] *= en;
}



/*!

THIS ROUTINE TAKES NUMERICAL RECIPES 1...n,. SO CALL WITH (data-1)

!	Calculates the fourier transform of 2*N real data points.
!	Replaces data with the positive frequency half of the
!	complex fourier transform. The real parts of the first
!	and last frequency components are returned in data(1)
!	and data(2) (i.e. for frequencies of zero and N/2). The
!	other spectral components are given as complex pairs
!	in data(3),data(4) etc. The inverse transform is obtained
!	with ISIGN=-1, and dividing the data or result by N.
!	Calls routine four1(data,n,isign) for FFT.
!
*/

void 
rick_realft_nr (rdata, n, isign)
SH_RICK_PREC *rdata;
int n;
int isign;
{
  SH_RICK_PREC c1,c2,h1r,h1i,h2r,h2i;
  SH_RICK_HIGH_PREC theta,wi,wpi,wpr,wr,wtemp;
  int i,n2p3,ilim,i1,i2,i3,i4,n2;
  static int negunity = -1,unity = 1;

  theta = RICK_PI/(SH_RICK_HIGH_PREC)(n); 
  
  wr = 1.0;
  wi = 0.0;
  c1 = 0.5;
  
  /* offsets */
  n2 = 2*n;
  n2p3 = n2+3;    

  if (isign == 1) {
    c2 = -0.5;
    rick_four1_nr(rdata,n,unity);	/*four1 also wants 1..n */
    rdata[n2+1] = rdata[1];	
    rdata[n2+2] = rdata[2];	
  }

  else {
    c2 = 0.5;
    theta = -theta;
    rdata[n2+1] =  rdata[2];	/* rdata indices changed */
    rdata[n2+2] =  0.0;	
    rdata[2]=0.0;	
  }
  wtemp = sin(0.5 * theta);
  wpr = -2.0 * wtemp * wtemp;
  wpi = sin(theta);

  ilim = n/2 + 1;
  for (i=1;i <= ilim;i++) {		
    i1 = 2*i-1;
    i2 = i1+1;	
    i3 = n2p3 - i2;	
    i4 = i3+1;
    h1r =  c1*(rdata[i1] + rdata[i3]);
    h1i =  c1*(rdata[i2] - rdata[i4]);
    h2r = -c2*(rdata[i2] + rdata[i4]);
    h2i =  c2*(rdata[i1] - rdata[i3]);
    rdata[i1]= h1r+wr*h2r-wi*h2i;
    rdata[i2]= h1i+wr*h2i+wi*h2r;
    rdata[i3]= h1r-wr*h2r+wi*h2i;
    rdata[i4]=-h1i+wr*h2i+wi*h2r; /*end of index changes */
    wtemp=wr;
    wr=wr*wpr-wi*wpi+wr;
    wi=wi*wpr+wtemp*wpi+wi;
  }
  if (isign == 1) {
    rdata[2]=rdata[n2+1];	
  }else {
    rick_four1_nr(rdata,n,negunity);  /*Again, sending the rdata[0..n] array to four1 (which works in [1..n+1]) requires passing rdata-1 */
  }
}


/* 

CALL THIS ROUTINE 1...N FASHION, IE. WITH (RDATA-1) FROM REGULAR 
C


*/
void 
rick_four1_nr (rdata, nn, isign)
SH_RICK_PREC *rdata;
int nn;
int isign;
{
  //
  //	FFT routine from Numerical Recipes. Replaces data by
  //	its discrete fourier transform if isign=1, or by
  //	NN times its inverse transform if isign=-1. Array
  //	data is made up of NN complex numbers (2*NN pairs)
  //	and NN must be a power of 2. Spectral components
  //	are complex, and ordered from frequency zero to
  //	+-NN/2 to -1 in the standard fashion.
  //
  // local
  SH_RICK_HIGH_PREC  tempr,tempi;
  // this should be SH_RICK_HIGH_PREC precision locally, regardless
  SH_RICK_HIGH_PREC  wr,wi,wpr,wpi,wtemp,theta,temp2;
  
  int n,m,i,j,mmax,istep;
    
  n=2*nn;
  j=1;
  for(i=1;i <= n;i += 2){
    if(j > i){
      tempr = rdata[j];
      tempi = rdata[j+1];
      rdata[j] = rdata[i];
      rdata[j+1] = rdata[i+1];
      rdata[i]  = tempr;
      rdata[i+1]=tempi;
    }
    m=n/2;
    while((m >= 2) && (j > m)){
      j=j-m;
      m/=2;
    }
    j += m;
  }
  mmax=2;
  while(n > mmax){
    istep=2*mmax;
    theta = 6.28318530717959/(isign*mmax);
    temp2 = sin(0.5 * theta);
    wpr = -2.0 * temp2 * temp2;
    wpi = sin(theta);
    wr = 1.0;
    wi = 0.0;
    for(m=1;m <= mmax;m += 2){
      for(i=m;i <= n;i += istep){
	j = i + mmax;
	tempr=wr*rdata[j]-wi*rdata[j+1];
	tempi=wr*rdata[j+1]+wi*rdata[j];
	rdata[j]=rdata[i]-tempr;
	rdata[j+1]=rdata[i+1]-tempi;
	rdata[i]=rdata[i]+tempr;
	rdata[i+1]=rdata[i+1]+tempi;
      }
      wtemp=wr;
      wr=wr*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}


