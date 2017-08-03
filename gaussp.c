#include "hc.h"

/* 
given an N of points, return the Gauss point locations and weights in latitude
from -1 ... 1

Thorsten Becker (twb@ig.utexas.edu)

usage: 

gaussp N

*/

int main(int argc, char **argv)
{
  int nlat,lmax,i;
  HC_PREC *x, *w;
  const HC_PREC fac = 180./M_PI;
  /* 
     command line parameters
  */
  if(argc != 2){
    fprintf(stderr,"%s lmax\nprints Gauss points in latitude and weights\n",argv[0]);
    exit(-1);
  }
  sscanf(argv[1],"%i",&lmax);
  nlat = lmax+1;
  fprintf(stderr,"%s: L_max: %i corresponds to %i points in latitude\n",
	  argv[0],lmax,nlat);
  
  hc_vecalloc(&x,nlat,"gaussp");
  hc_vecalloc(&w,nlat,"gaussp");
  
  /* spaced in cos(theta) space */
  rick_gauleg(-1.0,1.0,x,w,nlat);
  fprintf(stderr,"%s: output is latitude[deg] cos(theta) weight\n",argv[0]);

  for(i=0;i<nlat;i++)
    fprintf(stderr,"%12.8f %14.10e %14.10e\n",90-acos(x[i])*fac,x[i],w[i]);
  
  free(x);free(w);
  return 0;
}
