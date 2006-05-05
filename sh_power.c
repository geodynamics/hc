#include "hc.h"

/* 

read in a spherical harmonic expansion in  Dahlen and Tromp normalization 
and computer power per degree and unit area

Thorsten Becker (twb@usc.edu)


usage: cat ab.sh_ana | sh_power [short_format, 0]

if short format is 0, will expect 

  type lmax shps ilayer nset zlabel ivec
  A-00 B-00
  A-10 B-10
  ...

format. if short format is set, will expect 

  lmax 
  A-00 B-00
  A-10 B-10
  ...

format



$Id: sh_power.c,v 1.5 2004/12/20 05:18:12 becker Exp $

*/

int main(int argc, char **argv)
{
  int type,lmax,shps,ilayer,nset,ivec,i,l;
  hc_boolean verbose = TRUE, short_format = FALSE ,binary = FALSE;
  float *power = NULL;
  HC_PREC fac[3] = {1.,1.,1.},zlabel;
  struct sh_lms *exp;
  if(argc>1){
    sscanf(argv[1],"%i",&i);
    if(i)
      short_format = TRUE;
  }
  fprintf(stderr,"%s: awaiting spherical harmonics expansion (%s) from stdin\n",
	  argv[0],short_format ? "short format" : "long format");
  while(sh_read_parameters(&type,&lmax,&shps,&ilayer,&nset,
			   &zlabel,&ivec,stdin,short_format,
			   binary,verbose)){
    fprintf(stderr,"%s: computing power per degree and unit area for lmax %i ivec: %i at z: %g\n",
	    argv[0],lmax,ivec,zlabel);
    /* 
       input and init 
    */
    sh_allocate_and_init(&exp,shps,lmax, type,ivec,verbose);
    sh_read_coefficients(exp,shps,-1,stdin,binary,fac,verbose);
    /* get space */
    hc_svecrealloc(&power,exp->lmaxp1 * shps,"sh_power");
    /* compute the powers */
    for(i=0;i<shps;i++)
      sh_compute_power_per_degree((exp+i),(power+i*exp->lmaxp1));
    for(l=0;l<=exp->lmax;l++){
      fprintf(stdout,"%5i ",l);
      for(i=0;i<shps;i++)
	fprintf(stdout,"%15.7e",power[l+i*exp->lmaxp1]);
      fprintf(stdout,"\n");
    }
    free(exp);
  }
  return 0;
}
