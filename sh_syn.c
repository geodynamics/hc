#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "hc.h"

/* 

read in spherical harmonics coefficients (stdin) and expand to spatial basis (stdout)

Thorsten Becker (twb@usc.edu)


$Id: sh_syn.c,v 1.6 2006/01/22 01:11:34 becker Exp $

*/

int main(int argc, char **argv)
{
  int type,lmax,shps,ilayer,nset,ivec,i;
  hc_boolean verbose = TRUE, short_format = FALSE ,short_format_ivec = FALSE ,binary = FALSE;
  float *data;
  HC_PREC fac[3] = {1.,1.,1.},zlabel;
  double *dbl_dummy;
  struct sh_lms *exp;
  if(argc > 1){
    if((strcmp(argv[1],"-h")==0)||(strcmp(argv[1],"--help")==0)||(strcmp(argv[1],"-help")==0))
      argc = -1000;
    else{
      sscanf(argv[1],"%i",&i);
      if(i)
	short_format = TRUE;
    }
  } 
  if(argc > 2){
    sscanf(argv[1],"%i",&i);
    if(i)
      short_format_ivec = TRUE;
  }
  if((argc > 3)||(argc<0)){
    fprintf(stderr,"usage: %s [short_format, %i] [short_ivec,%i]\n",
	    argv[0],short_format,short_format_ivec);
    exit(-1);
  }
  fprintf(stderr,"%s: waiting to read spherical harmonic coefficients from stdin\n",
	  argv[0]);
  while(sh_read_parameters_from_file(&type,&lmax,&shps,&ilayer,&nset,
				     &zlabel,&ivec,stdin,short_format,
				     binary,verbose)){
    if(short_format_ivec){
      ivec = 1;
      shps = 2;
    }
    fprintf(stderr,"%s: converting lmax %i ivec: %i at z: %g\n",
	    argv[0],lmax,ivec,zlabel);

    /* input and init */
    sh_allocate_and_init(&exp,shps,lmax,type,ivec,verbose);
    sh_read_coefficients_from_file(exp,shps,-1,stdin,binary,fac,verbose);
    /* expansion */
    hc_svecalloc(&data,exp[0].npoints * shps,"sh_shsyn");
    sh_compute_spatial(exp,ivec,FALSE,&dbl_dummy,data,verbose);
    /* output */
    sh_print_spatial_data_to_file(exp,shps,data,(nset>1)?(TRUE):(FALSE),zlabel,stdout);
    free(exp);free(data);
  }

  return 0;
}
