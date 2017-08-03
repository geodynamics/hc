#include "hc.h"
#include "gmt.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* 

read in two sets of spherical harmonics coefficients (stdin) and
compute linear correlation coefficients up to L = min(L1,L2), or
limited by input parameter

USAGE:

for total correlation

cat geoid.ab itg-hc-geoid.ab | sh_corr -1 0

for correlation per degree

cat geoid.ab itg-hc-geoid.ab | sh_corr -1 0 1

for correlation between l = 4 and 9

cat geoid.ab itg-hc-geoid.ab | sh_corr 9 0 0 4  


Thorsten Becker (twb@ig.utexas.edu)


*/

int main(int argc, char **argv)
{
  int type,lmax[2],llim,shps,ilayer,nset,ivec,i,l,lmin;
  /* 
     switches 
  */
  hc_boolean verbose = TRUE, short_format = FALSE , binary = FALSE, per_degree = FALSE;
  struct sh_lms *exp1,*exp2;
  HC_PREC zlabel,fac[3] = {1.,1.,1.};
  llim = -1;			/* max l */
  lmin = 1;			/* min l */
  if(argc > 1){
    if((strcmp(argv[1],"-h")==0)||(strcmp(argv[1],"--help")==0)||(strcmp(argv[1],"-help")==0))
      argc = -1000;
    else{
      sscanf(argv[1],"%i",&llim);
      if(llim == 0){
	fprintf(stderr,"%s: error: llim should not be zero\n",argv[0]);
	exit(-1);
      }
    }
  }
  if((argc > 5)|| (argc < 0)){
    fprintf(stderr,"usage: cat sh.1 sh.2 | %s [llim, %i] [short_format, %i] [per_degree, %i] [lmin, %i]\n",
	    argv[0],llim,short_format,per_degree,lmin);
    fprintf(stderr,"reads two spherical harmonic expansions from stdin and computes\n");
    fprintf(stderr,"linear correlation coefficients for L=min(L1,L2) if llim == -1\n");
    fprintf(stderr,"if llim > 0, will limit the maximum degree to min(llim,L1,L2)\n");
    fprintf(stderr,"short_format:\n\t0: expects regular format with long header\n");
    fprintf(stderr,"\t1: expects short format with only llim in header\n");
    fprintf(stderr,"per_degree: if set, will print out per degree, if not total correlation\n\n");
    fprintf(stderr,"lmin: will compute total correlation starting at lmin\n");
    exit(-1);
  }
  if(argc > 2){
    sscanf(argv[2],"%i",&i);
    if(i)
      short_format = TRUE;
  }
  if(argc > 3){
    sscanf(argv[3],"%i",&i);
    if(i)
      per_degree = TRUE;
  }
  if(argc > 4)
    sscanf(argv[4],"%i",&lmin);
  
  if(verbose)

    fprintf(stderr,"%s: waiting to read two spherical harmonic coefficients %s from stdin (use %s -h for help)\n",
	    argv[0],(short_format)?("(short format)"):("(long format)"),argv[0]);
  /* allocate two expansions */
  i=0;
  while(sh_read_parameters_from_stream(&type,(lmax+i),&shps,&ilayer,&nset,
				       &zlabel,&ivec,stdin,short_format,
				       binary,verbose)){
    if(ivec){
      fprintf(stderr,"%s: error, vector fields not implemented yet\n",argv[0]);
      exit(-1);
    }
    if(verbose)
      fprintf(stderr,"%s: expansion %i lmax %i type %i shps %i ...",
	      argv[0],i+1,lmax[i],type,shps);
    /* input and init */
    if(i==0){
      sh_allocate_and_init(&exp1,shps,lmax[i],type,ivec,verbose,0);
      sh_read_coefficients_from_stream(exp1,shps,-1,stdin,binary,fac,verbose);
    }else{
      sh_allocate_and_init(&exp2,shps,lmax[i],type,ivec,verbose,0);
      sh_read_coefficients_from_stream(exp2,shps,-1,stdin,binary,fac,verbose);
    }
    if(verbose)
      fprintf(stderr,"ok.\n");
    i++;
  }
  if(i!=2){
    fprintf(stderr,"%s: read error, expecting two and only two expansions\n",
	    argv[0]);
    exit(-1);
  }
  /* check bounds */
  if(llim < 0){
    llim = (lmax[0] < lmax[1])?(lmax[0]):(lmax[1]);
  }else{
    llim = (llim < lmax[0])?(llim):(lmax[0]);
    llim = (llim < lmax[1])?(llim):(lmax[1]);
  }
  if(lmin>llim){
    fprintf(stderr,"%s: error: lmin (%i) needs to be <=  llim (%i)\n",
	    argv[0],lmin,llim);
    exit(-1);
  }
  /*  */
  if(per_degree){
    if(verbose)
      fprintf(stderr,"%s: computing linear correlation per degree from %i to %i \n",argv[0],lmin,llim);
    for(l=lmin;l<=llim;l++)
      fprintf(stdout,"%5i %14.7e\n",l,(double)sh_correlation_per_degree(exp1,exp2,l,l));
  }else{
    if(verbose)
      fprintf(stderr,"%s: computing total linear correlation from %i to %i\n",argv[0],lmin,llim);
    fprintf(stdout,"%14.7e\n",(double)sh_correlation_per_degree(exp1,exp2,lmin,llim));
  }

  
  free(exp1);free(exp2);

  return 0;
}
