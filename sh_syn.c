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
  int type,lmax,shps,ilayer,nset,ivec,i,j,npoints,nphi,ntheta;
  /* 
     switches 
  */
  hc_boolean verbose = TRUE, short_format = FALSE ,short_format_ivec = FALSE ,binary = FALSE;
  double regular_basis = 0;
  /*  */
  float *data,*theta,*phi;
  /* spacing for irregular output */
  double dphi,x,y;
  HC_PREC fac[3] = {1.,1.,1.},zlabel;
  SH_RICK_PREC *dummy;
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
    sscanf(argv[2],"%i",&i);
    if(i)
      short_format_ivec = TRUE;
  }
  if(argc > 3)
    sscanf(argv[3],"%lf",&regular_basis);

  if((argc > 4)||(argc<0)){
    fprintf(stderr,"usage: %s [short_format, %i] [short_ivec,%i] [regular_basis,%g]\n",
	    argv[0],short_format,short_format_ivec,regular_basis);
    fprintf(stderr,"short_format:\n\t0: expects regular format with long header\n");
    fprintf(stderr,"\t1: expects short format with only lmax in header\n\n");
    fprintf(stderr,"short_ivec:\n\t0: for short format, expect AB for scalar expansion\n");
    fprintf(stderr,"\t1: for short format, expect poloidal toroidal AP BP AT BT for vector expansion\n\n");
    fprintf(stderr,"regular_basis:\n\t0: use Gauss latitudes and FFT divided longitudes dependening on lmax\n");
    fprintf(stderr,"\t>0 use even spacing with regular_basis degree inrecrments on the globe\n\n");
    fprintf(stderr,"The output format will depend on the type of SH input.\n");
    fprintf(stderr,"\tfor scalara: lon lat scalar if a single SH is read in, else lon lat zlabel scalar.\n");
    fprintf(stderr,"\tfor vectors: lon lat v_theta v_phi if a single SH is read in, else lon lat zlabel v_theta v_phi.\n\n\n");
    exit(-1);
  }
  if(verbose)
    fprintf(stderr,"%s: waiting to read spherical harmonic coefficients from stdin\n",
	    argv[0]);
  while(sh_read_parameters_from_file(&type,&lmax,&shps,&ilayer,&nset,
				     &zlabel,&ivec,stdin,short_format,
				     binary,verbose)){
    if(short_format_ivec){
      ivec = 1;
      shps = 2;
    }
    if(verbose)
      fprintf(stderr,"%s: converting lmax %i ivec: %i at z: %g\n",
	      argv[0],lmax,ivec,zlabel);

    /* input and init */
    sh_allocate_and_init(&exp,shps,lmax,type,ivec,verbose,((regular_basis>0)?(1):(0)));
    sh_read_coefficients_from_file(exp,shps,-1,stdin,binary,fac,verbose);
    if(regular_basis > 0){
      /* 
	 irregular basis output 
      */
      dphi = regular_basis;
      if(verbose)
	fprintf(stderr,"sh_syn: using regular spaced grid with %g deg spacing\n",dphi);
      /*  */
      dphi = DEG2RAD(dphi);
      nphi = (TWOPI-dphi)/dphi + 1;
      ntheta = nphi/2;
      npoints = nphi * ntheta;
      /*  */
      hc_svecalloc(&phi,nphi,"sh_shsyn");
      hc_svecalloc(&theta,ntheta,"sh_shsyn");
      for(x=0,i=0;i < nphi;i++,x+=dphi)
	phi[i] = x;
      for(y=dphi/2,j=0;j<ntheta;y+=dphi,j++)
	theta[j] = y;
      hc_svecalloc(&data,npoints * shps,"sh_shsyn data");
      /* compute the expansion */
      sh_compute_spatial_irreg(exp,ivec,FALSE,&dummy,
			       theta,ntheta,phi,nphi,data,
			       verbose,FALSE);
      /* output */
      sh_print_irreg_spatial_data_to_file(exp,shps,data,
					  (nset>1)?(TRUE):(FALSE),
					  zlabel, theta,ntheta,
					  phi,nphi,stdout);
      
      

    }else{			/* use the built in spatial basis (Gaussian) */
      /* expansion */
      hc_svecalloc(&data,exp[0].npoints * shps,"sh_syn");
      sh_compute_spatial(exp,ivec,FALSE,&dummy,data,verbose);
      /* output */
      sh_print_spatial_data_to_file(exp,shps,data,
				    (nset>1)?(TRUE):(FALSE),
				    zlabel,stdout);
    }


    free(exp);free(data);
  }

  return 0;
}
