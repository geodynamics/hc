#include "hc.h"
#include "gmt.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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
  hc_boolean verbose = TRUE, short_format = FALSE ,short_format_ivec = FALSE ,binary = FALSE,
    regular_basis = FALSE;
  double w,e,s,n,dx,dy;
  /*  */
  float *data,*theta,*phi;
  /* spacing for irregular output */
  double dphi,x,y,dtheta;
  HC_PREC fac[3] = {1.,1.,1.},zlabel;
  SH_RICK_PREC *dummy;
  struct sh_lms *exp;
  dx = 1.0;
  w=0;e=360.;s=-90;n=90;
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
  if(argc > 3){
    sscanf(argv[3],"%lf",&w);
    regular_basis = TRUE;
  }
  if(argc > 4)
    sscanf(argv[4],"%lf",&e);
  if(argc > 5)
    sscanf(argv[5],"%lf",&s);
  if(argc > 6)
    sscanf(argv[6],"%lf",&n);
   if(argc > 7)
    sscanf(argv[7],"%lf",&dx);
  if(argc > 8)
    sscanf(argv[8],"%lf",&dy);
  else
    dy = dx;
  if((argc > 9)|| (argc < 0)){
    fprintf(stderr,"usage: %s [short_format, %i] [short_ivec, %i] [w, %g] [e, %g] [s, %g] [n, %g] [dx, %g] [dy, dx] (in that order)\n",
	    argv[0],short_format,short_format_ivec,w,e,s,n,dx);
    fprintf(stderr,"short_format:\n\t0: expects regular format with long header\n");
    fprintf(stderr,"\t1: expects short format with only lmax in header\n\n");
    fprintf(stderr,"short_ivec:\n\t0: for short format, expect AB for scalar expansion\n");
    fprintf(stderr,"\t1: for short format, expect poloidal toroidal AP BP AT BT for vector expansion\n\n");
    fprintf(stderr,"w,e,...\n\tif none of those are set, will use Gauss latitudes and FFT divided longitudes dependening on lmax\n");
    fprintf(stderr,"\tif set, will switch to regular spaced output with -Rw/e/s/n -Idx/dy type output\n\n");
    fprintf(stderr,"The output format will depend on the type of SH input.\n");
    fprintf(stderr,"\tfor scalara: lon lat scalar if a single SH is read in, else lon lat zlabel scalar.\n");
    fprintf(stderr,"\tfor vectors: lon lat v_theta v_phi if a single SH is read in, else lon lat zlabel v_theta v_phi.\n\n\n");
    exit(-1);
  }
  if(verbose)
    fprintf(stderr,"%s: waiting to read spherical harmonic coefficients from stdin (use %s -h for help)\n",
	    argv[0],argv[0]);
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
    sh_allocate_and_init(&exp,shps,lmax,type,ivec,verbose,((regular_basis)?(1):(0)));
    sh_read_coefficients_from_file(exp,shps,-1,stdin,binary,fac,verbose);
    if(regular_basis){
      /* 
	 irregular basis output 
      */
      if(verbose)
	fprintf(stderr,"sh_syn: using regular spaced grid with -R%g/%g/%g/%g -I%g/%g spacing\n",
		w,e,s,n,dx,dy);
      if((w > e)||(s>n)||(s<-90)||(s>90)||(n<-90)||(n>90)){
	fprintf(stderr,"%s: range error\n",argv[0]);
	exit(-1);
      }

      if((ivec) && (s == -90)&&(n == 90)){
	s += dy/2;
	n -= dy/2;
	fprintf(stderr,"sh_syn: vector fields: adjusting to -R%g/%g/%g/%g\n",
		w,e,s,n);
      }
      /*  */
      dphi = DEG2RAD(dx);
      nphi = DEG2RAD(e-w)/dphi + 1;
      dtheta = DEG2RAD(dy);
      ntheta = DEG2RAD(n-s)/dtheta + 1;
      npoints = nphi * ntheta;

      /*  */
      hc_svecalloc(&phi,nphi,"sh_shsyn");
      hc_svecalloc(&theta,ntheta,"sh_shsyn");
      for(x=LON2PHI(w),i=0;i < nphi;i++,x += dphi)
	phi[i] = x;
      for(y = LAT2THETA(n),j=0;j < ntheta;y += dtheta,j++)
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
