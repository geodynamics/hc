#include "hc.h"
#include "gmt.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* 

read in spherical harmonics coefficients (stdin) and expand to spatial basis (stdout)

Thorsten Becker (twb@ig.utexas.edu)


$Id: sh_syn.c,v 1.6 2006/01/22 01:11:34 becker Exp $

*/

int main(int argc, char **argv)
{
  int type,lmax,shps,ilayer,nset,ivec,i,j,npoints,nphi,ntheta;
  /* 
     switches 
  */
  hc_boolean verbose = TRUE, short_format = FALSE ,short_format_ivec = FALSE ,binary = FALSE;
  int regular_basis = 0;
  HC_PREC w,e,s,n,dx,dy;
  /*  */
  FILE *in;
  HC_PREC *data,*theta,*phi;
  /* spacing for reg_ular output */
  HC_PREC dphi,x,y,dtheta;
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
    sscanf(argv[3],HC_FLT_FORMAT,&w);
    if(w == 999)
      regular_basis = -1;
    else
      regular_basis = 1;
  }
  if(argc > 4)
    sscanf(argv[4],HC_FLT_FORMAT,&e);
  if(argc > 5)
    sscanf(argv[5],HC_FLT_FORMAT,&s);
  if(argc > 6)
    sscanf(argv[6],HC_FLT_FORMAT,&n);
   if(argc > 7)
    sscanf(argv[7],HC_FLT_FORMAT,&dx);
  if(argc > 8)
    sscanf(argv[8],HC_FLT_FORMAT,&dy);
  else
    dy = dx;
  if((argc > 9)|| (argc < 0)){
    fprintf(stderr,"usage: %s [short_format, %i] [short_ivec, %i] [w, %g] [e, %g] [s, %g] [n, %g] [dx, %g] [dy, dx] (in that order)\n",
	    argv[0],short_format,short_format_ivec,(double)w,(double)e,(double)s,(double)n,(double)dx);
    fprintf(stderr,"short_format:\n\t0: expects regular format with long header\n");
    fprintf(stderr,"\t1: expects short format with only lmax in header\n\n");
    fprintf(stderr,"short_ivec:\n\t0: for short format, expect AB for scalar expansion\n");
    fprintf(stderr,"\t1: for short format, expect poloidal toroidal AP BP AT BT for vector expansion\n\n");
    fprintf(stderr,"w,e,...\n\tif none of those are set, will use Gauss latitudes and FFT divided longitudes dependening on lmax\n");
    fprintf(stderr,"\tif w is set to anything but 999, will switch to regular spaced output with -Rw/e/s/n -Idx/dy type output\n");
    fprintf(stderr,"\tif w is set to 999, will read lon lat in deg from \"tmp.lonlat\", and expand on those locations\n\n");
    fprintf(stderr,"The output format will depend on the type of SH input.\n");
    fprintf(stderr,"\tfor scalara: lon lat scalar if a single SH is read in, else lon lat zlabel scalar.\n");
    fprintf(stderr,"\tfor vectors: lon lat v_theta v_phi if a single SH is read in, else lon lat zlabel v_theta v_phi.\n\n\n");
    exit(-1);
  }
  if(verbose)
    fprintf(stderr,"%s: waiting to read spherical harmonic coefficients from stdin (use %s -h for help)\n",
	    argv[0],argv[0]);
  while(sh_read_parameters_from_stream(&type,&lmax,&shps,&ilayer,&nset,
				       &zlabel,&ivec,stdin,short_format,
				       binary,verbose)){
    if(short_format_ivec){
      ivec = 1;
      shps = 2;
    }
    if(verbose)
      fprintf(stderr,"%s: converting lmax %i ivec: %i at z: %g\n",
	      argv[0],lmax,ivec,(double)zlabel);

    /* input and init */
    sh_allocate_and_init(&exp,shps,lmax,type,ivec,verbose,((regular_basis != 0)?(1):(0)));
    sh_read_coefficients_from_stream(exp,shps,-1,stdin,binary,fac,verbose);
    if(regular_basis == 1){
      /* 
	 regular basis output on regular grid
      */
      if(verbose)
	fprintf(stderr,"sh_syn: using regular spaced grid with -R%g/%g/%g/%g -I%g/%g spacing\n",
		(double)w,(double)e,(double)s,(double)n,(double)dx,(double)dy);
      if((w > e)||(s>n)||(s<-90)||(s>90)||(n<-90)||(n>90)){
	fprintf(stderr,"%s: range error\n",argv[0]);
	exit(-1);
      }

      if((ivec) && (s == -90)&&(n == 90)){
	s += dy/2;
	n -= dy/2;
	fprintf(stderr,"sh_syn: vector fields: adjusting to -R%g/%g/%g/%g\n",
		(double)w,(double)e,(double)s,(double)n);
      }
      /*  */
      dphi = DEG2RAD(dx);
      nphi = ((e-w)/dx) + 1;
      dtheta = DEG2RAD(dy);
      ntheta = ((n-s)/dy) + 1;
      npoints = nphi * ntheta;

      /*  */
      hc_vecalloc(&phi,nphi,"sh_shsyn");
      hc_vecalloc(&theta,ntheta,"sh_shsyn");
      for(x=LON2PHI(w),i=0;i < nphi;i++,x += dphi)
	phi[i] = x;
      for(y = LAT2THETA(n),j=0;j < ntheta;y += dtheta,j++)
	theta[j] = y;
      hc_vecalloc(&data,npoints * shps,"sh_shsyn data");
      /* compute the expansion */
      sh_compute_spatial_reg(exp,ivec,FALSE,&dummy,
			     theta,ntheta,phi,nphi,data,
			     verbose,FALSE);
      /* output */
      sh_print_reg_spatial_data_to_stream(exp,shps,data,
					  (nset>1)?(TRUE):(FALSE),
					  zlabel, theta,ntheta,
					  phi,nphi,stdout);
    }else if(regular_basis == -1){
      /* output on locations input lon lat file */
      if(verbose)
	fprintf(stderr,"sh_syn: reading locations lon lat from stdin for expansion\n");
      npoints = 0;
      hc_vecalloc(&phi,1,"sh_syn");
      hc_vecalloc(&theta,1,"sh_syn");
      in = fopen("tmp.lonlat","r");
      if(!in){
	fprintf(stderr,"sh_syn: error, could not open tmp.lonlat for reading lon lat locations\n");
	exit(-1);
      }
      while(fscanf(in,HC_TWO_FLT_FORMAT,&dphi,&dtheta)==2){
	phi[npoints] = LON2PHI(dphi);
	theta[npoints] = LAT2THETA(dtheta);
	npoints++;
	hc_vecrealloc(&phi,npoints+1,"sh_syn");
	hc_vecrealloc(&theta,npoints+1,"sh_syn");
      }
      if(verbose)
	fprintf(stderr,"sh_syn: read %i locations lon lat from tmp.lonlat for expansion\n",npoints);
      fclose(in);
      hc_vecalloc(&data,npoints * shps,"sh_shsyn data");
      sh_compute_spatial_irreg(exp,ivec,theta,phi,npoints,data,verbose);
      sh_print_irreg_spatial_data_to_stream(exp,shps,data,(nset>1)?(TRUE):(FALSE),
					    zlabel,theta,phi,npoints,stdout);
     }else{			/* use the built in spatial basis (Gaussian) */
      /* expansion */
      hc_vecalloc(&data,exp[0].npoints * shps,"sh_syn");
      sh_compute_spatial(exp,ivec,FALSE,&dummy,data,verbose);
      /* output */
      sh_print_spatial_data_to_stream(exp,shps,data,(nset>1)?(TRUE):(FALSE),
				      zlabel,stdout);
    }
    free(exp);free(data);
  }

  return 0;
}
