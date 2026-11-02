#include "hc.h"

/* 

this directory contains a few subroutines to test using HealPix and
Rick's spherical harmonics subroutines. It is intended to illustrate
how those packages can be used interchangably using the wrappers this
package provides

Thorsten Becker (twb@ig.utexas.edu)


usage: 

hctest type mode [lmax, 63]

type: 0=Healpix, 1=Rick
mode: 0=coord out, 1=expand 

$Id: sh_test.c,v 1.6 2004/12/20 05:09:42 becker Exp becker $

*/

int main(int argc, char **argv)
{
  FILE *out,*in;
  int nside,i,mode,type,nset,lmax,ivec;
  struct sh_lms_model model;
  float *ndata;
  hc_boolean verbose = TRUE, coord_out = FALSE;
  char filename[HC_CHAR_LENGTH];
  ndata = NULL;

  /* defaults */
  lmax = 63;
  ivec= 0;

  if(argc < 3){
    fprintf(stderr,"usage:\nhctest type[0=Healpix, 1=Rick] mode[0=coord out,1=expand] [lmax, %i] [ivec, %i\n\n",
	    lmax,ivec);
    exit(-1);
  }
  /* read in arguments */
  sscanf(argv[1],"%i",&type);
  sscanf(argv[2],"%i",&mode);
  if(argc>3)
    sscanf(argv[3],"%i",&lmax);
  if(argc>4)
    sscanf(argv[4],"%i",&ivec);
  if(mode == 0)			/* select mode */
    coord_out = TRUE;
  else
    coord_out = FALSE;
  /* number of sets */
  nset = 1;
  fprintf(stderr,"%s: type: %i mode: %i lmax: %i nset: %i\n",
	  argv[0],type,mode,lmax,nset);
  /* 
     parameters for Healpix 
  */
  nside = lmax/2+2;		/* lmax/2 should be OK */
  /* 

  initialize an expansion of lmax, allocating space for spectral
  storage

  */
  sh_init_model(&model,lmax,type,nset,ivec,nside,SH_HEALPIX_RING,
		verbose);
  if(coord_out){
    /* 
       print the spatial basis coordinates on which the 
       spherical harmonics routines operate
    */
    sprintf(filename,"hc.%i.dat",type);
    out = ggrd_open(filename,"w",argv[0]);
    sh_print_model_spatial_basis(&model,out,verbose);
    fclose(out);
    if(verbose)
      fprintf(stderr,"%s: written %i x %i = %i coordinates type %i to %s\n",
	      argv[0],model.exp[0].npoints,
	      model.nset,model.tnpoints,type,filename);
    /* that's it */
    exit(0);
  }
  /* 
     read in data at pixel locations , allocates odata array
  */
  sprintf(filename,"test.%i.data",type);
  in = ggrd_open(filename,"r",argv[0]);
  sh_read_model_spatial_data(&model,&model.data,in,verbose);
  fclose(in);
  
  /* 
     compute the original expansion 
  */
  sh_compute_model_spectral(&model,model.data,verbose);
  /* 
     output of expansion 
  */
  sh_print_model_coefficients(&model,stdout,FALSE,verbose);
  /* 
     compute the spatial basis of that expansion, 
     also allocates ndata array
  */
  sh_compute_model_spatial(&model,&ndata,verbose);
  /* 
     output of spatial basis of that expansion 
  */
  sprintf(filename,"testn.%i.data",type);
  out = ggrd_open(filename,"w",argv[0]);
  sh_print_model_spatial_data(&model,ndata,out,verbose);
  fclose(out);
  for(i=0;i<3;i++){
    /* 
       re-expand 
    */
    fprintf(stderr,"expansion: %i: misfit: %.6e\n",
	    i+1,hc_svec_rms_diff(model.data,ndata,
			      model.tnpoints));
    sh_compute_model_spectral(&model,ndata,verbose);
    hc_a_equals_b_svector(model.data,ndata,model.tnpoints);
    sh_compute_model_spatial(&model,&model.data,verbose);
  }

 
  /* free the spectal data */
  sh_free_model(&model);
  /* free the spatial data */
  free(ndata);

  return 0;
}
