#include "hc.h"

/* 

read in data in lon lat z format, spaced as required by the particular
type of spherical harmonic expansion desired and performs a spherical
harmonics analysis. input spatial data is read from stdin, output
coefficients in Dahlen and Tromp normalization to stdout

Thorsten Becker (twb@usc.edu)


usage: cat data.lonlatz | sh_ana l_max ivec

       l_max: max order of expansion. if negative, will print out the spatial 
              locations needed on input
       
       ivec:  0: expand scalar field
              1: expand vector field

where data.lonlatz has the data at the required locations as put out by sh_ana -lmax


for scalars, the input format is 

lon[deg] lat[deg] scalar

for vectors

lon[deg] lat[deg] v_theta v_phi

where theta and phi are the vector components in South and East direction, respectively. 


$Id: sh_ana.c,v 1.6 2006/01/22 01:11:34 becker Exp $

*/

int main(int argc, char **argv)
{
  int type = SH_RICK,lmax,shps, nset=1,ivec=0,ilayer=0;
  struct sh_lms *exp;
  hc_boolean verbose = TRUE, use_3d = FALSE, short_format = FALSE,
    binary = FALSE, print_spatial_base = FALSE;
  float *data, zlabel = 0,*flt_dummy;
  SH_RICK_PREC *dummy;
  HC_PREC fac[3] = {1.,1.,1.};

  /* 
     command line parameters
  */
  if(argc > 1){
    if((strcmp(argv[1],"-h")==0)||(strcmp(argv[1],"--help")==0)||(strcmp(argv[1],"-help")==0))
      argc = -1000;
    else{			/* max order of expansion */
      sscanf(argv[1],"%i",&lmax);
      if(lmax < 0){
	print_spatial_base = TRUE;
	lmax = -lmax;
      }
    }
  }
  if(argc > 2)
    sscanf(argv[2],"%i",&ivec);
  if(argc > 3)
    sscanf(argv[3],"%i",&type);
  if((argc > 4)||(argc<=1)){
    fprintf(stderr,"usage: %s l_max [ivec, %i] [type, %i]\n",
	    argv[0],ivec,type);
    fprintf(stderr,"       l_max: max order of expansion. if negative, will print out the spatial\n");
    fprintf(stderr,"                 locations needed on input\n");
    fprintf(stderr,"              for Rick, lmax needs to be 2**n-1\n\n");
    fprintf(stderr,"       ivec:  0: expand scalar field (input: lon lat scalar)\n");
    fprintf(stderr,"              1: expand vector field (input: lon lat v_t v_p\n\n");
    fprintf(stderr,"       type:  %i: use Healpix's routines\n",SH_HEALPIX);
    fprintf(stderr,"              %i: use Rick's routines\n",SH_RICK);
    exit(-1);
  }
  if(print_spatial_base)
    fprintf(stderr,"%s: printing spatial base for lmax: %i\n",
	    argv[0],lmax);
  else
    fprintf(stderr,"%s: expanding to lmax: %i, expecting %s\n",
	    argv[0],lmax,(ivec)?("lon lat vt vp"):("lon lat x"));
  /* 
     select numbers of expansions, scalar or pol/tor for vector field
  */
  shps = (ivec)?(2):(1);
  /* intialize expansion first */
  sh_allocate_and_init(&exp,shps*nset,lmax,type,ivec,verbose,FALSE);
  /* make room for data */
  hc_svecalloc(&data,shps * exp->npoints,"sh_ana");
  if(print_spatial_base){
    /* 
       print out spatial basis 
    */
    sh_compute_spatial_basis(exp,stdout,use_3d,zlabel,&flt_dummy,0,verbose);
  }else{
    /* 
       perform spherical harmonic analysis
    */
    /* read in data from stdin */
    sh_read_spatial_data_from_file(exp,stdin,use_3d,shps,data,&zlabel);
    /* 
       perform spherical harmonic expansion 
    */
    sh_compute_spectral(data,ivec,FALSE,&dummy,
			exp,verbose);
    /* print parameters of expansion */
    sh_print_parameters_to_file(exp,shps,ilayer,nset,zlabel,
				stdout,short_format,binary,
				verbose);
    /* print coefficients */
    sh_print_coefficients_to_file(exp,shps,stdout,fac,binary, 
				  verbose);
  }
  fprintf(stderr,"%s: printing to stdout, done\n",argv[0]);
  free(exp);free(data);
  return 0;
}
