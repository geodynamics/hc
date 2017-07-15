#include "hc.h"

/* 

read in data in lon lat z format, spaced as required by the particular
type of spherical harmonic expansion desired and performs a spherical
harmonics analysis. input spatial data is read from stdin, output
coefficients in Dahlen and Tromp normalization to stdout

Thorsten Becker (twb@ig.utexas.edu)


usage: 

  cat data.lonlatz | sh_ana l_max ivec type short_format

       l_max: max order of expansion. if negative, will print out the spatial 
              locations needed on input
       
       ivec:  0: expand scalar field
              1: expand vector field
	      file: if filename other than 0, 1, or vec_t.grd will read from Netcdf/GMT grid file
	      vec_t.grd: will read theta and phi components of vector field from vec_t.grd and vec_p.grd 

       type:  0: Rick's spherical harmonics
              1: Healpix (defunct)

       short_format: header format

              0: long output for hc
	      1: short output, only lmax
       

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
  int type = SH_RICK,lmax,shps, nset=1,ivec=0,ilayer=0,i,failed;
  struct sh_lms *exp;
  hc_boolean verbose = TRUE, use_3d = FALSE, short_format = FALSE,read_grd =FALSE,
    binary = FALSE, print_spatial_base = FALSE;
  HC_PREC *data, zlabel = 0,*flt_dummy;
  struct ggrd_gt ggrd[2];
  SH_RICK_PREC *dummy;
  HC_PREC fac[3] = {1.,1.,1.};
  char cdummy;

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
  if(argc > 2){
    if((strcmp(argv[2],"0")==0)||(strcmp(argv[2],"1")==0)) /* regular input */
      sscanf(argv[2],"%i",&ivec);
    else{
      /* GMT grd file input */
#ifndef USE_GMT3      
      failed = ggrd_grdtrack_init_general(FALSE,argv[2],&cdummy,"-fg",ggrd,TRUE,FALSE,FALSE);
#else
      failed = ggrd_grdtrack_init_general(FALSE,argv[2],&cdummy,"-Lg",ggrd,TRUE,FALSE,FALSE);
#endif
      if(failed){
	fprintf(stderr,"%s: error opening netcdf grd %s file\n",argv[0],argv[2]);
	exit(-1);
      }
      if(strcmp(argv[2],"vec_t.grd")==0){ /* vectors */
#ifndef USE_GMT3      
	failed = ggrd_grdtrack_init_general(FALSE,"vec_p.grd",&cdummy,"-fg",(ggrd+1),TRUE,FALSE,FALSE);
#else
	failed = ggrd_grdtrack_init_general(FALSE,"vec_p.grd",&cdummy,"-Lg",(ggrd+1),TRUE,FALSE,FALSE);
#endif
	if(failed){
	  fprintf(stderr,"%s: error opening second netcdf grd file %s\n",argv[0],"vec_p.grd");
	  exit(-1);
	}
	ivec = 1;
      }else{
	ivec = 0;
      }
      read_grd = TRUE;
    }
  }
  if(argc > 3)
    sscanf(argv[3],"%i",&type);
  if(argc > 4){
    sscanf(argv[4],"%i",&i);
    short_format = (hc_boolean)i;
  }
  if((argc > 5)||(argc<=1)){
    fprintf(stderr,"usage: %s l_max [ivec, %i] [type, %i] [short_format, %i]\n",
	    argv[0],ivec,type,short_format);
    fprintf(stderr,"        l_max: max order of expansion. if negative, will print out the spatial\n");
    fprintf(stderr,"                  locations needed on input\n");
    fprintf(stderr,"               for Rick SH format, lmax needs to be 2**n-1\n\n");
    fprintf(stderr,"        ivec:  0: expand scalar field (input: lon lat scalar)\n");
    fprintf(stderr,"               1: expand vector field (input: lon lat v_t v_p\n");
    fprintf(stderr,"               vec_t.grd: expand vector field in vec_t.grd (theta) and vec_p.grd (phi)\n");
    fprintf(stderr,"               file.grd: expand scalar in file.grd, where .grd indicates global Netcdf files\n\n");
    fprintf(stderr,"        type   %i: use Rick's routines internally (output format DT convection)\n",SH_RICK);
    fprintf(stderr,"               %i: use Healpix's routines internally (output format DT convection)\n\n",SH_HEALPIX);
    fprintf(stderr,"short_format:  0: use long header files format as in HC\n");
    fprintf(stderr,"               1: use short header file format (only lmax)\n\n");
    fprintf(stderr,"Note that integration accuracy will depend on the choice of l_max,\n");
    fprintf(stderr,"because only l_max-1 Gauss points are used.\n\n");
    exit(-1);
  }
  if(print_spatial_base)
    fprintf(stderr,"%s: printing spatial base for lmax: %i\n",
	    argv[0],lmax);
  else{
    if(read_grd)
     fprintf(stderr,"%s: expanding to lmax: %i, reading from %s %s\n",
	     argv[0],lmax,argv[2],(ivec)?("vec_p.grd"):(""));

    else
      fprintf(stderr,"%s: expanding to lmax: %i, expecting %s\n",
	      argv[0],lmax,(ivec)?("lon lat vt vp"):("lon lat x"));

  }
  /* 
     select numbers of expansions, scalar or pol/tor for vector field
  */
  shps = (ivec)?(2):(1);
  /* intialize expansion first */
  sh_allocate_and_init(&exp,shps*nset,lmax,type,ivec,verbose,FALSE);
  /* make room for data */
  hc_vecalloc(&data,shps * exp->npoints,"sh_ana");
  if(print_spatial_base){
    /* 
       print out spatial basis 
    */
    sh_compute_spatial_basis(exp,stdout,use_3d,zlabel,&flt_dummy,0,verbose);
  }else{
    /* 
       perform spherical harmonic analysis
    */
    if(read_grd){
      sh_read_spatial_data_from_grd(exp,ggrd,use_3d,shps,data,&zlabel);
    }else{
      /* read in data from stdin */
      sh_read_spatial_data_from_stream(exp,stdin,use_3d,shps,data,&zlabel);
    }
    /* 
       perform spherical harmonic expansion 
    */
    sh_compute_spectral(data,ivec,FALSE,&dummy,
			exp,verbose);
    /* print parameters of expansion */
    sh_print_parameters_to_stream(exp,shps,ilayer,nset,zlabel,
				  stdout,short_format,binary,
				  verbose);
    /* print coefficients */
    sh_print_coefficients_to_stream(exp,shps,stdout,fac,binary, 
				    verbose);
  }
  fprintf(stderr,"%s: printing to stdout, done\n",argv[0]);
  free(exp);free(data);
  return 0;
}
