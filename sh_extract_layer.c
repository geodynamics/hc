#include "hc.h"
/* 


extract part of a spherical harmonics file that has several layers,
such as the hc geoid.ab output for the case of all layers (hc -ag)

*/

int main(int argc, char **argv)
{
  int ilayer,i,shps,nset,ivec,lmax,type;
  FILE *in;
  struct sh_lms *exp=NULL;
  HC_PREC unitya[3] = {1.0,1.0,1.0},*zdepth;
  hc_boolean binary = FALSE, verbose = TRUE, 
    short_format_in = FALSE,
    short_format_output = FALSE;
  /* 
     deal with parameters
  */
  hc_vecalloc(&zdepth,1,"");
  ilayer = 0;
  switch(argc){
  case 3:
    sscanf(argv[2],"%i",&ilayer);
    break;
  case 4:
    sscanf(argv[2],"%i",&ilayer);
    sscanf(argv[3],"%i",&i);
    short_format_output = (i)?(TRUE):(FALSE);
    break;
  case 5:
    sscanf(argv[2],"%i",&ilayer);
    sscanf(argv[3],"%i",&i);
    short_format_output = (i)?(TRUE):(FALSE);
    sscanf(argv[4],"%i",&i);
    short_format_in = (i)?(TRUE):(FALSE);
    break;

  default:
    fprintf(stderr,"%s: usage\n\n%s value.ab layer [short_format_out, %i] [short_format_in, %i]\n\n",
	    argv[0],argv[0],short_format_output,short_format_in);
    fprintf(stderr,"extracts one SH layer (e.g. for use in sh_syn) from a (long format, hc) spherical harmonic file value.ab\n");
    fprintf(stderr,"layer: 1...nset\n");
    fprintf(stderr,"\tif ilayer= 1..nset, will print one layer\n");
    fprintf(stderr,"\t          -1, will select nset\n");
    fprintf(stderr,"if short_format_out=1, will use short format output, else long\n");
    fprintf(stderr,"if short_format_in= 1, will use short format input,  else long\n");
    exit(-1);
    break;
  }
  /* 
     read in spherical harmonics
  */
  in = ggrd_open(argv[1],"r","sh_extract_layer");
  /* start loop */
  i = 0;
  sh_read_parameters_from_stream(&type,&lmax, &shps,&i,&nset,(zdepth+i),
				 &ivec,in,short_format_in,binary,verbose);
  if(verbose)
    fprintf(stderr,"sh_extract_layer: detected %i layers, vec: %i, type %i SH file, shps: %i, layer %i at depth %g\n",
	    nset,ivec,type,shps,i+1,(double)zdepth[i]);
  sh_allocate_and_init(&exp,nset*shps,lmax,type,ivec,verbose,FALSE);
  hc_vecrealloc(&zdepth,nset,"sh_extract_layer: zdepth mem error");
  /* which layer to select */
  if(ilayer == -1)
    ilayer = nset - 1;
  else
    ilayer--;
  /* check bounds */
  if(ilayer < 0){
    fprintf(stderr,"sh_extract_layer: ilayer should be given between 1 and %i for file %s\n",
	    nset,argv[1]);
    exit(-1);
  }
  for(;i<nset;i++){
    if(i != 0)
      sh_read_parameters_from_stream(&type,&lmax,&shps,&i,&nset,(zdepth+i),
				     &ivec,in,short_format_in,binary,verbose);
    sh_read_coefficients_from_stream((exp+i),shps,-1,in,binary,unitya,
				     verbose);
    if(i == ilayer){		/* output */
      if(verbose)
	fprintf(stderr,"%s: printing SH from %s at layer %i out of %i to stdout (depth: %g)\n",
		argv[0],argv[1],i+1,nset,(double)zdepth[i]);
      /* 
	 output will remove the layer number information 
      */
      sh_print_parameters_to_stream((exp+i),shps,ilayer,1,zdepth[i],
				    stdout,short_format_output,FALSE,verbose);
      sh_print_coefficients_to_stream((exp+i),shps,stdout,unitya,FALSE,verbose);
    }
  }
  fclose(in);

  sh_free_expansion(exp,nset);

  return 0;
}
