#include "hc.h"
/* 


extract part of a spherical harmonics solution of a HC run


$Id: hc_extract_sh_layer.c,v 1.9 2006/01/22 01:11:34 becker Exp becker $


*/

int main(int argc, char **argv)
{
  int ilayer,nsol,i,mode,shps=1,nset=1,loop,i1,i2,shps_read;
  FILE *in;
  struct sh_lms *sol=NULL;
  struct hcs *model;
  HC_PREC fac[3] = {1.0,1.0,1.0};
  hc_boolean binary = TRUE, verbose = TRUE, short_format = FALSE;

  
  hc_struc_init(&model);
  /* 
     deal with parameters
  */
  ilayer = 0;
  mode = 1;
  switch(argc){
  case 3:
    sscanf(argv[2],"%i",&ilayer);
    break;
  case 4:
    sscanf(argv[2],"%i",&ilayer);
    sscanf(argv[3],"%i",&mode);
    break;
  case 5:
    sscanf(argv[2],"%i",&ilayer);
    sscanf(argv[3],"%i",&mode);
    sscanf(argv[4],"%i",&i);
    if(i)
      short_format = TRUE;
    else
      short_format = FALSE;
    break;
  default:
    fprintf(stderr,"%s: usage\n%s sol.file layer [mode,%i] [short_format, %i]\n\n",
	    argv[0],argv[0],mode,short_format);
    fprintf(stderr,"extracts spherical harmonic solution x (vel or str) from HC run\n");
    fprintf(stderr,"layer: 1...nset\n");
    fprintf(stderr,"\tif ilayer= 1..nset, will print one layer\n");
    fprintf(stderr,"\t          -1, will select nset\n");
    fprintf(stderr,"\t          -2, will print all layers\n");
    fprintf(stderr,"mode: 1...6\n");
    fprintf(stderr,"\tif mode = 1, will print x_r \n");
    fprintf(stderr,"\t          2, will print x_pol x_tor \n");
    fprintf(stderr,"\t          3, will print x_r x_pol x_tor\n");
    fprintf(stderr,"\t          4, will print the depth levels of all layers\n");
    fprintf(stderr,"\t          5, will print x_pol\n");
    fprintf(stderr,"\t          6, will print x_tor\n");
    exit(-1);
    break;
  }
  if(mode == 4)
    ilayer = -2;
  /* 
     read in solution
  */
  in = hc_fopen(argv[1],"r","hc_extract_sh_layer",argv[0]);
  shps_read = hc_read_sh_solution(model,&sol,in,binary,verbose);
  fclose(in);
  nsol = model->nradp2 * shps_read;
  /* 
     deal with selection
  */
  loop = 0;
  if(ilayer == -1)
    ilayer = model->nradp2;
  if(ilayer == -2){
    ilayer = model->nradp2;
    loop =1;
  }
  if((ilayer<1)||(ilayer > model->nradp2)){
    fprintf(stderr,"%s: ilayer (%i) out of range, use 1 ... %i\n",
	    argv[0],ilayer,model->nradp2);
    exit(-1);
  }
  if(loop){
    i1=0;i2=model->nradp2-1;
    if(short_format)
      printf("%i\n",model->nradp2);
  }else{
    i1=ilayer-1;i2 = i1;
  }
  /* 
     detect number of expansions 
  */
  if((mode == 1)||(mode == 5)||(mode == 6))
    shps = 1;
  else if(mode == 2)
    shps = 2;
  else if(mode == 3)
    shps = 3;

  if(shps > shps_read){
    fprintf(stderr,"%s: solution file only had %i expansions, mode %i requests %i\n",
	    argv[0],shps_read,mode,shps);
    exit(-1);
    
  }
  for(ilayer=i1;ilayer <= i2;ilayer++){
    /* 
       output 
    */
    if(mode != 4){
      /* SH header */
      if(short_format && loop)
	fprintf(stdout,"%g\n",(double)HC_Z_DEPTH(model->r[ilayer]));
      sh_print_parameters_to_stream((sol+ilayer*shps_read),shps,
				    ilayer,nset,(HC_PREC)(HC_Z_DEPTH(model->r[ilayer])),
				    stdout,short_format,FALSE,verbose);
    }
    switch(mode){
    case 1:
      /*  */
      if(verbose)
	fprintf(stderr,"%s: printing x_r SHE at layer %i (depth: %g)\n",
		argv[0],ilayer,(double)HC_Z_DEPTH(model->r[ilayer]));
      sh_print_coefficients_to_stream((sol+ilayer*shps_read),shps,stdout,fac,FALSE,verbose);
      break;
    case 2:
      /*  */
      if(verbose)
	fprintf(stderr,"%s: printing x_pol x_tor SHE at layer %i (depth: %g)\n",
		argv[0],ilayer,(double)HC_Z_DEPTH(model->r[ilayer]));
      sh_print_coefficients_to_stream((sol+ilayer*shps_read+1),shps,stdout,fac,FALSE,verbose);
      break;
    case 3:
      /* mode == 3 */
      if(verbose)
	fprintf(stderr,"%s: printing x_r x_pol x_tor SHE at layer %i (depth: %g)\n",
		argv[0],ilayer,(double)HC_Z_DEPTH(model->r[ilayer]));
      sh_print_coefficients_to_stream((sol+ilayer*shps_read),shps,stdout,fac,FALSE,verbose);
      break;
    case 4:
      fprintf(stdout,"%5i %11g\n",ilayer,(double)HC_Z_DEPTH(model->r[ilayer]));
      break;
    case 5:
      /*  */
      if(verbose)
	fprintf(stderr,"%s: printing x_pol SHE at layer %i (depth: %g)\n",
		argv[0],ilayer,(double)HC_Z_DEPTH(model->r[ilayer]));
      sh_print_coefficients_to_stream((sol+ilayer*shps_read+1),shps,stdout,fac,FALSE,verbose);
      break;
    case 6:
      /*  */
      if(verbose)
	fprintf(stderr,"%s: printing x_tor SHE at layer %i (depth: %g)\n",
		argv[0],ilayer,(double)HC_Z_DEPTH(model->r[ilayer]));
      sh_print_coefficients_to_stream((sol+ilayer*shps_read+2),shps,stdout,fac,FALSE,verbose);
      break;
 
    default:
      fprintf(stderr,"%s: error, mode %i undefined\n",argv[0],mode);
      exit(-1);
      break;
    }
  }
  /* clear and exit */
  sh_free_expansion(sol,nsol);

  return 0;
}
