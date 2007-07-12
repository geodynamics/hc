#include "hc.h"

/* 

routines dealing with input, pass sol initialized as NULL


$Id: hc_input.c,v 1.8 2004/12/20 05:18:12 becker Exp $

*/
void hc_read_sh_solution(struct hcs *hc, struct sh_lms **sol, FILE *in, 
			 hc_boolean binary, hc_boolean verbose)
{
  int nset,ilayer,shps,lmax,type,ivec,nsol,i,os,n;
  HC_PREC zlabel,unity[3]={1.,1.,1.};
   /* 

   read all layes as spherical harmonics assuming real Dahlen & Tromp
     (physical) normalization
     
  */
  n = os = 0;
  while(sh_read_parameters_from_file(&type,&lmax,&shps,&ilayer,&nset,
				     &zlabel,&ivec,in,FALSE,binary,
				     verbose)){
    hc->sh_type = type;
    if((shps != 3)||(ilayer != n)){
      fprintf(stderr,"hc_read_sh_solution: error: shps %i ilayer %i n %i\n",
	      shps,ilayer,n);
      exit(-1);
    }
    if(n == 0){			/* initialize */
      /* solution expansions */
      nsol = shps * nset;
      *sol = (struct sh_lms *)realloc(*sol,nsol * sizeof(struct sh_lms));
      if(!(*sol))
	HC_MEMERROR("hc_read_sh_solution: sol");
      hc->sh_type = type;
      for(i=0;i < nsol;i++)	/* init expansions on irregular grid */
	sh_init_expansion((*sol+i),lmax,hc->sh_type,1,verbose,FALSE);
      hc->nrad = nset - 2;
      hc_vecrealloc(&hc->r,nset,"hc_read_sh_solution");
    }
    /* assign depth */
    hc->r[ilayer] = HC_ND_RADIUS(zlabel);
    /* 
       read coefficients
    */
    sh_read_coefficients_from_file((*sol+os),shps,-1,in,
				   binary,unity,verbose);
    if(verbose)
      fprintf(stderr,"hc_read_sh_solution: z: %8.3f |r|: %12.5e |pol|: %12.5e |tor|: %12.5e\n",
	      HC_Z_DEPTH(hc->r[ilayer]),
	      sqrt(sh_total_power((*sol+os))),sqrt(sh_total_power((*sol+os+1))),
	      sqrt(sh_total_power((*sol+os+2))));

    n++;
    os += shps;
  }
  if(n != nset)
    HC_ERROR("hc_read_sh_solution","read error");
  if(verbose)
    fprintf(stderr,"hc_read_sh_solution: read %i solution layer\n",nset);
}
