#include "hc.h"

/* 

routines dealing with input, pass sol initialized as NULL

returns shps, the number of expansions. should be 3 for velocities and
tractions, and 1 for density anomalies

$Id: hc_input.c,v 1.8 2004/12/20 05:18:12 becker Exp $

*/
int hc_read_sh_solution(struct hcs *hc, 
			 struct sh_lms **sol, FILE *in, 
			 hc_boolean binary, hc_boolean verbose)
{
  int nset,ilayer,shps,lmax,type,ivec,nsol,i,os,n;
  HC_PREC zlabel,unity[3]={1.,1.,1.};
   /* 

   read all layes as spherical harmonics assuming real Dahlen & Tromp
     (physical) normalization
     
  */
  n = os = 0;
  while(sh_read_parameters_from_stream(&type,&lmax,&shps,&ilayer,
				       &nset,&zlabel,&ivec,in,
				       FALSE,binary,verbose)){
    hc->sh_type = type;
    if(ilayer != n){
      fprintf(stderr,"hc_read_sh_solution: error: ilayer %i n %i\n",
	      ilayer,n);
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
      hc->nradp2 = hc->nrad + 2;
      hc_vecrealloc(&hc->r,nset,"hc_read_sh_solution");
    }
    /* assign depth */
    hc->r[ilayer] = HC_ND_RADIUS(zlabel);
    /* 
       read coefficients
    */
    sh_read_coefficients_from_stream((*sol+os),shps,-1,in,
				     binary,unity,verbose);
    if(verbose){
      if(shps == 3)
	fprintf(stderr,"hc_read_sh_solution: z: %8.3f |r|: %12.5e |pol|: %12.5e |tor|: %12.5e\n",
		(double)HC_Z_DEPTH(hc->r[ilayer]),
		(double)sqrt(sh_total_power((*sol+os))),
		(double)sqrt(sh_total_power((*sol+os+1))),
		(double)sqrt(sh_total_power((*sol+os+2))));
      else
	fprintf(stderr,"hc_read_sh_solution: z: %8.3f |scalar|: %12.5e\n",
		(double)HC_Z_DEPTH(hc->r[ilayer]),
		(double)sqrt(sh_total_power((*sol+os))));
    }
    n++;
    os += shps;
  }
  if(n != nset)
    HC_ERROR("hc_read_sh_solution","read error");
  if(verbose){
    fprintf(stderr,"hc_read_sh_solution: read %i solution layer\n",nset);
    if(shps != 3)
      fprintf(stderr,"hc_read_sh_solution: WARNING: only scalar field read\n");
  }
  return shps;
}
