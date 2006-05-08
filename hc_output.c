/* 


output routines of Hager & Connell code

$Id: hc_output.c,v 1.11 2006/01/22 01:11:34 becker Exp becker $


*/
#include "hc.h"
/* 


print the spherical harmonics version of a solution set


*/
void hc_print_spectral_solution(struct hcs *hc,struct sh_lms *sol,
				FILE *out,int sol_mode, 
				hc_boolean binary, 
				hc_boolean verbose)
{
  int nradp2,i,os;
  static int ntype = 3;			/* three sets of solutions, r/pol/tor */
  HC_PREC fac[3];
  if(!hc->spectral_solution_computed)
    HC_ERROR("hc_print_spectral_solution","spectral solution not computed");
  /* 
     number of solution sets of ntype solutions 
  */
  nradp2 = hc->nrad+2;
  hc_compute_solution_scaling_factors(hc,sol_mode,fac);
  for(i=os=0;i < nradp2;i++,os += ntype){
    /* 
       write parameters, convert radius to depth in [km]  
    */
    sh_print_parameters((sol+os),ntype,i,nradp2,
			HC_Z_DEPTH(hc->r[i]),
			out,FALSE,binary,verbose);
    /* 
       
       write the set of coefficients in D&T convention
       
    */
    sh_print_coefficients((sol+os),ntype,out,fac,binary,verbose);
    if(verbose >= 2)
      fprintf(stderr,"hc_print_spectral_solution: z: %8.3f |r|: %11.3e |pol|: %11.3e |tor|: %11.3e (scale: %g cm/yr)\n",
	      HC_Z_DEPTH(hc->r[i]),sqrt(sh_total_power((sol+os))),
	      sqrt(sh_total_power((sol+os+1))),
	      sqrt(sh_total_power((sol+os+2))),
	      fac[0]/8.99321605918731);
  }
  if(verbose)
    fprintf(stderr,"hc_print_spectral_solution: wrote solution at %i levels\n",
	    nradp2);
}
/* 

print the spatial solution in 

lon lat vr vt vp 

format to nrad+2 files named filename.i.pre, where i is 0...nrad+1,
and pre is dat or bin, depending on ASCII or binary output.

will also write the corresponding depth layers to dfilename

*/
void hc_print_spatial_solution(struct hcs *hc, struct sh_lms *sol,
			       float *sol_x, char *name, 
			       char *dfilename, 
			       int sol_mode, hc_boolean binary, 
			       hc_boolean verbose)
{
  int nradp2,i,j,k,os[3],los,np,np2,np3;
  FILE *file_dummy=NULL,*out,*dout;
  float flt_dummy=0,*xy=NULL,value[3];
  HC_PREC fac[3];
  char filename[300];
  if(!hc->spatial_solution_computed)
    HC_ERROR("hc_print_spatial_solution","spectral solution not computed");
  /* number of solution sets of ntype solutions */
  nradp2 = hc->nrad+2;
  /* number of lateral points */
  np = sol[0].npoints;
  np2 = np*2;
  np3 = np*3;
  if(!np)
    HC_ERROR("hc_print_spatial_solution","npoints is zero");
  /* 
     compute the lateral coordinates
  */
  sh_compute_spatial_basis(sol, file_dummy, FALSE,flt_dummy, &xy,
			   1,verbose);
  /* 
     compute the scaling factors 
  */
  hc_compute_solution_scaling_factors(hc,sol_mode,fac);
  /* depth file */
  dout = hc_open(dfilename,"w","hc_print_spatial_solution");
  if(verbose >= 2)
    fprintf(stderr,"hc_print_spatial_solution: writing depth levels to %s\n",
	    dfilename);
  for(i=0;i < nradp2;i++){
    /* write depth in [km] to dout file */
    fprintf(dout,"%g\n",HC_Z_DEPTH(hc->r[i]));
    for(k=0;k < 3;k++)		/* pointers */
      os[k] = i * np3 + k*np;
    /* 

    format:


    lon lat vr vt vp

    
    */

    if(binary){
      /* binary output */
      sprintf(filename,"%s.%i.bin",name,i+1);
      out = hc_open(filename,"w","hc_print_spatial_solution");
      for(j=los=0;j < np;j++,los+=2){ /* loop through all points in layer */
	fwrite((xy+los),sizeof(float),2,out);
	for(k=0;k<3;k++)
	  value[k] = sol_x[os[k]] * fac[k];
	fwrite(value,sizeof(float),3,out);
	os[0]++;os[1]++;os[2]++;
      }     
      fclose(out);
    }else{
      /* ASCII output */
      sprintf(filename,"%s.%i.dat",name,i+1);
      out = hc_open(filename,"w","hc_print_spatial_solution");
      for(j=los=0;j < np;j++,los+=2){ /* loop through all points in layer */
	for(k=0;k<3;k++)
	  value[k] = sol_x[os[k]] * fac[k];
	fprintf(out,"%11g %11g\t%12.5e %12.5e %12.5e\n",
		xy[los],xy[los+1],value[0],value[1],value[2]);
	os[0]++;os[1]++;os[2]++;
      }
      fclose(out);
    }
    if(verbose >= 2)
      fprintf(stderr,"hc_print_spatial_solution: layer %3i: RMS: r: %12.5e t: %12.5e p: %12.5e file: %s\n",
	      i+1,hc_svec_rms((sol_x+i*np3),np),hc_svec_rms((sol_x+i*np3+np),np),
	      hc_svec_rms((sol_x+i*np3+np2),np),
	      filename);
  }
  fclose(dout);
  if(verbose)
    fprintf(stderr,"hc_print_spatial_solution: wrote solution at %i levels\n",
	    nradp2);
  free(xy);
}

/* 

print the depth layers solution

*/
void hc_print_depth_layers(struct hcs *hc, FILE *out, 
			   hc_boolean verbose)
{
  int nradp2,i;
  /* number of solution sets of ntype solutions */
  nradp2 = hc->nrad + 2;
  for(i=0;i < nradp2;i++)
    fprintf(out,"%g\n",HC_Z_DEPTH(hc->r[i]));
}


/* 

print a [3][3] matrix

*/
void hc_print_3x3(HC_PREC a[3][3], FILE *out)
  {
  int i,j;
  for(i=0;i<3;i++){
    for(j=0;j<3;j++)
      fprintf(out,"%11.4e ",a[i][j]);
    fprintf(out,"\n");
  }
}
/* 

print a [6][4] solution matrix

*/
void hc_print_sm(HC_PREC a[6][4], FILE *out)
{
  int i,j;
  for(i=0;i < 6;i++){
    for(j=0;j<4;j++)
      fprintf(out,"%11.4e ",a[i][j]);
    fprintf(out,"\n");
  }
}

void hc_print_vector(HC_PREC *a, int n,FILE *out)
{
  int i;
  for(i=0;i<n;i++)
    fprintf(out,"%11.4e ",a[i]);
  fprintf(out,"\n");
}
void hc_print_vector_label(HC_PREC *a, int n,FILE *out,
			   char *label)
{
  int i;
  fprintf(out,"%s: ",label);
  for(i=0;i<n;i++)
    fprintf(out,"%11.4e ",a[i]);
  fprintf(out,"\n");
}
void hc_print_matrix_label(HC_PREC *a, int m,
			   int n,FILE *out,char *label)
{
  int i,j;
  for(j=0;j<m;j++){
    fprintf(out,"%s: ",label);
    for(i=0;i<n;i++)
      fprintf(out,"%11.4e ",a[j*n+i]);
    fprintf(out,"\n");
  }
}


void hc_print_vector_row(HC_PREC *a, int n,FILE *out)
{
  int i;
  for(i=0;i<n;i++)
    fprintf(out,"%11.4e\n",a[i]);
}

/* 
   compute the r, theta, phi fac[3] scaling factors for the solution output 
*/
void hc_compute_solution_scaling_factors(struct hcs *hc,int sol_mode,HC_PREC *fac)
{

 switch(sol_mode){
  case HC_VEL:
    fac[0]=fac[1]=fac[2] = 1.0/hc->vel_scale; /* go to cm/yr  */
    break;
  case HC_STRESS:
    fac[0]=fac[1]=fac[2] = 1.0; /* go to ??? */
    break;
  case HC_GEOID:
    fac[0]=fac[1]=fac[2] = 1.0; /* go to ??? */
    break;
  default:
    HC_ERROR("hc_print_spectral_solution","mode undefined");
    break;
  }

}
/* 
   
output of poloidal solution up to l_max

 */
void hc_print_poloidal_solution(struct sh_lms *pol_sol,
				struct hcs *hc,
				int l_max, char *filename,
				hc_boolean verbose)
{
  int l,m,i,j,a_or_b,ll,nl,os,alim;
  FILE *out;
  HC_PREC value[2];
  /* 
     output of poloidal solution vectors in internal convention
  */
  if(verbose)
    fprintf(stderr,"hc_print_poloidal_solution: printing poloidal colution vector to %s\n",
	    filename);
  /* find max output degree */
  ll = HC_MIN(l_max,pol_sol[0].lmax);
  /* number of output layers */
  nl = hc->nrad + 2;
  
  out = hc_open(filename,"w","hc_print_poloidal_solution");
  for(l=1;l <= ll;l++){
    for(m=0;m <= l;m++){
      alim = (m==0)?(1):(2);
      for(a_or_b=0;a_or_b < alim;a_or_b++){
	for(i=os=0;i < nl;i++,os+=6){
	  fprintf(out,"%3i %3i %1i %3i %8.5f ",l,m,a_or_b,i+1,
		  hc->r[i]);
	  for(j=0;j < 6;j++){
	    sh_get_coeff((pol_sol+os+j),l,m,a_or_b,FALSE,value);
	    fprintf(out,"%11.4e ",value[0]);
	  } /* end u_1 .. u_4 nu_1 nu_2 loop */
	  fprintf(out,"\n");
	} /* end layer loop */
      } /* and A/B coefficient loop */
    }	/* end m loop */
  } /* end l loop */
  fclose(out);
}

/* 
   print toroidal solution vector (kernel), not expansion
*/
void hc_print_toroidal_solution(double *tvec,int lmax,
				struct hcs *hc,int l_max_out, 
				char *filename,
				hc_boolean verbose)
{
  FILE *out;
  int ll,l,i,nl,lmaxp1,os,os2;
  ll = HC_MIN(l_max_out,lmax); /* output lmax */
  nl = hc->nrad + 2;		/* number of layers */
  lmaxp1 = lmax + 1;		/* max expansion */
  os2 = lmaxp1 * nl;
  /* 
     kernel output
  */
  if(verbose)
    fprintf(stderr,"hc_print_toiroidal_solution: writing toroidal solutions 1 and 2 as f(l,r) to %s\n",
	    filename);
  out = hc_open(filename,"w","hc_toroidal_solution");
  for(l=1;l <= ll;l++){
    for(os=i=0;i < nl;i++,os+=lmaxp1)
      fprintf(out,"%3i %16.7e %16.7e %16.7e\n",
	      l,hc->r[i],tvec[os+l],tvec[os2+os+l]);
    fprintf(out,"\n");
  }
  fclose(out);
}
