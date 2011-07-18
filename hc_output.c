/* 


output routines of Hager & Connell code

$Id: hc_output.c,v 1.11 2006/01/22 01:11:34 becker Exp becker $


*/
#include "hc.h"
/* 


print the spherical harmonics version of a solution set


*/
void hc_print_spectral_solution(struct hcs *hc,struct sh_lms *sol,FILE *out,int sol_mode, 
				hc_boolean binary, hc_boolean verbose)
{
  int i,os;
  static int ntype = 3;			/* three sets of solutions, r/pol/tor */
  HC_PREC fac[3];
  if(!hc->spectral_solution_computed)
    HC_ERROR("hc_print_spectral_solution","spectral solution not computed");
  /* 
     number of solution sets of ntype solutions 
  */
  for(i=os=0;i < hc->nradp2;i++,os += ntype){
    /* 
       
    scale to cm/yr, or MPa  for stress solutions
    
    */
    hc_compute_solution_scaling_factors(hc,sol_mode,hc->r[i],hc->dvisc[i],fac);
    /* 
       write parameters, convert radius to depth in [km]  
    */
    sh_print_parameters_to_file((sol+os),ntype,i,hc->nradp2,
				HC_Z_DEPTH(hc->r[i]),
				out,FALSE,binary,verbose);
    /* 
       
       write the set of coefficients in D&T convention
       
    */
    sh_print_coefficients_to_file((sol+os),ntype,out,fac,
				  binary,verbose);
    if(verbose >= 2){
      switch(sol_mode){
      case HC_VEL:
	fprintf(stderr,"hc_print_spectral_solution: z: %8.3f vel: |r|: %11.3e |pol|: %11.3e |tor|: %11.3e (scale: %g cm/yr)\n",
		HC_Z_DEPTH(hc->r[i]),sqrt(sh_total_power((sol+os))),
		sqrt(sh_total_power((sol+os+1))),
		sqrt(sh_total_power((sol+os+2))),
		fac[0]/11.1194926644559);
	break;
      case HC_RTRACTIONS:
	fprintf(stderr,"hc_print_spectral_solution: z: %8.3f rtrac: |r|: %11.3e |pol|: %11.3e |tor|: %11.3e (scale: %g MPa)\n",
		HC_Z_DEPTH(hc->r[i]),sqrt(sh_total_power((sol+os))),
		sqrt(sh_total_power((sol+os+1))),
		sqrt(sh_total_power((sol+os+2))),
		fac[0]/(0.553073278428428/hc->r[i]));
	break;
      case HC_HTRACTIONS:
	fprintf(stderr,"hc_print_spectral_solution: z: %8.3f htrac: |r|: %11.3e |pol|: %11.3e |tor|: %11.3e (scale: %g MPa)\n",
		HC_Z_DEPTH(hc->r[i]),sqrt(sh_total_power((sol+os))),
		sqrt(sh_total_power((sol+os+1))),
		sqrt(sh_total_power((sol+os+2))),
		fac[0]/(0.553073278428428/hc->r[i]));
	break;
      default:
	fprintf(stderr,"hc_print_spectral_solution: sol mode %i undefined\n",sol_mode);
	exit(-1);
	break;
      }
    }
  }
  if(verbose)
    fprintf(stderr,"hc_print_spectral_solution: wrote solution at %i levels\n",
	    hc->nradp2);
}

/* 

print a single scalar field to file

*/

void hc_print_sh_scalar_field(struct sh_lms *sh, FILE *out, hc_boolean short_format,
			      hc_boolean binary, hc_boolean verbose)
{
  HC_CPREC fac[1] = {1.0};
  sh_print_parameters_to_file(sh,1,0,1,0.0,out,short_format,binary,verbose); /* parameters in long format */
  sh_print_coefficients_to_file(sh,1,out,fac,binary,verbose); /* coefficients */
}


/* 

print the spatial solution in 

lon lat vr vt vp 

format to nrad+2 files named filename.i.pre, where i is 0...nrad+1,
and pre is dat or bin, depending on ASCII or binary output.

will also write the corresponding depth layers to dfilename

*/
void hc_print_spatial_solution(struct hcs *hc, 
			       struct sh_lms *sol,
			       float *sol_x, char *name, 
			       char *dfilename, 
			       int sol_mode, hc_boolean binary, 
			       hc_boolean verbose)
{
  int i,j,k,os[3],los,np,np2,np3;
  FILE *file_dummy=NULL,*out,*dout;
  float flt_dummy=0,*xy=NULL,value[3];
  HC_PREC fac[3];
  char filename[300];
  if(!hc->spatial_solution_computed)
    HC_ERROR("hc_print_spatial_solution","spectral solution not computed");
  /* number of solution sets of ntype solutions */

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

  /* depth file */
  dout = ggrd_open(dfilename,"w","hc_print_spatial_solution");
  if(verbose >= 2)
    fprintf(stderr,"hc_print_spatial_solution: writing depth levels to %s\n",
	    dfilename);
  for(i=0;i < hc->nradp2;i++){
    /* 

    compute the scaling factors, those do depend on radius
    in the case of the stresses, so leave inside loop!

    */
    hc_compute_solution_scaling_factors(hc,sol_mode,hc->r[i],
					hc->dvisc[i],fac);

    /* write depth in [km] to dout file */
    fprintf(dout,"%g\n",HC_Z_DEPTH(hc->r[i]));
    for(k=0;k < 3;k++)		/* pointers */
      os[k] = i * np3 + k*np;
    /* 

    format:


    lon lat vr vt vp   OR

    lon lat srr srt srp 

    
    */

    if(binary){
      /* binary output */
      sprintf(filename,"%s.%i.bin",name,i+1);
      out = ggrd_open(filename,"w","hc_print_spatial_solution");
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
      out = ggrd_open(filename,"w","hc_print_spatial_solution");
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
	    hc->nradp2);
  free(xy);
}

/* 

print the depth layers solution

*/
void hc_print_depth_layers(struct hcs *hc, FILE *out, 
			   hc_boolean verbose)
{
  int i;
  /* number of solution sets of ntype solutions */
  for(i=0;i < hc->nradp2;i++)
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
   compute the r, theta, phi fac[3] scaling factors for the solution
   output
*/
void hc_compute_solution_scaling_factors(struct hcs *hc,
					 int sol_mode,
					 HC_PREC radius,
					 HC_PREC viscosity,
					 HC_PREC *fac)
{

 switch(sol_mode){
  case HC_VEL:
    fac[0]=fac[1]=fac[2] = hc->vel_scale; /* go to cm/yr  */
    break;
  case HC_RTRACTIONS:		/* radial tractions */
    fac[0]=fac[1]=fac[2] = hc->stress_scale/radius; /* go to MPa */
    break;
  case HC_HTRACTIONS:		/* horizontal tractions, are actually
				   given as strain-rates  */
    fac[0]=fac[1]=fac[2] = 2.0*viscosity*hc->stress_scale/radius; /* go to MPa */
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
				hc_boolean convert_to_dt, /* convert spherical harmonic coefficients 
							     to Dahlen & Tromp format */
				hc_boolean verbose)
{
  int l,m,i,j,a_or_b,ll,nl,os,alim;
  FILE *out;
  HC_PREC value[2];
  /* 
     output of poloidal solution vectors 
  */
  if(verbose)
    fprintf(stderr,"hc_print_poloidal_solution: printing poloidal solution vector %s to %s\n",
	    (convert_to_dt)?("(physical convention"):("(internal convention)"),filename);
  /* find max output degree */
  ll = HC_MIN(l_max,pol_sol[0].lmax);
  /* number of output layers */
  nl = hc->nrad + 2;
  
  out = ggrd_open(filename,"w","hc_print_poloidal_solution");
  for(l=1;l <= ll;l++){
    for(m=0;m <= l;m++){
      alim = (m==0)?(1):(2);
      for(a_or_b=0;a_or_b < alim;a_or_b++){
	for(i=os=0;i < nl;i++,os+=6){
	  fprintf(out,"%3i %3i %1i %3i %8.5f ",l,m,a_or_b,i+1,hc->r[i]);
	  for(j=0;j < 6;j++){
	    sh_get_coeff((pol_sol+os+j),l,m,a_or_b,convert_to_dt,value);
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
  out = ggrd_open(filename,"w","hc_toroidal_solution");
  for(l=1;l <= ll;l++){
    for(os=i=0;i < nl;i++,os+=lmaxp1)
      fprintf(out,"%3i %16.7e %16.7e %16.7e\n",
	      l,hc->r[i],tvec[os+l],tvec[os2+os+l]);
    fprintf(out,"\n");
  }
  fclose(out);
}
/* 

print a simple VTK file given already expanded input
(called from hc_print_spatial)


*/

void hc_print_vtk(FILE *out,float *xloc,float *xvec,
		  int npoints_orig,int nlay,
		  hc_boolean binary,int shps_d,
		  float *xscalar,int nlon, int nlat)
{
  int i,ilay,ndata,poff,j,ndata_d,nele_lay,nele_x,nele_y,
    npe,npe1,ncon[12],k,nleft,nlon_m1,npoints,
    nele_lay_reg,nele_lay_pole,tl,tr,nlay_m1,
    nele_brick = 8,nele_tri = 6;
  
  hc_boolean little_endian;
  float xtmp[3],r,spole[3],npole[3];
  /* determine machine type */
  little_endian = hc_is_little_endian();
  /*  */
  npoints = npoints_orig + 2;
  ndata =   npoints_orig*3;
  ndata_d = npoints_orig * shps_d;

  nlay_m1 = nlay - 1;
  /* print VTK */
  fprintf(out,"# vtk DataFile Version 4.0\n");
  fprintf(out,"generated by hc_print_vtk\n");
  if(binary)
    fprintf(out,"BINARY\n");
  else
    fprintf(out,"ASCII\n");
  fprintf(out,"DATASET UNSTRUCTURED_GRID\n");
  /* 
     nodes locations
  */
  fprintf(out,"POINTS %i float\n",npoints * nlay);
  for(ilay=0;ilay < nlay;ilay++){ /* bottom up */
    for(i=0;i < npoints_orig;i++){	  /* S to N, W to E */
      poff = ilay * ndata + i*3;
      if(binary)
	hc_print_be_float((xloc+poff),3,out,little_endian);
      else
	fprintf(out,"%.6e %.6e %.6e\n",
		xloc[poff],xloc[poff+1],xloc[poff+2]);
    }
    /* 
       south and north poles, add two, to go to npoints per layer 
    */
    poff = ilay * ndata;
    r = sqrt(xloc[poff]*xloc[poff] + 
	     xloc[poff+1]*xloc[poff+1] + 
	     xloc[poff+2]*xloc[poff+2]);
    /* south pole */
    xtmp[0] = xtmp[1] = 0.0;xtmp[2]=-r; 
    if(binary)
      hc_print_be_float((xtmp),3,out,little_endian);
    else 
      fprintf(out,"%.6e %.6e %.6e\n",xtmp[0],xtmp[1],xtmp[2]);
    /* north pole */
    xtmp[2] = r;		
    if(binary)
      hc_print_be_float((xtmp),3,out,little_endian);
    else 
      fprintf(out,"%.6e %.6e %.6e\n",xtmp[0],xtmp[1],xtmp[2]);
  }
  /*  */
  nele_x = nlon;nlon_m1=nlon-1; /* elements in longitude */
  nele_y = (nlat - 1);	  /* proper elements and pole connections */
  /* top row node names */
  tl = (nlat-1)*nlon;tr = tl + nlon;

  nele_lay_reg = nele_x * nele_y; /* regular elements per layer */
  nele_lay_pole =  2 * nele_x;
  nele_lay = nele_lay_reg + nele_lay_pole; /* total per layer */
  /* 
     element connectivity 
  */
  fprintf(out,"CELLS %i %i\n",nlay_m1 * nele_lay,
	  nlay_m1 * (nele_lay_reg * (1+nele_brick) 
		  + nele_x * 2 * (1+nele_tri)));
  for(ilay = 0; ilay < nlay_m1; ilay++){
    /* 
       loop bottom up, until one layer below top
    */
    npe = nele_brick;		/* real nodes per element */
    npe1 = npe+1;
    ncon[0] = npe;		/* counter */
    for(i=0;i < nele_y;i++){
      for(j=0;j < nele_x;j++){
	nleft = ilay * npoints + i * nlon + j;
	if(j == nlon_m1){	/* at edge, wrap around */
	  ncon[4] = nleft+nlon;ncon[3] = nleft-nlon_m1+nlon;
	  ncon[1] = nleft;     ncon[2] = nleft-nlon_m1;
	}else{			/* regular */
	  ncon[4] = nleft+nlon;ncon[3] = nleft+nlon+1;
	  ncon[1] = nleft;     ncon[2] = nleft+1;
	}
	/* top level */
	for(k=0;k < 4;k++)
	  ncon[5+k] = ncon[k+1] + npoints;
	if(binary){
	  hc_print_be_int(ncon,npe1,out,little_endian);
	}else{
	  for(k=0;k < npe1;k++)
	    fprintf(out,"%i ",ncon[k]);
	  fprintf(out,"\n");
	}
      }
    }
    /* 
       south and north polar rows 
    */
    npe = nele_tri;npe1 = npe+1;
    ncon[0] = npe;
    for(j=0;j < 2;j++){
      for(k=0;k < nele_x;k++){
	if(j == 0){
	  /* south pole */
	  nleft = ilay * npoints + k; 
	  ncon[1] = ilay * npoints + npoints_orig;
	  if(k == nlon_m1){	/* at edge, wrap around */
	    ncon[2] = nleft-nlon_m1;ncon[3]=nleft;
	  }else{
	    ncon[2] = nleft+1;ncon[3]=nleft;
	  }
	}else{		/* north pole */
	  nleft = ilay * npoints + (nele_y-1) * nlon + k; 
	  if(k == nlon_m1){	/* at edge, wrap around */
	    ncon[1] = nleft;ncon[2] = nleft-nlon_m1;
	  }else{
	    ncon[1] = nleft;ncon[2] = nleft+1;
	  }
	  ncon[3] = ilay * npoints + npoints_orig+1;
	}
	for(i=0;i < 3;i++)
	  ncon[4+i] = ncon[i+1] + npoints;
	if(binary){
	  hc_print_be_int(ncon,npe1,out,little_endian);
	}else{
	  for(i=0;i < npe1;i++)
	    fprintf(out,"%i ",ncon[i]);
	  fprintf(out,"\n");
	}
      }
    }
  }
  /* 
     print cell types 
  */
  fprintf(out,"CELL_TYPES %i\n",nlay_m1 * nele_lay);
  for(ilay = 0; ilay < nlay_m1; ilay++){
    if(binary){
      /* binary */
      ncon[0] = 12;		/* VTK quad */
      for(i=0;i < nele_lay_reg;i++)
	hc_print_be_int(ncon,1,out,little_endian);
      ncon[0] = 13;		/* VTK triangle */
      for(i=0;i < nele_lay_pole;i++)
	hc_print_be_int(ncon,1,out,little_endian);
    }else{
      /* ascicc */
      for(i=0;i < nele_lay_reg;i++){
	fprintf(out,"%i ",12);	/* vtk quad */
	if(i % 80 == 0)
	  fprintf(out,"\n");
      }
      for(i=0;i < nele_lay_pole;i++){
	fprintf(out,"%i ",13);	/* vtk triagnle */
	if(i % 80 == 0)
	  fprintf(out,"\n");
      }
      fprintf(out,"\n");
    }
  }
  fprintf(out,"POINT_DATA %i\n",npoints*nlay);
  if(shps_d){
    for(j=0;j < shps_d;j++){
      fprintf(out,"SCALARS scalar%i float 1\n",j+1);
      fprintf(out,"LOOKUP_TABLE default\n");
      for(ilay=0;ilay < nlay;ilay++){
	spole[0] = npole[0] = 0.0; /* pole avg */
	for(i=0,poff = j * nlay * ndata_d + ilay * ndata_d;
	    i < npoints_orig;i++,poff++){ /* all points */
	  if(binary)
	    hc_print_be_float((xscalar+poff),1,out,little_endian);
	  else{
	    fprintf(out,"%.6e ",xscalar[poff]);
	    if(i%20 == 0)fprintf(out,"\n");
	  }
	  if(i < nlon)
	    spole[0] += xscalar[poff];
	  if((i >= tl) && (i < tr))
	    npole[0] += xscalar[poff];
	}
	spole[0] /= (float)nlon;
	npole[0] /= (float)nlon;
	if(!binary){		/* ascii */
	  fprintf(out,"\n");
	  fprintf(out,"%.6e %.6e\n",spole[0],npole[0]);
	}else{			/* binary */
	  hc_print_be_float(spole,1,out,little_endian);
	  hc_print_be_float(npole,1,out,little_endian);
	}
      }
    }
  }
  fprintf(out,"VECTORS velocity float\n");
  for(ilay=0;ilay < nlay;ilay++){
    spole[0] = spole[1] = spole[2] = 
      npole[0] = npole[1] = npole[2] = 0.0; /* pole avg */
    for(i=0;i < npoints_orig;i++){
      poff = ilay * ndata + i*3;
      if(i < nlon){
	for(k=0;k<3;k++)
	  spole[k] += xvec[poff+k];
      }
      if((i >= tl) && (i < tr)){
	for(k=0;k<3;k++)
	  npole[k] += xvec[poff+k];
      }
      if(binary)		/* binary */
	hc_print_be_float((xvec+poff),3,out,little_endian);
      else			/* ascii */
	fprintf(out,"%.6e %.6e %.6e\n",xvec[poff],xvec[poff+1],
		xvec[poff+2]);
    }
    for(k=0;k<3;k++){
      spole[k] /= (float)nlon;
      npole[k] /= (float)nlon;
    }
    if(binary){
      /* binary */
      hc_print_be_float(spole,3,out,little_endian);  
      hc_print_be_float(npole,3,out,little_endian);  
    }else{
      /* ascii */
      fprintf(out,"%.6e %.6e %.6e\n",spole[0],spole[1],spole[2]);
      fprintf(out,"%.6e %.6e %.6e\n",npole[0],npole[1],npole[2]);
    }
  }
  
}
/* 
   print big endian binary to file, no matter what hardware

*/
void hc_print_be_float(float *x, int n, FILE *out, 
		       hc_boolean little_endian)
{
  int i;
  float *xcopy;
  static size_t len = sizeof(float);
  if(little_endian){
    /* need to flip the byte order */
    hc_svecalloc(&xcopy,n,"hc_print_be_float");
    memcpy(xcopy,x,len*n);
    for(i=0;i < n;i++)
      hc_flip_byte_order((void *)(xcopy+i),len);
    fwrite(xcopy,len,n,out);
    free(xcopy);
  }else{
    /* can write as is */
    fwrite(x,len,n,out);
  }
}
void hc_print_be_int(int *x, int n, FILE *out, 
		     hc_boolean little_endian)
{
  int i, *xcopy;
  static size_t len = sizeof(int);
  if(little_endian){
    /* need to flip the byte order */
    hc_ivecalloc(&xcopy,n,"hc_print_be_int");
    memcpy(xcopy,x,len*n);
    for(i=0;i < n;i++)
      hc_flip_byte_order((void *)(xcopy+i),len);
    fwrite(xcopy,len,n,out);
    free(xcopy);
  }else{
    /* can write as is */
    fwrite(x,len,n,out);
  }
}
/* 
   check if we're on a little endian machine
 */
hc_boolean hc_is_little_endian(void)
{
  static const unsigned long a = 1;
  return *(const unsigned char *)&a;
}

/* 

   flip endianness of x

*/
void hc_flip_byte_order(void *x, size_t len)
{
  void *copy;
  copy = (void *)malloc(len);
  if(!copy){
    fprintf(stderr,"flip_byte_order: memerror with len: %i\n",(int)len);
    exit(-1);
  }
  memcpy(copy,x,len);
  hc_flipit(x,copy,len);
  free(copy);
}
/* 

   actually flip the big endianness this should not be called with
   (i,i,size i)

*/
void hc_flipit(void *d, void *s, size_t len)
{
  unsigned char *dest = d;
  unsigned char *src  = s;
  src += len - 1;
  for (; len; len--)
    *dest++ = *src--;
}
/* print the density anomaly field interpolated to the nodal radii */
void hc_print_dens_anom(struct hcs *hc, 
			FILE *out, hc_boolean binary,
			hc_boolean verbose)
{
  int i,i1,i2;
  HC_PREC f1,f2;
  HC_PREC fac[3] = {1.,1.,1.};
  struct sh_lms *exp;
  sh_allocate_and_init(&exp,3,hc->dens_anom[0].lmax,hc->sh_type,0,FALSE,FALSE);
  for(i=0;i < hc->nradp2;i++){
    /* interpolate density depth to velocity node layer depth */
    hc_linear_interpolate(hc->rden,hc->inho,hc->r[i],&i1,&i2,&f1,&f2);
    
    sh_copy_lms((hc->dens_anom+i1),(exp+0));sh_scale_expansion((exp+0),f1); 
    sh_copy_lms((hc->dens_anom+i2),(exp+1));sh_scale_expansion((exp+1),f2);
    sh_c_is_a_plus_b_coeff((exp+2),(exp+0),(exp+1)); /* c = a+b */

    /* print to file */
    sh_print_parameters_to_file((exp+2),1,i,hc->nradp2,HC_Z_DEPTH(hc->r[i]),out,FALSE,binary,verbose);
    sh_print_coefficients_to_file((exp+2),1,out,fac,binary,verbose);
    if(verbose>2)fprintf(stderr,"hc_print_dens_anom: z: %8.3f (f1: %6.3f f2: %6.3f) %3i/%3i pow: %10.3e %10.3e %10.3e\n",
			 HC_Z_DEPTH(hc->r[i]),f1,f2,i+1,hc->nradp2,
			 sqrt(sh_total_power((hc->dens_anom+i1))),sqrt(sh_total_power((hc->dens_anom+i2))),
			 sqrt(sh_total_power((exp+2))));
  }
  
  sh_free_expansion(exp,3);
}
