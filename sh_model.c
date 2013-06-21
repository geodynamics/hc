/* 


routines to handle models with several spherical harmonics expansions

the lower level routines that operate on single (or triple) 
sets of expansions are in sh_exp.c

$Id: sh_model.c,v 1.11 2006/01/22 01:11:34 becker Exp $

*/
#include "hc.h"
/* 
   
initialize a model structure that may hold several spherical harmonics
expansions

input:

model structure (holds all the expansions)

input:

lmax: max order of expansions
type: type of expansion (SH_RICK, SH_HEALPIX,...)
nset: number of spherical harmonic sets
ivec: 0: scalar expansion, one SH expansion per set
      1: u-r u-theta u-phi vector field, three expansions per set


for SH_HEALPIX:
nside: order of skymap
order: sorting order of skymap


*/
void sh_init_model(struct sh_lms_model *model,int lmax,int type,
		   int nset, int ivec, int nside, int order, 
		   hc_boolean verbose)
{
  int i;
  if(nset <= 0){
    fprintf(stderr,"sh_init_model: error: nset (%i) out of bounds\n",
	    nset);
    exit(-1);
  }
  model->nset = nset;
  if((ivec<0)||(ivec>1)){
    fprintf(stderr,"sh_init_model: error: ivec (%i) out of bounds\n",
	    ivec);
    exit(-1);
  }
  model->ivec = ivec;		/* is vector expansion? */
  if(model->ivec==1)/* SH expansions per set */
    model->shps = 3;		/* spectrally:
				   u_r poloidal toroidal
				   spatially:
				   u_r u_theta u_phi 
				*/
  else
    model->shps = 1;
  /* 
     initialize expansions 
  */
  model->nexp = model->shps * model->nset; 
  model->exp = (struct sh_lms *)
    calloc(model->nexp,sizeof(struct sh_lms));
  if(!model->exp)
    HC_MEMERROR("sh_init_model");
  model->tnpoints = 0;
  for(i=0;i < model->nexp;i++){	/* initialize expansions and add 
				   up total number of points 
				   in spatial domain
				*/
    /* use irregular grid */
    sh_init_expansion((model->exp+i),lmax,type,model->ivec,
		      verbose,FALSE);
    model->tnpoints += model->exp[i].npoints;
  }
  /* logic  flag for spatial data */
  model->spatial_init = FALSE;
  /* should we attempt to precompute and store legendre factors? */
  if((type == SH_RICK)||(type == SH_HEALPIX))
    model->save_plm = TRUE;
  else
    model->save_plm = FALSE;
  /* will be allocated later */
  model->plm = NULL;
  for(i=1;i < model->nexp;i++) 
    /* check if all expansisons are the
       same type so that we need only
       one Plm array 
    */
    if(model->exp[i].tn_plm != model->exp[0].tn_plm)
      HC_ERROR("init_sh_model","tnplm mismatch");
  /* 
     layer indicator attributes 
  */
  model->z = (HC_PREC *)calloc(model->nset,sizeof(HC_PREC));
  if(!model->z)
    HC_MEMERROR("sh_init_model: z");
  /* 
     data pointer, initialize as NULL
  */
  model->data = NULL;
  model->initialized = TRUE;
}
/* 
   free a model structure 
*/
void sh_free_model(struct sh_lms_model *model)
{
  free(model->z);
  sh_free_expansion(model->exp,model->nexp);
  if(model->save_plm)
    free(model->plm);
}
/* 

write the coefficients of a spherical harmonic expansion to out stream

output will be real spherical harmonics coefficients as in Dahlen and
Tromp p. 859

for vec = 0, then output is in format

   lmax i z_i nset 1 type_parameters
   
   A00 B00
   A10 B10
   A11 B11
   A20 ...

and so on, for each layer i (1..nset) of nset sets. 
if ntype=3, the coefficients will be 

   lmax i z_i nset 3 type_parameters

   A00_s B00_s   A00_p B00_p    A00_t B00_t
   A10_s B10_s   A10_p B10_p    A10_t B10_t
   ....

where s, p, and t are the scalar (radial), poloidal, and toroidal expansions, 
respectively

type_parameters is:

type parameters

where type can be SH_HEALPIX (0) or SH_RICK(1) and then the parameters 
depend on the type of expansion

z_i is the label for this set

*/
void sh_print_model_coefficients(struct sh_lms_model *model, 
				 FILE *out,hc_boolean binary,
				 hc_boolean verbose)
{
  int i;
  HC_CPREC *fac;
  hc_vecalloc(&fac,model->shps,"sh_print_model_coefficients");
  for(i=0;i<model->shps;i++)
    fac[i] = 1.0;
  for(i=0;i < model->nexp;i++)	/* check all expansions */
    if(!model->exp[i].spectral_init)
      HC_ERROR("sh_print_model_coefficients","coefficients not initialized expansion");
  for(i = 0; i < model->nset;i++){	/* loop through sets */
    /* output of parameters */
    sh_print_parameters_to_stream((model->exp+i*model->shps),
				model->shps,i,model->nset,
				model->z[i],out,FALSE,binary,verbose);
    /* output of coefficient(s) (1 (ivec=0) or 3 (ivec=1)) */
    sh_print_coefficients_to_stream((model->exp+i*model->shps), 
				    model->shps,out,fac,binary,verbose);
  } /* end set loop */
  free(fac);
}
/* 
   
write the coordinates of the spatial basis of a model including
possibly several spherical harmonics expansion. if nset==1, the format is 

lon lat

if nset != 1, the format is 

lon lat z

if z has been read in, else it's zero


*/
void sh_print_model_spatial_basis(struct sh_lms_model *model, 
				  FILE *out, 
				  hc_boolean verbose)
{
  int i;
  HC_PREC **flt_dummy=NULL;
  if(verbose)
    fprintf(stderr,"sh_print_model_spatial_basis: printing spatial basis for nset %i expansions\n",
	    model->nset);
  for(i=0;i < model->nset;i++)
    sh_compute_spatial_basis((model->exp+i*model->shps),out,
			     (model->nset==1)?(FALSE):(TRUE),
			     model->z[i],flt_dummy,0,verbose);
}
/*

   read in lon lat data tripels at all spatial basis point location of
   the spherical harmonic expansions

   if model.nset == 1:
   
   reads lon lat data from stream "in"

   if model.nset != 1:

   read lon lat z data from stream "in"


   data is either a scalar, or u_r u_theta u_phi for 
   model.ivec == 1
   
   input:  model parameter, including nset and ivec
   output: data, will be dimnesionalized exp->npoints * exp->nset

   data has to be initialized, eg. as NULL
   
*/
void sh_read_model_spatial_data(struct sh_lms_model *model, 
				HC_PREC **data,FILE *in,
				hc_boolean verbose)
{
  int i,j;
  /* 
     make room for the spatial data
  */
  fprintf(stderr,"sh_read_model_spatial_data: reading %i sets, expecting %s locations and %s data\n",
	  model->nset,(model->nset == 1)?("lon lat"):("lon lat z"),
	  (model->ivec==1)?("u_r u_theta u_phi"):("scalar"));
  /* 
     allocate space for all the data points
  */
  hc_vecrealloc(data, model->tnpoints,
		"sh_read_model_spatial_data");
  /* 
     check all sets, number of points have to be the same for each
     expansion
  */
  for(i=0;i<model->nset;i++)
    for(j=((i==0)?(1):(0));j<model->shps;j++){
      if(model->exp[i*model->shps+j].npoints != 
	 model->exp[0].npoints){
	fprintf(stderr,"sh_read_model_spatial_data: error: set %i expansion %i out of %i: npoints: %i npoints(0): %i\n",
		i+1,j+1,model->shps,
		model->exp[i*model->shps+j].npoints,
		model->exp[0].npoints);
	exit(-1);
      }
    }
  for(i=0;i < model->nset;i++)
    /* read in the data for this layer */
    sh_read_spatial_data_from_stream((model->exp+i*model->shps),
				     in,(model->nset!=1)?(TRUE):(FALSE),
				     model->shps,
				     (*data+model->shps*model->exp[i*model->shps].npoints),
				     (model->z+i));
  model->spatial_init = TRUE;
}

/* 
   compute the spherical harmonic expansion of a model with
   several sets. this calls the spherical harmonics routines
   
*/
void sh_compute_model_spectral(struct sh_lms_model *model,
			       HC_PREC *data,
			       hc_boolean verbose)
{
  int i;
  const int unity = 1,zero = 0;
  if(!model->spatial_init)
    HC_ERROR("sh_compute_model_spectral","spatial set not initialized");
  if(((model->ivec)&&(model->shps != 3))||
     ((!model->ivec)&&(model->shps != 1))){
    fprintf(stderr,"sh_compute_model_spectral: error: ivec: %i nshp: %i\n",
	    model->ivec,model->shps);
    exit(-1);
  }
  for(i=0;i < model->nset;i++){	/* loop through sets */
    if(model->ivec == 1){
      /* 
	 first call for theta, phi components 
	 to assign poloidal and toroidal part
      */
      sh_compute_spectral((data+model->exp[i*model->shps+0].npoints), 
			  unity,model->save_plm,&model->plm,
			  (model->exp+i*model->shps+1), 
			  verbose);
      /* then the u_r scalar components, which is the first 
	 expansion */
      sh_compute_spectral(data,zero,model->save_plm,&model->plm,
			  (model->exp+i*model->shps+0), 
			  verbose);
    }else{
      /* scalar expansion */
      sh_compute_spectral(data,zero,model->save_plm,&model->plm,
			  (model->exp+i),verbose);
    }
  }
}
/* 
   
expand the coefficients of a model into the spatial space
pass data initialized at least at NULL

*/
void sh_compute_model_spatial(struct sh_lms_model *model,
			      HC_PREC **data,hc_boolean verbose)
{
  int i;
  const int unity = 1, zero = 0;
  /* 
     resize data
  */
  if(verbose)
    fprintf(stderr,"sh_compute_model_spatial: computing %i spatial points, type %i\n",
	    model->tnpoints,model->exp[0].type);
  /* 
     make room for data 
  */
  hc_vecrealloc(data,model->tnpoints,"sh_compute_model_spatial");
  if(((model->ivec)&&(model->shps != 3))||
     ((!model->ivec)&&(model->shps != 1))){
    fprintf(stderr,"sh_compute_model_spatial: error: ivec: %i nshp: %i\n",
	    model->ivec,model->shps);
    exit(-1);
  }
  for(i=0;i < model->nset;i++){
    if(model->ivec == 1){
      /* 
	 first call for theta, phi components from poloidal and
	 toroidal expansions
      */
      sh_compute_spatial((model->exp+i*model->shps+1), 
			 unity,model->save_plm,&model->plm,
			 (*data+model->exp[i*model->shps+0].npoints),
			 verbose);
      /* then the u_r scalar expansion */
      sh_compute_spatial((model->exp+i*model->shps+0),
			 zero,model->save_plm,&model->plm,
			 *data,verbose);
    }else{
      /* scalar expansion */
      sh_compute_spatial((model->exp+i*model->shps),zero,
			 model->save_plm,&model->plm,
			 *data,verbose);
    }
  }
}
/* 
  
print a whole model of spatial data with different levels

*/
void sh_print_model_spatial_data(struct sh_lms_model *model,
				 HC_PREC *data, FILE *out,
				 hc_boolean verbose)
{
  int i,j,os;
  os = 0;
  for(i=0;i < model->nset;i++){
    /* print out data for the set */
    sh_print_spatial_data_to_stream((model->exp+i*model->shps),
				    (model->ivec)?(3):(1),(data+os),
				    (model->nset == 1)?(FALSE):(TRUE),
				    model->z[i],out);
    if(model->ivec)
      for(j=0;j<3;j++)
	os += model->exp[i*model->shps+j].npoints;
    else
      os += model->exp[i*model->shps+0].npoints;
  }
}
