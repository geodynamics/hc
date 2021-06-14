/* 


lower level routines to handle spherical harmonics expansions

the higher level routines that operate on single (or triple) 
sets of expansions are in sh_model.c

$Id: sh_exp.c,v 1.15 2006/03/20 05:32:48 becker Exp becker $

*/
#include "hc.h"
/* 

allocates and initializes spherical harmonics structure

*/
void 
sh_allocate_and_init (exp, n, lmax, type, ivec, verbose, regular)
struct sh_lms **exp;
int n;
int lmax;
int type;
int ivec;
hc_boolean verbose;
hc_boolean regular;
{
  int i;
  /* init as zeroes! (but this won't necessarily set the logic flags!) */
  *exp = (struct sh_lms *)calloc(n,sizeof(struct sh_lms));
  if(!(*exp))
    HC_MEMERROR("sh_allocate_and_init");
  for(i=0;i < n;i++){
    sh_init_expansion((*exp+i),lmax,type,ivec,verbose,regular);
  }
}

/* 

compute the parameters needed for a single expansion of order lmax and
allocate initial array for spectral coefficients. the type of array
depends on the type of expansion. 

if you want vector harmonics at any point during a computation, set
ivec to unity initially

coefficients are initialized as zero

if regular is set, will not use Gauss points

*/
void 
sh_init_expansion (exp, lmax, type, ivec, verbose, regular)
struct sh_lms *exp;
int lmax;
int type;
int ivec;
hc_boolean verbose;
hc_boolean regular;
{
  /* 
     initialize logic flags 
  */
  exp->spectral_init = FALSE;
  /* type of expansion, e.g. SH_HEALPIX or SH_RICK, this will be checked later */
  exp->type = type;
#ifdef HC_DEBUG
  /* 
     l and m bounds 
  */
  if(lmax < 1){
    fprintf(stderr,"sh_init_expansion: error: lmax out of bounds: %i\n",
	    lmax);
    exit(-1);
  }
#endif
  exp->lmax = lmax;
  /* same as above, plus one */
  exp->lmaxp1 = exp->lmax+1;
  /* 
     size of one set of coefficients (l,m) if stored by wasting space
     (assuming we are using complex numbers)
  */
  exp->lmbig = exp->lmaxp1 * exp->lmaxp1;
  /* 
     size of coefficients when stored compactly like 
     (l+1)*l/2 + m, times two for A and B
  */
  exp->lmsmall2 = (exp->lmaxp1)*(exp->lmaxp1+1); /* for A and B */
  
  /* 
   
  */
  exp->plm_computed = FALSE;

  /* 

  allocate the spectral (coefficients) storage and initialize possibly
  other arrays

  */
  switch(exp->type){
#ifdef HC_USE_HEALPIX

  case SH_HEALPIX:			/* SH_HEALPIX part */
    if(regular)
      HC_ERROR("regular init not implemented for healpix");
    /* 
       
       get single precision complex array which holds A and B 
       
    */
    exp->n_lm = exp->lmbig;
    hc_scmplx_vecalloc(&exp->alm_c,exp->n_lm,"init_expansion");
    sh_clear_alm(exp);
    /* 
       init the Healpix parameters and determine the number of points
       in the spatial domain
    */
    heal_init_parameters(&exp->heal,(lmax/2)+2,SH_HEALPIX_RING,
			 ivec,exp->lmax,&exp->npoints,
			 &exp->n_plm,&exp->tn_plm);

    break;
#endif
  case SH_RICK:			/* SH_RICK PART */
    exp->rick.was_called = FALSE; /* this used to work via a calloc
				     call, but because not int, make
				     sure to work */
    exp->rick.computed_legendre =
      exp->rick.initialized = 
      exp->rick.vector_sh_fac_init =
      exp->rick.sin_cos_saved = FALSE;
    /* 
       make room for the coefficients A and B in compact storage
    */
    exp->n_lm = exp->lmsmall2;
    /* 
       use single precision vector 
    */
    rick_vecalloc(&exp->alm,exp->n_lm,"sh_init_expansion");
    sh_clear_alm(exp);		/* set to zero */
    /* 
       
    init the parameters for Rick subroutines

    */
#ifdef NO_RICK_FORTRAN
    rick_init(exp->lmax,ivec,&exp->npoints,
	      &exp->n_plm,&exp->tn_plm,&exp->rick,regular);
#else
    /* f90 version */
    rick_f90_init(&exp->lmax,&ivec,&exp->npoints,
		  &exp->n_plm,&exp->tn_plm,regular);
#endif
    break;
#ifdef HC_USE_SPHEREPACK
  case SH_SPHEREPACK_GAUSS:
  case SH_SPHEREPACK_EVEN:
    HC_ERROR("init_expansion","Spherepack not implemented");
    /* 
       make room for coefficients
    */

    /* 
       initialize 
    */
    
    break;
#endif
  default:
    sh_exp_type_error("sh_init_expansion",exp);
    break;
  }
}
/* 
   free an expansion structure with n elements
*/
void 
sh_free_expansion (exp, n)
struct sh_lms *exp;
int n;
{
  int i;
  for(i=0;i<n;i++){
    switch(exp[i].type){
#ifdef HC_USE_HEALPIX

    case SH_HEALPIX:			/* SH_HEALPIX part */
      heal_free_structure(&exp[i].heal);
      free(exp[i].alm_c);
      break;
#endif
    case SH_RICK:
      free(exp[i].alm);
      break;
#ifdef HC_USE_SPHEREPACK
    case SH_SPHEREPACK_GAUSS:
    case SH_SPHEREPACK_EVEN:
      HC_ERROR("free_expansion","Spherepack not implemented");
      break;

#endif
    default:
      sh_exp_type_error("sh_free_expansion",(exp+i));
      break;
    }
  }
}
/* 
   zero out all coefficients
*/
void 
sh_clear_alm (exp)
struct sh_lms *exp;
{
  int i;
  switch(exp->type){
#ifdef HC_USE_HEALPIX
  case SH_HEALPIX:
    /* init with zeroes */
    for(i=0;i<exp->n_lm;i++)
      exp->alm_c[i].dr = exp->alm_c[i].di = 0.0;
    break;
#endif
  case SH_RICK:
    for(i=0;i<exp->n_lm;i++)
      exp->alm[i] = 0.0;
    break;
#ifdef HC_USE_SPHEREPACK
  case SH_SPHEREPACK_GAUSS:
  case SH_SPHEREPACK_EVEN:
    HC_ERROR("clear_alm","Spherepack not implemented");
    break;
#endif
  default:
    sh_exp_type_error("sh_clear_alm",exp);
    break;
  }
}
/* 

compute | \Psi | ^2 = \sum (2l+1) \sigma^2

kindof RMS^2

*/
HC_CPREC 
sh_total_power (exp)
struct sh_lms *exp;
{
  HC_PREC *power;
  double sum;
  int l;
  hc_vecalloc(&power,exp->lmaxp1,"sh_total_power");
  sh_compute_power_per_degree(exp,power);
  for(sum=0.0,l=0;l<=exp->lmax;l++)
    sum += (double)(2.0*(HC_CPREC)l+1.0) * (double)power[l];
  free(power);
  return (HC_CPREC)sum;
}
/* 

compute the power (sigma^2) per degree and unit area

power[lmaxp1]

*/
void 
sh_compute_power_per_degree (exp, power)
struct sh_lms *exp;
HC_PREC *power;
{
  int l,m;
  HC_CPREC value[2];
  hc_boolean need_b;
  for(l=0;l<=exp->lmax;l++){
    power[l] = 0.0;
    for(m=0;m<=l;m++){
      need_b = (hc_boolean) ((m == 0) ? (0) : (2));
      sh_get_coeff(exp,l,m,need_b,TRUE,value); /* convert to DT
						  normalization  */
      power[l] += value[0] * value[0];
      if(need_b) 
	power[l] += value[1] * value[1];
    } /* end m loop */
    power[l] /= 2.0*((HC_CPREC)l)+1.0;
  } /* end l loop */
}
/* compute total correlation up to llim */
HC_PREC 
sh_correlation (exp1, exp2, llim)
struct sh_lms *exp1;
struct sh_lms *exp2;
int llim;
{
  return sh_correlation_per_degree(exp1,exp2,1,llim);
}


HC_PREC 
sh_correlation_per_degree (exp1, exp2, lmin, lmax)
struct sh_lms *exp1;
struct sh_lms *exp2;
int lmin;
int lmax;
{
  int l,m;
  HC_CPREC sum[3],tmp,atmp,btmp,ctmp,value1[2],value2[2];
  double dtmp;
  hc_boolean need_b;

  sum[0]=sum[1]=sum[2]=0.0;

  if((lmax > exp1->lmax)||(lmax > exp2->lmax)||(lmax < 1)||(lmin < 1)){
    fprintf(stderr,"sh_compute_correlation_per_degree: error: L1 %i L2 %i lmin %i lmax %i\n",
	    exp1->lmax,exp2->lmax,lmin,lmax);
    exit(-1);
  }
  for(l=lmin;l <= lmax;l++){

    for(m=0;m<=l;m++){

      need_b = (hc_boolean) ((m == 0) ? (0) : (2));
      sh_get_coeff(exp1,l,m,need_b,TRUE,value1); /* convert to DT normalization  */
      sh_get_coeff(exp2,l,m,need_b,TRUE,value2); /* convert to DT normalization  */

      atmp = value1[0];
      ctmp = value2[0];
      sum[0] += atmp * ctmp;
      sum[1] += atmp * atmp;
      sum[2] += ctmp * ctmp;

      if(need_b){
	btmp = value1[1];
	dtmp = value2[1];
	sum[0] += btmp* dtmp;
	sum[1] += btmp * btmp;
	sum[2] += dtmp * dtmp;
      }
    } /* end m loop */
    
  } /* end l loop */
  tmp = sqrt(sum[1]*sum[2]);
  return sum[0]/tmp;
}
void 
sh_single_par_and_exp_to_file (exp, name, binary, verbose)
struct sh_lms *exp;
char *name;
hc_boolean binary;
hc_boolean verbose;
{
  FILE *out;
  out = fopen(name,"w");
  if(!out){
    fprintf(stderr,"sh_single_par_and_exp_to_file: ERROR: problem openeing %s\n",name);
    exit(-1);
  }
  sh_single_par_and_exp_to_stream(exp,out,binary,verbose);
  fclose(out);
  if(verbose)
    fprintf(stderr,"sh_single_par_and_exp_to_file: written to %s\n",name);
}



void 
sh_single_par_and_exp_to_stream (exp, out, binary, verbose)
struct sh_lms *exp;
FILE *out;
hc_boolean binary;
hc_boolean verbose;
{
  HC_PREC fac[1]={1.0};
  const hc_boolean short_format = FALSE;
  sh_print_parameters_to_stream(exp,1,0,1,0,out,short_format,binary,verbose);
  sh_print_coefficients_to_stream(exp,1,out,fac,binary,verbose);
}


/* 

print one line with all parameters needed to identify a spherical
harmonics expansion for a scalar (shps == 1), poloidal/toroidal (shps
== 2), or a vector field (shps == 3)

exp[shps]

ilayer: the layer index (0...nset-1) for nset layers

zlabel: label for layer ilayer


if short_format is selected, will only print

lmax

*/
void 
sh_print_parameters_to_stream (exp, shps, ilayer, nset, zlabel, out, short_format, binary, verbose)
struct sh_lms *exp;
int shps;
int ilayer;
int nset;
HC_CPREC zlabel;
FILE *out;
hc_boolean short_format;
hc_boolean binary;
hc_boolean verbose;
{
  HC_PREC fz;
  /* 
     print
     
     lmax i z[i] nset shps expansion_type
     
  */
  if(binary){
    fz = (HC_PREC)zlabel;
    fwrite(&exp[0].lmax,sizeof(int),1,out);
    if(!short_format){
      fwrite(&ilayer,sizeof(int),1,out);
      hc_print_float(&fz,1,out);
      fwrite(&nset,sizeof(int),1,out);
      fwrite(&shps,sizeof(int),1,out);
      fwrite(&exp[0].type,sizeof(int),1,out);
     }
  }else{
    if(!short_format)
      fprintf(out,"%6i %6i %.8e %6i %2i %2i ",
	      exp[0].lmax,ilayer,
	      (double)zlabel,nset,
	      shps,exp[0].type);
    else
      fprintf(out,"%6i ",
	      exp[0].lmax);
  }
  /* 
     additional parameters? 
  */
  switch(exp[0].type){
  case SH_RICK:
    break;
#ifdef HC_USE_HEALPIX
  case SH_HEALPIX:
    break;
#endif
#ifdef HC_USE_SPHEREPACK
  case SH_SPHEREPACK_GAUSS:
  case SH_SPHEREPACK_EVEN:
    break;
#endif
  default:
    sh_exp_type_error("sh_print_parameters",exp);
    break;
  }
  if(!binary)
    fprintf(out,"\n");		/* finish header line */
}
/* 
   
   read parameters needed to initialize an expansion
   can be used to read several layers of a model
   will return TRUE, if success
   
   default (long format) input is 

   type lmax shps ilayer nset zlabel ivec


   short format is just 

   lmax

   explanations:

   type: type of spherical harmonics expansion
   lmax: max degree
   shps: 1 or 3 for different sets of expansions, scalar or scalar + vector
   ilayer: 0..nset
   nset: number of sets
   zlabel: float label of this set
   ivec: scalar/vector flag. 0 for shps==1, 1 else
   
   
   
*/
hc_boolean 
sh_read_parameters_from_stream (type, lmax, shps, ilayer, nset, zlabel, ivec, in, short_format, binary, verbose)
int *type;
int *lmax;
int *shps;
int *ilayer;
int *nset;
HC_CPREC *zlabel;
int *ivec;
FILE *in;
hc_boolean short_format;
hc_boolean binary;
hc_boolean verbose;
{
  int input1[2],input2[3];
  HC_PREC fz;
  double dtmp;
  /* 
     read
     
     lmax i+1 z[i] nset shps expansion_type
     
  */
  if(binary){
    if(short_format){
      if(fread(input1,sizeof(int),1,in) != 1)
	return FALSE;
      *lmax = input1[0]; 
    }else{
      if(fread(input1,sizeof(int),2,in)+
	 hc_read_float(&fz,1,in) +
	 fread(input2,sizeof(int),3,in) != 6)
	return FALSE;
      *lmax = input1[0];    
      *ilayer=input1[1];
      *zlabel = (HC_CPREC) fz; 
      *nset = input2[0];
      *shps = input2[1];
      *type = input2[2];
    }
  }else{
    if(short_format){
      if(fscanf(in,"%i",lmax)!=1){
	return FALSE;
      }
    }else{
      if(fscanf(in,"%i %i %lf %i %i %i",
		lmax,ilayer,&dtmp,nset,shps,type)!=6){
	return FALSE;
      }
      *zlabel = (HC_PREC)dtmp;
    }
  }
  if(short_format){
    /* default settings for short format */
    *zlabel = 0.0; 
    *ilayer = 0; 
    *nset=1;
    *shps = 1;     
    *type = HC_DEFAULT_INTERNAL_FORMAT;
  }
  if(*shps == 1)
    *ivec = 0;
  else
    *ivec = 1;
  /* 
     additional parameters? 
  */
  switch(*type){
  case SH_RICK:
    break;
#ifdef HC_USE_HEALPIX
  case SH_HEALPIX:
    break;
#endif
#ifdef HC_USE_SPHEREPACK
  case SH_SPHEREPACK_EVEN:
  case SH_SPHEREPACK_GAUSS:
    break;
#endif
  default:
    fprintf(stderr,"sh_read_parameters: type %i undefined\n",
	    *type);
    exit(-1);
    break;
  }
  return TRUE;
}

/* 
   write the coefficients of a spherical harmonic expansion to out stream
   
   output will be real spherical harmonics coefficients as in Dahlen and
   Tromp p. 859
   
   for shps = 1
   
   A00 B00
   A10 B10
   A11 B11
   A20 ...

   for shps = 2

   A00_p B00_p    A00_t B00_t
   A10_p B10_p    A10_t B10_t
   ....
  

   for shps = 3

   A00_s B00_s   A00_p B00_p    A00_t B00_t
   A10_s B10_s   A10_p B10_p    A10_t B10_t
   ....
   
   where s, p, and t are the scalar (radial), poloidal, and toroidal expansions, 
   respectively


   pass shps as unity, if scalar, as 2 if pol/tor part of field, and
   as 3 if u_r pol tor components with three expansions for velocity
   field

   basically, all values of shps can be passed, as long as exp[shps]
   
   fac[3] scales the coefficients

*/
void 
sh_print_coefficients_to_stream (exp, shps, out, fac, binary, verbose)
struct sh_lms *exp;
int shps;
FILE *out;
HC_CPREC *fac;
hc_boolean binary;
hc_boolean verbose;
{
  int j,l,m;
  HC_PREC value[2];
  HC_PREC fvalue[2];
  /* 
     test  other expansions this set 
  */
  for(j=0;j < shps;j++){ /* check the lmax */
    if(exp[j].lmax != exp[0].lmax){
      fprintf(stderr,"sh_print_coefficients: error: lmax(%i):%i != lmax(0):%i\n",
	      j+1,exp[j].lmax,exp[0].lmax);
      exit(-1);
    }
    if(exp[j].type != exp[0].type ){
      fprintf(stderr,"sh_print_coefficients: error: type(%i):%i != type(0):%i\n",
	      j+1,exp[j].type,exp[0].type);
      exit(-1);
    }
  } /* end test */
  if(binary){
    for(l=0;l <= exp[0].lmax;l++)
      for(m=0;m <= l;m++)
	for(j=0;j < shps;j++){
	  /* 
	     output is in physical convention, convert from whatever
	     we are using internally
	  */
	  sh_get_coeff((exp+j),l,m,2,TRUE,value);
	  fvalue[0] = value[0]*fac[j];
	  fvalue[1] = value[1]*fac[j];
	  hc_print_float(fvalue, 2, out);
	}
  }else{
    for(l=0;l <= exp[0].lmax;l++){
      for(m=0;m <= l;m++){
	for(j=0;j < shps;j++){
	  /* output in physical convention, convert from internal
	     convention */
	  sh_get_coeff((exp+j),l,m,2,TRUE,value);
	  fprintf(out,"%15.7e %15.7e\t",
		  (double)(value[0]*fac[j]),
		  (double)(value[1]*fac[j]));
	}
	fprintf(out,"\n");
      } /* end m loop */
    }	/* end l loop */
    //fprintf(out,"\n");
  }
}
/* 

read in spherical harmonic coefficients in real, physics convention of
Dahlen and Tromp p. 859, and convert to whatever internal format may
be used

shps: number of spherical harmonics sets

scale coefficients with factor fac[shps]

lmax: -1: use lmax from expansion
       else: read up to lmax


*/
void 
sh_read_coefficients_from_stream (exp, shps, lmax, in, binary, fac, verbose)
struct sh_lms *exp;
int shps;
int lmax;
FILE *in;
hc_boolean binary;
HC_CPREC *fac;
hc_boolean verbose;
{
  int j,k,l,m,lmax_loc;
  HC_CPREC value[2]={0,0};
  HC_PREC fvalue[2]={0,0};
  if(lmax < 0)
    lmax_loc = exp[0].lmax;
  else
    lmax_loc = lmax;
  /* 
     test other expansions of this set 
  */
  for(j=1;j < shps;j++){ /* check the lmax */
    if(exp[j].lmax != exp[0].lmax){
      fprintf(stderr,"sh_read_coefficients_from_stream: error: lmax(%i):%i != lmax(0):%i\n",
	      j+1,exp[j].lmax,exp[0].lmax);
      exit(-1);
    }
    if(exp[j].type != exp[0].type ){
      fprintf(stderr,"sh_read_coefficients_from_stream: error: type(%i):%i != type(0):%i\n",
	      j+1,exp[j].type,exp[0].type);
      exit(-1);
    }
  } /* end test */
  if(binary){
    for(l=0;l <= lmax_loc;l++)
      for(m=0;m <= l;m++)
	for(j=0;j < shps;j++){
	  if(hc_read_float(fvalue,2,in)!=2){
	    fprintf(stderr,"sh_read_coefficients_from_stream: read error: set %i l %i m %i\n",
		    j+1,l,m);
	    exit(-1);
	  }
	  for(k=0;k<2;k++)
	    value[k] = (HC_CPREC)fvalue[k];
	  /* read in real, Dahlen & Tromp normalized coefficients and
	     convert to whatever format we are using internally */
	  sh_write_coeff((exp+j),l,m,(m==0)?(0):(2),TRUE,value);
	}
  }else{
    for(l=0;l <= lmax_loc;l++)
      for(m=0;m <= l;m++)
	for(j=0;j < shps;j++){
	  if(fscanf(in,HC_TWO_FLT_FORMAT,value,(value+1))!=2){
	    fprintf(stderr,"sh_read_coefficients_from_stream: read error: set %i l %i m %i, last val: %g %g\n",
		    j+1,l,m,(double)value[0],(double)value[1]);
	    exit(-1);
	  }
	  /* read in real, Dahlen & Tromp normalized coefficients and
	     convert to whatever format we are using internally */
	  sh_write_coeff((exp+j),l,m,(m==0)?(0):(2),TRUE,value);
	}
  }
  /* fill up rest with zeroes, if we are limiting to lmax_loc */
  if(lmax_loc < exp[0].lmax){
    value[0] = value[1] = 0.0;
    for(l=lmax_loc+1;l <= exp[0].lmax;l++)
      for(m=0;m <= l;m++)
	for(j=0;j < shps;j++)
	  sh_write_coeff((exp+j),l,m,(m==0)?(0):(2),TRUE,value);
  }



  for(j=0;j < shps;j++){
    //sh_print_nonzero_coeff((exp+j),stderr);
    sh_scale_expansion((exp+j),fac[j]);
    //sh_print_nonzero_coeff((exp+j),stderr);
    exp[j].spectral_init = TRUE;
  }
}
/* print a raw set of coefficients to out if nonzero, for
   debugging  */
void 
sh_print_nonzero_coeff (exp, out)
struct sh_lms *exp;
FILE *out;
{
  int l,m;
  HC_CPREC value[2];
  for(l=0;l <= exp->lmax;l++){
    for(m=0;m <= l;m++){
      sh_get_coeff(exp,l,m,2,FALSE,value);
      if(fabs(value[0])+fabs(value[1]) > 1e-8)
	fprintf(out,"%5i %5i %15.7e %15.7e\n",l,m,
		(double)value[0],(double)value[1]);
    }
  }
}

/* 

given an initializez expansion structure, read the corresponding spatial data in 

lon lat data

format from FILE *in

the output is in the data array, which has to be passed 

data[shps * exp->npoints]

if use_3d is set, will read in 

lon lat z data

instead

*/
void 
sh_read_spatial_data_from_stream (exp, in, use_3d, shps, data, z)
struct sh_lms *exp;
FILE *in;
my_boolean use_3d;
int shps;
HC_PREC *data;
HC_PREC *z;
{
  sh_read_spatial_data(exp,in,use_3d,shps,data,z);
}
/* 

generic function

*/
void 
sh_read_spatial_data (exp, in, use_3d, shps, data, z)
struct sh_lms *exp;
FILE *in;
my_boolean use_3d;
int shps;
HC_PREC *data;
HC_PREC *z;

{
  HC_PREC lon,lat,xp[3];
  int j,k;
  /* 
     read in data for each layer 
  */
  for(j=0;j < exp->npoints;j++){
    /* 
       get expected coordinates to check if the input is OK
    */
    switch(exp->type){
#ifdef HC_USE_HEALPIX
    case SH_HEALPIX:
      switch(exp->heal.ordering){
      case SH_HEALPIX_RING:  
	pix2ang_ring((long)exp->heal.nside,(long)j,
		     (xp+HC_THETA),(xp+HC_PHI));
	break;
      case SH_HEALPIX_NEST:  
	pix2ang_nest((long)exp->heal.nside,(long)j,
		     (xp+HC_THETA),(xp+HC_PHI));
	break;
      default:
	fprintf(stderr,"sh_read_spatial_data: error: ordering %i undefined\n",
		exp->heal.ordering);
	exit(-1);
	break;
      }
      break;			/* end Healpix branch */
#endif
    case SH_RICK:
      /* for Rick's type routine */
#ifdef NO_RICK_FORTRAN
      rick_pix2ang(j,exp->lmax,(xp+HC_THETA),(xp+HC_PHI),
		   &exp->rick);
#else
      rick_f90_pix2ang(&j,&exp->lmax,(SH_RICK_PREC)(xp+HC_THETA),
		       (SH_RICK_PREC)(xp+HC_PHI));
#endif
      break;
#ifdef HC_USE_SPHEREPACK
  case SH_SPHEREPACK_GAUSS:
  case SH_SPHEREPACK_EVEN:
    break;
#endif
    default:
      sh_exp_type_error("sh_read_model_spatial_data",exp);
      break;
    }	/* end type branch */
    /* 
       read coordinates 
    */
    if(!use_3d){
      /* 
	 read in lon lat  
      */
      if(fscanf(in,HC_TWO_FLT_FORMAT,&lon,&lat) != 2){
	fprintf(stderr,"sh_read_spatial_data: error: lon lat format: pixel %i: read error\n",
		(int)j);
	exit(-1);
      }
    }else{
      /* 
	 read in lon lat z[i] 
      */
      if(fscanf(in,HC_THREE_FLT_FORMAT,&lon,&lat,z) != 3){
	fprintf(stderr,"sh_read_spatial_data: error: lon lat z format: pixel %i: read error\n",
		(int)j);
	exit(-1);
      }
    }
    /* 
       read data 
    */
    for(k=0;k < shps;k++){
      if(fscanf(in,HC_FLT_FORMAT,(data+k*exp[0].npoints+j))!=1){
	fprintf(stderr,"sh_read_spatial_data: error: scalar format: pixel %i: read error\n",
		(int)j);
	exit(-1);
      }
    }
    /* 
       adjust longitude range 
    */
    if(lon < 0)
      lon += 360.0;
    /* 
       check if location is OK (we don't know about z defaults) 
    */
    if(((fabs(PHI2LON(xp[HC_PHI])-lon) > 1e-3)&&(fabs(fabs(PHI2LON(xp[HC_PHI])-lon)-360) > 1e-3))||
       (fabs(THETA2LAT(xp[HC_THETA])-lat) > 1e-3)){
      fprintf(stderr,"sh_read_model_spatial_data: error: pixel %i coordinate mismatch:\n",
	      (int)j);
      fprintf(stderr,"sh_read_model_spatial_data: orig: %g, %g file: %g, %g\n",
	      (double)PHI2LON(xp[HC_PHI]),(double)THETA2LAT(xp[HC_THETA]),
	      (double)lon,(double)lat);
      exit(-1);
    }
  }	/* end points in layer loop */
}

/* 
   
   print the spatial basis coordinates of a single spherical
   harmonics expansion 

   for show_z_label is FALSE, the format is

   lon lat 

   for all the npoints of an expansion

   for show_z_label is TRUE, the format is

   lon lat z

   as passed 

   if out_mode == 0, will write to file *out, 
   
   else will use x[] and store the values

*/
void 
sh_compute_spatial_basis (exp, out, use_3d, z, x, out_mode, verbose)
struct sh_lms *exp;
FILE *out;
hc_boolean use_3d;
HC_PREC z;
HC_PREC **x;
int out_mode;
hc_boolean verbose;
{
  int j,os,inc;
  HC_PREC xp[3];
  if(out_mode)			/* make room for storing x,y,z */
    hc_vecrealloc(x,exp->npoints*(2+((use_3d)?(1):(0))),"sh_compute_spatial_basis");
  inc = (use_3d)?(3):(2);
  for(j=os=0;j < exp->npoints;j++,os+=inc){
    /* 
       get coordinates 
    */
    switch(exp->type){
#ifdef HC_USE_HEALPIX
    case SH_HEALPIX:
      switch(exp->heal.ordering){
      case SH_HEALPIX_RING:  
	pix2ang_ring((long)exp->heal.nside,(long)j,
		     (xp+HC_THETA),(xp+HC_PHI));
	break;
      case SH_HEALPIX_NEST:  
	pix2ang_nest((long)exp->heal.nside,(long)j,
		     (xp+HC_THETA),(xp+HC_PHI));
	break;
      default:

	fprintf(stderr,"sh_compute_spatial_basis: error: ordering %i undefined\n",
		exp->heal.ordering);
	exit(-1);
	break;
      }
      break;			/* end Healpix part */
#endif
    case SH_RICK:
      /* compute the coordinates for Rick's type of expansion */
#ifdef NO_RICK_FORTRAN
      rick_pix2ang(j,exp->lmax,(xp+HC_THETA),(xp+HC_PHI),
		   &exp->rick);
#else
      rick_f90_pix2ang(&j,&exp->lmax,(xp+HC_THETA),(xp+HC_PHI));
#endif
      break;
#ifdef HC_USE_SPHEREPACK
  case SH_SPHEREPACK_GAUSS:
  case SH_SPHEREPACK_EVEN:
    break;
#endif
    default:
      sh_exp_type_error("sh_compute_spatial_basis",exp);
      break;
    }	/* end type of expansion branch */
    if(out_mode){		/* save */
     if(!use_3d){
       /* save in lon lat format */
       *(*x+os)   = (HC_PREC) PHI2LON(xp[HC_PHI]);
       *(*x+os+1) = (HC_PREC) THETA2LAT(xp[HC_THETA]);
     }else{
       /* save in lon lat z format */
       *(*x+os)   = (HC_PREC) PHI2LON(xp[HC_PHI]);
       *(*x+os+1) = (HC_PREC) THETA2LAT(xp[HC_THETA]);
       *(*x+os+2) = (HC_PREC) z;
     }
    }else{			/* write to file */
      if(!use_3d){
	/* write in lon lat format */
	fprintf(out,"%14.7f %14.7f\n",(double)PHI2LON(xp[HC_PHI]),
		(double)THETA2LAT(xp[HC_THETA]));
      }else{
	/* write in lon lat z format */
	fprintf(out,"%14.7f %14.7f %g\n",(double)PHI2LON(xp[HC_PHI]),
		(double)THETA2LAT(xp[HC_THETA]),(double)z);
      }
    }
  }
}
/* 


compute the coefficients of a spherical harmonic expansion 

depending on the implementation, this will save the Plm Legendre
factors

input:  data[exp->npoints * (1 + ivec)]
        if ivec == 0, will expand a scalar
	if ivec == 1, will expand the theta, phi part of a vector

	ivec: switch for scalar/vector expansion

	IF IVEC IS EVER TO BE SET TO UNITY AND THE PLM FUNCTIONS
	ARE SAVED, CALL WITH IVEC == 1 THE FIRST TIME TO INITIALIZE 
	VECTOR PLM FACTORS

	save_plm: if TRUE, will save the Legendre factors

input/output:

        plm: Legendre factors, has to be a pointer that is 
	     initialized, e.g. as NULL

	exp: spherical harmonics expansion(s), 
	     has parameters like
             lmax, and will hold the ALM arrays


	     THIS IS ONE EXPANSION FOR IVEC == 0 AND TWO FOR
	     IVEC = 1

*/
void 
sh_compute_spectral (data, ivec, save_plm, plm, exp, verbose)
HC_PREC *data;
int ivec;
hc_boolean save_plm;
SH_RICK_PREC **plm;
struct sh_lms *exp;
hc_boolean verbose;
{
  if(save_plm){
    /* 
       this routine will only compute the plm once and performs 
       some sanity checks
    */
    sh_compute_plm((exp+0),ivec,plm,verbose); 
  }
  switch(exp[0].type){
#ifdef HC_USE_HEALPIX
  case SH_HEALPIX:
    HC_ERROR("sh_compute_spectral","healpix: not properly implemented yet");
    if(ivec){
      /* vector part */
      HC_ERROR("sh_compute_spectral","healpix: ivec==1 not implemented yet");
    }else{
      /* 
	 scalar part 
      */
      if(save_plm){
	/* using precomputed Plm */
	heal_map2alm_sc_pre(&exp[0].heal.nside,&exp[0].lmax,
			    &exp[0].lmax,data,exp[0].alm_c, 
			    &exp[0].heal.cos_theta_cut,
			    exp[0].heal.weights,*plm);
      }else{
	/* compute Plm now */
	heal_map2alm_sc(&exp[0].heal.nside,&exp[0].lmax,
			&exp[0].lmax,data,exp[0].alm_c,
			&exp[0].heal.cos_theta_cut,
			exp[0].heal.weights);
      }
    }
    break;
#endif
  case SH_RICK:
    /* 
       compute either scalar expansion of data[n] array
       poloidal and toroidal fields of data[2*n].

       if ivec == 0, exp[1].alm will not be referenced
    */
#ifdef NO_RICK_FORTRAN
    if(save_plm)
      rick_shd2c_pre(data,(data+exp[0].npoints),exp[0].lmax,
		     *plm,(*plm+exp->n_plm),
		     ivec,exp[0].alm,exp[1].alm,&exp->rick);
    
    else
      rick_shd2c(data,(data+exp[0].npoints),exp[0].lmax,ivec,
		 exp[0].alm,exp[1].alm,&exp->rick);
#else
    if(save_plm)
      rick_f90_shd2c_pre(data,(data+exp[0].npoints),&exp[0].lmax,
		     *plm,(*plm+exp->n_plm),
		     &ivec,exp[0].alm,exp[1].alm);
    
    else
      rick_f90_shd2c(data,(data+exp[0].npoints),&exp[0].lmax,&ivec,
		 exp[0].alm,exp[1].alm);
#endif
    break;
#ifdef HC_USE_SPHEREPACK
  case SH_SPHEREPACK_GAUSS:
  case SH_SPHEREPACK_EVEN:
    break;
#endif
  default:
    sh_exp_type_error("sh_compute_spectral",exp);
    break;
  }
  exp[0].spectral_init = TRUE;
  if(ivec)
    exp[1].spectral_init = TRUE;
}
/* 


compute the spatial correspondence of a spherical harmonic expansion
depending on the implementation, this will save the Plm Legendre
factors

spatial output will be in data in the ordering required by the spherical
haronics method

input: exp (with parameters and alm) [1+ivec]

ivec: use vectors or not 

input/output: data
 
the DATA array has to be exp->npoints * (1 + ivec) = exp->tnpoints

*/
void 
sh_compute_spatial (exp, ivec, save_plm, plm, data, verbose)
struct sh_lms *exp;
int ivec;
hc_boolean save_plm;
SH_RICK_PREC **plm;
HC_PREC *data;
hc_boolean verbose;
{
  if((!exp[0].spectral_init)||(ivec && !exp[1].spectral_init)){
    fprintf(stderr,"sh_compute_spatial: coefficients set not initialized, ivec: %i\n",
	    ivec);
    exit(-1);
  }
  /* 
     make sure that the Legendre factors are computed 
  */
  if(save_plm){
    /* 
       this routine will only compute the plm once and performs 
       some sanity checks
    */
    sh_compute_plm(exp,ivec,plm,verbose); 
  }
  switch(exp->type){
#ifdef HC_USE_HEALPIX
  case SH_HEALPIX:
    if(ivec){
      HC_ERROR("sh_compute_spatial","healpix: ivec==1 not implemented yet");
    }else{
      /* 
	 scalar
      */
      if(save_plm){/* using precomputed Plm */
	heal_alm2map_sc_pre(&exp[0].heal.nside,&exp[0].lmax,
			    &exp[0].lmax,exp[0].alm_c,data,
			    *plm);
      }else{/* compute Plm now */
	heal_alm2map_sc(&exp[0].heal.nside,&exp[0].lmax,
			&exp[0].lmax,exp[0].alm_c,data);
      }
    }
    break;			/* end Healpix */
#endif
  case SH_RICK:
#ifdef NO_RICK_FORTRAN
    if(save_plm)
      rick_shc2d_pre(exp[0].alm,exp[1].alm,exp[0].lmax,
		     *plm,(*plm+exp->n_plm),
		     ivec,data,(data+exp[0].npoints),
		     &exp->rick);
    else
      rick_shc2d(exp[0].alm,exp[1].alm,exp[0].lmax,ivec,
		 data,(data+exp[0].npoints),&exp->rick);
#else
    if(save_plm)
      rick_f90_shc2d_pre(exp[0].alm,exp[1].alm,&exp[0].lmax,
		     *plm,(*plm+exp->n_plm),
		     &ivec,data,(data+exp[0].npoints));
    else
      rick_f90_shc2d(exp[0].alm,exp[1].alm,&exp[0].lmax,&ivec,
		 data,(data+exp[0].npoints));
#endif
    break;
#ifdef HC_USE_SPHEREPACK
  case SH_SPHEREPACK_GAUSS:
  case SH_SPHEREPACK_EVEN:
    break;
#endif
  default:
    sh_exp_type_error("sh_compute_spatial",exp);
    break;
  }
}

/* 

compute a spatial expansion on an regular grid given in theta and
phi arrays of npoints length

cos(theta) and phi are cos(colatitude) and longitude in radians

*/
void 
sh_compute_spatial_reg (exp, ivec, save_plm, plm, theta, ntheta, phi, nphi, data, verbose, save_sincos_fac)
struct sh_lms *exp;
int ivec;
hc_boolean save_plm;
SH_RICK_PREC **plm;
HC_PREC *theta;
int ntheta;
HC_PREC *phi;
int nphi;
HC_PREC *data;
hc_boolean verbose;
hc_boolean save_sincos_fac;
{
  int npoints;
  npoints = nphi * ntheta;
  if((!exp[0].spectral_init)||(ivec && !exp[1].spectral_init)){
    fprintf(stderr,"sh_compute_spatial_reg: coefficients set not initialized, ivec: %i\n",
	    ivec);
    exit(-1);
  }
  /* 
     make sure that the Legendre factors are computed 
  */
  if(save_plm){
    /* 
       this routine will only compute the plm once and performs 
       some sanity checks
    */
    sh_compute_plm_reg(exp,ivec,plm,verbose,theta,ntheta); 
  }
  switch(exp->type){
#ifdef HC_USE_HEALPIX
  case SH_HEALPIX:
    HC_ERROR("sh_compute_spatial_reg","healpix: not implemented yet");
#endif
  case SH_RICK:
#ifdef NO_RICK_FORTRAN
    if(save_plm)
      rick_shc2d_pre_reg(exp[0].alm,exp[1].alm,exp[0].lmax,
			   *plm,(*plm+exp->n_plm),
			   ivec,data,(data+npoints),
			   &exp->rick,theta,ntheta,phi,nphi,
			   save_sincos_fac);
    else
      rick_shc2d_reg(exp[0].alm,exp[1].alm,exp[0].lmax,ivec,
		       data,(data+npoints),&exp->rick,
		       theta,ntheta,phi,nphi,save_sincos_fac);
#else
    HC_ERROR("sh_compute_spatial_reg","Rick fortran not implemented");
#endif
    break;
#ifdef HC_USE_SPHEREPACK
  case SH_SPHEREPACK_GAUSS:
  case SH_SPHEREPACK_EVEN:
    HC_ERROR("sh_compute_spatial_reg","Spherepack not implemented");
    break;
#endif
  default:
    sh_exp_type_error("sh_compute_spatial",exp);
    break;
  }
}
void 
sh_compute_spatial_irreg (exp, ivec, theta, phi, npoints, data, verbose)
struct sh_lms *exp;
int ivec;
HC_PREC *theta;
HC_PREC *phi;
int npoints;
HC_PREC *data;
hc_boolean verbose;
{
  if((!exp[0].spectral_init)||(ivec && !exp[1].spectral_init)){
    fprintf(stderr,"sh_compute_spatial_irreg: coefficients set not initialized, ivec: %i\n",
	    ivec);
    exit(-1);
  }
  switch(exp->type){
#ifdef HC_USE_HEALPIX
  case SH_HEALPIX:
    HC_ERROR("sh_compute_spatial_irreg","healpix: not implemented yet");
#endif
  case SH_RICK:
#ifdef NO_RICK_FORTRAN
    rick_shc2d_irreg(exp[0].alm,exp[1].alm,exp[0].lmax,ivec,
		   data,(data+npoints),&exp->rick,
		   theta,phi,npoints);
#else
    HC_ERROR("sh_compute_spatial_irreg","Rick fortran not implemented");
#endif
    break;
#ifdef HC_USE_SPHEREPACK
  case SH_SPHEREPACK_GAUSS:
  case SH_SPHEREPACK_EVEN:
    HC_ERROR("sh_compute_spatial_irreg","Spherepack not implemented");
    break;
#endif
  default:
    sh_exp_type_error("sh_compute_spatial",exp);
    break;
  }
}

/* 
   
print an error message and exit if a spherical harmonics type is not
defined

*/
void 
sh_exp_type_error (subroutine, exp)
char *subroutine;
struct sh_lms *exp;
{
  fprintf(stderr,"%s: error: spherical harmonics type %i undefined\n",
	  subroutine,exp->type);
  exit(-1);
}
/* 

print the Plm factors

*/
void 
sh_print_plm (plm, n_plm, ivec, type, out)
SH_RICK_PREC *plm;
int n_plm;
int ivec;
int type;
FILE *out;
{
  int i,j,jlim;
  switch(type){
#ifdef HC_USE_HEALPIX
  case SH_HEALPIX:
    jlim=(ivec)?(3):(1);
    for(i=0;i < n_plm;i++){	/* number of points loop */
      for(j=0;j < jlim;j++)	/* scalar or scalar + pol? */
	fprintf(out,"%16.7e ",(double)plm[j*n_plm+i]);
      fprintf(out,"\n");
    }
    break;
#endif
  case SH_RICK:
    jlim=(ivec)?(2):(1);
    for(i=0;i < n_plm;i++){	/* number of points loop */
      for(j=0;j < jlim;j++)	/* scalar or scalar + pol? */
	fprintf(out,"%16.7e ",(double)plm[j*n_plm+i]);
      fprintf(out,"\n");
    }
    break;
#ifdef HC_USE_SPHEREPACK
  case SH_SPHEREPACK_GAUSS:
  case SH_SPHEREPACK_EVEN:
    break;
#endif
  default:
    fprintf(stderr,"sh_print_plm: expansion type %i undefined\n",
	    type);
    exit(-1);
  }
}
/* 

print the spatial equivalent of a spherical harmonic expansion with
npoints spatial points in the format

coordinates data

for 
use_3d: coordinates = lon lat z
else  : coordinates = lon lat

shps is the number of scalars that are passed in the data[shps * npoints] array

*/
void 
sh_print_spatial_data_to_stream (exp, shps, data, use_3d, z, out)
struct sh_lms *exp;
int shps;
HC_PREC *data;
hc_boolean use_3d;
HC_PREC z;
FILE *out;
{
  int j,k;
  HC_PREC lon,lat;
  for(j=0;j < exp[0].npoints;j++){
    /* 
       get coordinates
    */
    sh_get_coordinates(exp,j,&lon,&lat);
    /* print coordinates */
    if(!use_3d){
      /* print lon lat  */
      fprintf(out,"%14.7f %14.7f\t",(double)lon,(double)lat);
    }else{
      /* print lon lat z[i] */
      fprintf(out,"%14.7f %14.7f %14.7f\t",(double)lon,(double)lat,(double)z);
    }
    for(k=0;k < shps;k++)		/* loop through all scalars */
      fprintf(out,"%14.7e ",(double)data[j+exp[0].npoints*k]);
    fprintf(out,"\n");
  }	/* end points in lateral space loop */
}

void 
sh_get_coordinates (exp, i, lon, lat)
struct sh_lms *exp;
int i;
HC_PREC *lon;
HC_PREC *lat;
{
  HC_PREC xp[3];
  switch(exp->type){
#ifdef HC_USE_HEALPIX
    
  case SH_HEALPIX:
    switch(exp[0].heal.ordering){
    case SH_HEALPIX_RING:  
      pix2ang_ring((long)exp[0].heal.nside,
		   (long)i,(xp+HC_THETA),(xp+HC_PHI));
      break;
    case SH_HEALPIX_NEST:  
      pix2ang_nest((long)exp[0].heal.nside,
		   (long)i,(xp+HC_THETA),(xp+HC_PHI));
      break;
    default:
      fprintf(stderr,"print_sh_spatial_data: error: ordering %i undefined\n",
	      exp[0].heal.ordering);
      exit(-1);
      break;
    }
    break;			/* end Healpix branch */
#endif
  case SH_RICK:
    /* compute location */
#ifdef NO_RICK_FORTRAN
    rick_pix2ang(i,exp[0].lmax,(xp+HC_THETA),(xp+HC_PHI),
		 &exp->rick);
#else
    rick_f90_pix2ang(&i,&exp[0].lmax,(xp+HC_THETA),(xp+HC_PHI));
#endif
    break;
#ifdef HC_USE_SPHEREPACK
  case SH_SPHEREPACK_GAUSS:
  case SH_SPHEREPACK_EVEN:
    break;
#endif
  default:
    sh_exp_type_error("sh_print_spatial_data",exp);
    break;
  }	/* end type branch */
  *lon = PHI2LON(xp[HC_PHI]);
  *lat = THETA2LAT(xp[HC_THETA]);
}

/* 

regular grid version

*/
void 
sh_print_reg_spatial_data_to_stream (exp, shps, data, use_3d, z, theta, ntheta, phi, nphi, out)
struct sh_lms *exp;
int shps;
HC_PREC *data;
hc_boolean use_3d;
HC_PREC z;
HC_PREC *theta;
int ntheta;
HC_PREC *phi;
int nphi;
FILE *out;
{
  int i,j,k,l,npoints;
  HC_PREC lon,lat;
  npoints = nphi * ntheta;
  /* 
     get coordinates
  */
  switch(exp->type){
#ifdef HC_USE_HEALPIX
    
  case SH_HEALPIX:
    HC_ERROR("sh_print_reg_spatial_data_to_stream","healpix not implemented");
    break;		
#endif
  case SH_RICK:
    for(i=l=0;i<ntheta;i++){
      lat = THETA2LAT(theta[i]);
      for(j=0;j<nphi;j++,l++){
	lon = PHI2LON(phi[j]);
	/* print coordinates */
	if(!use_3d){
	  /* print lon lat  */
	  fprintf(out,"%12.3f %12.3f\t",(double)lon,(double)lat);
	}else{
	  /* print lon lat z[i] */
	  fprintf(out,"%12.3f %12.3f %14.7f\t",(double)lon,(double)lat,(double)z);
	}
	for(k=0;k<shps;k++)		/* loop through all scalars */
	  fprintf(out,"%14.7e ",(double)data[l+npoints*k]);
	fprintf(out,"\n");
      }
    }

    break;
#ifdef HC_USE_SPHEREPACK
  case SH_SPHEREPACK_GAUSS:
  case SH_SPHEREPACK_EVEN:
    HC_ERROR("sh_print_reg_spatial_data_to_stream","spherepack not implemented");
    break;
#endif
  default:
      sh_exp_type_error("sh_print_spatial_data",exp);
      break;
  }	/* end type branch */
}
/* 

irregular, arbitrary version

*/
void 
sh_print_irreg_spatial_data_to_stream (exp, shps, data, use_3d, z, theta, phi, npoints, out)
struct sh_lms *exp;
int shps;
HC_PREC *data;
hc_boolean use_3d;
HC_PREC z;
HC_PREC *theta;
HC_PREC *phi;
int npoints;
FILE *out;
{
  int i,k;
  HC_PREC lon,lat;
  /* 
     get coordinates
  */
  switch(exp->type){
#ifdef HC_USE_HEALPIX
  case SH_HEALPIX:
    HC_ERROR("sh_print_irreg_spatial_data_to_stream","healpix not implemented");
    break;		
#endif
  case SH_RICK:
    for(i=0;i<npoints;i++){
      lat = THETA2LAT(theta[i]);
      lon = PHI2LON(phi[i]);
      /* print coordinates */
      if(!use_3d){
	/* print lon lat  */
	fprintf(out,"%14.7f %14.7f\t",(double)lon,(double)lat);
      }else{
	/* print lon lat z[i] */
	fprintf(out,"%14.7f %14.7f %14.7f\t",(double)lon,(double)lat,(double)z);
      }
      for(k=0;k < shps;k++)		/* loop through all scalars */
	fprintf(out,"%14.7e ",(double)data[i+npoints*k]);
      fprintf(out,"\n");
    }
    break;
#ifdef HC_USE_SPHEREPACK
  case SH_SPHEREPACK_GAUSS:
  case SH_SPHEREPACK_EVEN:
    HC_ERROR("sh_print_irreg_spatial_data_to_stream","spherepack not implemented");
    break;
#endif
  default:
      sh_exp_type_error("sh_print_spatial_data",exp);
      break;
  }	/* end type branch */
}
/* 

compute the associated Legendre functions for all (l,m) at 
all latidutinal lcoations once and only once

input:
exp: holds the expansion parameters
ivec_global: if 1, will construct vector arrays, else only for scalar

output:

plm: will be re-allocated, has to be passed at least as NULL

*/
void 
sh_compute_plm (exp, ivec, plm, verbose)
struct sh_lms *exp;
int ivec;
SH_RICK_PREC **plm;
hc_boolean verbose;
{
  if(!exp->plm_computed){
    if((!exp->lmax)||(!exp->n_plm)||(!exp->tn_plm)){
      fprintf(stderr,"sh_compute_plm: error, expansion not initialized?\n");
      fprintf(stderr,"sh_compute_plm: lmax: %i n_plm: %i tn_plm: %i\n",
	      exp->lmax,exp->n_plm,exp->tn_plm);
      exit(-1);
    }
    /* 
       allocate 
    */
    rick_vecalloc(plm,exp->tn_plm,"sh_compute_plm");
    /* 
       compute the Legendre polynomials 
    */
    switch(exp->type){
#ifdef HC_USE_HEALPIX

    case SH_HEALPIX:
      if(verbose)
	fprintf(stderr,"sh_compute_plm: healpix: computing Plm for lmax %i\n",
		exp->lmax);
      heal_plmgen(*plm,&exp->heal.nside,&exp->lmax,&ivec);
      break;
#endif
    case SH_RICK:
      if(verbose)
	fprintf(stderr,"sh_compute_plm: Rick: computing all Plm for lmax %i\n",
		exp->lmax);
#ifdef NO_RICK_FORTRAN
      rick_compute_allplm(exp->lmax,ivec,*plm,
			  (*plm+exp->n_plm),&exp->rick);
#else
      rick_f90_compute_allplm(&exp->lmax,&ivec,*plm,
			      (*plm+exp->n_plm));
#endif
      break;
#ifdef HC_USE_SPHEREPACK
  case SH_SPHEREPACK_GAUSS:
  case SH_SPHEREPACK_EVEN:
    break;
#endif
    default:
      sh_exp_type_error("compute_plm",exp);
      break;
    }
    exp->plm_computed = TRUE;
    exp->old_lmax = exp->lmax;
    exp->old_ivec = ivec;
    exp->old_tnplm = exp->tn_plm;
  }else{
    /* 
       simple checks
       
       first lmax
    */
    if(exp->old_lmax != exp->lmax){
      fprintf(stderr,"sh_compute_plm: error: lmax initially %i, now %i\n",
	      exp->old_lmax,exp->lmax);
      exit(-1);
    }
    /* check if ivec was initialized if ever used  */
    if(ivec > exp->old_ivec){
      fprintf(stderr,"sh_compute_plm: error: plm are to be saved but routine was initialized\n");
      fprintf(stderr,"sh_compute_plm: error: with ivec: %i and now we want vectors, ivec: %i\n",
	      exp->old_ivec,ivec);
      exit(-1);
    }
    if(exp->tn_plm != exp->old_tnplm){
      fprintf(stderr,"sh_compute_plm: error: tn_plm initially %i, now %i\n",
	      exp->old_tnplm,exp->tn_plm);
      exit(-1);
    }
  }
}

/* 

compute the associated Legendre functions for all (l,m) at all
latidutinal lcoations once and only once for all regularly placed latitudes in theta

input:
exp: holds the expansion parameters
ivec_global: if 1, will construct vector arrays, else only for scalar

output:

plm: will be re-allocated, has to be passed at least as NULL

*/
void 
sh_compute_plm_reg (exp, ivec, plm, verbose, theta, npoints)
struct sh_lms *exp;
int ivec;
SH_RICK_PREC **plm;
hc_boolean verbose;
HC_PREC *theta;
int npoints;
{
  /*  */
  exp->tn_plm_irr = (1+ivec) * exp->lmsmall2 * npoints;

  if((!exp->plm_computed_irr)||(exp->tn_plm_irr != exp->old_tnplm_irr)){
    if((!exp->lmax)||(!exp->n_plm)||(!exp->tn_plm)){
      fprintf(stderr,"sh_compute_plm_reg: error, expansion not initialized?\n");
      fprintf(stderr,"sh_compute_plm_reg: lmax: %i n_plm: %i tn_plm: %i\n",
	      exp->lmax,exp->n_plm,exp->tn_plm);
      exit(-1);
    }
    /* 
       allocate 
    */
    exp->old_tnplm_irr =  exp->tn_plm_irr;
    rick_vecrealloc(plm,exp->old_tnplm_irr,"sh_compute_plm_reg");
    /* 
       compute the Legendre polynomials 
    */
    switch(exp->type){
#ifdef HC_USE_HEALPIX
      HC_ERROR("compute_plm_reg","healpix not implemented");
      break;
#endif
    case SH_RICK:
      if(verbose)
	fprintf(stderr,"sh_compute_plm_reg: Rick: computing all Plm for lmax %i and %i points\n",
		exp->lmax,npoints);
#ifdef NO_RICK_FORTRAN
      rick_compute_allplm_reg(exp->lmax,ivec,*plm,(*plm+exp->old_tnplm_irr),&exp->rick,
				theta,npoints);
#else
      HC_ERROR("compute_plm_reg","rick fortran not implemented");
#endif
      break;
#ifdef HC_USE_SPHEREPACK
  case SH_SPHEREPACK_GAUSS:
  case SH_SPHEREPACK_EVEN:
      HC_ERROR("compute_plm_reg","spherepack implemented");
    break;
#endif
    default:
      sh_exp_type_error("compute_plm",exp);
      break;
    }
    exp->plm_computed_irr = TRUE;
    exp->old_lmax_irr = exp->lmax;
    exp->old_ivec_irr = ivec;
  }else{
    /* 
       simple checks
       
       first lmax
    */
    if(exp->old_lmax_irr != exp->lmax){
      fprintf(stderr,"sh_compute_plm_reg: error: lmax initially %i, now %i\n",
	      exp->old_lmax_irr,exp->lmax);
      exit(-1);
    }
    /* check if ivec was initialized if ever used  */
    if(ivec > exp->old_ivec_irr){
      fprintf(stderr,"sh_compute_plm_reg: error: plm are to be saved but routine was initialized\n");
      fprintf(stderr,"sh_compute_plm_reg: error: with ivec: %i and now we want vectors, ivec: %i\n",
	      exp->old_ivec_irr,ivec);
      exit(-1);
    }
  }
}
/* 
   
returns the l,m A (use_b=0) or B (use_b=1) or A and B (use_b=2)
coefficient(s) of an expansion. if phys_normalization is set, 
will change the normalization so that the output is in the 
real spherical harmonics as in Dahlen and Tromp. 

else, will leave in the original convention, what ever it is 
stored in

for m==0, B will be returned as zero

*/
void  sh_get_coeff (exp, l, m, use_b, phys_norm, value)
     struct sh_lms *exp;
     int l;
     int m;
     int use_b;
     hc_boolean phys_norm;
     HC_CPREC *value;
{
#ifdef HC_USE_HEALPIX
  const HC_CPREC sqrt2 = SQRT_TWO;
#endif
  HC_CPREC s1;
  int index;
#ifdef HC_DEBUG
  /* do some checks */
  if((m > l)||(l<0)||(m<0)){
    fprintf(stderr,"sh_get_coeff: error: attempting to read l=%i and m=%i\n",
	    l,m);
    exit(-1);
  }
  if(l > exp->lmax){
    fprintf(stderr,"sh_get_coeff: error: attempting to read l=%i to lmax=%i expansion\n",
	    l,exp->lmax);
    exit(-1);
  }
  if((use_b < 0)||(use_b>2)){
    fprintf(stderr,"sh_get_coeff: error: attempting to read use_b=%i expansion (should be 0,1, or 2)\n",
	    use_b);
    exit(-1);
  }
#endif
  switch(exp->type){
#ifdef HC_USE_HEALPIX

  case SH_HEALPIX:
    /* 
       healpix is stored in complex number format, convert to real
       expansion coefficients
    */
    if(phys_norm){		/* convert */
      s1=(m==0)?(1.0):(sqrt2);
      if(use_b == 0)
	*value = (HC_CPREC) s1 * exp->alm_c[heal_index(&exp->heal,l,m)].dr;
      else if(use_b == 1)
	*value =   (m != 0)?((HC_CPREC)-s1 * exp->alm_c[heal_index(&exp->heal,l,m)].di):(0.0);
      else{			/* a and b needed */
	value[0] = (HC_CPREC) s1 * 
	  exp->alm_c[(index=heal_index(&exp->heal,l,m))].dr;
	value[1] = (m != 0)?((HC_CPREC)-s1 * exp->alm_c[index].di):(0.0);
      }
    }else{			/* leave as is  */
      if(use_b == 0)
	*value = (HC_CPREC)exp->alm_c[heal_index(&exp->heal,l,m)].dr;
      else if(use_b == 1){
	*value = (m != 0)?((HC_CPREC)exp->alm_c[heal_index(&exp->heal,l,m)].di):(0.0);
      }else{
	value[0] =  (HC_CPREC)exp->alm_c[(index=heal_index(&exp->heal,l,m))].dr;
	value[1] =  (m != 0)?((HC_CPREC)exp->alm_c[index].di):(0.0);
      }
    }
    break;
#endif
  case SH_RICK:
    if(use_b < 2){		/* A or B */
      if(use_b && (m==0)){
	*value = 0.0;
      }else{
	if(phys_norm){
	  *value = (HC_CPREC)exp->alm[LM_INDEX(l,m,use_b)] * 
	    SH_RICK_FACTOR(l,m);
	}else{
	  *value = (HC_CPREC)exp->alm[LM_INDEX(l,m,use_b)];
	}
      }
    }else{ 			/* both  */
      index = LM_INDEX(l,m,0);
      if(phys_norm){
	s1 = SH_RICK_FACTOR(l,m);
	value[0] = (HC_CPREC)exp->alm[index  ]* s1;
	value[1] = (m != 0)?((HC_CPREC)exp->alm[index+1]* s1):(0.0);
      }else{
	value[0] = (HC_CPREC)exp->alm[index];
	value[1] = (m != 0)?((HC_CPREC)exp->alm[index+1]):(0.0);
      }
    }
    break;
#ifdef HC_USE_SPHEREPACK
  case SH_SPHEREPACK_GAUSS:
  case SH_SPHEREPACK_EVEN:
    break;
#endif
  default:
    sh_exp_type_error("sh_read_coefficients_from_stream",exp);
    break;
  }
}
/*  
    
write value "value" to the expansion,. if phys_norm is set, will
convert from real spherical harmonics as in Dahlen and Tromp to
different internal convention

if use_b = 0: write A coefficient from value[0]
if use_b = 1: write B coefficient from value[0]
if use_b = 2: write A and B coefficient from value[0] and value[1]


*/
void sh_write_coeff (exp, l, m, use_b, phys_norm, value)
     struct sh_lms *exp;
     int l;
     int m;
     int use_b;
     hc_boolean phys_norm;
     HC_CPREC *value;
{
#ifdef HC_USE_HEALPIX
  const HC_CPREC sqrt2 = SQRT_TWO;
#endif
  HC_CPREC s1;
  int index;
#ifdef HC_DEBUG
  /* do some checks */
  if((m > l)||(l<0)||(m<0)){
    fprintf(stderr,"sh_write_coeff: error: attempting to write l=%i and m=%i\n",
	    l,m);
    exit(-1);
  }
  if(l > exp->lmax){
    fprintf(stderr,"sh_write_coeff: error: attempting to write l=%i to lmax=%i expansion\n",
	    l,exp->lmax);
    exit(-1);
  }
  if((use_b < 0)||(use_b>2)){
    fprintf(stderr,"sh_write_coeff: error: attempting to write use_b=%i expansion (should be 0,1, or 2)\n",
	    use_b);
    exit(-1);
  }
#endif
  if((m==0) && (use_b > 0)){
    fprintf(stderr,"sh_write_coeff: error: can't assign B coefficient for m == 0 (l: %i)\n",
	    l);
    exit(-1);
  }
  switch(exp->type){
#ifdef HC_USE_HEALPIX

  case SH_HEALPIX:
    /* 
       healpix is stored in complex number format, convert to real
       expansion coefficients
    */
    index = heal_index(&exp->heal,l,m);
    if(phys_norm){		/* convert */
      s1=(m==0)?(1.0):(sqrt2);
      if(use_b == 0)
	exp->alm_c[index].dr =   *value/s1;
      else if(use_b == 1)
	exp->alm_c[index].di = -(*value)/s1;
      else{
	exp->alm_c[index].dr =   value[0]/s1;
	exp->alm_c[index].di = - value[1]/s1;
      }
    }else{			/* use as is */
       if(use_b == 0)
	 exp->alm_c[index].dr = *value;
       else if(use_b == 1)
	 exp->alm_c[index].di = *value;
      else{
	exp->alm_c[index].dr = value[0];
	exp->alm_c[index].di = value[1];
      }
    }
    break;
#endif
  case SH_RICK:
    if(phys_norm){		/* convert */
      if(use_b < 2){		/* A or B */
	exp->alm[LM_INDEX(l,m,use_b)] = 
	  *value / SH_RICK_FACTOR(l,m);
      }else{			/* both */
	exp->alm[(index=LM_INDEX(l,m,0))] = 
	  value[0] / (s1=SH_RICK_FACTOR(l,m));
	exp->alm[index+1] = value[1] / s1;
      }
    }else{			/* as is */
      if(use_b < 2){		/* A or B */
	exp->alm[LM_INDEX(l,m,use_b)] = *value;
      }else{			/* both */
	exp->alm[(index = LM_INDEX(l,m,0))]   = value[0];
	exp->alm[index+1] = value[1];
      }
    } 
    break;
#ifdef HC_USE_SPHEREPACK
  case SH_SPHEREPACK_GAUSS:
  case SH_SPHEREPACK_EVEN:
    break;
#endif
  default:
    sh_exp_type_error("sh_write_coefficients",exp);
    break;
  }
}

void 
sh_add_coeff (exp, l, m, use_b, phys_norm, value)
struct sh_lms *exp;
int l;
int m;
int use_b;
hc_boolean phys_norm;
HC_CPREC *value;
{
#ifdef HC_USE_HEALPIX
  const HC_CPREC sqrt2 = SQRT_TWO;
#endif
  HC_CPREC s1;
  int index;
  if((m==0) && (use_b > 0)){
    fprintf(stderr,"sh_add_coeff: error: can't assign B coefficient for m == 0 (l: %i)\n",
	    l);
    exit(-1);
  }
  switch(exp->type){
#ifdef HC_USE_HEALPIX
  case SH_HEALPIX:
    /* 
       healpix is stored in complex number format, convert to real
       expansion coefficients
    */
    index = heal_index(&exp->heal,l,m);
    if(phys_norm){		/* convert */
      s1=(m==0)?(1.0):(sqrt2);
      if(use_b == 0)
	exp->alm_c[index].dr +=   *value/s1;
      else if(use_b == 1)
	exp->alm_c[index].di += -(*value)/s1;
      else{
	exp->alm_c[index].dr +=   value[0]/s1;
	exp->alm_c[index].di += - value[1]/s1;
      }
    }else{			/* use as is */
       if(use_b == 0)
	 exp->alm_c[index].dr += *value;
       else if(use_b == 1)
	 exp->alm_c[index].di += *value;
      else{
	exp->alm_c[index].dr += value[0];
	exp->alm_c[index].di += value[1];
      }
    }
    break;
#endif
  case SH_RICK:
    if(phys_norm){		/* convert */
      if(use_b < 2){		/* A or B */
	exp->alm[LM_INDEX(l,m,use_b)] += 
	  *value / SH_RICK_FACTOR(l,m);
      }else{			/* both */
	exp->alm[(index=LM_INDEX(l,m,0))] += 
	  value[0] / (s1=SH_RICK_FACTOR(l,m));
	exp->alm[index+1] += value[1] / s1;
      }
    }else{			/* as is */
      if(use_b < 2){		/* A or B */
	exp->alm[LM_INDEX(l,m,use_b)] += *value;
      }else{			/* both */
	exp->alm[(index = LM_INDEX(l,m,0))]   += value[0];
	exp->alm[index+1] += value[1];
      }
    } 
    break;
#ifdef HC_USE_SPHEREPACK
  case SH_SPHEREPACK_GAUSS:
  case SH_SPHEREPACK_EVEN:
    break;
#endif
  default:
    sh_exp_type_error("sh_add_coefficients",exp);
    break;
  }
}

/* 
   copy a whole expansion structure 
   a --> b, i.e.
   b = a


*/
void 
sh_copy_lms (a, b)
struct sh_lms *a;
struct sh_lms *b;
{

  b->type = a->type;
  b->lmax = a->lmax;
  b->spectral_init = a->spectral_init;
  b->lmaxp1 = a->lmaxp1;
  b->lmbig = a->lmbig;
  b->lmsmall2 = a->lmsmall2;
  b->n_lm = a->n_lm;
  b->n_plm = a->n_plm;
  b->tn_plm = a->tn_plm;
  b->tn_plm_irr = a->tn_plm_irr;;
  b->npoints = a->npoints;
  b->plm_computed = a->plm_computed;
  b->plm_computed_irr = a->plm_computed_irr;
  b->old_lmax = a->old_lmax;
  b->old_ivec = a->old_ivec;
  b->old_tnplm = a->old_tnplm;
  b->old_tnplm_irr = a->old_tnplm_irr;
  b->old_lmax_irr  = a->old_lmax_irr;
  b->old_ivec_irr  = a->old_ivec_irr;
  /* this still needs to be fixed */
  b->rick.nlat =  a->rick.nlat;
  b->rick.nlon =  a->rick.nlon;
  sh_aexp_equals_bexp_coeff(b,a);
  
}
/* 
   copy the coefficients of one expansion to another 

   a = b

*/
void 
sh_aexp_equals_bexp_coeff (a, b)
struct sh_lms *a;
struct sh_lms *b;
{
  int i;
  if(a->type != b->type){
    fprintf(stderr,"sh_aexp_equals_bexp_coeff: error: type mix (%i vs. %i) not implemented yet\n",
	    a->type,b->type);
    exit(-1);
  }
  /* 
     make sure lmax(b) <= lmax(a). if lmax(b) < lmax(a), the rest 
     of the coefficients in a will be set to zero

  */
  if(a->lmax < b->lmax){
    fprintf(stderr,"sh_aexp_equals_bexp_coeff: error: lmax(a): %i < lmax(b): %i\n",
	    a->lmax,b->lmax);
    exit(-1);
  }
  switch(a->type){
#ifdef HC_USE_HEALPIX

  case SH_HEALPIX:			/* nplm should reflect the lmax */
    for(i=0;i < b->n_lm;i++){
      a->alm_c[i].dr = b->alm_c[i].dr;
      a->alm_c[i].di = b->alm_c[i].di;
    }
    for(i=b->n_lm;i < a->n_lm;i++)
      a->alm_c[i].dr = a->alm_c[i].di = 0.0;
    break;
#endif
  case SH_RICK:
    for(i=0;i < b->n_lm;i++)
      a->alm[i] = b->alm[i];
    for(i=b->n_lm;i < a->n_lm;i++)
      a->alm[i] = 0.0;
    break;
#ifdef HC_USE_SPHEREPACK
  case SH_SPHEREPACK_GAUSS:
  case SH_SPHEREPACK_EVEN:
    break;
#endif
  default:
    sh_exp_type_error("sh_aexp_equals_bexp_coeff",a);
    break;
  }
}
/* sum two expansions, c = a + b */

void 
sh_c_is_a_plus_b_coeff (c, a, b)
struct sh_lms *c;
struct sh_lms *a;
struct sh_lms *b;
{
  int i;
  if((a->type != b->type)||(b->type != c->type)){
    fprintf(stderr,"sh_c_is_a_plus_b_coeff: error: type mix (%i vs. %i vs. %i) not implemented yet\n",
	    a->type,b->type,c->type);
    exit(-1);
  }
  /* 
     make sure lmax(b) = lmax(a) = lmax(c)

  */
  if((a->lmax != b->lmax)||(b->lmax != c->lmax)){
    fprintf(stderr,"sh_c_is_a_plus_b_coeff: error: a!=b || b != c lmax %i %i %i\n",
	    a->lmax,b->lmax,c->lmax);
    exit(-1);
  }
  switch(a->type){
#ifdef HC_USE_HEALPIX

  case SH_HEALPIX:			/* nplm should reflect the lmax */
    for(i=0;i < b->n_lm;i++){
      c->alm_c[i].dr = a->alm_c[i].dr + b->alm_c[i].dr;
      c->alm_c[i].dr = a->alm_c[i].di + b->alm_c[i].di;
    }
    break;
#endif
  case SH_RICK:
    for(i=0;i < b->n_lm;i++)
      c->alm[i] = a->alm[i] + b->alm[i];
    break;
#ifdef HC_USE_SPHEREPACK
  case SH_SPHEREPACK_GAUSS:
  case SH_SPHEREPACK_EVEN:
    break;
#endif
  default:
    sh_exp_type_error("sh_c_is_a_plus_b_coeff",a);
    break;
  }
}
/* 

given a spherical harmonic expansion, scale it with factor fac[0...lmax] 
which only depends on l


*/
void 
sh_scale_expansion_l_factor (exp, lfac)
struct sh_lms *exp;
HC_CPREC *lfac;
{
  int l,m,index;
  HC_CPREC fac;
  switch(exp->type){
#ifdef HC_USE_HEALPIX

  case SH_HEALPIX:
    for(l=0;l <= exp->lmax;l++){
      fac = lfac[l];
      for(m=0;m <= l;m++){
	exp->alm_c[(index = heal_index(&exp->heal,l,m))].dr * fac;
	exp->alm_c[index].di *= fac;
      }
    }
    break;
#endif
  case SH_RICK:
    for(l=0;l <= exp->lmax;l++){
      fac = lfac[l];
      for(m=0;m <= l;m++){
	exp->alm[(index = LM_INDEX(l,m,0))] *= fac;
	exp->alm[index+1] *= fac;
      }
    }
    break;
#ifdef HC_USE_SPHEREPACK
  case SH_SPHEREPACK_GAUSS:
  case SH_SPHEREPACK_EVEN:
    break;
#endif
  default:
    sh_exp_type_error("sh_scale_expansion_l_factor",exp);
    break;
 }
}


/* 
   scale all coefficients 
*/
void 
sh_scale_expansion (exp, fac)
struct sh_lms *exp;
HC_CPREC fac;
{
  int l,m,index;
  switch(exp->type){
#ifdef HC_USE_HEALPIX

  case SH_HEALPIX:
    for(l=0;l <= exp->lmax;l++)
      for(m=0;m <= l;m++){
	exp->alm_c[(index = heal_index(&exp->heal,l,m))].dr * fac;
	exp->alm_c[index].di *= fac;
      }
    break;
#endif
  case SH_RICK:
    for(l=0;l <= exp->lmax;l++)
      for(m=0;m <= l;m++){
	exp->alm[(index = LM_INDEX(l,m,0))] *= fac;
	exp->alm[index+1] *= fac;
      }
    break;
#ifdef HC_USE_SPHEREPACK
  case SH_SPHEREPACK_GAUSS:
  case SH_SPHEREPACK_EVEN:
    break;
#endif
  default:
    sh_exp_type_error("sh_scale_expansion",exp);
    break;
 }
}
