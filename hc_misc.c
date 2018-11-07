#include "hc.h"

/* 

miscellaneous functions for allocating arrays, copying arrays, opening
file streams, the like


$Id: hc_misc.c,v 1.5  2004/12/20 05:18:12 becker Exp becker $

*/



/* high precision vector allocation */
void hc_hvecalloc(HC_HIGH_PREC **x,int n,char *message)
{
  *x = (HC_HIGH_PREC *)malloc(sizeof(HC_HIGH_PREC)*(size_t)n);
  if(! (*x))
    HC_MEMERROR(message);
}
/* double */
void hc_dvecalloc(double **x,int n,char *message)
{
  *x = (double *)malloc(sizeof(double)*(size_t)n);
  if(! (*x))
    HC_MEMERROR(message);
}
/* single prec vector allocation */
void hc_svecalloc(float **x,int n,char *message)
{
  *x = (float *)malloc(sizeof(float)*(size_t)n);
  if(! (*x))
    HC_MEMERROR(message);
}
/* integer vector allocation */
void hc_ivecalloc(int **x,int n,char *message)
{
  *x = (int *)malloc(sizeof(int)*(size_t)n);
  if(! (*x))
    HC_MEMERROR(message);
}


/* general floating point vector allocation */
void hc_vecalloc(HC_PREC **x,int n,char *message)
{
  *x = (HC_PREC *)malloc(sizeof(HC_PREC)*(size_t)n);
  if(! (*x))
    HC_MEMERROR(message);
}
/* single prec complex vector allocation */
void hc_scmplx_vecalloc(struct hc_scmplx **x,int n,char *message)
{
  *x = (struct hc_scmplx *)malloc(sizeof(struct hc_scmplx)*(size_t)n);
  if(! (*x))
    HC_MEMERROR(message);
}
/* single vector reallocation */
void hc_svecrealloc(float **x,int n,char *message)
{
  *x = (float *)realloc(*x,sizeof(float)*(size_t)n);
  if(!(*x))
    HC_MEMERROR(message);
}
/* HC_HIGH_PREC vector reallocation */
void hc_dvecrealloc(HC_HIGH_PREC **x,int n,char *message)
{
  *x = (HC_HIGH_PREC *)realloc(*x,sizeof(HC_HIGH_PREC)*(size_t)n);
  if(!(*x))
    HC_MEMERROR(message);
}
/* general version */
void hc_vecrealloc(HC_PREC **x,int n,char *message)
{
  *x = (HC_PREC *)realloc(*x,sizeof(HC_PREC)*(size_t)n);
  if(!(*x))
    HC_MEMERROR(message);
}
/* 
   sqrt(sum(squared diff)) between [n] vectors a and b
*/
float hc_vec_rms_diff(HC_PREC *a,HC_PREC *b,int n)
{
  int i;
  HC_PREC sum=0.0,tmp;
  for(i=0;i<n;i++){
    tmp = a[i] - b[i];
    sum += tmp*tmp;
  }
  return sqrt(sum/n);
}
float hc_vec_rms(HC_PREC *a,int n)
{
  int i;
  HC_PREC sum=0.0;
  for(i=0;i<n;i++){
    sum += a[i] * a[i];
  }
  return sqrt(sum/n);
}
/* a[n] = b[n], single precision version */
void hc_a_equals_b_svector(float *a,float *b,int n)
{
  memcpy(a,b,sizeof(float)*n);
}
/* general version */
void hc_a_equals_b_vector(HC_PREC *a,HC_PREC *b,int n)
{
  memcpy(a,b,sizeof(HC_PREC)*n);
}
/* compute the mean of a single precision vector */
float hc_mean_svec(float *x,int n)
{
  float sum=0.0;
  int i;
  for(i=0;i<n;i++){
    sum += x[i];
  }
  sum /= (float)n;
  return sum;
}
/* compute the mean of a vector */
HC_PREC hc_mean_vec(HC_PREC *x,int n)
{
  HC_PREC sum=0.0;
  int i;
  for(i=0;i<n;i++){
    sum += x[i];
  }
  sum /= (HC_PREC)n;
  return sum;
}

/* zero a HC_HIGH_PREC precision vector */
void hc_zero_dvector(HC_HIGH_PREC *x, int n)
{
  int i;
  for(i=0;i<n;i++)
    x[i] = 0.0;
}
/* zero a vector of type logical */
void hc_zero_lvector(hc_boolean *x, int n)
{
  int i;
  for(i=0;i<n;i++)
    x[i] = FALSE;
}
/* 

assign floating point formats to a string as used by sscanf 
of fscanf

if append is TRUE, will add a format string, else will create 
anew 

*/
void hc_get_flt_frmt_string(char *string, int n, 
			    hc_boolean append)
{
  static hc_boolean init=FALSE;	/* that's OK, multiple instances calling are fine */
  static char type_s[3];
  int i;
  if(!init){
    if(sizeof(HC_PREC) == sizeof(float)){
      sprintf(type_s,"f");
    }else if (sizeof(HC_PREC) == sizeof(double)){
      sprintf(type_s,"lf");
    }else if (sizeof(HC_PREC) == sizeof(long double)){
      sprintf(type_s,"lf");
    }else{
      fprintf(stderr,"hc_get_flt_frmt_string: assignment error\n");
      exit(-1);
    }
    init=TRUE;
  }
  if(!append)
    sprintf(string,"%%%s",type_s);
  for(i=1;i<n;i++)
    sprintf(string,"%s %%%s",string,type_s);
}
//
// deal with boolean values/switches
char *hc_name_boolean(hc_boolean value)
{
  if(value)
    return("ON");
  else
    return("OFF");
}

hc_boolean hc_toggle_boolean(hc_boolean *variable)
{
  if(*variable){
    *variable=FALSE;
    return(FALSE);
  }else{
    *variable=TRUE;
    return(TRUE);
  }
}
//
// check, if we can read a value for the option flag in a chain of command line
// arguments
//
void hc_advance_argument(int *i,int argc, char **argv)
{
  if(argc <= *i + 1){// no arguments left
    fprintf(stderr,"%s: input parameters: error: option \"%s\" needs a value\n",
	    argv[0],argv[*i]);
    exit(-1);
  }
  *i += 1;
}

/* 
   compute the correlation between two scalar fields
   mode 0 : full correlation
        1 : up to 20, and between 4...9
        2 : up to 20, and between 4...9, and between 2...4

 */
void hc_compute_correlation(struct sh_lms *g1,struct sh_lms *g2,
			    HC_PREC *c,int mode,hc_boolean verbose)
{
  int lmaxg;
  lmaxg = HC_MIN(g1->lmax,g1->lmax);

  switch(mode){
  case 0:			/* 1...LMAX */
    if(verbose)
      fprintf(stderr,"hc_compute_correlation: computing 1...%i\n",lmaxg);
    c[0] = sh_correlation_per_degree(g1,g2,1,lmaxg);    
    break;
  case 1:			/* 1...20 and 4..9 correlations */
    lmaxg = HC_MIN(20,lmaxg);
    if(verbose)
      fprintf(stderr,"hc_compute_correlation: computing 1...%i and 4..9 correlations\n",lmaxg);
    c[0] = sh_correlation_per_degree(g1,g2,1,lmaxg);
    c[1] = sh_correlation_per_degree(g1,g2,4,9);
    break;
  case 2:			/* 1...20, 4..9, 2..4 correlations */
    lmaxg = HC_MIN(20,lmaxg);
    if(verbose)
      fprintf(stderr,"hc_compute_correlation: computing 1...%i and 4..9 correlations\n",lmaxg);
    c[0] = sh_correlation_per_degree(g1,g2,1,lmaxg);
    c[1] = sh_correlation_per_degree(g1,g2,4,9);
    c[2] = sh_correlation_per_degree(g1,g2,2,4);
    break;
  default:
    fprintf(stderr,"sh_compute_correlation: mode %i undefined\n",mode);
    exit(-1);
  }
}

/* 
   convert polar vector in r,theta,phi format to cartesian 
   vector x 

*/

void lonlatpv2cv(HC_PREC lon, float lat, 
		 HC_PREC *polar_vec,HC_PREC *cart_vec)
{
  HC_PREC theta,phi;
  theta = LAT2THETA((HC_HIGH_PREC)lat);
  phi   = LON2PHI((HC_HIGH_PREC)lon);
  thetaphipv2cv(theta,phi,polar_vec,cart_vec);
}
/* theta, phi version */
void thetaphipv2cv(HC_PREC theta, float phi, 
		   HC_PREC *polar_vec,HC_PREC *cart_vec)
{
  HC_HIGH_PREC polar_base[9];
  calc_polar_base_at_theta_phi(theta,phi,polar_base);
  lonlatpv2cv_with_base(polar_vec,polar_base,cart_vec);
}

void lonlatpv2cv_with_base(HC_PREC *polar_vec,
			   HC_HIGH_PREC *polar_base,
			   HC_PREC *cart_vec)
{
  int i;
  // convert vector
  for(i=0;i<3;i++){
    cart_vec[i]  = polar_base[i]   * polar_vec[0]; /* r,theta,phi */
    cart_vec[i] += polar_base[3+i] * polar_vec[1];
    cart_vec[i] += polar_base[6+i] * polar_vec[2];
  }
}
void calc_polar_base_at_theta_phi(HC_PREC theta, HC_PREC phi, 
				  HC_HIGH_PREC *polar_base)
{
  HC_HIGH_PREC cp,sp,ct,st;
  // base vecs
  ct=cos((HC_HIGH_PREC)theta);cp=cos((HC_HIGH_PREC)phi);
  st=sin((HC_HIGH_PREC)theta);sp=sin((HC_HIGH_PREC)phi);
  //
  polar_base[0]= st * cp;
  polar_base[1]= st * sp;
  polar_base[2]= ct;
  //
  polar_base[3]= ct * cp;
  polar_base[4]= ct * sp;
  polar_base[5]= -st;
  //
  polar_base[6]= -sp;
  polar_base[7]=  cp;
  polar_base[8]= 0.0;
}
/* 

   given a sorted vector y (y0<y1<...<yn-1) with n elements, return
   the weights f1 for element i1 and weight f2 for element i2 such
   that the interpolation is at y1

 */
void hc_linear_interpolate(HC_PREC *y, int n, HC_PREC y1,
			   int *i1, int *i2,
			   HC_PREC *f1, HC_PREC *f2)
{
  int n1;
  n1 = n-1;
  *i2 = 0;
  while((*i2 < n1) && (y[*i2] < y1))
    *i2 += 1;
  if(*i2 == 0)
    *i2 = 1;
  *i1 = *i2 - 1;
  *f2 = (y1 - y[*i1])/(y[*i2]-y[*i1]);
  *f1 = 1.0 - *f2;
}
