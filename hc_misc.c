#include "hc.h"

/* 

miscellaneous functions for allocating arrays, copying arrays, opening
file streams, the like


$Id: hc_misc.c,v 1.5  2004/12/20 05:18:12 becker Exp becker $

*/



/* double vector allocation */
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
/* double vector reallocation */
void hc_dvecrealloc(double **x,int n,char *message)
{
  *x = (double *)realloc(*x,sizeof(double)*(size_t)n);
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
float hc_svec_rms_diff(float *a,float *b,int n)
{
  int i;
  float sum=0.0,tmp;
  for(i=0;i<n;i++){
    tmp = a[i] - b[i];
    sum += tmp*tmp;
  }
  return sqrt(sum/n);
}
float hc_svec_rms(float *a,int n)
{
  int i;
  float sum=0.0;
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

/* zero a double precision vector */
void hc_zero_dvector(double *x, int n)
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
  static char type[2];
  int i;
  if(!init){
    if(sizeof(HC_PREC) == sizeof(float)){
      sprintf(type,"f");
    }else if (sizeof(HC_PREC) == sizeof(double)){
      sprintf(type,"lf");
    }else{
      fprintf(stderr,"hc_get_flt_frmt_string: assignment error\n");
      exit(-1);
    }
    init=TRUE;
  }
  if(!append)
    sprintf(string,"%%%s",type);
  for(i=1;i<n;i++)
    sprintf(string,"%s %%%s",string,type);
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
