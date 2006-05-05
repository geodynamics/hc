#include "hc.h"

/* 

miscellaneous functions for allocating arrays, copying arrays, opening
file streams, the like


$Id: hc_misc.c,v 1.5 2004/12/20 05:18:12 becker Exp becker $

*/



/* 
   open a file safely and give an error message if there was
   a problem
*/
FILE *hc_open(char *name, char *mode, char *program)
{
  FILE *in;
  if((in=fopen(name,mode)) == NULL){
    fprintf(stderr,"%s: error: can not open file %s for mode %s access\n",
	    program,name,mode);
    exit(-1);
  }
  return in;
}
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
void hc_scmplx_vecalloc(struct scmplx **x,int n,char *message)
{
  *x = (struct scmplx *)malloc(sizeof(struct scmplx)*(size_t)n);
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
  static hc_boolean init=FALSE;
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
/*

  calculate mean and standard deviation of x_i
  if hypoth is set, calc mean and stddev of sqrt(x_i^2 + y_i^2)
  if weighted is set, uses weights[n] for weighting the mean and so on

  this way of calculating the stddev is inaccurate but fast

*/
void hc_calc_mean_and_stddev(HC_PREC *x, HC_PREC *y,int n,HC_PREC *mean,
			     HC_PREC *stddev,HC_PREC *rms, 
			     hc_boolean hypoth, 
			     hc_boolean weighted,HC_PREC *weight)
{
  HC_PREC sum1=0.0,sum2=0.0,tmp,ws;
  int i;
  if(n <= 1){
    fprintf(stderr,"hc_calc_mean_and_stddev: error: n: %i\n",n);
    exit(-1);
  }
  ws=0.0;
  if(hypoth){// sqrt(x^2+y+2)
    for(i=0;i<n;i++){
      if(weighted){
	tmp = hypot(x[i],y[i]) * weight[i];
	ws += weight[i];
      }else{
	tmp = hypot(x[i],y[i]);
	ws += 1.0;
      }
      sum1 += tmp;sum2 += tmp * tmp;
    }
  }else{
    for(i=0;i<n;i++){
      if(weighted){
	tmp  = x[i] * weight[i];
	ws += weight[i];
      }else{
	tmp = x[i];
	ws += 1.0;
      }
      sum1 += tmp;sum2 += tmp*tmp;
    }
  }
  // standard deviation
  tmp = (ws * sum2 - sum1 * sum1) / (ws*(ws-1.0));
  if(tmp > 0)
    *stddev = sqrt(tmp);
  else
    *stddev = 0.0;
  *rms  = sqrt(sum2 / ws);// RMS 
  *mean = sum1 / ws;      // mean 
}
/*

  sort array by order, returns index

  GIVE INPUT INDEX SHIFTED
  (ie.    x-1)
  output will be 0...n-1

  based on numerical recipes

  
*/


#define HC_INDEXX_SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define HC_INDEXX_M 7
#define HC_INDEXX_NSTACK 5000

void hc_indexx(int n,HC_PREC *arr, int *indx)
{
  int i,indxt,ir=n,itemp,j,k,l=1;
  int jstack=0,*istack;
  HC_PREC a;
  //
  istack=(int *)malloc(sizeof(int)*(HC_INDEXX_NSTACK+1));
  if(!istack)
    HC_MEMERROR("indexx");
  //
  for (j=1;j<=n;j++) 
    indx[j]=j;
  for (;;) {
    if (ir-l < HC_INDEXX_M) {
      for (j=l+1;j<=ir;j++) {
	indxt=indx[j];
	a=arr[indxt];
	for (i=j-1;i>=1;i--) {
	  if (arr[indx[i]] <= a) break;
	  indx[i+1]=indx[i];
	}
	indx[i+1]=indxt;
      }
      if (jstack == 0) break;
      ir=istack[jstack--];
      l=istack[jstack--];
    } else {
      k=(l+ir) >> 1;
      HC_INDEXX_SWAP(indx[k],indx[l+1]);
      if (arr[indx[l+1]] > arr[indx[ir]]) {
	HC_INDEXX_SWAP(indx[l+1],indx[ir]);
      }
      if (arr[indx[l]] > arr[indx[ir]]) {
	HC_INDEXX_SWAP(indx[l],indx[ir]);
      }
      if (arr[indx[l+1]] > arr[indx[l]]) {
	HC_INDEXX_SWAP(indx[l+1],indx[l]);
      }
      i=l+1;
      j=ir;
      indxt=indx[l];
      a=arr[indxt];
      for (;;) {
	do i++; while (arr[indx[i]] < a);
	do j--; while (arr[indx[j]] > a);
	if (j < i) break;
	HC_INDEXX_SWAP(indx[i],indx[j]);
      }
      indx[l]=indx[j];
      indx[j]=indxt;
      jstack += 2;
      if (jstack > HC_INDEXX_NSTACK) {
	fprintf(stderr,
		"indexx: HC_INDEXX_NSTACK too small");
	exit(-1);
      }
      if (ir-i+1 >= j-l) {
	istack[jstack]=ir;
	istack[jstack-1]=i;
	ir=j-1;
      } else {
	istack[jstack]=j-1;
	istack[jstack-1]=l;
	l=i;
      }
    }
  }
  free(istack);
  // go back to normal numbering
  for(i=1;i<=n;i++)
    indx[i]--;
}
#undef HC_INDEXX_M
#undef HC_INDEXX_NSTACK
#undef HC_INDEXX_SWAP

