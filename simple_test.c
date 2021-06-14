#include <stdio.h>
#include <math.h>


#define vshd2c vshd2c_
#define gauleg gauleg_

extern void vshd2c(int *,double *,double *,double *,double *,int *,
		   double *,double *);
extern void gauleg(double *,double *,double*,double *,int*);


void 
main ()
{
  int lmax,nlat,nlon,lmsize,i,j,index;
  double *datax,*datay,*cpol,*ctor,*z,*w,unity,negunity;

  unity = 1.0;negunity=-1.0;
  /*  */
  lmax = 63;
  /* sizes */
  nlat = lmax+1;
  nlon = nlat * 2;
  lmsize = (lmax+1)*(lmax+2)/2;
  /* allocate */
  hc_dvecalloc(&cpol,lmsize*2,"");
  hc_dvecalloc(&ctor,lmsize*2,"");
  hc_dvecalloc(&datax,nlon*nlat,"");
  hc_dvecalloc(&datay,nlon*nlat,"");

  hc_dvecalloc(&z,nlat,"");
  hc_dvecalloc(&w,nlon*nlat,"");
  
  /* gauss points */
  gauleg(&negunity,&unity,z,w,&nlat);
  /* read in data */
  for(index=0;index<nlon*nlat;index++){
    if(fscanf(stdin,"%*lf %*lf %*lf %lf %lf",
	      (datay+index),(datax+index))!=2){
      fprintf(stderr,"read error data\n");
      exit(-1);
    }
  }
  /* get coefficients */
  vshd2c(&lmax,z,w,cpol,ctor,&lmax,datax,datay);
  /* 
     output of pol/tor components, flip sign since Gauleq is now
     called -1 to 1
  */
  fprintf(stderr,"%i\n",lmax);
  for(i=0;i<lmsize*2;i++)
    fprintf(stderr,"%i %15.7e %15.7e\t%15.7e %15.7e\t%15.7e %15.7e\n",i,
	    0.,0.,-cpol[i*2],-cpol[i*2+1],-ctor[i*2],-ctor[i*2+1]);
}


