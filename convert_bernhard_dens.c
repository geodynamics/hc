#include "hc.h"
/* 

convert Bernhard's density files

*/
int main(int argc, char **argv)
{
  int i,j,k,lmax,nlat,nlon,nr;
  SH_RICK_PREC *z,*w,r,dx,lon,dr,val;
  FILE *in;
  if(argc != 3){
    fprintf(stderr,"%s den-file lmax\n",argv[0]);
    exit(-1);
  }
  sscanf(argv[2],"%i",&lmax);
  nlat = lmax+1;nlon = nlat * 2;
  z = (SH_RICK_PREC *)malloc(sizeof(SH_RICK_PREC)*nlat);
  w = (SH_RICK_PREC *)malloc(sizeof(SH_RICK_PREC)*nlat);
  rick_gauleg(1.0,-1.0,z,w,nlat);
  for(i=0;i < nlat;i++)		/* convert theta colat to lat[deg] */
    z[i] = THETA2LAT(acos(z[i]));
  fprintf(stderr,"%s: reading from %s, assuming lmax %i and %i points in latitude\n",
	  argv[0],argv[1],lmax,nlat);
  
  if((in = fopen(argv[1],"r"))==NULL){
      fprintf(stderr,"%s: error, cannot open %s\n",argv[0],argv[1]);
      exit(-1);
  }
  fscanf(in,SH_RICK_FLT_FMT,&dx);
  fscanf(in,"%i",&nr);
  i = 0;r = 0.55;dr = (0.99-0.55)/(SH_RICK_PREC)(nr-1);
  j = 0;			/* latitude */
  k = 0;lon = 0.;		/* longitude */
  while(fscanf(in,SH_RICK_FLT_FMT,&val)==1){
    fprintf(stdout,"%11g %11g %11g %11g\n",lon,z[j],HC_Z_DEPTH(r),val);
    k++;lon+=dx;
    if(k == nlon){
      k=0;lon=0;j++;
    }
    if(j == nlat){
      j=0;i++;r+=dr;
    }
  }
  fprintf(stderr,"done: r: %i/%i\n",i,nr);
  fclose(in);
  free(z);free(w);
  return 0;
}

