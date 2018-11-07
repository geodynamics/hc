#include "hc.h"
/* 

print Gauss latitudes

 */
int main(int argc, char **argv)
{
  int i,nlat;
  SH_RICK_PREC *z,*w;
  if(argc != 2){
    fprintf(stderr,"%s nlat - print Gauss latitudes\n",argv[0]);
    exit(-1);
  }
  sscanf(argv[1],"%i",&nlat);
  z = (SH_RICK_PREC *)malloc(sizeof(SH_RICK_PREC)*nlat);
  w = (SH_RICK_PREC *)malloc(sizeof(SH_RICK_PREC)*nlat);
  /*  */
  rick_gauleg(-1.0,1.0,z,w,nlat);
  for(i=0;i < nlat;i++)		/* convert theta colat to lat[deg] */
    fprintf(stdout,"%.8e\n",THETA2LAT(acos(z[i])));
  
  free(z);free(w);
  return 0;
}

