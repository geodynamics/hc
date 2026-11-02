/* 


read in two grid files and compute the correlation of
points within a certain radius around each location


*/
#include "ggrd_grdtrack_util.h"

int main(int argc, char **argv)
{
  struct ggrd_gt g1[1],g2[1];
  char char_dummy;
  float x,y,x1,y1,dx,dy,rad_km;
  int n;
  double *z1,*z2;
  z1=(double *)malloc(sizeof(double));
  z2=(double *)malloc(sizeof(double));

  rad_km = 500;
  if(argc < 3){
    fprintf(stderr,"%s: usage:\n\t%s file1.grd file2.grd [radius, %g km]\n\n",
	    argv[0],argv[0],rad_km);
    exit(-1);
  }
  if(argc > 3 )
    sscanf(argv[3],"%f",&rad_km);

  fprintf(stderr,"%s: assuming grids %s and %s are global, computing for radius %g km\n",
	  argv[0],argv[1],argv[2],rad_km);
  
  if(ggrd_grdtrack_init_general(FALSE,argv[1],&char_dummy,
				"-Lx",g1,TRUE,FALSE,FALSE)){
    fprintf(stderr,"%s: error reading %s\n",argv[0],argv[1]);
    exit(-1);
  }
  if(ggrd_grdtrack_init_general(FALSE,argv[2],&char_dummy,
				"-Lx",g2,TRUE,FALSE,FALSE)){
    fprintf(stderr,"%s: error reading %s\n",argv[0],argv[2]);
    exit(-1);
  }

  /* 

  compute total correlation and best-fit slopes

  */
  dy = 2;
  for(n=0,y1 = 90-dy/2+1e-7,y = -90+dy/2;y <= y1;y += dy){
    dx = dy/cos(y/GGRD_PIF);	/* adjust for sphericity */
    for(x=0,x1=360-dx+1e-7;x<=x1;x+=dx,n++){
      ggrd_grdtrack_interpolate_xy((double)x,(double)y,g1,(z1+n),FALSE);
      ggrd_grdtrack_interpolate_xy((double)x,(double)y,g2,(z2+n),FALSE);
      
      fprintf(stderr,"%g %g %g %g %.6e\n",x,y,z1[n],z2[n],fabs(z1[n]-z2[n]));
      z1=(double *)realloc(z1,(n+2)*sizeof(double));
      z2=(double *)realloc(z2,(n+2)*sizeof(double));
    }
  }



}
