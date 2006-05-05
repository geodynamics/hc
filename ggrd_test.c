#include "hc.h"

int main(int argc, char **argv)
{
  struct ggrd_vel *v;
  HC_PREC xloc[3],time,vr[4],vphi[4],vtheta[4];
  static int order = 3;
  hc_boolean calc_derivatives ;
  double lon,lat,z;
  /* 
     initialize velocity structure
  */
  v = (struct ggrd_vel *)calloc(1,sizeof(struct ggrd_vel));
  ggrd_init_vstruc(v);
  /* 
     read in velocities 
  */
  ggrd_read_vel_grids(v,1.0,TRUE,TRUE,"");

  time = 0.0;
  calc_derivatives = FALSE;

  while(fscanf(stdin,"%lf %lf %lf",&lon,&lat,&z)==3){
    xloc[HC_R] = HC_ND_RADIUS(z);
    xloc[HC_THETA] = LAT2THETA(lat);
    xloc[HC_PHI] = LON2PHI(lon);
    /* 
       interpolate
    */
    ggrd_find_vel_and_der(xloc,time,v,order,calc_derivatives,
			  TRUE,vr,vtheta,vphi);
    fprintf(stdout,"%11g %11g %11g\n",vr[0],vtheta[0],vphi[0]);
  }

  return 0;
}
