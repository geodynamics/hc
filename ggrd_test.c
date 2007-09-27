#include "hc.h"

int main(int argc, char **argv)
{
  struct ggrd_vel *v;
  HC_PREC xloc[3],time,vr[4],vphi[4],vtheta[4],dtrange;
  static int order = 3;
  hc_boolean calc_derivatives ;
  double lon,lat,age;
  /* 
     initialize velocity structure
  */
  v = (struct ggrd_vel *)calloc(1,sizeof(struct ggrd_vel));
  ggrd_init_vstruc(v);
  v->history = TRUE;		/* expect history */
  v->use_age = TRUE;		/* expect seafloor age files */
  v->age_bandlim = 900;
  /* 
     read in velocities 
  */
  if(ggrd_read_vel_grids(v,1.0,FALSE,TRUE,"/home/walter/becker/data/plates/past/clb/hall/")){
    fprintf(stderr,"error opening grids\n");
    exit(-1);
  }

  if(argc>1)
    sscanf(argv[1],"%lf",&time);
  else
    time = 0.0;
  dtrange = 1.0;			/* transition width, in Ma */

  fprintf(stderr,"%s: using time %g\n",argv[0],time);

  calc_derivatives = FALSE;

  xloc[HC_R] = HC_ND_RADIUS(0.0);
  for(lat=-89;lat<=89;lat+=2)
    for(lon=0;lon<=358;lon+=2){
  //lon=270;lat=-15;{
      xloc[HC_THETA] = LAT2THETA(lat);
      xloc[HC_PHI] = LON2PHI(lon);
      /* 
	 interpolate
      */
      if(ggrd_find_vel_and_der(xloc,time,dtrange,v,order,calc_derivatives,
			       TRUE,vr,vtheta,vphi))
	exit(-1);
      if(interpolate_seafloor_ages(xloc[HC_THETA], 
				   xloc[HC_PHI],time,v, &age))
	exit(-1);

      fprintf(stdout,"%11g %11g\t%11g %11g %11g\t%11g\n",lon,lat,vphi[0],-vtheta[0],vr[0],age);
    }

  return 0;
}
