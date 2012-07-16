#include "hc.h"

int main(int argc, char **argv)
{
  struct ggrd_master *ggrd;
  double time,vr[4],vphi[4],vtheta[4],dtrange;
  double xloc[3];
  static int order = 3;
  hc_boolean calc_derivatives ;
  double lon,lat,age;
  int mode;

  mode = 1;


  switch(mode){
  case 1:
    if(argc>1)
      ggrd_grdinfo(argv[1]);
    break;
  case 2:
    /* 
       initialize velocity structure
    */
    ggrd = (struct ggrd_master *)calloc(1,sizeof(struct ggrd_master));
    ggrd_init_master(ggrd);
    ggrd->v.history = TRUE;		/* expect history */
    ggrd->age_control = TRUE;		/* expect seafloor age files */
    ggrd->age_bandlim = 900;
    /* 
       read in velocities 
    */
    if(ggrd_read_vel_grids(ggrd,1.0,FALSE,TRUE,
			   "/home/walter/becker/data/plates/past/clb/hall/",FALSE)){
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
	if(ggrd_find_vel_and_der(xloc,time,dtrange,ggrd,
				 order,calc_derivatives,
				 TRUE,vr,vtheta,vphi))
	  exit(-1);
	if(interpolate_seafloor_ages(xloc[HC_THETA], 
				     xloc[HC_PHI],time,ggrd, &age))
	  exit(-1);
	
	fprintf(stdout,"%11g %11g\t%11g %11g %11g\t%11g\n",lon,lat,vphi[0],-vtheta[0],vr[0],age);
      }
    break;			/* end mode 2 */
  }
  return 0;
}
