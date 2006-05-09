//
//     max order of interpolation - 1
//     ie. if 3, this means will do fourth order on 0,1,2,GGRD_DEF_ORDER arrays
#define GGRD_MAX_ORDER 5
  //     order of derivatives needed in general for the weighting routine
#define GGRD_MAX_IORDER 1

#define GGRD_MAX_ORDERP1 (GGRD_MAX_ORDER+1)
#define GGRD_MAX_IORDERP1 (GGRD_MAX_IORDER+1)

#ifndef GGRD_CPREC			/* 
				   precision for most C functions
				*/
#define GGRD_CPREC double
#define GGRD_EPS 5e-15

#define HC_FLT_FORMAT "%lf"
#endif

/* string length */
#define GGRD_STRLEN 300

/* errors */
#define GGRD_PE(x) {fprintf(stderr,"ggrd: %s\n",x);}

/* radius of CMB */
#define GGRD_RCMB_ND 0.546225

#define GGRD_PI 3.1415926535897932384626433832795

/* 
   
modes

*/
#define GGRD_NORMAL 0
#define GGRD_ONLY_VEL_STATS 1


#include "ggrd_struc.h"

/* filenames */

#define GGRD_VSFILE "vrms.dat"	/* vel stat file */
#define GGRD_THFILE "vtimes.dat" /* times for velocities */
#define GGRD_DFILE "vdepth.dat"	/* depth layer file */
