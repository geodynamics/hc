/* 

structure for GGRD stuff (scalar and velocity interpolation)


*/
#ifdef USE_GMT4
#include "gmt.h"
#endif
/* 
   
plate tectonic stages interpolation structure

*/
#ifndef GGRD_STRUC_INIT

struct ggrd_t{
  GGRD_CPREC *vtimes;		/* times at which velocities 
				   are specified. this will hold
				   t_left t_mid t_right
				   ....

				*/
  
  int nvtimes,nvtimes3;		/* number of times where specified,
				   and that number times 3 */

  GGRD_CPREC tmin,tmax;		/* range of times */
  unsigned char init;
};

/* 

several GMT grid file structure 

*/
struct ggrd_gt{
  /* 
     grd info 
  */
  struct GRD_HEADER *grd;
  struct GMT_EDGEINFO *edgeinfo;
  /* 
     data 
  */
  float *f,*fmaxlim,bandlim;
  int mm;
  
  float *z;			/* depth levels */
  int nz;
  unsigned char zlevels_are_negative;

  unsigned char init,
    is_three;			/* is it a 3-D set? */
#ifdef USE_GMT4
  struct GMT_BCR bcr;
#endif
};

/*


velocity structure for interpolation

*/
struct ggrd_vel{
  
  GGRD_CPREC *vr,*vt,*vp;	/* velocity field */
  int n[5];		/* dimensions in r, theta, and 
				   phi directions */
  int ntnp,nrntnp;		/*  */
  GGRD_CPREC *rlevels;		/* levels where velocities are 
				   specified */
  GGRD_CPREC dtheta,dphi;		/* spacing in theta and phi */
  GGRD_CPREC velscale,rcmb;
  struct ggrd_t thist;
  unsigned char init,		/* initialized? */
    history,			/* time-dependent? */
    use_age,			/* use an additional age file */
    read_gmt;		/*  read GMT grd files or binary format?*/
  int amode;
  struct ggrd_gt *ages;		/* for seafloor ages */
  int nage;			/* ntime for velo + 1 */
  GGRD_CPREC *age_time;		/* times  */
  float age_bandlim;		/* bandlim for age to decide on
				   continent 
				*/
};
#define GGRD_STRUC_INIT
#endif
