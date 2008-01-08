/* 

structure for GGRD stuff (scalar and velocity interpolation)


*/
#ifndef __GMT_INCLUDED__
#include "gmt.h"
#include "gmt_bcr.h"
#define __GMT_INCLUDED__
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

  GGRD_CPREC xllimit,xrlimit;
  GGRD_CPREC f1_loc,f2_loc,time_old;
  int ileft_loc,iright_loc,ntlim;
  unsigned char called;

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

  double west,east,south,north;

#ifdef USE_GMT4
  struct GMT_BCR loc_bcr[1];
#else
  struct BCR loc_bcr[1];
#endif
};

/* velocity interpolation structure */
struct ggrd_vip{
  int ider[1+3*GGRD_MAX_IORDER],istencil[3],
    ixtracer[3],old_order,orderp1,isshift[3];
  unsigned char init,reduce_r_stencil,z_warned,w_warned;

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
  unsigned char rl_warned,vd_init,vd_reduce_r_stencil;	/*  */
  int amode;
  struct ggrd_gt *ages;		/* for seafloor ages */
  int nage;			/* ntime for velo + 1 */
  GGRD_CPREC *age_time;		/* times  */
  float age_bandlim;		/* bandlim for age to decide on
				   continent 
				*/
  struct ggrd_vip vd;		/* velocity interpolation structure */
  /* seafloor stuff */
  unsigned short sf_init;
  GGRD_CPREC  sf_old_age,sf_old_f1,sf_old_f2;
  int sf_old_left,sf_old_right,sf_ntlim;
};
#define GGRD_STRUC_INIT
#endif
