/* 

structure for GGRD stuff (scalar and velocity interpolation)


*/
#ifndef __GGRD_STRUC_INIT__


#ifndef __GMT_INCLUDED__
#include "gmt.h"
#define __GMT_INCLUDED__
#endif


#include "prem.h"

/* 
   
plate tectonic stages interpolation structure

*/
struct ggrd_t{
  char file[1000];		/* filename */
  GGRD_CPREC *vtimes;		/* times at which velocities or materials 
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

  GGRD_CPREC vstage_transition;

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



struct ggrd_temp_init{
  /* 
     for temperature init from grd files 
  */
  int init,scale_with_prem;
  int override_tbc,limit_trange;
  double scale,offset;
  char gfile[1000];
  char dfile[1000];
  struct ggrd_gt d[1];		/* grid structure */
  struct prem_model prem; 	/* PREM model */
};

struct ggrd_master{		/* master structure */

  int mat_control,mat_control_init;
  int vel_control,vel_control_init;
  
  char mat_file[1000];
  char vel_file[1000];
  
  /* grid structures */
  struct ggrd_gt *mat;		/* material grids */

  /* different for velocities */
  struct ggrd_vel *ggrd_v;	/* velocity grids */
  struct ggrd_t time_hist;	/* time history structure */

  /* temperature init */
  struct ggrd_temp_init temp_init;
};



#define __GGRD_STRUC_INIT__
#endif
