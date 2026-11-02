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
  char file[GGRD_CHAR_LENGTH];		/* filename */
  GGRD_CPREC *vtimes;		/* times at which velocities or materials 
				   are specified. this will hold
				   
				   t_left t_mid t_right

				   ....


				*/
  
  int nvtimes,nvtimes3;		/* number of times where specified,
				   and that number times 3 */

  GGRD_CPREC tmin,tmax;		/* range of times */
  ggrd_boolean init;

  GGRD_CPREC xllimit,xrlimit;
  GGRD_CPREC f1_loc,f2_loc,time_old;
  int ileft_loc,iright_loc,ntlim;

  GGRD_CPREC vstage_transition;

  ggrd_boolean interpol_time_lin;
  
  GGRD_CPREC *tl;		/* for linear interpolation */
  

  ggrd_boolean called;

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

  ggrd_boolean zlevels_are_negative;

  ggrd_boolean init,geographic_in,
    is_three;			/* is it a 3-D set? */

  double west,east,south,north;

#ifndef USE_GMT3
  struct GMT_BCR loc_bcr[1];
#else
  struct BCR loc_bcr[1];
#endif
};

/* velocity interpolation structure */
struct ggrd_vip{
  int ider[1+3*GGRD_MAX_IORDER],istencil[3],
    ixtracer[3],old_order,orderp1,isshift[3];
  ggrd_boolean init,reduce_r_stencil,z_warned,w_warned;

};
/*


structure for 3-D velocity interpolation

*/
struct ggrd_vel{
  GGRD_CPREC *vr,*vt,*vp;	/* velocity field */
  int n[5];		/* dimensions in r, theta, and 
				   phi directions */
  int ntnp,nrntnp;		/*  */
  GGRD_CPREC *rlevels;		/* levels where velocities are 
				   specified */
  GGRD_CPREC dtheta,dphi;	/* spacing in theta and phi */
  GGRD_CPREC velscale,rcmb;
  ggrd_boolean init,		/* initialized? */
    history,			/* time-dependent? */
    read_gmt;		/*  read GMT grd files or binary format?*/
  ggrd_boolean rl_warned,vd_init,vd_reduce_r_stencil;	/*  */

  struct ggrd_vip vd;		/* velocity interpolation structure */
};



struct ggrd_temp_init{
  /* 
     for temperature init from grd files 
  */
  int init,scale_with_prem;
  int override_tbc,limit_trange;
  double scale,offset;
  char gfile[GGRD_CHAR_LENGTH];
  char dfile[GGRD_CHAR_LENGTH];
  struct ggrd_gt d[1];		/* grid structure */
  struct prem_model prem; 	/* PREM model */
};

struct ggrd_master{		/* master structure */
  /* citcom use flags */
  int mat_control,mat_control_init;
  int ray_control,ray_control_init;
  int vtop_control,vtop_control_init;
  int age_control,age_control_init;
  
  char mat_file[GGRD_CHAR_LENGTH];
  char ray_file[GGRD_CHAR_LENGTH];
  char vtop_dir[GGRD_CHAR_LENGTH];
  char age_dir[GGRD_CHAR_LENGTH];
  
  /* grid structures */
  struct ggrd_gt *mat;		/* material grids */
  /* rayleigh number grids */
  struct ggrd_gt *ray;

  /* surface velocity */
  struct ggrd_gt *svp,*svt;	/* phi/theta surface velocities */
  /* age stuff */
  struct ggrd_gt *ages;
  int nage,amode;
  GGRD_CPREC *age_time;		/* times  */
  float age_bandlim;		/* bandlim for age to decide on
				   continent  */

  unsigned short sf_init;	/* seafloor stuff */
  GGRD_CPREC  sf_old_age,sf_old_f1,sf_old_f2;
  int sf_old_left,sf_old_right,sf_ntlim;

  struct ggrd_vel v;	/* 3D velocity grid structure */

  /* time history */
  struct ggrd_t time_hist;	/* time history structure */
  /* temperature init */
  struct ggrd_temp_init temp;int use_temp;
  /* composition init  */
  struct ggrd_temp_init comp;int use_comp;
};



#define __GGRD_STRUC_INIT__
#endif
