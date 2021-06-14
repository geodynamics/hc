/* 
   
   header file for Hager & O'Connell experimental code

   Thorsten (twb@usc.edu)

   $Id: hc.h,v 1.12 2006/05/01 17:46:18 becker Exp becker $

*/
#ifndef __HC_READ_HEADER_FILE__
#include "hc_filenames.h"
#include <math.h>



/* 

general variable type defines 

*/
#ifndef hc_boolean
#define hc_boolean unsigned short
#endif

#ifndef HC_PRECISION			/* not defined */
/* 
   default precision for most C functions
*/
#define HC_PRECISION 16		/* 8,16, or 32 denoting float, double,
				   or quad only for the purposes of
				   switching */
#endif

#ifndef __HC_PREC_DEFINED__

#if HC_PRECISION == 32	/* quad prec */

#define HC_PREC long double
#define HC_HIGH_PREC long double
#define HC_FLT_FORMAT "%Lf"
#define HC_TWO_FLT_FORMAT "%Lf %Lf"
#define HC_THREE_FLT_FORMAT "%Lf %Lf %Lf"
#define HC_EPS_PREC 3e-25

#elif HC_PRECISION == 16		/* double prec */

#define HC_PREC double
#define HC_HIGH_PREC double	/* some computations at higher prec */
#define HC_FLT_FORMAT "%lf"
#define HC_TWO_FLT_FORMAT "%lf %lf"
#define HC_THREE_FLT_FORMAT "%lf %lf %lf"
#define HC_EPS_PREC 3e-16

#else  /* single prec */

#define HC_PREC float
#define HC_HIGH_PREC double
#define HC_FLT_FORMAT "%f"
#define HC_TWO_FLT_FORMAT "%f %f"
#define HC_THREE_FLT_FORMAT "%f %f %f"
#define HC_EPS_PREC 3e-7

#endif

#define HC_BIN_PREC float	/* precision for binary I/O */

#ifndef HC_CPREC
#define HC_CPREC HC_PREC
#endif
#define __HC_PREC_DEFINED__

#endif

#ifndef GGRD_CHAR_LENGTH
#define HC_CHAR_LENGTH 1000		/* length of char arrays */
#else
#define HC_CHAR_LENGTH GGRD_CHAR_LENGTH
#endif



#ifndef __HC_DEF_COMPLEX__
/* complex variables */
struct hc_scmplx{			/* single precision */
  float dr,di;
};
struct hc_dcmplx{			/* double precision */
  double dr,di;
};
#define __HC_DEF_COMPLEX__
#endif

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

#define HC_PI M_PI

#define HC_TWOPI 6.2831853071795864769252867665590


#define HC_SCALAR 0
#define HC_VECTOR 1

/* 

PREM

*/
#include "prem.h"

/* 

Other Earth model, and modeling constants

*/
#include "hc_constants.h"


/*

spherical harmonics 
   
*/
#include "sh.h"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
/* do we need this for 64 bit? OS-X didn't like it */
//#include <malloc.h>

#include <limits.h>


/* 

use short format for spherical harmonic models

*/
#define HC_DEFAULT_INTERNAL_FORMAT SH_RICK /* by default, use Rick's
					      spherical harmonic
					      routines internally  */


/* 

for H & C solutions

*/
struct hc_sm{
  HC_PREC u[6][4];
};
/* 
   for polsol and solve
*/
struct hc_ps{
  /* scaling factors will only be computed once */
  int ncalled;
  HC_HIGH_PREC rho_scale;
  HC_HIGH_PREC alpha, beta, geoid_factor;
  hc_boolean rho_init, 
    prop_params_init, 	/* parameters for propagator computation */
    abg_init  ,		/* alpha, beta factors */
    prop_mats_init;	/* will be true only if save_prop_mats is 
			   requested  */
  int solver_kludge_l;		/* for treating higher orders
				   differently  */

  /* for solve */
  hc_boolean tor_init, pol_init;
};
/* 


parameter structure to allow for settings that are specific to the
implementation and higher level than the calls to the subroutines


 */
struct hc_parameters{
  hc_boolean compressible;	/* compressibility following Panasyuk
				   & Steinberger */
  hc_boolean free_slip;		/* surface mechanical boundary condition */
  hc_boolean no_slip;		/* 
				   zero surface velocities
				*/
  hc_boolean platebc;		/* plate velocity BC */
  HC_PREC dens_anom_scale;	/* default density anomaly scaling to
				   go from PREM percent traveltime
				   anomalies to density anomalies */
  hc_boolean scale_dens_anom_with_prem; /* scale input anomalies with
					   PREM, or average mantle
					   density */
  hc_boolean verbose;		/* debugging output? (0,1,2,3,4...) */
  hc_boolean sol_binary_out;	/* binary or ASCII output of SH expansion */
  hc_boolean print_spatial;	/* print the spatial solution */
  hc_boolean compute_geoid; 	/* compute and print the geoid */
  hc_boolean print_density_field;	 /* print the scaled density field */
  int solution_mode;	/* velocity or stress */

  int pvel_mode;		/* plate velocity mode */
  HC_PREC pvel_time;		/* time to use */
  hc_boolean remove_nr;		/* remove net rotation from plate models? */
  int solver_kludge_l;		/* for CMB BC tricks */

  hc_boolean solver_mode;	
  hc_boolean visc_init_mode;

  hc_boolean vscan_rlv; /* restrict the v[4] values to
			   something "reasonable" in
			   terms of relative layer
			   strength */

  HC_PREC vscan_dv;		/* spacing for viscosity scan */
  int vscan_n;		/* how many layers to use */

  HC_PREC elayer[HC_VSCAN_NLAYER_MAX]; /* use full four layers for array */

  hc_boolean read_short_dens_sh, read_short_pvel_sh; /* short SH format for density or plate velocity files
							files? */
  int dd_dens_scale; /* 0: constant scaling
			1: read the density/velocity scaling from file 
			2: polynomial applies
		     */
  HC_PREC *rdf,*sdf;
  int ndf;
  struct sh_lms *ref_geoid,*ref_dtopo;
  
  HC_PREC rlayer[HC_VSCAN_NLAYER_MAX-1];		/* for four
							   layer
							   approaches
							   (first
							   layer
							   radius is
							   CMB
							   radius) */

  hc_boolean print_pt_sol;	/* output of p[6] and t[2] vectors */
  hc_boolean print_kernel_only;

  
  char visc_filename[HC_CHAR_LENGTH];	/* name of viscosity profile file */
  char pvel_filename[HC_CHAR_LENGTH];	/* name of plate velocities file */
  char dens_filename[HC_CHAR_LENGTH];	/* name of density model file */
  char prem_model_filename[HC_CHAR_LENGTH];	/* PREM model filename */
  char dens_scaling_filename[HC_CHAR_LENGTH];	/*  */
  char ref_geoid_file[HC_CHAR_LENGTH]; /* reference geoid */
  char ref_dtopo_file[HC_CHAR_LENGTH]; /* reference dynamic topography */
};

/* plate velocity structure */
struct pvels{
  int n;
  HC_PREC *t;
  struct  sh_lms *p;
};
/* 


general HAGER & O'CONNELL flow code structure


*/
struct hcs{

  int lmax;			/* max order of fields */
  int sh_type;			/* type of spherical harmonics to use, e.g.
				   HEALPIX or RICK

				*/
  int nrad;			/* number of output radial 
				   layers 

				   this is determined by the density
				   model

				*/
  HC_PREC *r;			/* non-dim radii */
  HC_PREC *dvisc;		/* viscosity at the density layers */
  /* 
     viscosity structure
  */
  int nvis;		/* number of viscosity layers */
  HC_PREC *visc;	/* values of normalized viscosities */
  HC_PREC *rvisc; 	/* radii of normalized viscosities */

  /* 
     density contribution
  */
  int inho;		/* number of density layers */
  HC_PREC *dfact; 	/* density factors, derived from layer thickness */
  HC_PREC *rden; 	/* radii of density layers [normalized] */
  HC_PREC rho_top_kg;	/* top density in kg/m^3 */
  struct sh_lms *dens_anom; /* 
			       expansions of density
			       anomalies has to be [inho] (if
			       those change, that's OK, no
			       need to call with
			       dens_fac_changed == TRUE)
			    */
  struct sh_lms *dens_anom_orig;
  HC_PREC dens_scale;		/* scale for density file */

 

  hc_boolean compressible; /* 
			      if TRUE, will use PREM
			      densities, else, average
			      mantle density
			   */
  int npb;			/* number of phase boundaries */
  HC_PREC *rpb;	/* radius and F factor for phase
		   boundaries */
  HC_PREC *fpb;
  my_boolean free_slip;  /* 
			    include plate velocities?
			    possible, if free_slip is
			    FALSE
			 */

  struct pvels pvel; /* 
			   
			   poloidal and toroidal part of plate motions
			   
			   (typically only one set, meaning two expansions)
			   
			*/

  struct hc_ps psp;
  hc_boolean save_solution; /* memory intensive speedup in poloidal
				solution by saving propagator matrices 
				this will also keep the toroidal solution 
				for speedup in case the velocities have not 
				changed between calls to hc_solve
			     */

  /* 
     SOLUTION RELATED 
  */
  /* poloidal solution */
  struct sh_lms *pol_sol;
  HC_PREC *rho,*rho_zero;	/* 
				   density factors 
				*/
  HC_HIGH_PREC *rprops,*pvisc,*props,*ppots, /* propagator related */
    *den;
  /* 
     propagator related factors as well 
  */
  int nprops,nradp2,nvisp1,inho2;
  hc_boolean *qwrite;
  struct sh_lms *tor_sol; /* 
			  toroidal solution
		       */
  hc_boolean initialized,const_init,visc_init,dens_init,pvel_init;	/* logic flags */
  hc_boolean orig_danom_saved;
  /* sqrt(l(l+1)) and 1/lfac factors */
  HC_PREC *lfac,*ilfac;
  int lfac_init;
  /* PREM related */
  struct prem_model prem[1];	/* PREM model constants */
  hc_boolean prem_init; 
  /* 
     constants
  */
  HC_PREC timesc;			/* timescale in [yr] */
  HC_PREC visnor;			/* visocsity for nomalization [Pas] */
  HC_PREC g;			/* gravitational constant */
  HC_PREC gacc;			/* gravitational acceleration  [m/s^2]*/
  HC_PREC re;			/* radius of Earth [m] */
  HC_PREC avg_den_core, avg_den_mantle; /* average densities of 
					 core and mantle */
  HC_PREC secyr;		/* seconds per year */
  HC_PREC vel_scale;		/* 
				   velocity scale to change from input
				   [cm/yr] to internal
				*/
  HC_PREC stress_scale;		/* to go to MPa */
  HC_PREC r_cmb;		/* radius of CMB */
  /* Legendre functions */
  SH_RICK_PREC *plm;

  /* more logic flags */
  my_boolean spectral_solution_computed, spatial_solution_computed;
};


/* 

solution modes

*/
#define HC_VEL 0		/* compute velocities 
				   v_r v_t v_p
				*/
#define HC_RTRACTIONS 1		/* compute tractions 
				   trr trt trp  */
#define HC_HTRACTIONS 2		/* compute tractions ttt, ttp, tpp */

/* 

init and assignment modes

*/
#define HC_INIT_E_FROM_FILE 0	/* viscosity */
#define HC_INIT_E_FOUR_LAYERS 1

#define HC_INIT_D_FROM_FILE 0	/* density */
#define HC_RESCALE_D 1

/* different depth dependent initialization */
#define HC_INIT_DD_FROM_FILE 0	/* depth dependent density scaling */
#define HC_INIT_DD_FOUR_LAYERS 1	/*  */

/* 

density scaling modes

*/

#define HC_DD_CONSTANT 0
#define HC_DD_READ_FROM_FILE 1
#define HC_DD_POLYNOMIAL 2


/* plate initialization mode */
#define HC_INIT_P_FROM_FILE 0	/* single plate expansion */
#define HC_INIT_P_FROM_DIRECTORY 1 /* directory with a range of
				      expansions at different times */

#define HC_SOLVER_MODE_DEFAULT 0
#define HC_SOLVER_MODE_VISC_SCAN 1
#define HC_SOLVER_MODE_DYNTOPO_INVERT 2

/* 
   
declarations

*/
#ifndef GENERATE_PROTO
/* local declarations */
#include "hc_auto_proto.h"
#endif
/* 
   
macro defintions

*/
/* ordering for spherical coordinates */
#define HC_R 0
#define HC_THETA 1
#define HC_PHI 2
#define HC_TPPROD 3
#define HC_NRNTNP 4

#define HC_X 0
#define HC_Y 1
#define HC_Z 2

/* for solution */
#define HC_RAD 0
#define HC_POL 1
#define HC_TOR 2
/* 
   boolean 
*/
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif


/* in-line functions */
#define HC_MEMERROR(x) {fprintf(stderr,"%s: memory allocation error, exiting\n",x);exit(-1);}
#define HC_ERROR(x,y) {fprintf(stderr,"%s: error: %s, exiting\n",x,y);exit(-1);}

#define HC_MIN(x,y) (( (x) < (y)) ? (x) : (y))

#define HC_DIFFERENT(x,y) ((fabs((x)-(y)) > 1e-7)?(1):(0))

/* 
   trigonometry stuff 
*/
#ifndef SQRT_TWO 
#ifdef M_SQRT2
#define SQRT_TWO M_SQRT2
#else
#define SQRT_TWO  1.41421356237309504880168872420970 
#endif
#endif
#ifndef PIOVERONEEIGHTY
#define PIOVERONEEIGHTY  0.017453292519943295769236907684886 
#endif
#ifndef ONEEIGHTYOVERPI
#define ONEEIGHTYOVERPI  57.295779513082320876798154814105
#endif
#ifndef DEG2RAD
#define DEG2RAD(x) ((x) * PIOVERONEEIGHTY)
#endif
#ifndef RAD2DEG
#define RAD2DEG(x) ((x) * ONEEIGHTYOVERPI)
#endif

#ifndef HC_FLT_MAX 
#define HC_FLT_MAX 1e20
#endif
#ifndef HC_FLT_MIN
#define HC_FLT_MIN -1e20
#endif 

#ifndef THETA2LAT
#define THETA2LAT(x) ( (90.0 - RAD2DEG(x)) )
#endif

#ifndef LAT2THETA
#define LAT2THETA(x) ( (90.0 - (x))*PIOVERONEEIGHTY )
#endif

#ifndef PHI2LON
#define PHI2LON(x) ( RAD2DEG(x) )
#endif

#ifndef LON2PHI
#define LON2PHI(x) ( DEG2RAD(x) )
#endif
/* 


*/
/* 
   convert depth (>0) in km to non-dimensionalizd radius
*/
#define HC_ND_RADIUS(x) (1.0-((x)/HC_RE_KM))
/* 
   the other way around
*/
#define HC_Z_DEPTH(x) ((HC_RE_KM * (1.0-(x))))


#define __HC_READ_HEADER_FILE__
#endif
