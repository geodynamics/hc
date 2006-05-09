/* 
   
   header file for Hager & O'Connell experimental code

   Thorsten (twb@usc.edu)

   $Id: hc.h,v 1.12 2006/05/01 17:46:18 becker Exp becker $

*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "hc_filenames.h"
/* 

general variable type defines 

*/
#ifndef hc_boolean
#define hc_boolean unsigned short
#endif
#ifndef HC_PREC			/* 
					   precision for most C functions
					*/
#define HC_PREC double
#define HC_FLT_FORMAT "%lf"
#define HC_TWO_FLT_FORMAT "%lf %lf"
#define HC_EPS_PREC 5e-15
#endif

#ifndef HC_CPREC
#define HC_CPREC HC_PREC
#endif

#define HC_CHAR_LENGTH 300		/* length of char arrays */

#define HC_BIN_PREC float	/* precision for binary I/O */

#ifndef DEFINED_MY_COMPLEX
/* complex variables */
struct scmplx{			/* single precision */
  float dr,di;
};
struct dcmplx{			/* double precision */
  double dr,di;
};
#define DEFINED_MY_COMPLEX
#endif



/* 

PREM

*/
#include "prem.h"
/* 
   dealing with velocity grids 

*/
#include "ggrd.h"
#include "ggrd_grdtrack_util.h"
/*

spherical haronics 
   
*/
#include "sh.h"

/* 

for H & C solutions

*/
struct hc_sm{
  HC_PREC u[6][4];
};

struct hc_parameters{
  hc_boolean compressible;	/* compressibility following Panasyuk
				   & Steinberger */
  hc_boolean free_slip;		/* surface mechanical boundary condition */
  hc_boolean vel_bc_zero;		/* 
				   if false, plate velocities, else no
				   slip if free_slip is false as well
				*/
  HC_PREC dens_anom_scale;	/* default density anomaly scaling to
				   go from PREM percent traveltime
				   anomalies to density anomalies */
  hc_boolean verbose;		/* debugging output? (0,1,2,3,4...) */
  hc_boolean sol_binary_out;	/* binary or ASCII output of SH expansion */
  hc_boolean print_spatial;	/* print the spatial solution */
  int solution_mode;	/* velocity, stress, or geoid */
  hc_boolean print_pt_sol;	/* output of p[6] and t[2] vectors */
  char visc_filename[HC_CHAR_LENGTH];	/* name of viscosity profile file */
  char pvel_filename[HC_CHAR_LENGTH];	/* name of plate velocities file */
  char dens_filename[HC_CHAR_LENGTH];	/* name of density model file */
  char prem_model_filename[HC_CHAR_LENGTH];	/* PREM model filename */
};

/* 


general HAGER & O'CONNELL flow code structure


*/
struct hcs{

  int lmax;			/* max order of fields */
  int sh_type;			/* type of spherical harmonics to use, e.g.
				   HEALPIX or RICK

				*/
  int nrad;			/* number of output radial layers */
  HC_PREC *r;			/* non-dim radii */
  /* 
     viscosity structure
  */
  int nvis;		/* number of viscosity
			   layers */
  HC_PREC *visc;	/* values of normalized viscosities */
  HC_PREC *rvisc; 	/* radii of normalized viscosities */

  /* 
     density contribution
  */
  int inho;		/* number of density layers */
  HC_PREC *dfact; 	/* density factors, derived from layer thickness */
  HC_PREC *rden; 	/* radii of density layers [normalized] */
  
  struct sh_lms *dens_anom; /* 
			    expansions of density
			    anomalies has to be [inho] (if
			    those change, that's OK, no
			    need to call with
			    dens_fac_changed == TRUE)
			 */
  HC_PREC dens_scale[1];		/* scale for density file */

 

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

  struct sh_lms pvel[2]; /* 
			 
		      poloidal and toroidal part of plate motions
		      (only one expansion)
		      
		      */


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
  struct sh_lms *pol_sol,*geoid;
  double *rho,*rho_zero;	/* 
				   density factors 
				*/
  double *rprops,*pvisc,*props,*ppots, /* propagator related */
    *den;
  /* 
     propagator related factors as well 
  */
  int nprops,nradp2,nvisp1,inho2;
  hc_boolean *qwrite;
  struct sh_lms *tor_sol; /* 
			  toroidal solution
		       */
  hc_boolean initialized;	/* logic flag */
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
  HC_PREC r_cmb;		/* radius of CMB */
  /* Legendre functions */
  double *plm;

  /* more logic flags */
  my_boolean spectral_solution_computed, spatial_solution_computed;
};


/* 

solution modes

*/
#define HC_GEOID 1		/* compute geoid */
#define HC_VEL 0		/* compute velocities 
				   v_r v_t v_p
				*/
#define HC_STRESS -1		/* compute tractions 
				   trr trt trp  */
/* 

init and assignment modes

*/
#define HC_INIT_FROM_FILE 0


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
#define SQRT_TWO 1.41421356237309504880168872420970 
#endif
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
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

#ifndef FLT_MAX 
#define FLT_MAX 1e20
#endif
#ifndef FLT_MIN
#define FLT_MIN -1e20
#endif 

#define THETA2LAT(x) ( (90.0 - (x)*ONEEIGHTYOVERPI) )
#define LAT2THETA(x) ( (90.0 - (x))*PIOVERONEEIGHTY )
#define PHI2LON(x) ( RAD2DEG(x) )
#define LON2PHI(x) ( DEG2RAD(x) )
/* 

other constants

*/
// now taken from earth model
#define HC_RE_KM 6371.0		/* radius(Earth) in [km] */

//#define HC_RCMB_ND 0.546225    /* non-dim radius, ~10km above
//				  CMB  */

#define HC_TIMESCALE_YR 1e6	/* timescale [yr] */
/* 
   convert depth (>0) in km to non-dimensionalizd radius
*/
#define HC_ND_RADIUS(x) (1.0-((x)/HC_RE_KM))
/* 
   the other way around
*/
#define HC_Z_DEPTH(x) ((HC_RE_KM * (1.0-(x))))

