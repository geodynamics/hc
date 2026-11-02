#ifndef  __SH_HEADER_READ__

#ifdef HC_USE_HEALPIX
/* 
   load definitions for complex variable structures, Healpix structure
   and declarations of Healpix related functions
*/
#include "myhealpix.h"
#endif
/* 
   
same for Rick's spherical harmonic routines

*/
#define SH_RICK_PRECISION  HC_PRECISION

#include "sh_rick.h"
/* 

shana

 */
#include "sh_shana.h"

/* 
   spherical harmonics types 
*/
#define SH_RICK 0		/* Rick's Gauss quadrature/FFT routines */
#define SH_HEALPIX 1		/* Healpix package */
#define SH_SHANA 2		/* grid based spherical harmonics, similar to SH_RICK 
				   but no FFT/Gauss points
				*/


/* 

 my spherical harmonics structures


 
*/
/* 

   individual expansion
   
*/
struct sh_lms{			
  int type;			/* this is the type of expansion 

				HEALPIX: ab will be in imaginary convention

				*/
  /* have the spatial or the spherical versions been initialized? */
  my_boolean spectral_init;
  /* l_max  */
  int lmax;
  /*  bounds from above + 1 */
  int lmaxp1;
  /*  
      two ways of storing the (l,m) arrays 
  */
  int lmbig; 			/* (lmax+1)**2 */
  int lmsmall2; 		/* (exp->lmax+1)*(exp->lmax+2)for A and B */
  /* 
     number of A,B entries, and total number of entries for a
     spherical harmonics expansion set (depends on ivec)
  */
  int n_lm;
  /* 
     number of entries for the Legendre function array 
  */
  int n_plm,tn_plm,tn_plm_irr;
  /* 
     number of points in each layer the spatial domain 
  */
  int npoints;
  /* 

  */
  hc_boolean plm_computed,plm_computed_irr;
  int  old_lmax,old_ivec,old_tnplm,old_tnplm_irr,
    old_lmax_irr,old_ivec_irr;

  /*

    holds the coefficients:
    
  */
#ifdef HC_USE_HEALPIX
  /* 
     for Healpix 
  */
  struct scmplx *alm_c;		/* single prec complex */
  /* 
     for HEALPIX
  */
  struct healpix_parameters heal;
#endif
  /* 
     for Rick type 
  */
  SH_RICK_PREC *alm;
  struct rick_module rick;
};
/* 

spherical harmonics model structure, this holds several spherical
harmonic expansions

*/
struct sh_lms_model{
  /* 
     number of sets 
  */
  int nset;
  /* 
     Legendre polynomial flags
  */
  my_boolean save_plm;
  /* layer indicators if nset != 1 */
  HC_PREC *z;
  /* scalar only? ivec=0 or velocities? ivec=1 */
  int ivec;
  /* number of expansions per set */
  int shps;
  /* expansions  */
  struct sh_lms *exp;
  int nexp;			/* number of expansions */
  SH_RICK_PREC *plm;			/* precomputed Legendre 
					   functions */
  
  /* 
     spatial data points
  */
  int tnpoints;			/* number of the total datapoints */
  /* data structure */
  HC_PREC *data;
  my_boolean spatial_init, initialized; 
};
/* 

 compute the (l,m) index of a tighly packed array of lmsize2
 size. arrays are C style (0...lmsize2-1).  

 lmsize2 = (lmax+1)*(lmax+2), which holds all A and B coefficients
 of an expansion of maximum order lmax
 
 this assumes that A and B coefficients are stored next to each
 other. pass a_or_b as 0 or 1 for A or B coefficients, respectively
 
 0 <= l <= lmax, 0 <= m <= l, 0 <= a_or_b <= 1
 

*/
#define LM_INDEX(l,m,a_or_b) ((((l)+1)*(l)/2+(m))*2+(a_or_b))

#define __SH_HEADER_READ__
#endif
