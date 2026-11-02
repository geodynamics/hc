 /* 

factor by which to multiply Rick's internal AB coefficients to
physical convention, real spherical haronics coefficients as in Dahlen
and Tromp


*/
#ifndef __READ_RICK_SH_HEADER__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef RICK_TWOPI
#define RICK_TWOPI 6.283185307179586476925286766559005768394
#endif

#ifdef M_PI
#define RICK_PI M_PI
#else
#define RICK_PI 3.1415926535897932384626433832795
#endif

#define SH_RICK_TWO_SQRT_PI 3.5449077018110320545963349666823 /* sqrt(4 pi) */

// 
// how to convert from rick's spherical harmonic convention to dahlen and tromp?
// C_DT(l,m) = SH_RICK_FACTOR(l, m) * C_RICK(l,m)
//
#define SH_RICK_FACTOR(l, m) (pow(-1.0,(SH_RICK_PREC)(m))*SH_RICK_TWO_SQRT_PI)

//


#if SH_RICK_PRECISION == 32

#define SH_RICK_HIGH_PREC long double
#define SH_RICK_PREC long double
#define rick_vecalloc hc_vecalloc
#define rick_vecrealloc hc_vecrealloc
#define SH_RICK_FLT_FMT "%Lf"

#elif SH_RICK_PRECISION == 16	/* double */

#define SH_RICK_HIGH_PREC double
#define SH_RICK_PREC double
#define rick_vecalloc hc_dvecalloc
#define rick_vecrealloc hc_dvecrealloc
#define SH_RICK_FLT_FMT "%lf"

#else  /* single */

#define SH_RICK_HIGH_PREC double
#define SH_RICK_PREC float
#define rick_vecalloc hc_svecalloc
#define rick_vecrealloc hc_svecrealloc
#define SH_RICK_FLT_FMT "%f"

#endif

#ifndef my_boolean
#define my_boolean unsigned short 
#endif
#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE 1
#endif


struct rick_module{
  // read the defines for single/double precision
  // other stuff needed by more than one subroutine
  // Gauss points: cos(theta), weights, and actual theta
  SH_RICK_PREC  *gauss_z, *gauss_w, *gauss_theta;
  //
  //
  SH_RICK_PREC *cfac,*sfac;

  SH_RICK_PREC *lfac, *ilfac;
  // those are for Legendre polynomials (fac1  fac2 only for ivec=1)
  // make those double precision
  //
  SH_RICK_PREC  *plm_f1,*plm_f2,*plm_fac1,*plm_fac2,*plm_srt;
  // this is for vector harmonics, only for ivec=1
  SH_RICK_PREC  *sin_theta,*ell_factor;
  // spacing in longitudes
  SH_RICK_PREC dphi;
  // int (bounds and such)
  int nlat,nlon,lmsize,lmsize2,nlonm1;
  // logic flags
  my_boolean initialized,computed_legendre,
    vector_sh_fac_init,sin_cos_saved;
  // init
  my_boolean was_called;


  /* more book keeping stuff */
  int old_lmax,old_ivec,old_npoints,old_nplm,old_tnplm;

};

/* 


C declarations for F90 functions from Rick's 
spherical harmonics routines


*/

#ifdef LINUX_SUBROUTINE_CONVENTION
#define rick_f90_shd2c rick_f90_shd2c_
#define rick_f90_shd2c_pre rick_f90_shd2c_pre_
#define rick_f90_shc2d rick_f90_shc2d_
#define rick_f90_shc2d_pre rick_f90_shc2d_pre_
#define rick_f90_init rick_f90_init_
#define rick_f90_pix2ang rick_f90_pix2ang_
#define rick_f90_free_module rick_f90_free_module_
#define rick_f90_index rick_f90_index_
#define rick_f90_compute_allplm rick_f90_compute_allplm_
#define rick_f90_ab2cs  rick_f90_ab2cs_
#define rick_f90_cs2ab  rick_f90_cs2ab_
#define rick_f90_realft rick_f90_realft_
#define rick_f90_four1 rick_f90_four1_
#endif

/* 

f90 versions

*/

extern void rick_f90_shd2c(float *,float *,int *,int *, 
			   SH_RICK_PREC *,SH_RICK_PREC *);
extern void rick_f90_shd2c_pre(float *, float *,int *,SH_RICK_PREC *,
			   SH_RICK_PREC *, int *,SH_RICK_PREC *,
			   SH_RICK_PREC *);
extern void rick_f90_shc2d(SH_RICK_PREC *,SH_RICK_PREC *,int *,int *, 
		       float *, float *);
extern void rick_f90_shc2d_pre(SH_RICK_PREC *,SH_RICK_PREC *,
			   int *,SH_RICK_PREC *,
			   SH_RICK_PREC *,int *,float *,float *);
extern void rick_f90_init(int *,int *,int *,int *, int *);
extern void rick_f90_pix2ang(int *, int *, SH_RICK_PREC *, SH_RICK_PREC *);
extern void rick_f90_free_module(int *);
extern void rick_f90_index(int *,int *,int *,int *);
extern void rick_f90_compute_allplm(int *,int *,SH_RICK_PREC *,SH_RICK_PREC *);

/* FFT stuff */
extern void rick_f90_cs2ab( SH_RICK_PREC *, int *);
extern void rick_f90_ab2cs( SH_RICK_PREC *, int *);
extern void rick_f90_realft(SH_RICK_PREC *,int *,int *);
extern void rick_f90_four1(SH_RICK_PREC *,int *,int *);

#define  __READ_RICK_SH_HEADER__
#endif
