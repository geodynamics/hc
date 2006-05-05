#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RICK_PI 3.1415926535897932384626433832795
#define RICK_TWOPI 6.2831853071795864769252867665590

#ifndef rick_boolean
#define rick_boolean unsigned short 
#endif
/* main rick spherical harmonics structure  */
#define TRUE 1
#define FALSE 0

struct rick_module{
  //
  // other stuff needed by more than one subroutine
  // Gauss points: cos(theta), weights, and actual theta
  double *gauss_z,*gauss_w,*gauss_theta;
  double *lfac,ilfac;
  //
  // those are for Legendre polynomials (fac1 & fac2 only for ivec=1)
  // make those double precision
  //
  double *plm_f1,*plm_f2,*plm_fac1,*plm_fac2,*plm_srt;
  // this is for vector harmonics, only for ivec=1
  double *sin_theta,*ell_factor;
  // spacing in longitudes
  double dphi;
  // integer (bounds and such)
  int nlat,nlon,lmsize,lmsize2,nlonm1;
  // logic flags
  rick_boolean initialized,computed_legendre,
    vector_sh_fac_init;
};


#ifndef GENERATE_PROTO
#include "auto_proto.h"
#endif
