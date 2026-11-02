/* 

header file for shana approach to spherical harmonics 

*/

struct shana_module{

  HC_CPREC *lfac, *ilfac;
  //
  double  *plm_f1,*plm_f2,*plm_fac1,*plm_fac2,*plm_srt;
  // this is for vector harmonics, only for ivec=1
  HC_CPREC  *sin_theta,*ell_factor;
  // spacing in longitudes
  double dphi,dtheta;
  // int (bounds and such)
  int nlat,nlon,lmsize,lmsize2,nlonm1;
  // logic flags
  hc_boolean initialized,computed_legendre,
    vector_sh_fac_init;
  // init
  hc_boolean was_called;

  int old_ivec,old_lmax,old_npoints,old_tnplm;
};
