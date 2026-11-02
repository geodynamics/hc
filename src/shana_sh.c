#include "hc.h"

/* 

//
// compute Legendre function (l,m) evaluated on nlat points 
// in latitude and their derviatives with respect to theta, if 
// ivec is set to 1
// 

*/

void shana_compute_allplm(int lmax,int ivec,double *plm,
			  double *dplm, struct shana_module *shana) 
{
  int i,os;
  os = 0;
  for (i=0;i < shana->nlat;i++) { 
    shana_plmbar1((plm),(dplm),ivec,lmax,shana->gauss_z[i],shana);
    os += shana->lmsize;
  }
}


/* //
// detemine the colatidude and longitude of PIxel index
// where index goes from 0 ... (nlon * nlat)-1
// */
void 
shana_pix2ang (index, lmax, theta, phi, shana)
int index;
int lmax;
double *theta;
double *phi;
struct shana_module *shana; {
  int  i,j;
  if(!shana->initialized){
    fprintf(stderr,"shana_pix2ang: error: shana module not initialized\n");
    exit(-1);
  }
  
  j = index;
  i=0;
  while(j > shana->nlonm1) { /* */
    j -= shana->nlon;
    i++;
  }
  *theta = shana=>dtheta * (HC_CPREC)i;
  *phi =   shana->dphi   * (HC_CPREC)j;
}


void 
shana_shc2d (cslm, dslm, lmax, ivec, rdatax, rdatay, shana)
HC_CPREC *cslm;
HC_CPREC *dslm;
int lmax;
int ivec;
HC_CPREC *rdatax;
HC_CPREC *rdatay;
struct shana_module *shana;
{
  //
  // Transforms spherical harmonic coefficients of a scalar (ivec = 0)
  // or a vector (ivec=1) field
  // to data points on a grid. Reverse transform of subroutine
  // shd2c.f so long as the same degree is used and the points
  // in latitude are Gaussian integration points. Maximum degree
  // must be 2**n-1 in order for the FFT to work; this determines
  // the grid spacing in longitude. nlat = lmax+1, nlon = 2*nlat
  //
  // INPUT:
  //
  //		lmax	spherical harmonic degree used in expansion
  //		cslm	(lmsize2) spherical harmonic coefficients
  //             dslm    (lmsize2)
  //                     if ivec, will hold poloidal and toroidal coeff,else
  //                     dslm will not be referenced
  //
  //             ivec   0: scalar 1: vector field
  // OUTPUT:
  //		rdatax	(nlon * nlat) values on grid points
  //             rdatay  (nlon * nlat). x and y will be theta and phi for 
  //                     ivec = 1, else rdatay will not be referenced
  //
  // if Ivec is set, assume velocities to be expanded instead
  //
  // input: lmax,ivec
  // local
  double *plm,*dplm;
  if(!shana->initialized){
    fprintf(stderr,"shana_shc2d: error: initialize first\n");
    exit(-1);
  }
  // compute the Plm first
  if(lmax != shana->nlat-1){
    fprintf(stderr,"shana_shc2d: error: lmax mismatch: nlat: %i lmax: %i\n",
	    shana->nlat,lmax);
    exit(-1);
  }
  /* allocate memory */
  hc_dvecalloc(&plm,shana->nlat*shana->lmsize,"shana_shc2d: mem 1");
  if(ivec)
    hc_dvecalloc(&dplm,shana->nlat*shana->lmsize,"shana_shc2d: mem 2");
  //
  // compute the Plm first
  shana_compute_allplm(lmax,ivec,plm,dplm,shana);

  //
  // call the precomputed subroutine
  shana_shc2d_pre(cslm,dslm,lmax,plm,dplm,ivec,rdatax,rdatay,shana);

  /* free legendre functions if not needed */
  free(plm);
  if(ivec)
    free(dplm);
}



/* //
// the actual routine to go from spectral to spatial: added structure shana
// */
void 
shana_shc2d_pre (cslm, dslm, lmax, plm, dplm, ivec, rdatax, rdatay, shana)
HC_CPREC *cslm;
HC_CPREC *dslm;
int lmax;
double *plm;
double *dplm;
int ivec;
float *rdatax;
float *rdatay;
struct shana_module *shana;
{
  /* //
  // Legendre functions are precomputed
  // */
  HC_CPREC  *valuex, *valuey;
  double  dpdt,dpdp;
  static int negunity = -1;
  int  i,j,m,m2,j2,ios1,l,lmaxp1,lmaxp1t2,oplm,nlon2,lm1;
  if(!shana->initialized){
    fprintf(stderr,"shana_shc2d_pre: error: initialize modules first\n");
    exit(-1);
  }
  // check bounds 
  lmaxp1 = lmax + 1;                // this is nlat
  lmaxp1t2 = 2 * lmaxp1;               // this is nlon
  if((shana->nlat != lmaxp1)||(shana->nlon != lmaxp1t2)){
    fprintf(stderr,"shana_shc2d_pre: dimension mismatch: lmax: %i nlon: %i nlat:%i\n",
	    lmax,shana->nlon,shana->nlat);
    exit(-1);
  }
  if(ivec){
    if(!shana->vector_sh_fac_init){
      fprintf(stderr,"shana_shc2d_pre: error: vector harmonics factors not initialized\n");
      exit(-1);
    }
  }
  nlon2 = shana->nlon + 2;
  /* allocate value arrays */
  shana_vecalloc(&valuex,nlon2,"shana_shc2d_pre 1");
  if(ivec)
    shana_vecalloc(&valuey,nlon2,"shana_shc2d_pre 2");
  
  for (i=0;i < shana->nlat; i++) {	
    /*  
	loop through latitudes 
    */
    oplm = i * shana->lmsize;        /*  offset for Plm array. (changed indices TWB) */
    ios1 = i * shana->nlon;        /* offset for data array */
    if(!ivec){          
      /* 
	 scalar
      */
      for(j=0;j < nlon2;j++)
	valuex[j] = 0.0;               /* init with 0.0es */
      l = 0; m = -1;
      for (j=j2=0; j < shana->lmsize; j++,j2+=2) { /* loop through l,m */
	m++;
	if (m > l) {
	  m=0;
	  l++;
	}
	m2 = 2*m;
	//fprintf(stderr,"%11i %11i %11g %11g\n",l,m,cslm[j2],cslm[j2+1]);
	/*  add up contributions from all l,m  */
	valuex[m2]   +=  /* cos term */
	  plm[oplm+j] * (double)cslm[j2]; /* A coeff */
	valuex[m2+1] +=  /* sin term */
	  plm[oplm+j] * (double)cslm[j2+1];   /* B coeff */
      }

      /* compute inverse FFT  */
#ifdef NO_SHANA_FORTRAN      
      shana_cs2ab(valuex,shana->nlon);
      shana_realft_nr((valuex-1),shana->nlat,negunity);	
#else
      shana_f90_cs2ab(valuex,&shana->nlon);	
      shana_f90_realft(valuex,&shana->nlat,&negunity);	
#endif

      for (j=0; j < shana->nlon; j++) { /* can't vectorize */
	rdatax[ios1 + j] = valuex[j]/(HC_CPREC)(shana->nlat);
      }
      /* end scalar part */
    } else {
      /* 
	 vector harmonics
      */
      for(j=0;j < nlon2;j++){
	valuex[j] = valuey[j] = 0.0;
      }
      l = 1; 
      m = -1;       
      /* start at l = 1 */
      for (j=1,j2=2; j < shana->lmsize; j++,j2+=2) { 
	/* loop through l,m */
	m++;
	if (m > l) {
	  m=0;
	  l++;
	}
	m2  = 2*m;
	/* derivative factors */
	lm1 = l - 1;
	dpdt = dplm[oplm+j] * (double)shana->ell_factor[lm1]; /* d_theta(P_lm) factor */
	dpdp  = ((double)m) * plm[oplm+j]/ (double)shana->sin_theta[i];
	dpdp *= (double)shana->ell_factor[lm1]; /* d_phi (P_lm) factor */
	/* add up contributions from all l,m 

	u_theta

	*/
	valuex[m2] +=   /* cos term */
	  cslm[j2]   * dpdt + dslm[j2+1] * dpdp;
	valuex[m2+1] +=   /* sin term */
	  cslm[j2+1] * dpdt - dslm[j2]   * dpdp;
	/* 
	   u_phi
	*/
	valuey[m2] +=  // cos term
	  cslm[j2+1] * dpdp  - dslm[j2]   * dpdt;
	valuey[m2+1] +=   // sin term
	  - cslm[j2] * dpdp  - dslm[j2+1] * dpdt;
      }	/* end l,m loop */
        /* do inverse FFTs */
#ifdef NO_SHANA_FORTRAN
      shana_cs2ab(valuex,shana->nlon);
      shana_cs2ab(valuey,shana->nlon);
      shana_realft_nr((valuex-1),shana->nlat,negunity);
      shana_realft_nr((valuey-1),shana->nlat,negunity);
#else
      shana_f90_cs2ab(valuex,&shana->nlon);
      shana_f90_cs2ab(valuey,&shana->nlon);
      shana_f90_realft(valuex,&shana->nlat,&negunity);
      shana_f90_realft(valuey,&shana->nlat,&negunity);
#endif
      /* assign to output array */
      for (j=0; j < shana->nlon; j++) {   
	rdatax[ios1 + j] = valuex[j]/(HC_CPREC)(shana->nlat);
	rdatay[ios1 + j] = valuey[j]/(HC_CPREC)(shana->nlat);
      }
    }
  } /* end latitude loop */

  /* free temporary arrays */
  free(valuex);
  if(ivec)
    free(valuey);
  
}

void 
shana_shd2c (rdatax, rdatay, lmax, ivec, cslm, dslm, shana)
HC_CPREC *rdatax;
HC_CPREC *rdatay;
int lmax;
int ivec;
HC_CPREC *cslm;
HC_CPREC *dslm;
struct shana_module *shana;
{
  //
  //	Calculates spherical harmonic coefficients cslm(l,m) of
  //	a scalar (ivec = 0) or vector (ivec=1) function on a sphere. 
  //
  //     if ivec == 1, then cslm will be poloidal, and dslm toroidal
  //                   rdata will be nlon*nlat
  //
  //     Coefficients stored with
  //	cosine term and sine term alternating, starting at l=0
  //	and m=0 with m ranging from 0 to l before l is incremented.
  //
  //     Harmonics normalized with mean square = 1.0. Expansion is:
  //     SUM( Plm(cos(theta))*(cp(l,m)*cos(m*phi)+sp(l,m)*sin(m*phi))
  //
  //     Uses FFT in longitude and gaussian integration of spectral
  //     coefficients (for fixed m) times associated Legendre
  //     functions (of same order m) to get expansion coefficients.
  //
  //     Uses Num Rec routine for Gaussian points and weights, which should
  //     be initialized first by a call to shana_init. 
  //
  //     Uses recursive routine to generate associated Legendre functions.
  //     
  //     INPUT:
  //
  //     lmax    maximum possible degree of expansion (2**n-1). There
  //             are nlon = 2*lmax+2 data points in longitude for each latitude
  //             which has nlat = lmax+1 points
  //
  //     rdatax,rdatay   data((nlon * nlat)) arrays for theta and phi
  //                     components
  //
  //     ic_pd: should be either 0.0_cp or 1.0_cp
  //            if set to 1.0_cp, will expand velocities instead
  //            in this case, data should hold the theta and phi components
  //            of a vector field and cslm will hold the poloidal and toroidal
  //            on output components, respectively
  //
  //     OUTPUT:
  //
  //     cslm,dslm    coefficients, (2*lmsize)
  //
  //     dslm and rdatay will only be referenced when ivec = 1
  //
  //
  // local
  double *plm,*dplm;
  /* allocate memory */
  hc_dvecalloc(&plm,shana->nlat*shana->lmsize,"shana_shd2c: mem 1");
  hc_dvecalloc(&dplm,shana->nlat*shana->lmsize,"shana_shd2c: mem 2");
  // check
  if(!shana->initialized){
    fprintf(stderr,"shana_shd2c: error: initialize first\n");
    exit(-1);
  }
  if(lmax != shana->nlat-1){
    fprintf(stderr,"shana_shd2c: error: lmax mismatch: nlat: %i lmax: %i\n",
	    shana->nlat,lmax);
    exit(-1);
  }
  
  // compute the Plm first
  shana_compute_allplm(lmax,ivec,plm,dplm,shana);
  //
  // call the precomputed version
  shana_shd2c_pre(rdatax,rdatay,lmax,plm,dplm,ivec,cslm,dslm,shana);
  /* free legendre functions if not needed */
  free(plm);
  free(dplm);
}
//
// the actual routine to go from spatial to spectral, 
// for comments, see above
//
void 
shana_shd2c_pre (rdatax, rdatay, lmax, plm, dplm, ivec, cslm, dslm, shana)
HC_CPREC *rdatax;
HC_CPREC *rdatay;
int lmax;
double *plm;
double *dplm;
int ivec;
HC_CPREC *cslm;
HC_CPREC *dslm;
struct shana_module *shana;
{
  // local
  HC_CPREC *valuex, *valuey;
  double dfact,dpdt,dpdp;
  static int unity = 1;
  //
  int  lmaxp1,lmaxp1t2,i,j,l,m,ios1,m2,j2,oplm,nlon2,lm1;
  // check
  if(!shana->initialized){
    fprintf(stderr,"shana_shd2c_pre: error: initialize first\n");
    exit(-1);
  }
  // check some more and compute bounds
  lmaxp1 = lmax + 1;
  lmaxp1t2 = lmaxp1 * 2;
  nlon2 = shana->nlon + 2;
  if((lmaxp1 != shana->nlat)||(shana->nlon != lmaxp1t2)||
     ((lmax+1)*(lmax+2)/2 != shana->lmsize)){
    fprintf(stderr,"shana_shd2c_pre: dimension error, lmax %i\n",lmax);
    fprintf(stderr,"shana_shd2c_pre: nlon %i nlat %i\n",shana->nlon,shana->nlat);
    fprintf(stderr,"shana_shd2c_pre: lmsize %i\n",shana->lmsize);
    exit(-1);
  }
  /* allocate */
  shana_vecalloc(&valuex,nlon2,"shana_shd2c_pre 1");
  if(ivec)
    shana_vecalloc(&valuey,nlon2,"shana_shd2c_pre 2");
  
  //
  // initialize the coefficients
  //
  for(i=0;i < shana->lmsize2;i++)
    cslm[i] = 0.0;
  if(ivec){
    if(! shana->vector_sh_fac_init){
      fprintf(stderr,"shana_shd2c_pre: error: vector harmonics factors not initialized\n");
      exit(-1);
    }
    for(i=0;i < shana->lmsize2;i++)
      dslm[i] = 0.0;
  }
  for(i=0;i < shana->nlat;i++){
    //
    // loop through latitudes
    //
    ios1 = i * shana->nlon;          // offset for data array
    oplm = i * shana->lmsize;        // offset for Plm array
    //
    if(!ivec){
      //
      // scalar expansion
      //
      for(j=0;j < shana->nlon;j++){      
	valuex[j] = rdatax[ios1 + j];
      }
      //
      // compute the FFT 
      //

#ifdef NO_SHANA_FORTRAN
      shana_realft_nr((valuex-1),shana->nlat,unity);
      shana_ab2cs(valuex,shana->nlon);
#else
      shana_f90_realft(valuex,&shana->nlat,&unity);
      shana_f90_ab2cs(valuex,&shana->nlon);
#endif
      // sum up for integration
      l = 0;m = -1;
      for(j=j2=0;j < shana->lmsize;j++,j2+=2){
	m++;
	if( m >  l ) {
	  l++;m=0;
	}
	// we incorporate the Gauss integration weight and Plm factors here
	if (m == 0) {
	  dfact = ((double)shana->gauss_w[i] * plm[oplm+j])/2.0;
	}else{
	  dfact = ((double)shana->gauss_w[i] * plm[oplm+j])/4.0;
	}
	m2 = m * 2;
	cslm[j2]   += valuex[m2]   * dfact; // A coefficient
	cslm[j2+1] += valuex[m2+1] * dfact; // B coefficient
      }	/* end lmsize loop */
      /* end scalar */
    }else{
      //
      // vector field expansion
      //
      for(j=0;j < shana->nlon;j++){    
	valuex[j] = rdatax[ios1 + j]; // theta
	valuey[j] = rdatay[ios1 + j]; // phi
      }
      // perform the FFTs on both components
#ifdef NO_SHANA_FORTRAN
      shana_realft_nr((valuex-1),shana->nlat,unity);
      shana_realft_nr((valuey-1),shana->nlat,unity);
      shana_ab2cs(valuex,shana->nlon);
      shana_ab2cs(valuey,shana->nlon);
#else
      shana_f90_realft(valuex,&shana->nlat,&unity);
      shana_f90_realft(valuey,&shana->nlat,&unity);
      shana_f90_ab2cs(valuex,&shana->nlon);
      shana_f90_ab2cs(valuey,&shana->nlon);
#endif
      //
      l=1;m=-1;                // there's no l=0 term
      //
      for(j = 1,j2 = 2;j < shana->lmsize;j++,j2+=2){
	m++;
	if( m >l ) {
	  l++;m=0;
	}
	lm1 = l - 1;
	if (m == 0){ // ell_factor is 1/sqrt(l(l+1))
	  dfact = shana->gauss_w[i] * shana->ell_factor[lm1]/2.0;
	}else{
	  dfact = shana->gauss_w[i] * shana->ell_factor[lm1]/4.0;
	}
	//
	// some more factors
	//
	// d_theta(P_lm) factor
	dpdt = dplm[oplm+j];
	// d_phi (P_lm) factor
	dpdp = ((double)m) * plm[oplm+j]/(double)shana->sin_theta[i];
	//
	m2 = m * 2;
	//           print *,m,l,dpdt*dfact,dpdp*dfact
	/* poloidal */
	cslm[j2]   += // poloidal A
	  (dpdt * valuex[m2]   - dpdp * valuey[m2+1])*dfact;
	cslm[j2+1] += //   poloidal B
	  (dpdt * valuex[m2+1] + dpdp * valuey[m2]  )*dfact;
	/* toroidal */
	dslm[j2]   += // toroidal A 
	  (-dpdp * valuex[m2+1] - dpdt * valuey[m2]  )*dfact;
	dslm[j2+1] += // toroidal B 
	  ( dpdp * valuex[m2]   - dpdt * valuey[m2+1])*dfact;
      }
    }                      // end vector field
  }                        // end latitude loop

  free(valuex);
  if(ivec)
    free(valuey);
}


//
// initialize all necessary arrays for Shana type expansion
//
// if ivec == 1, will initialize for velocities/polarizations
//
void 
shana_init (lmax, ivec, npoints, nplm, tnplm, shana)
int lmax;
int ivec;
int *npoints;
int *nplm;
int *tnplm;
struct shana_module *shana;
{

  // input: lmax,ivec			
  // output: npoints,nplm,tnplm
  // local 
  HC_CPREC xtemp;
  int i,l;

  if(!shana->was_called){
    if(lmax == 0){
      fprintf(stderr,"shana_init: error: lmax is zero: %i\n",lmax);
      exit(-1);
    }
    //
    // test if lmax is 2**n-1
    //
    xtemp = log((HC_CPREC)(lmax+1))/log(2.0);
    if(fabs((int)(xtemp+.5)-xtemp) > 1e-7){
      fprintf(stderr,"shana_init: error: lmax has to be 2**n-1\n");
      fprintf(stderr,"shana_init: maybe use lmax = %i\n",
	      (int)pow(2,(int)(xtemp)+1)-1);
      fprintf(stderr,"shana_init: instead of %i\n",lmax);
      exit(-1);
    }
    //
    // number of longitudinal and latitudinal points
    //
    shana->nlat = lmax + 1;
    shana->nlon = 2 * shana->nlat;
    //
    // number of points in one layer
    //
    *npoints = shana->nlat * shana->nlon; 
    shana->old_npoints = *npoints;
    //
    //
    // for coordinate computations
    //
    shana->dphi = TWOPI / (HC_CPREC)(shana->nlon);
    shana->nlonm1 = shana->nlon - 1;
    //
    // size of tighly packed arrays with l,m indices
    shana->lmsize  = (lmax+1)*(lmax+2)/2;
    shana->lmsize2 = shana->lmsize * 2;          //for A and B
    //
    //
    // size of the Plm array 
    *nplm = shana->lmsize * shana->nlat;
    *tnplm = *nplm * (1+ivec);           // for all layers
    shana->old_tnplm = *tnplm;
    shana->old_nplm = *nplm;
    //
    // initialize the Gauss points, at which the latitudes are 
    // evaluated
    //
    shana_vecalloc(&shana->gauss_z,shana->nlat,"shana_init 1");
    shana_vecalloc(&shana->gauss_w,shana->nlat,"shana_init 2");
    shana_vecalloc(&shana->gauss_theta,shana->nlat,"shana_init 3");
    /* 
       gauss weighting 
    */
    shana_gauleg(-1.0,1.0,shana->gauss_z,shana->gauss_w,shana->nlat);
    //
    // theta values of the Gauss quadrature points
    //
    for(i=0;i < shana->nlat;i++)
      shana->gauss_theta[i] = acos(shana->gauss_z[i]);
    //
    // those will be used by plmbar to store some of the factors
    //
    hc_dvecalloc(&shana->plm_f1,shana->lmsize,"shana_init 4");
    hc_dvecalloc(&shana->plm_f2,shana->lmsize,"shana_init 5");
    hc_dvecalloc(&shana->plm_fac1,shana->lmsize,"shana_init 6");
    hc_dvecalloc(&shana->plm_fac2,shana->lmsize,"shana_init 7");
    hc_dvecalloc(&shana->plm_srt,shana->nlon,"shana_init 8");
    if(ivec){
      //
      // additional arrays for vector spherical harmonics
      // (perform the computations in HC_CPREC precision)
      //
      shana_vecalloc(&shana->ell_factor,shana->nlat,"shana init 9");
      shana_vecalloc(&shana->sin_theta,shana->nlat,"shana init 9");
      
      // 1/(l(l+1)) factors
      for(i=0,l=1;i < shana->nlat;i++,l++){
	// no l=0 term, obviously
	// go from l=1 to l=lmax+1
	shana->ell_factor[i] = 1.0/sqrt((HC_CPREC)(l*(l+1)));
	shana->sin_theta[i] = sqrt((1.0 - shana->gauss_z[i])*
				  (1.0+shana->gauss_z[i]));
	shana->vector_sh_fac_init = TRUE;
      }
    }else{
      shana->vector_sh_fac_init = FALSE;
    }
    //
    // logic flags
    //
    shana->computed_legendre = FALSE;
    shana->initialized = TRUE;
    shana->was_called = TRUE;
    /* 
       save initial call lmax and ivec settings
    */
    shana->old_lmax = lmax;
    shana->old_ivec = ivec;

    /* end initial branch */
  }else{
    if(lmax != shana->old_lmax){
      fprintf(stderr,"shana_init: error: was init with lmax %i, now: %i. (ivec: %i, now: %i)\n",
	      shana->old_lmax,lmax,shana->old_ivec,ivec);
      exit(-1);
    }
    if(ivec > shana->old_ivec){
      fprintf(stderr,"shana_init: error: original ivec %i, now %i\n",shana->old_ivec,ivec);
      exit(-1);
    }
    *npoints = shana->old_npoints;
    *nplm = shana->old_nplm;
    *tnplm = shana->old_tnplm;
  }
} /* end shana init */


//
// free all arrays that were allocate for the module   
//
void 
shana_free_module (shana, ivec)
struct shana_module *shana;
int ivec;
{
  // input: ivec
  
  free(shana->gauss_z);
  free(shana->gauss_w);
  free(shana->gauss_theta);
  if(shana->computed_legendre){
    // those are from the Legendre function action

    free(shana->plm_f1);free(shana->plm_f2);
    free(shana->plm_fac1);free(shana->plm_fac2);
    free(shana->plm_srt);
  }
  if(ivec){
    free(shana->ell_factor);free(shana->sin_theta);
  }
}
void 
shana_plmbar1 (p, dp, ivec, lmax, z, shana)
double *p;
double *dp;
int ivec;
int lmax;
HC_CPREC z;
struct shana_module *shana;
{
 
  double plm,pm1,pm2,pmm,sintsq,fnum,fden;
  //
  int i,l,m,k,kstart,l2,mstop,lmaxm1;

  if(!shana->initialized){
    fprintf(stderr,"shana_plmbar1: error: module not initialized, call shana_init first\n");
    exit(-1);
  }
  if ((lmax < 0) || (fabs(z) > 1.0)) {
    fprintf(stderr,"shana_plmbar1: error: bad arguments\n");
    exit(-1);
  }
  if(!shana->computed_legendre) {
    /* 
       need to initialize the legendre factors 
    */
    if(shana->nlon != (lmax+1)*2){
      fprintf(stderr,"shana_plmbar1: factor mismatch, lmax: %i vs nlon: %i (needs to be (lmax+1)*2)\n",
	      lmax,shana->nlon);
      exit(-1);
    }
    //
    // first call, set up some factors. the arrays were allocated in shana_init
    //
    for(k=0,i=1;k < shana->nlon;k++,i++){
      shana->plm_srt[k] = sqrt((double)(i));
    }
    // initialize plm factors
    for(i=0;i < shana->lmsize;i++){
      shana->plm_f1[i] = shana->plm_fac1[i] = 0.0;
      shana->plm_f2[i] = shana->plm_fac2[i] = 0.0;
    }
    //     --case for m > 0
    kstart = 0;
    for(m=1;m <= lmax;m++){
      //     --case for P(m,m) 
      kstart += m+1;
      if (m != lmax) {
	//     --case for P(m+1,m)
	k = kstart+m+1;		
	//     --case for P(l,m) with l > m+1
	if (m < (lmax-1)) {
	  for(l = m+2;l <= lmax;l++){
	    l2 = l * 2;	
	    k = k+l;
	    shana->plm_f1[k] = shana->plm_srt[l2] * shana->plm_srt[l2-2]/
	      (shana->plm_srt[l+m-1] * shana->plm_srt[l-m-1]);
	    shana->plm_f2[k]=(shana->plm_srt[l2] * shana->plm_srt[l-m-2]*shana->plm_srt[l+m-2])/
	      (shana->plm_srt[l2-4] * shana->plm_srt[l+m-1] * shana->plm_srt[l-m-1]);
	  }
	}
      }
    }
    if(ivec){
      //
      // for derivative of Plm with resp. to theta
      // (if we forget to call Plm with ivec=1, but use
      // ivec=1 later, we will notice as all should be zero)
      //
      k=2;			
      for(l=2;l<=lmax;l++){
	k++;
	mstop = l - 1;
	for(m=1;m <= mstop;m++){
	  k++;
	  shana->plm_fac1[k] = shana->plm_srt[l-m-1] * shana->plm_srt[l+m];
	  shana->plm_fac2[k] = shana->plm_srt[l+m-1] * shana->plm_srt[l-m];
	  if(m == 1){
	    shana->plm_fac2[k] = shana->plm_fac2[k] * shana->plm_srt[1];
	  }
	}
	k++;
      }
    } /* end ivec==1 */
    shana->old_lmax = lmax;
    shana->old_ivec = ivec;
    shana->computed_legendre = TRUE;
    /* 
       end first call 
    */
  }else{
    /* 
       subsequent call 
    */
    // test if lmax has changed
    if(lmax != shana->old_lmax){
      fprintf(stderr,"shana_plmbar1: error: factors were computed for lmax %in",shana->old_lmax);
      fprintf(stderr,"shana_plmbar1: error: now, lmax is %i\n",lmax);
      exit(-1);
    }
    if(ivec > shana->old_ivec){
      fprintf(stderr,"shana_plmbar1: error: init with %i, now ivec %i\n",shana->old_ivec,ivec);
      exit(-1);
    }
  }
  /* 

  what follows will be executed for each z

  */
  //     --start calculation of Plm etc.
  //     --case for P(l,0) 
  
  pm2  = 1.0;
  p[0] = 1.0;                   // (0,0)
  if(ivec)
    dp[0] = 0.0;// else, don't refer to this array
  if (lmax == 0){
    fprintf(stderr,"lmax is zero. what the hell?//\n");
    exit(-1);
  }
  pm1  = z;
  p[1] = shana->plm_srt[2] * pm1;             // (1,0)
  k=1;
  for(l = 2;l<=lmax;l++){               // now all (l,0)
    k  += l;
    l2 = l * 2;
    plm = ((double)(l2-1) * z * pm1 - (double)(l-1) * pm2)/((double)l);
    p[k] = shana->plm_srt[l2] * plm;
    pm2 = pm1;
    pm1 = plm;
  }
  //       --case for m > 0
  pmm = 1.0;
  sintsq = (1.0 - z) * (1.0 + z);
  fnum = -1.0;
  fden =  0.0;
  kstart = 0;
  lmaxm1 = lmax - 1;
  for(m = 1;m<=lmax;m++){
    //     --case for P(m,m) 
    kstart += m+1;
    fnum += 2.0;
    fden += 2.0;
    pmm = pmm * sintsq * fnum / fden;
    pm2 = sqrt((HC_CPREC)(4*m+2)*pmm);
    p[kstart] = pm2;
    if (m != lmax) {
      //     --case for P(m+1,m)
      pm1 = z * shana->plm_srt[2*m+2]*pm2;
      k = kstart + m + 1;
      p[k] = pm1;
      //     --case for P(l,m) with l > m+1
      if (m < lmaxm1) {
	for(l = m+2;l <= lmax;l++){
	  k += l;
	  plm = z * shana->plm_f1[k] * pm1 - 
	    shana->plm_f2[k] * pm2;
	  p[k] = plm;
	  pm2 = pm1;
	  pm1 = plm;
	}
      }
    }
  }
  if(ivec){
    // 
    // derivatives
    //     ---derivatives of P(z) wrt theta, where z=cos(theta)
    //     
    dp[1] = -p[2];
    dp[2] =  p[1];
    k = 2;
    for(l=2;l <= lmax;l++){
      k++;
      //     treat m=0 and m=l separately
      dp[k] =  -shana->plm_srt[l-1] * shana->plm_srt[l] / shana->plm_srt[1] * p[k+1];
      dp[k+l] = shana->plm_srt[l-1] / shana->plm_srt[1] * p[k+l-1];
      mstop = l-1;
      for(m=1;m <= mstop;m++){
	k++;
	dp[k] = shana->plm_fac2[k] * p[k-1] - 
	  shana->plm_fac1[k] * p[k+1];
	dp[k] *= 0.5;
      }
      k++;
    }
  }
}



