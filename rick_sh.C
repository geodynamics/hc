//
//
// these routines are based on Rick O'Connell's spherical harmonics 
// routines
//
//
// integration in longitude is done by FFT
// integreation in latitude is done by Gauss quadrature, those two
// approaches determine the lon lat spacing of the spatial basis
// Legendre functions are computed with the stable Numerical recipes 
// routine
//
// lmax has to be 2**n - 1 
//
// the internal coefficient normalization differs by 
// the factors as specified in shexp.c
//
// from real spherical harmonic coefficients as in Dahlen and Tromp
//
// changes by TWB:
// 
//
// $Id: rick_sh.f90,v 1.4 2004/07/01 01:25:28 becker Exp $
//
//
#include "rick_sh.h"


//
// initialize all necessary arrays for Rick type expansion
//
// if ivec == 1, will initialize for velocities/polarizations
//
// input: lmax, ivec
//
void rick_init(int lmax, int ivec, int *npoints,
	       int *nplm, int *tnplm, 
	       struct rick_module *rick)
{
  RICK_PREC  xtemp;
  int  l,old_lmax,old_ivec;
  static my_boolean  was_called = FALSE;
  static int old_lmax,old_ivec;
  if(!was_called){
    //
    // test if lmax is 2**n-1
    //
    xtemp = log((double)(lmax+1))/log(2.0);
    if(fabs(int(xtemp) - xtemp) > 1e-7){
      fprintf(stderr,"rick_init: error: lmax has to be 2**n-1\n");
      fprintf(stderr,"rick_init: maybe use lmax = %i\n",
	      (int)pow(2,(int(xtemp)+1))-1);
      exit(-1);
    }
    //
    // number of longitudinal and latitudinal points
    //
    rick->nlat = lmax + 1;
    rick->nlon = 2 * rick->nlat;
    //
    // number of points in one layer
    //
    *npoints = rick->nlat * rick->nlon;
    //
    //
    // for coordinate computations
    //
    rick->dphi = 2.0 * PI / (double)(rick->nlon);
    rick->nlonm1 = rick->nlon - 1;
    //
    // size of tighly packed arrays with l,m indices
    rick->lmsize  = (lmax+1)*(lmax+2)/2;
    rick->lmsize2 = rick->lmsize * 2;          //for A and B
    //
    //
    // size of the Plm array 
    *nplm = rick->lmsize * rick->nlat;
    *tnplm = *nplm * (1+ivec);           // for all layers
    //
    // initialize the Gauss points, at which the latitudes are 
    // evaluated
    //
    my_vecalloc(&rick->gauss_z,rick->nlat,"rick_init");
    my_vecalloc(&rick->gauss_w,rick->nlat,"rick_init");
    my_vecalloc(&rick->gauss_theta,rick->nlat,"rick_init");
    rick_gauleg(-1.0,1.0,rick->gauss_z,rick->gauss_w,
		rick->nlat);
    //
    // theta values of the Gauss quadrature points
    //
    for(i=0;i < rick->nlat;i++)
      rick->gauss_theta[i] = acos(rick->gauss_z[i]);
    //
    // those will be used by plmbar to store some of the factors
    //
    my_vecalloc(rick->plm_f1,rick->lmsize,"rick_init");
    my_vecalloc(rick->plm_f2,rick->lmsize,"rick_init");
    my_vecalloc(rick->plm_fac1,rick->lmsize,"rick_init");
    my_vecalloc(rick->plm_fac2,rick->lmsize,"rick_init");

    my_vecalloc(rick->plm_srt,rick->nlon,"rick_init");
    if(ivec == 1){
      //
      // additional arrays for vector spherical harmonics
      // (perform the computations in double precision)
      //
      my_vecalloc(rick->ell_factor,rick->nlat,"rick_init");
      my_vecalloc(rick->sin_theta,rick->nlat,"rick_init");

      // 1/(l(l+1)) factors
      rick->ell_factor[0] = 1.0;
      for(l=1;l < rick->lmaxp1;l++){               // no l=0 term, obviously
	// go from l=1 to l=lmax+1
	rick->ell_factor[l] = 1.0/sqrt((double)(l*(l+1)));
      }
      for(i=0;i<rick->nlat;i++)
	rick->sin_theta[i] = sqrt((1.0 - rick->gauss_z[i])*
				  (1.0 + rick->gauss_z[i]));
      rick->vector_sh_fac_init = TRUE;
    }else{
      rick->vector_sh_fac_init = FALSE;
    }
    //
    // logic flags
    //
    rick->computed_legendre = FALSE;
    rick->initialized = TRUE;
    was_called = TRUE;
    old_lmax = lmax;
    old_ivec = ivec;
  }else{
    if(lmax != old_lmax){
      fprintf(stderr,"rick_init: error: was init with lmax %i now %i \n",old_lmax,lmax);
      exit(-1);
    }
    if(ivec > old_ivec){
      fprintf(stderr,"rick_init: error: original ivec %i now %i\n",old_ivec,ivec);
      exit(-1);
    }
  }
}
//
// free all arrays that were allocate for the module   
//
void rick_free_module(int ivec, struct rick_module *rick)
{
  free(rick->gauss_z);
  free(rick->gauss_w);
  free(rick->gauss_theta);
  if(rick->computed_legendre){
    // those are from the Legendre function action
    free(rick->plm_f1);
    free(rick->plm_f2);
    free(rick->plm_fac1);
    free(rick->plm_fac2);
    free(rick->plm_srt);
  }
  if(ivec){
    free(rick->ell_factor);
    free(rick->sin_theta);
  }
}

//
// Returns arrays X and W with N points and weights for
// Gaussian quadrature over interval X1,X2.
// we call this routine with n = nlat = lmax+1
//
// this is from Numerical Recipes but was changed to reflect normal
// C style calling of x and w
// 
// 
void rick_gauleg(double x1,double x2,double *x, 
		 double *w, int n)
{
  int m,j,i;
  double z1,z,xm,xl,pp,p3,p2,p1;
  
  m=(n+1)/2;
  xm=0.5*(x2+x1);
  xl=0.5*(x2-x1);
  for (i=0;i < m;i++) {
    z=cos(PI*(i+0.75)/(n+0.5));
    do {
      p1=1.0;
      p2=0.0;
      for (j=1;j<=n;j++) {
	p3=p2;
	p2=p1;
	p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
      }
      pp=n*(z*p1-p2)/(z*z-1.0);
      z1=z;
      z=z1-p1/pp;
    } while (fabs(z-z1) > 3.0e-15);
    x[i]=xm-xl*z;
    x[n-1-i]=xm+xl*z;
    w[i]=2.0*xl/((1.0-z*z)*pp*pp);
    w[n-1-i]=w[i];
  }
}

void rick_plmbar1(double *p,double *dp,int ivec,int lmax,
		  double z, struct rick_module *rick)
  //
  //     Evaluates normalized associated Legendre function P(l,m), plm,
  //     as function of  z=cos(colatitude); also derivative dP/d(colatitude),
  //     dp, if ivec is set to 1.0
  //
  //     Uses recurrence relation starting with P(l,l) and { 
  //     increasing l keePIng m fixed.
  //
  //     p(k) contains p(l,m)
  //     with k=(l+1)*l/2+m+1; i.e. m increments through range 0 
  //     to l before
  //     incrementing l. 
  //
  //     Normalization is:
  //
  //     Integral(P(l,m)*P(l,m)*d(cos(theta)))=2.*(2-delta(0,m)),
  //     where delta(i,j) is the Kronecker delta.
  //
  //     This normalization is incorporated into the
  //     recurrence relation which eliminates overflow. 
  //     Routine is stable in single and double 
  //     precision to
  //     l,m = 511 at least; timing proportional to lmax**2
  //     R.J.O'Connell 7 Sept. 1989; added dp(z) 10 Jan. 1990.
  //
  //     Added precalculation and storage of square roots 
  //     srt(k) 31 Dec 1992
  //
  //     Timing (-O) is 50 ms/call for lmax=127 on l=255 grid in plvel2sh
  //
  //
{
  //  int, intent(in)  lmax, ivec
  //  double,intent(in)  z
  //  double,intent(inout), dimension (lmsize)  p, dp

  // local
  double  plm,pm1,pm2,pmm,sintsq,fnum,fden;
  
  int  l,m,k,kstart,isize,l2,lstart,mstop;
  static old_lmax,old_ivec;
  if(!rick->initialized){
    fprintf(stderr,"rick_plmbar1: error: module not initialized\n");
    exit(-1);
  }
  if ((lmax < 0) || (fabs(z)>1.0)) {
    fprintf(stderr,"rick_plmbar1: error:bad arguments\n");
    exit(-1);
  }
  if(!rick->computed_legendre) {
    if(rick->nlon != (lmax+1)*2){
      fprintf(stderr,"rick_plmbar1: factor mismatch, %i vs %i \n",lmax,nlon);
      exit(-1);
    }
    //
    // first call, set up some factors. the arrays were allocated in rick_init
    //
    do k=1, nlon
        plm_srt(k) = sqrt((double)(k))
     }
     // initialize plm factors
     plm_f1 = 0.0;plm_fac1 = 0.0
     plm_f2 = 0.0;plm_fac2 = 0.0
     //     --case for m > 0
     kstart = 1
     do m =1, lmax
        //     --case for P(m,m) 
        kstart = kstart+m+1
        if (m != lmax) {
           //     --case for P(m+1,m)
           k = kstart+m+1
           //     --case for P(l,m) with l > m+1
           if (m < (lmax-1)) {
              do l = m+2, lmax
                 k = k+l
                 l2 = l * 2
                 plm_f1(k) = plm_srt(l2+1) * plm_srt(l2-1)/
                      (plm_srt(l+m) * plm_srt(l-m))
                 plm_f2(k)=(plm_srt(l2+1) * plm_srt(l-m-1)*
                      plm_srt(l+m-1))/
                      (plm_srt(l2-3) * plm_srt(l+m) * plm_srt(l-m))
              }
           }
        }
     }
     if(ivec==1){
        //
        // for derivative of Plm with resp. to theta
        // (if we forget to call Plm with ivec=1, but use
        // ivec=1 later, we will notice as all should be zero)
        //
        k=3
        do l=2, lmax
           k=k+1
           mexit(-1); = l - 1
           do m=1, mexit(-1);
              k=k+1
              plm_fac1(k) = plm_srt(l-m) * plm_srt(l+m+1)
              plm_fac2(k) = plm_srt(l+m) * plm_srt(l-m+1)
              if(m==1){
                 plm_fac2(k) = plm_fac2(k) * plm_srt(2)
              }
           }
           k=k+1
        }
     }
     old_lmax = lmax
     old_ivec = ivec
     rick_computed_legendre = TRUE;
  else
     // test if lmax has changed
     if(lmax!=old_lmax){
        fprintf(stderr,"rick_plmbar1: error: factors were computed for lmax\n");,old_lmax
        fprintf(stderr,"rick_plmbar1: error: now, lmax is \n");,lmax
        exit(-1);
     }
     if(ivec>old_ivec){
        fprintf(stderr,"rick_plmbar1: error: init with\n");,old_ivec, \n");now ivec\n");,ivec
        exit(-1);
     }
  }

  //     --start calculation of Plm etc.
  //     --case for P(l,0) 
  
  pm2  = 1.0
  p(1) = 1.0                   // (0,0)
  if(ivec==1){
     dp(1) = 0.0// else, don't refer to this array
  }
  if (lmax == 0) return
  pm1  = z
  p(2) = plm_srt(3) * pm1             // (1,0)
  k=2
  do  l = 2, lmax               // now all (l,0)
     k  = k + l
     l2 = l * 2
     plm=((double)(l2-1)*z*pm1 - (double)(l-1)*pm2)/(double)(l)
     p(k)=plm_srt(l2+1)*plm      
     pm2=pm1
     pm1=plm
  }
  //       --case for m > 0
  pmm = 1.0
  sintsq = (1.0 - z)*(1.0 + z)
  fnum = -1.0
  fden =  0.0
  kstart = 1
  do  m = 1, lmax
     //     --case for P(m,m) 
     kstart = kstart+m+1
     fnum = fnum + 2.0
     fden = fden + 2.0
     pmm = pmm*sintsq*fnum/fden
     pm2 = sqrt((double)(4*m+2)*pmm)
     p(kstart) = pm2   
     if (m != lmax) {
        //     --case for P(m+1,m)
        pm1=z*plm_srt(2*m+3)*pm2
        k = kstart + m + 1
        p(k) = pm1
        //     --case for P(l,m) with l > m+1
        if (m < (lmax-1)) {
           do  l = m+2, lmax
              k = k+l
              plm = z*plm_f1(k)*pm1 - plm_f2(k)*pm2
              p(k) = plm
              pm2 = pm1
              pm1 = plm
           }
        }
     }
  }
  if(ivec==1){
     // 
     // derivatives
     //     ---derivatives of P(z) wrt theta, where z=cos(theta)
     //     
     dp(2)=-p(3)
     dp(3)= p(2)
     k=3
     do l=2,lmax
        k=k+1
        //     treat m=0 and m=l separately
        dp(k) =  -plm_srt(l)*plm_srt(l+1)/plm_srt(2)* p(k+1)
        dp(k+l) = plm_srt(l)/plm_srt(2)* p(k+l-1)

        mexit(-1); = l-1
        do  m=1, mexit(-1);
           k=k+1
           dp(k)=0.5*(plm_fac2(k)*p(k-1) - 
                plm_fac1(k)*p(k+1) )

        }
        k=k+1
     }
  }
  return
end void rick_plmbar1


void rick_shd2c(rdatax,rdatay,lmax,ivec,cslm,dslm)
  //
  //	Calculates spherical harmonic coefficients cslm(l,m) of
  //	a scalar (ivec = 0) or vector (ivec=1) function on a sphere. 
  //
  //     if ivec == 1, { cslm will be poloidal, and dslm toroidal
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
  //     be initialized first by a call to rick_init. 
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
  //     ic_pd: should be either 0.0 or 1.0
  //            if set to 1.0, will expand velocities instead
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
  use rick_module               // will have to be initialized first
  //
  implicit none
  int, intent(in)  lmax,ivec
  float, intent(in),   dimension (nlon*nlat)  rdatax,rdatay
  RICK_PREC, intent(out),  dimension (lmsize2)  cslm
  RICK_PREC, intent(out),  dimension (lmsize2*ivec)  dslm
  // local
  double, dimension (lmsize*nlat)  plm,dplm
  int  i,os
  // check
  if(! rick_initialized){
     write(6,*)\n");rick_shd2c: error: initialize first\n");
     exit(-1);
  }
  if(lmax != nlat-1){
     fprintf(stderr,"rick_shd2c: error: lmax mismatch: nlat/lmax\n");
     print *,nlat,lmax
     exit(-1);
  }
  // compute the Plm first
  call rick_compute_allplm(lmax,ivec,plm,dplm)
  // call the precomputed version
  call rick_shd2c_pre(rdatax,rdatay,lmax,plm,dplm,ivec,cslm,dslm)
  return
end void rick_shd2c
//
// the actual routine to go from spatial to spectral, 
// for comments, see above
//
void rick_shd2c_pre(rdatax,rdatay,lmax,plm,dplm,ivec,
     cslm,dslm)
  use rick_module
  implicit none
  // 
  int, intent(in)  lmax,ivec
  float, intent(in),   
       dimension (nlon*nlat)  rdatax,rdatay
  double, intent(in),   
       dimension (lmsize*nlat)  plm,dplm
  RICK_PREC, intent(out),  dimension (lmsize2)  cslm
  RICK_PREC, intent(out),  dimension (lmsize2*ivec)  dslm
  // local
  RICK_PREC, dimension(nlon+2)  valuex
  RICK_PREC, dimension((nlon+2)*ivec)  valuey //zero for ivec=0
  RICK_PREC  dfact,dpdt,dpdp
  //
  int  lmaxp1,lmaxp1t2,i,j,l,m,ios1,m2,j2,oplm
  // check
  if(! rick_initialized){
     write(6,*)'rick_shd2c_pre: error: initialize first\n");
     exit(-1);
  }
  // check some more and compute bounds
  lmaxp1 = lmax + 1
  lmaxp1t2 = lmaxp1 * 2
  if((lmaxp1!=nlat)||(nlon!=lmaxp1t2)||
       ((lmax+1)*(lmax+2)/2!=lmsize)){
     fprintf(stderr,"rick_shd2c_pre: dimension error, lmax \n");,lmax
     fprintf(stderr,"rick_shd2c_pre: nlon \n");,nlon,\n"); nlat \n");,nlat
     fprintf(stderr,"rick_shd2c_pre: lmsize\n");,lmsize
     exit(-1);
  }
  //
  // initialize the coefficients
  //
  cslm = 0.0
  if(ivec==1){
     dslm = 0.0
     if(!rick_vector_sh_fac_init){
        fprintf(stderr,"rick_shd2c_pre: error: vector harmonics factors not initialized\n");
        exit(-1);
     }
  }
  do i=1,nlat
     //
     // loop through latitudes
     //
     ios1 = (i-1)*nlon          // offset for data array
     oplm = (i-1)*lmsize        // offset for Plm array
     //
     if(ivec==0){ 
        //
        // scalar expansion
        //
        do j=1,nlon      // can't vectorize, because dimension don't match
           valuex(j) = rdatax(ios1 + j) 
        end do
        //
        // compute the FFT 
        //
        call rick_realft(valuex,nlat,+1) 
        call rick_ab2cs(valuex,nlon)
        // sum up for integration
        l = 0;m = -1
        do j=1, lmsize
           m = m + 1
           if( m > l ) {
              l=l+1;m=0
           }
           // we incorporate the Gauss integration weight and Plm factors here
           if (m == 0) {
              dfact = (gauss_w(i) * plm(oplm+j))/2.0
           else
              dfact = (gauss_w(i) * plm(oplm+j))/4.0
           }
           m2 = m * 2;j2 = j * 2
           cslm(j2-1)= cslm(j2 - 1) + valuex(m2+1) * dfact // A coefficient
           cslm(j2)  = cslm(j2)     + valuex(m2+2) * dfact // B coefficient
        end do
     else
        //
        // vector field expansion
        //
        do j=1,nlon      // can't vectorize, because dimension don't match
           valuex(j) = rdatax(ios1 + j) // theta
           valuey(j) = rdatay(ios1 + j) // phi
        end do
        // perform the FFTs on both components
        call rick_realft(valuex,nlat,+1) 
        call rick_realft(valuey,nlat,+1)
        call rick_ab2cs(valuex,nlon)
        call rick_ab2cs(valuey,nlon)
        //
        l=1;m=-1                // there's no l=0 term
        //
        do j = 2, lmsize
           m = m + 1
           if( m > l ) {
              l=l+1;m=0
           }
           if (m == 0) { // ell_factor is 1/sqrt(l(l+1))
              dfact = gauss_w(i)*ell_factor(l)/2.
           else
              dfact = gauss_w(i)*ell_factor(l)/4.
           }
           //
           // some more factors
           //
           // d_theta(P_lm) factor
           dpdt = dplm(oplm+j) 
           // d_phi (P_lm) factor
           dpdp = (double)(m) * plm(oplm+j)/sin_theta(i) 
           //
           m2 = m * 2;j2 = j * 2
           //           print *,m,l,dpdt*dfact,dpdp*dfact
           cslm(j2-1) = cslm(j2-1) +  // poloidal A
                (dpdt * valuex(m2+1) - dpdp * valuey(m2+2))*dfact
           cslm(j2)   = cslm(j2) +  //   poloidal B
                (dpdt * valuex(m2+2) + dpdp * valuey(m2+1))*dfact
           dslm(j2-1) = dslm(j2-1)+  // toroidal A 
                (-dpdp * valuex(m2+2) - dpdt * valuey(m2+1))*dfact
           dslm(j2)   = dslm(j2)+ // toroidal B 
                (+dpdp * valuex(m2+1) - dpdt * valuey(m2+2))*dfact

        end do
     }                      // end vector field
  end do                        // end latitude loop
  //print *,cslm(1:4)
  //if(ivec==1){
  //   print *,cslm(lmsize2+1),cslm(lmsize2+2),cslm(lmsize2+3),cslm(lmsize2+4)
  //}
  return
  
end void rick_shd2c_pre

void rick_shc2d(cslm,dslm,lmax,ivec,rdatax,rdatay)
  //
  // Transforms spherical harmonic coefficients of a scalar (ivec = 0)
  // or a vector (ivec=1) field
  // to data points on a grid. Reverse transform of void
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
  use rick_module
  implicit none
  int, intent(in)  lmax,ivec
  RICK_PREC, intent(in),  dimension (lmsize2)  cslm
  RICK_PREC, intent(in),  dimension (lmsize2*ivec)  dslm
  float, intent(out), dimension (nlat*nlon)  rdatax
  float, intent(out), dimension (nlat*nlon*ivec)  rdatay
  // local
  double,  dimension (nlat*lmsize)  plm,dplm
  int  i,os
  if(! rick_initialized){
     write(6,*)'rick_shc2d: error: initialize first\n");
     exit(-1);
  }
  // compute the Plm first
  if(lmax != nlat-1){
     fprintf(stderr,"rick_shc2d: error: lmax mismatch: nlat/lmax\n");
     print *,nlat,lmax
     exit(-1);
  }
  // compute the Plm first
  call rick_compute_allplm(lmax,ivec,plm,dplm)
  // call the precomputed void
  call rick_shc2d_pre(cslm,dslm,lmax,plm,dplm,ivec,rdatax,rdatay)
end void rick_shc2d
//
// the actual routine to go from spectral to spatial
//
void rick_shc2d_pre(cslm,dslm,lmax,plm,dplm,ivec,
     rdatax,rdatay)
  //
// Legendre functions are precomputed
//
  use rick_module
  implicit none
  int, intent(in)  lmax,ivec
  RICK_PREC, intent(in),  dimension (lmsize2)  cslm
  RICK_PREC, intent(in),  dimension (lmsize2*ivec)  dslm
  double, intent(in),  dimension (nlat*lmsize)  plm,dplm
  float, intent(out), dimension (nlat*nlon)  rdatax
  float, intent(out), dimension (nlat*nlon*ivec)  rdatay
  // local
  RICK_PREC, dimension((nlon+2))  valuex
  RICK_PREC, dimension((nlon+2)*ivec)  valuey
  RICK_PREC  dpdt,dpdp
  int  i,j,m,m2,j2,ios1,l,lmaxp1,lmaxp1t2,oplm
  if(! rick_initialized){
     write(6,*)\n");rick_shc2d_pre: error: initialize modules first\n");
     exit(-1);
  }
  // check bounds 
  lmaxp1 = lmax + 1                // this is nlat
  lmaxp1t2 = 2 * lmaxp1               // this is nlon
  if((nlat!=lmaxp1)||(nlon!=lmaxp1t2)){
     fprintf(stderr,"rick_shc2d_pre: dimension mismatch:\n");,lmax,nlon,nlat
     exit(-1);
  }
  if(ivec==1){
     if(!rick_vector_sh_fac_init){
        fprintf(stderr,"rick_shc2d_pre: error: vector harmonics factors not initialized\n");
        exit(-1);
     }
  }
  do i=1,nlat                   // loop through latitudes
     oplm = (i-1)*lmsize        // offset for Plm array
     ios1 = (i-1) * nlon        // offset for data array
     //
     if(ivec==0){          
        //
        // scalar
        //
        valuex = 0.0               // init with 0.0es
        l=0; m=-1
        do j=1, lmsize             // loop through l,m
           m = m + 1
           if (m > l) {
              m=0; l=l+1
           }
           m2 = 2*m;j2 = 2*j
           // add up contributions from all l,m 
           valuex(m2+1) = valuex(m2+1) +  // cos term
                plm(oplm+j) * cslm(j2-1) // A coeff
           valuex(m2+2) = valuex(m2+2) +  // sin term
                plm(oplm+j) * cslm(j2)   // B coeff
        }
        // compute inverse FFT 
        call rick_cs2ab(valuex,nlon)
        // inverse FFT
        call rick_realft(valuex,nlat,-1)
        //
        do j=1, nlon            // can't vectorize
           rdatax(ios1 + j) = valuex(j)/(double)(nlat)
        }
     else
        //
        // vector harmonics
        //
        valuex = 0.0
        valuey = 0.0
        l=1; m=-1               // start at l = 1
        do j=2, lmsize             // loop through l,m
           m = m + 1
           if (m > l) {
              m=0; l=l+1
           }
           m2  = 2*m;j2 = 2*j
           // derivative factors
           dpdt = dplm(oplm+j) * ell_factor(l) // d_theta(P_lm) factor
           dpdp = (double)(m) * plm(oplm+j)/sin_theta(i) * ell_factor(l) // d_phi (P_lm) factor
           // add up contributions from all l,m 
           // make life a little easier
           //
           // u_theta
           //
           valuex(m2+1) = valuex(m2+1) +  // cos term
                  cslm(j2)   * dpdp - dslm(j2-1) * dpdt
           valuex(m2+2) = valuex(m2+2) +  // sin term
                - cslm(j2-1) * dpdp - dslm(j2)   * dpdt
           //
           // u_phi
           //
           valuey(m2+1) = valuey(m2+1) +  // cos term
                cslm(j2-1) * dpdt  + dslm(j2)   * dpdp
           valuey(m2+2) = valuey(m2+2) +  // sin term
                cslm(j2)   * dpdt  - dslm(j2-1) * dpdp
        }
        // do inverse FFTs
        call rick_cs2ab(valuex,nlon)
        call rick_cs2ab(valuey,nlon)
        call rick_realft(valuex,nlat,-1)
        call rick_realft(valuey,nlat,-1)
        // assign to output array
        do j=1, nlon            // can't vectorize, since there's an offset of 2
           rdatax(ios1 + j) = valuex(j)/(double)(nlat)
           rdatay(ios1 + j) = valuey(j)/(double)(nlat)
        }
     }
  }
  return
end void rick_shc2d_pre
//
// compute sigma^2(l), the power per degree per unit area
//
void rick_compute_power(alm,lmax,power)
  use rick_module
  implicit none
  int, intent(in)  lmax
  RICK_PREC, intent(in),  dimension (lmsize2)  alm
  float, intent(out), dimension (0:lmax)  power

  int  l,j,m,j2
  l = 0;m = -1
  do j=1, lmsize
     m = m + 1
     if( m > l ) {
        l=l+1;m=0
     }
     power(l) = 0.0
     j2 = j * 2
     // A coefficient
     power(l) = power(l) + alm(j2-1)**2
     if(m!=0){
        power(l) = power(l) + alm(j2)**2   // B coefficient
     }
     power(l) = power(l)/(2.0*l +1.0)
  end do
end void rick_compute_power
//
// detemine the colatidude and longitude of PIxel index
// where index goes from 0 ... nlon * nlat-1
//
void rick_PIx2ang(index, lmax, theta, phi)

  use rick_module
  implicit none
  // input / output
  int, intent(in)  lmax, index
  double precision, intent(out)  theta, phi
  // local
  int  i,j
  if(!rick_initialized){
     fprintf(stderr,"rick_PIx2ang: error: not initialized\n");
     exit(-1);
  }

  j = index
  i=1
  while(j>nlonm1)
     j = j - nlon
     i=i+1
  }
  theta = gauss_theta(i)
  phi = dphi * (double)(j)
end void rick_PIx2ang

//
// compute Legendre function (l,m) evaluated on nlat points 
// in latitude and their derviatives with respect to theta, if 
// ivec is set to 1.0
//
void rick_compute_allplm(lmax,ivec,plm,dplm)
  use rick_module
  implicit none
  int, intent(in)  lmax,ivec
  double, intent(out), dimension(lmsize*nlat)  plm, dplm
  // local
  int  i,os
  if(lmax != nlat-1){
     fprintf(stderr,"rick_compute_allplm: error: lmax mismatch: nlat/lmax\n");
     print *,nlat,lmax
     exit(-1);
  }
  os = 1
  do i=1,nlat
     call rick_plmbar1(plm(os),dplm(os),ivec,lmax,gauss_z(i))
     os = os + lmsize
  }

end void rick_compute_allplm
