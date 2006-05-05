!
!
! these routines are based on Rick O'Connell's spherical harmonics 
! routines
!
!
! integration in longitude is done by FFT
! integreation in latitude is done by Gauss quadrature, those two
! approaches determine the lon lat spacing of the spatial basis
! Legendre functions are computed with the stable Numerical recipes 
! routine
!
! lmax has to be 2**n - 1 
!
! the internal coefficient normalization differs by 
! the factors as specified in shexp.c
!
! from real spherical harmonic coefficients as in Dahlen and Tromp
!
! changes by TWB:
! 
!
! $Id: rick_sh.f90,v 1.7 2006/01/22 01:11:34 becker Exp $
!
!
module rick_module
  ! read the defines for single/double precision
#include "sh_rick_ftrn.h"
  ! constants
  REAL(kind=cp), parameter, public :: pi = 3.141592653589793238462643383279502884197_cp
  REAL(kind=cp), parameter, public :: twopi = 6.283185307179586476925286766559005768394_cp
  ! other stuff needed by more than one subroutine
  ! Gauss points: cos(theta), weights, and actual theta
  real(kind=dprec), dimension(:),allocatable :: gauss_z,&
       gauss_w ,gauss_theta
  real(kind=dprec), dimension(:), allocatable :: lfac,ilfac
  ! those are for Legendre polynomials (fac1 & fac2 only for ivec=1)
  ! make those double precision
  !
  real(kind=dprec), dimension(:),allocatable :: plm_f1,plm_f2,&
       plm_fac1,plm_fac2,plm_srt
  ! this is for vector harmonics, only for ivec=1
  real(kind=cp), dimension(:),allocatable :: sin_theta,ell_factor
  ! spacing in longitudes
  double precision  :: dphi
  ! integer (bounds and such)
  integer nlat,nlon,lmsize,lmsize2,nlonm1
  ! logic flags
  logical rick_initialized,rick_computed_legendre,&
       rick_vector_sh_fac_init
end module rick_module

!
! initialize all necessary arrays for Rick type expansion
!
! if ivec == 1, will initialize for velocities/polarizations
!
subroutine rick_f90_init(lmax,ivec,npoints,nplm,tnplm)
  use rick_module
  implicit none
  
  integer, intent(in) :: lmax,ivec
  integer, intent(out) :: npoints,nplm,tnplm
  ! local 
  real(kind=cp) :: xtemp
  integer :: l,old_lmax,old_ivec,old_npoints,old_nplm,old_tnplm
  logical :: was_called
  data was_called /.false./
  save was_called, old_lmax, old_ivec, old_npoints,old_nplm,&
       old_tnplm
  if(.not.was_called)then
     !
     ! test if lmax is 2**n-1
     !
     xtemp = log(dble(lmax+1))/log(2.0_cp)
     if(abs(int(xtemp+.5)-xtemp).gt.1e-7)then
        print *,'rick_init: error: lmax has to be 2**n-1'
        print *,'rick_init: maybe use lmax = ',2**(int(xtemp)+1)-1
        print *,'rick_init: instead of ',lmax
        stop
     endif
     !
     ! number of longitudinal and latitudinal points
     !
     nlat = lmax + 1
     nlon = 2 * nlat
     !
     ! number of points in one layer
     !
     npoints = nlat * nlon 
     old_npoints = npoints
     !
     !
     ! for coordinate computations
     !
     dphi = 2.0_cp * pi / dble(nlon)
     nlonm1 = nlon - 1
     !
     ! size of tighly packed arrays with l,m indices
     lmsize  = (lmax+1)*(lmax+2)/2
     lmsize2 = lmsize * 2          !for A and B
     !
     !
     ! size of the Plm array 
     nplm = lmsize * nlat
     tnplm = nplm * (1+ivec)           ! for all layers
     old_tnplm = tnplm
     old_nplm = nplm
     !
     ! initialize the Gauss points, at which the latitudes are 
     ! evaluated
     !
     allocate(gauss_z(nlat),gauss_w(nlat),gauss_theta(nlat))
     call rick_f90_gauleg(-1.0_dprec,1.0_dprec,gauss_z,gauss_w,nlat)
     !
     ! theta values of the Gauss quadrature points
     !
     gauss_theta = acos(gauss_z)
     !
     ! those will be used by plmbar to store some of the factors
     !
     allocate(plm_f1(lmsize),plm_f2(lmsize),plm_fac1(lmsize),plm_fac2(lmsize))
     allocate(plm_srt(nlon))
     if(ivec.eq.1)then
        !
        ! additional arrays for vector spherical harmonics
        ! (perform the computations in double precision)
        !
        allocate(ell_factor(nlat),sin_theta(nlat))
        ! 1/(l(l+1)) factors
        do l=1,nlat               ! no l=0 term, obviously
           ! go from l=1 to l=lmax+1
           ell_factor(l) = 1.0d0/sqrt(dfloat(l*(l+1)))
        enddo
        sin_theta = sqrt((1.0d0 - gauss_z)*(1.0d0+gauss_z))
        rick_vector_sh_fac_init = .true.
     else
        rick_vector_sh_fac_init = .false.
     endif
     !
     ! logic flags
     !
     rick_computed_legendre = .false.
     rick_initialized = .true.
     was_called = .true.
     old_lmax = lmax
     old_ivec = ivec
  else
     if(lmax .ne. old_lmax)then
        print *,'rick_init: error: was init with lmax',old_lmax,' now',lmax
        stop
     endif
     if(ivec.gt.old_ivec)then
        print *,'rick_init: error: original ivec',old_ivec,' now ',ivec
        stop
     endif

     npoints = old_npoints
     nplm = old_nplm
     tnplm = old_tnplm
  endif
end subroutine rick_f90_init
!
! free all arrays that were allocate for the module   
!
subroutine rick_f90_free_module(ivec)
  use rick_module
  implicit none
  integer, intent(in) :: ivec
  deallocate(gauss_z,gauss_w,gauss_theta);
  if(rick_computed_legendre)then
     ! those are from the Legendre function action
     deallocate(plm_f1);deallocate(plm_f2)
     deallocate(plm_fac1);deallocate(plm_fac2)
     deallocate(plm_srt)
  endif
  if(ivec)then
     deallocate(ell_factor);deallocate(sin_theta)
  endif
end subroutine rick_f90_free_module

!
! Returns arrays X and W with N points and weights for
! Gaussian quadrature over interval X1,X2.
! we call this routine with n = nlat = lmax+1
!
! this is from Numerical Recipes
!       
!     
SUBROUTINE rick_f90_gauleg(x1,x2,x,w,n)
  use rick_module
  implicit none
  !
  integer, intent(in) :: n
  real(kind=dprec), intent(in) :: x1,x2
  real(kind=dprec), intent(out), dimension(n) :: x,w
  !
  ! local variables
  !
  INTEGER :: i,j,m
  real(kind=dprec) :: p1,p2,p3,pp,xl,xm,z,z1
  !
  m=(n+1)/2
  xm=0.5d0*(x2+x1)
  xl=0.5d0*(x2-x1)
  do i=1, m
     z = cos(pi * (i-.25d0)/(n+.5d0))
     z1 = -2.0d0
     do while(abs(z-z1).gt. 5.d-15)     
        p1 = 1.0d0
        p2 = 0.0d0
        do j = 1, n
           p3 = p2
           p2 = p1
           p1 = ((2.0d0*j - 1.0d0)*z*p2-(j- 1.0d0)*p3)/j
        enddo
        pp = n*(z*p1-p2)/(z*z-1.0d0)
        z1 = z
        z = z1-p1/pp
     enddo
     x(i) = xm-xl*z
     x(n+1-i) = xm+xl*z
     w(i) = 2.0d0*xl/((1.0d0-z*z)*pp*pp)
     w(n+1-i)=w(i)
  enddo
  return
END SUBROUTINE rick_f90_gauleg

subroutine rick_f90_plmbar1(p,dp,ivec,lmax,z)
  !
  !     Evaluates normalized associated Legendre function P(l,m), plm,
  !     as function of  z=cos(colatitude); also derivative dP/d(colatitude),
  !     dp, if ivec is set to 1.0_cp
  !
  !     Uses recurrence relation starting with P(l,l) and then 
  !     increasing l keeping m fixed.
  !
  !     p(k) contains p(l,m)
  !     with k=(l+1)*l/2+m+1; i.e. m increments through range 0 
  !     to l before
  !     incrementing l. 
  !
  !     Normalization is:
  !
  !     Integral(P(l,m)*P(l,m)*d(cos(theta)))=2.*(2-delta(0,m)),
  !     where delta(i,j) is the Kronecker delta.
  !
  !     This normalization is incorporated into the
  !     recurrence relation which eliminates overflow. 
  !     Routine is stable in single and double 
  !     precision to
  !     l,m = 511 at least; timing proportional to lmax**2
  !     R.J.O'Connell 7 Sept. 1989; added dp(z) 10 Jan. 1990.
  !
  !     Added precalculation and storage of square roots 
  !     srt(k) 31 Dec 1992
  !
  !     Timing (-O) is 50 ms/call for lmax=127 on l=255 grid in plvel2sh
  !
  !
  use rick_module
  implicit none
  integer, intent(in) :: lmax, ivec
  real(kind=dprec),intent(in) :: z
  real(kind=dprec),intent(inout), dimension (lmsize) :: p, dp

  ! local
  real(kind=dprec) :: plm,pm1,pm2,pmm,sintsq,fnum,fden
  
  integer :: l,m,k,kstart,isize,old_lmax,l2,lstart,mstop,old_ivec
  save old_lmax,old_ivec
  if(.not.rick_initialized)then
     print *,'rick_plmbar1: error: module not initialized'
     stop
  endif
  if (lmax.lt.0.or.abs(z).gt.1.0_cp) then
     print *,'bad arguments'
     stop
  endif
  if(.not.rick_computed_legendre) then
     if(nlon.ne.(lmax+1)*2)then
        print *,'rick_plmbar1: factor mismatch,',lmax,nlon
        stop
     endif
     !
     ! first call, set up some factors. the arrays were allocated in rick_init
     !
     do k=1, nlon
        plm_srt(k) = sqrt(dfloat(k))
     enddo
     ! initialize plm factors
     plm_f1 = 0.0_dprec;plm_fac1 = 0.0_dprec
     plm_f2 = 0.0_dprec;plm_fac2 = 0.0_dprec
     !     --case for m > 0
     kstart = 1
     do m =1, lmax
        !     --case for P(m,m) 
        kstart = kstart+m+1
        if (m .ne. lmax) then
           !     --case for P(m+1,m)
           k = kstart+m+1
           !     --case for P(l,m) with l > m+1
           if (m .lt. (lmax-1)) then
              do l = m+2, lmax
                 k = k+l
                 l2 = l * 2
                 plm_f1(k) = plm_srt(l2+1) * plm_srt(l2-1)/&
                      (plm_srt(l+m) * plm_srt(l-m))
                 plm_f2(k)=(plm_srt(l2+1) * plm_srt(l-m-1)*&
                      plm_srt(l+m-1))/&
                      (plm_srt(l2-3) * plm_srt(l+m) * plm_srt(l-m))
              enddo
           endif
        endif
     enddo
     if(ivec.eq.1)then
        !
        ! for derivative of Plm with resp. to theta
        ! (if we forget to call Plm with ivec=1, but use
        ! ivec=1 later, we will notice as all should be zero)
        !
        k=3
        do l=2, lmax
           k=k+1
           mstop = l - 1
           do m=1, mstop
              k=k+1
              plm_fac1(k) = plm_srt(l-m) * plm_srt(l+m+1)
              plm_fac2(k) = plm_srt(l+m) * plm_srt(l-m+1)
              if(m.eq.1)then
                 plm_fac2(k) = plm_fac2(k) * plm_srt(2)
              endif
           enddo
           k=k+1
        enddo
     endif
     old_lmax = lmax
     old_ivec = ivec
     rick_computed_legendre = .true.
  else
     ! test if lmax has changed
     if(lmax.ne.old_lmax)then
        print *,'rick_plmbar1: error: factors were computed for lmax',old_lmax
        print *,'rick_plmbar1: error: now, lmax is ',lmax
        stop
     endif
     if(ivec.gt.old_ivec)then
        print *,'rick_plmbar1: error: init with',old_ivec, 'now ivec',ivec
        stop
     endif
  endif

  !     --start calculation of Plm etc.
  !     --case for P(l,0) 
  
  pm2  = 1.0_dprec
  p(1) = 1.0_dprec                   ! (0,0)
  if(ivec.eq.1)then
     dp(1) = 0.0_dprec! else, don't refer to this array
  endif
  if (lmax .eq. 0) return
  pm1  = z
  p(2) = plm_srt(3) * pm1             ! (1,0)
  k=2
  do  l = 2, lmax               ! now all (l,0)
     k  = k + l
     l2 = l * 2
     plm=(dfloat(l2-1)*z*pm1 - dfloat(l-1)*pm2)/dfloat(l)
     p(k)=plm_srt(l2+1)*plm      
     pm2=pm1
     pm1=plm
  enddo
  !       --case for m > 0
  pmm = 1.0_dprec
  sintsq = (1.0_dprec - z)*(1.0_dprec + z)
  fnum = -1.0_dprec
  fden =  0.0_dprec
  kstart = 1
  do  m = 1, lmax
     !     --case for P(m,m) 
     kstart = kstart+m+1
     fnum = fnum + 2.0_dprec
     fden = fden + 2.0_dprec
     pmm = pmm*sintsq*fnum/fden
     pm2 = sqrt(dfloat(4*m+2)*pmm)
     p(kstart) = pm2   
     if (m .ne. lmax) then
        !     --case for P(m+1,m)
        pm1=z*plm_srt(2*m+3)*pm2
        k = kstart + m + 1
        p(k) = pm1
        !     --case for P(l,m) with l > m+1
        if (m .lt. (lmax-1)) then
           do  l = m+2, lmax
              k = k+l
              plm = z*plm_f1(k)*pm1 - plm_f2(k)*pm2
              p(k) = plm
              pm2 = pm1
              pm1 = plm
           enddo
        endif
     endif
  enddo
  if(ivec.eq.1)then
     ! 
     ! derivatives
     !     ---derivatives of P(z) wrt theta, where z=cos(theta)
     !     
     dp(2)=-p(3)
     dp(3)= p(2)
     k=3
     do l=2,lmax
        k=k+1
        !     treat m=0 and m=l separately
        dp(k) =  -plm_srt(l)*plm_srt(l+1)/plm_srt(2)* p(k+1)
        dp(k+l) = plm_srt(l)/plm_srt(2)* p(k+l-1)

        mstop = l-1
        do  m=1, mstop
           k=k+1
           dp(k)=0.5_dprec*(plm_fac2(k)*p(k-1) - &
                plm_fac1(k)*p(k+1) )

        enddo
        k=k+1
     enddo
  endif
  return
end subroutine rick_f90_plmbar1


subroutine rick_f90_shd2c(rdatax,rdatay,lmax,ivec,cslm,dslm)
  !
  !	Calculates spherical harmonic coefficients cslm(l,m) of
  !	a scalar (ivec = 0) or vector (ivec=1) function on a sphere. 
  !
  !     if ivec == 1, then cslm will be poloidal, and dslm toroidal
  !                   rdata will be nlon*nlat
  !
  !     Coefficients stored with
  !	cosine term and sine term alternating, starting at l=0
  !	and m=0 with m ranging from 0 to l before l is incremented.
  !
  !     Harmonics normalized with mean square = 1.0. Expansion is:
  !     SUM( Plm(cos(theta))*(cp(l,m)*cos(m*phi)+sp(l,m)*sin(m*phi))
  !
  !     Uses FFT in longitude and gaussian integration of spectral
  !     coefficients (for fixed m) times associated Legendre
  !     functions (of same order m) to get expansion coefficients.
  !
  !     Uses Num Rec routine for Gaussian points and weights, which should
  !     be initialized first by a call to rick_init. 
  !
  !     Uses recursive routine to generate associated Legendre functions.
  !     
  !     INPUT:
  !
  !     lmax    maximum possible degree of expansion (2**n-1). There
  !             are nlon = 2*lmax+2 data points in longitude for each latitude
  !             which has nlat = lmax+1 points
  !
  !     rdatax,rdatay   data((nlon * nlat)) arrays for theta and phi
  !                     components
  !
  !     ic_pd: should be either 0.0_cp or 1.0_cp
  !            if set to 1.0_cp, will expand velocities instead
  !            in this case, data should hold the theta and phi components
  !            of a vector field and cslm will hold the poloidal and toroidal
  !            on output components, respectively
  !
  !     OUTPUT:
  !
  !     cslm,dslm    coefficients, (2*lmsize)
  !
  !     dslm and rdatay will only be referenced when ivec = 1
  !
  !
  use rick_module               ! will have to be initialized first
  !
  implicit none
  integer, intent(in) :: lmax,ivec
  real(kind=sprec), intent(in),   dimension (nlon*nlat) :: rdatax,rdatay
  real(kind=cp), intent(out),  dimension (lmsize2) :: cslm
  real(kind=cp), intent(out),  dimension (lmsize2*ivec) :: dslm
  ! local
  real(kind=dprec), dimension (lmsize*nlat) :: plm,dplm
  integer :: i,os
  ! check
  if(.not. rick_initialized)then
     write(6,*)'rick_shd2c: error: initialize first'
     stop
  endif
  if(lmax .ne. nlat-1)then
     print *,'rick_shd2c: error: lmax mismatch: nlat/lmax'
     print *,nlat,lmax
     stop
  endif
  ! compute the Plm first
  call rick_f90_compute_allplm(lmax,ivec,plm,dplm)
  ! call the precomputed version
  call rick_f90_shd2c_pre(rdatax,rdatay,lmax,plm,&
       dplm,ivec,cslm,dslm)
  return
end subroutine rick_f90_shd2c
!
! the actual routine to go from spatial to spectral, 
! for comments, see above
!
subroutine rick_f90_shd2c_pre(rdatax,rdatay,lmax,plm,dplm,ivec,&
     cslm,dslm)
  use rick_module
  implicit none
  ! 
  integer, intent(in) :: lmax,ivec
  real(kind=sprec), intent(in),   &
       dimension (nlon*nlat) :: rdatax,rdatay
  real(kind=dprec), intent(in),   &
       dimension (lmsize*nlat) :: plm,dplm
  real(kind=cp), intent(out),  dimension (lmsize2) :: cslm
  real(kind=cp), intent(out),  dimension (lmsize2*ivec) :: dslm
  ! local
  real(kind=cp), dimension(nlon+2) :: valuex
  real(kind=cp), dimension((nlon+2)*ivec) :: valuey !zero for ivec=0
  real(kind=cp) :: dfact,dpdt,dpdp
  !
  integer :: lmaxp1,lmaxp1t2,i,j,l,m,ios1,m2,j2,oplm
  ! check
  if(.not. rick_initialized)then
     write(6,*)'rick_shd2c_pre: error: initialize first'
     stop
  endif
  ! check some more and compute bounds
  lmaxp1 = lmax + 1
  lmaxp1t2 = lmaxp1 * 2
  if((lmaxp1.ne.nlat).or.(nlon.ne.lmaxp1t2).or.&
       ((lmax+1)*(lmax+2)/2.ne.lmsize))then
     print *,'rick_shd2c_pre: dimension error, lmax ',lmax
     print *,'rick_shd2c_pre: nlon ',nlon,' nlat ',nlat
     print *,'rick_shd2c_pre: lmsize',lmsize
     stop
  endif
  !
  ! initialize the coefficients
  !
  cslm = 0.0_cp
  if(ivec.eq.1)then
     dslm = 0.0_cp
     if(.not.rick_vector_sh_fac_init)then
        print *,'rick_shd2c_pre: error: vector harmonics factors not initialized'
        stop
     endif
  endif
  do i=1,nlat
     !
     ! loop through latitudes
     !
     ios1 = (i-1)*nlon          ! offset for data array
     oplm = (i-1)*lmsize        ! offset for Plm array
     !
     if(ivec.eq.0)then 
        !
        ! scalar expansion
        !
        do j=1,nlon      ! can't vectorize, because dimension don't match
           valuex(j) = rdatax(ios1 + j) 
        end do
        !
        ! compute the FFT 
        !
        call rick_f90_realft(valuex,nlat,+1) 
        call rick_f90_ab2cs(valuex,nlon)
        ! sum up for integration
        l = 0;m = -1
        do j=1, lmsize
           m = m + 1
           if( m .gt. l ) then
              l=l+1;m=0
           endif
           ! we incorporate the Gauss integration weight and Plm factors here
           if (m .eq. 0) then
              dfact = (gauss_w(i) * plm(oplm+j))/2.0_dprec
           else
              dfact = (gauss_w(i) * plm(oplm+j))/4.0_dprec
           endif
           m2 = m * 2;j2 = j * 2
           cslm(j2-1)= cslm(j2 - 1) + valuex(m2+1) * dfact ! A coefficient
           cslm(j2)  = cslm(j2)     + valuex(m2+2) * dfact ! B coefficient
        end do
     else
        !
        ! vector field expansion
        !
        do j=1,nlon      ! can't vectorize, because dimension don't match
           valuex(j) = rdatax(ios1 + j) ! theta
           valuey(j) = rdatay(ios1 + j) ! phi
        end do
        ! perform the FFTs on both components
        call rick_f90_realft(valuex,nlat,+1) 
        call rick_f90_realft(valuey,nlat,+1)
        call rick_f90_ab2cs(valuex,nlon)
        call rick_f90_ab2cs(valuey,nlon)
        !
        l=1;m=-1                ! there's no l=0 term
        !
        do j = 2, lmsize
           m = m + 1
           if( m .gt. l ) then
              l=l+1;m=0
           endif
           if (m .eq. 0) then ! ell_factor is 1/sqrt(l(l+1))
              dfact = gauss_w(i)*ell_factor(l)/2._dprec
           else
              dfact = gauss_w(i)*ell_factor(l)/4._dprec
           endif
           !
           ! some more factors
           !
           ! d_theta(P_lm) factor
           dpdt = dplm(oplm+j) 
           ! d_phi (P_lm) factor
           dpdp = dfloat(m) * plm(oplm+j)/sin_theta(i) 
           !
           m2 = m * 2;j2 = j * 2
           !           print *,m,l,dpdt*dfact,dpdp*dfact
           cslm(j2-1) = cslm(j2-1) + & ! poloidal A
                (dpdt * valuex(m2+1) - dpdp * valuey(m2+2))*dfact
           cslm(j2)   = cslm(j2) + & !   poloidal B
                (dpdt * valuex(m2+2) + dpdp * valuey(m2+1))*dfact
           dslm(j2-1) = dslm(j2-1)+ & ! toroidal A 
                (-dpdp * valuex(m2+2) - dpdt * valuey(m2+1))*dfact
           dslm(j2)   = dslm(j2)+& ! toroidal B 
                (+dpdp * valuex(m2+1) - dpdt * valuey(m2+2))*dfact

        end do
     endif                      ! end vector field
  end do                        ! end latitude loop
  !print *,cslm(1:4)
  !if(ivec.eq.1)then
  !   print *,cslm(lmsize2+1),cslm(lmsize2+2),cslm(lmsize2+3),cslm(lmsize2+4)
  !Endif
  return
  
end subroutine rick_f90_shd2c_pre

subroutine rick_f90_shc2d(cslm,dslm,lmax,ivec,rdatax,rdatay)
  !
  ! Transforms spherical harmonic coefficients of a scalar (ivec = 0)
  ! or a vector (ivec=1) field
  ! to data points on a grid. Reverse transform of subroutine
  ! shd2c.f so long as the same degree is used and the points
  ! in latitude are Gaussian integration points. Maximum degree
  ! must be 2**n-1 in order for the FFT to work; this determines
  ! the grid spacing in longitude. nlat = lmax+1, nlon = 2*nlat
  !
  ! INPUT:
  !
  !		lmax	spherical harmonic degree used in expansion
  !		cslm	(lmsize2) spherical harmonic coefficients
  !             dslm    (lmsize2)
  !                     if ivec, will hold poloidal and toroidal coeff,else
  !                     dslm will not be referenced
  !
  !             ivec   0: scalar 1: vector field
  ! OUTPUT:
  !		rdatax	(nlon * nlat) values on grid points
  !             rdatay  (nlon * nlat). x and y will be theta and phi for 
  !                     ivec = 1, else rdatay will not be referenced
  !
  ! if Ivec is set, assume velocities to be expanded instead
  !
  use rick_module
  implicit none
  integer, intent(in) :: lmax,ivec
  real(kind=cp), intent(in),  dimension (lmsize2) :: cslm
  real(kind=cp), intent(in),  dimension (lmsize2*ivec) :: dslm
  real(kind=sprec), intent(out), dimension (nlat*nlon) :: rdatax
  real(kind=sprec), intent(out), dimension (nlat*nlon*ivec) :: rdatay
  ! local
  real(kind=dprec),  dimension (nlat*lmsize) :: plm,dplm
  integer :: i,os
  if(.not. rick_initialized)then
     write(6,*)'rick_shc2d: error: initialize first'
     stop
  endif
  ! compute the Plm first
  if(lmax .ne. nlat-1)then
     print *,'rick_shc2d: error: lmax mismatch: nlat/lmax'
     print *,nlat,lmax
     stop
  endif
  ! compute the Plm first
  call rick_f90_compute_allplm(lmax,ivec,plm,dplm)
  ! call the precomputed subroutine
  call rick_f90_shc2d_pre(cslm,dslm,lmax,plm,dplm,&
       ivec,rdatax,rdatay)
end subroutine rick_f90_shc2d
!
! the actual routine to go from spectral to spatial
!
subroutine rick_f90_shc2d_pre(cslm,dslm,lmax,plm,dplm,ivec,&
     rdatax,rdatay)
  !
! Legendre functions are precomputed
!
  use rick_module
  implicit none
  integer, intent(in) :: lmax,ivec
  real(kind=cp), intent(in),  dimension (lmsize2) :: cslm
  real(kind=cp), intent(in),  dimension (lmsize2*ivec) :: dslm
  real(kind=dprec), intent(in),  dimension (nlat*lmsize) :: plm,dplm
  real(kind=sprec), intent(out), dimension (nlat*nlon) :: rdatax
  real(kind=sprec), intent(out), dimension (nlat*nlon*ivec) :: rdatay
  ! local
  real(kind=cp), dimension((nlon+2)) :: valuex
  real(kind=cp), dimension((nlon+2)*ivec) :: valuey
  real(kind=cp) :: dpdt,dpdp
  integer :: i,j,m,m2,j2,ios1,l,lmaxp1,lmaxp1t2,oplm
  if(.not. rick_initialized)then
     write(6,*)'rick_shc2d_pre: error: initialize modules first'
     stop
  endif
  ! check bounds 
  lmaxp1 = lmax + 1                ! this is nlat
  lmaxp1t2 = 2 * lmaxp1               ! this is nlon
  if((nlat.ne.lmaxp1).or.(nlon.ne.lmaxp1t2))then
     print *,'rick_shc2d_pre: dimension mismatch:',lmax,nlon,nlat
     stop
  endif
  if(ivec.eq.1)then
     if(.not.rick_vector_sh_fac_init)then
        print *,'rick_shc2d_pre: error: vector harmonics factors not initialized'
        stop
     endif
  endif
  do i=1,nlat                   ! loop through latitudes
     oplm = (i-1)*lmsize        ! offset for Plm array
     ios1 = (i-1) * nlon        ! offset for data array
     !
     if(ivec.eq.0)then          
        !
        ! scalar
        !
        valuex = 0.0_cp               ! init with 0.0_cpes
        l=0; m=-1
        do j=1, lmsize             ! loop through l,m
           m = m + 1
           if (m .gt. l) then
              m=0; l=l+1
           endif
           m2 = 2*m;j2 = 2*j
           ! add up contributions from all l,m 
           valuex(m2+1) = valuex(m2+1) + & ! cos term
                plm(oplm+j) * cslm(j2-1) ! A coeff
           valuex(m2+2) = valuex(m2+2) + & ! sin term
                plm(oplm+j) * cslm(j2)   ! B coeff
        enddo
        ! compute inverse FFT 
        call rick_f90_cs2ab(valuex,nlon)
        ! inverse FFT
        call rick_f90_realft(valuex,nlat,-1)
        !
        do j=1, nlon            ! can't vectorize
           rdatax(ios1 + j) = valuex(j)/dfloat(nlat)
        enddo
     else
        !
        ! vector harmonics
        !
        valuex = 0.0_cp
        valuey = 0.0_cp
        l=1; m=-1               ! start at l = 1
        do j=2, lmsize             ! loop through l,m
           m = m + 1
           if (m .gt. l) then
              m=0; l=l+1
           endif
           m2  = 2*m;j2 = 2*j
           ! derivative factors
           dpdt = dplm(oplm+j) * ell_factor(l) ! d_theta(P_lm) factor
           dpdp = dfloat(m) * plm(oplm+j)/sin_theta(i) * ell_factor(l) ! d_phi (P_lm) factor
           ! add up contributions from all l,m 
           ! make life a little easier
           !
           ! u_theta
           !
           valuex(m2+1) = valuex(m2+1) + & ! cos term
                cslm(j2-1) * dpdt  + dslm(j2)   * dpdp
           valuex(m2+2) = valuex(m2+2) + & ! sin term
                cslm(j2)   * dpdt  - dslm(j2-1) * dpdp
           !
           ! u_phi
           !
           valuey(m2+1) = valuey(m2+1) + & ! cos term
                  cslm(j2)   * dpdp - dslm(j2-1) * dpdt
           valuey(m2+2) = valuey(m2+2) + & ! sin term
                - cslm(j2-1) * dpdp - dslm(j2)   * dpdt
        enddo
        ! do inverse FFTs
        call rick_f90_cs2ab(valuex,nlon)
        call rick_f90_cs2ab(valuey,nlon)
        call rick_f90_realft(valuex,nlat,-1)
        call rick_f90_realft(valuey,nlat,-1)
        ! assign to output array
        do j=1, nlon            ! can't vectorize, since there's an offset of 2
           rdatax(ios1 + j) = valuex(j)/dfloat(nlat)
           rdatay(ios1 + j) = valuey(j)/dfloat(nlat)
        enddo
     endif
  enddo
  return
end subroutine rick_f90_shc2d_pre

!
! detemine the colatidude and longitude of pixel index
! where index goes from 0 ... nlon * nlat-1
!
subroutine rick_f90_pix2ang(index, lmax, theta, phi)

  use rick_module
  implicit none
  ! input / output
  integer, intent(in) :: lmax, index
  double precision, intent(out) :: theta, phi
  ! local
  integer :: i,j
  if(.not.rick_initialized)then
     print *,'rick_pix2ang: error: not initialized'
     stop
  endif

  j = index
  i=1
  do while(j.gt.nlonm1)
     j = j - nlon
     i=i+1
  enddo
  theta = gauss_theta(i)
  phi = dphi * dble(j)
end subroutine rick_f90_pix2ang

!
! compute Legendre function (l,m) evaluated on nlat points 
! in latitude and their derviatives with respect to theta, if 
! ivec is set to 1.0_cp
!
subroutine rick_f90_compute_allplm(lmax,ivec,plm,dplm)
  use rick_module
  implicit none
  integer, intent(in) :: lmax,ivec
  real(kind=dprec), intent(out), dimension(lmsize*nlat) :: plm, dplm
  ! local
  integer :: i,os
  if(lmax .ne. nlat-1)then
     print *,'rick_compute_allplm: error: lmax mismatch: nlat/lmax'
     print *,nlat,lmax
     stop
  endif
  os = 1
  do i=1,nlat
     call rick_f90_plmbar1(plm(os),dplm(os),ivec,lmax,gauss_z(i))
     os = os + lmsize
  enddo

end subroutine rick_f90_compute_allplm
