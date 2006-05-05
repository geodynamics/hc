!
! fourier transform routines as used by rick_sh routines
! based on Rick O'Connell's subroutines, which are modified 
! Numerical Recipes
!
! $Id: rick_fft.f90,v 1.4 2006/01/22 02:11:12 becker Exp becker $
!

subroutine rick_f90_cs2ab(rdata,n)
  !
  !  Transforms spectral coefficients from cos-sin series to
  !    complex discrete fourier series. Function is real, and
  !    transformed by realft(rdata,n/2,1). Number of data points
  !    is n. Does not recover real component for frequency n/2.

  use rick_module
  implicit none
  integer :: n,i
  real(kind=cp), intent(inout),dimension(n) :: rdata
  real(kind=cp) :: en

  en = dfloat(n)
  rdata(1) = rdata(1) * en
  do  i=3,n
     rdata(i) = rdata(i)*en/2.0_cp
  enddo
  return
end subroutine rick_f90_cs2ab

subroutine rick_f90_ab2cs(rdata,n)

!  Changes coefficients of complex spectrum of a real function
!   transformed by realft.f to real coefficients of a series
!   of C*cos(m*x)+S*sin(mx). Coefficients are ordered as
!   C(0),S(0),C(1),S(1),C(2),...,C(n/2-1),S(n/2-1). This loses
!   the real part of spectrum for frequency n/2.
!   The number of data points is n, The call to realft is
!    call realft(rdata,n/2,1)
!
  use rick_module
  implicit none
  integer n,i
  real(kind=cp), intent(inout), dimension(n) :: rdata
  real(kind=cp) :: en

  en=dfloat(n)
  rdata(1)=rdata(1)/en
  rdata(2)= 0.0_cp
  do i=3,n
     rdata(i)=rdata(i)/en*2.0_cp
  enddo
  return
end subroutine rick_f90_ab2cs


subroutine rick_f90_realft(rdata,n,isign)
!
!	Calculates the fourier transform of 2*N real data points.
!	Replaces data with the positive frequency half of the
!	complex fourier transform. The real parts of the first
!	and last frequency components are returned in data(1)
!	and data(2) (i.e. for frequencies of zero and N/2). The
!	other spectral components are given as complex pairs
!	in data(3),data(4) etc. The inverse transform is obtained
!	with ISIGN=-1, and dividing the data or result by N.
!	Calls routine four1(data,n,isign) for FFT.
!
  use rick_module
  implicit none
  integer, intent(in) :: n,isign
  real(kind=cp), intent(inout), dimension(2*n+2) :: rdata
  ! local
  real(kind=cp) :: c1,c2,h1r,h1i,h2r,h2i
  double precision :: theta,wi,wpi,wpr,wr,wtemp
  integer :: i,n2p3,ilim,i1,i2,i3,i4,n2

  theta = pi / dble(n)          !this is different from Numrec

  wr = 1.0d0
  wi = 0.0d0
  c1 = 0.5_cp

  n2 =   2*n
  n2p3 = n2 + 3

  if (isign.eq.1) then
     c2=-0.5_cp
     call rick_f90_four1(rdata,n,+1)
     rdata(n2+1)=rdata(1)
     rdata(n2+2)=rdata(2)
  else
     c2=0.5_cp
     theta=-theta
     rdata(n2+1)=rdata(2)
     rdata(n2+2)=0.0_cp
     rdata(2)=0.0_cp
  endif
  wpr=-2.0d0 * sin(0.5d0 * theta)**2
  wpi=sin(theta)

  ilim = n/2+1
  do i=1,ilim
     i1=2*i-1                   ! 1,3,5,...
     i2=i1+1
     i3=n2p3-i2
     i4=i3+1
     h1r=c1*(rdata(i1)+rdata(i3))
     h1i=c1*(rdata(i2)-rdata(i4))
     h2r=-c2*(rdata(i2)+rdata(i4))
     h2i= c2*(rdata(i1)-rdata(i3))
     rdata(i1)= h1r+wr*h2r-wi*h2i
     rdata(i2)= h1i+wr*h2i+wi*h2r
     rdata(i3)= h1r-wr*h2r+wi*h2i
     rdata(i4)=-h1i+wr*h2i+wi*h2r
     wtemp=wr
     wr=wr*wpr-wi*wpi+wr
     wi=wi*wpr+wtemp*wpi+wi
  enddo
  if (isign.eq.1) then          
     rdata(2)=rdata(n2+1)
  else
     call rick_f90_four1(rdata,n,-1)
  endif
  return
end subroutine rick_f90_realft


subroutine rick_f90_four1(rdata,nn,isign)
!
!	FFT routine from Numerical Recipes. Replaces data by
!	its discrete fourier transform if isign=1, or by
!	NN times its inverse transform if isign=-1. Array
!	data is made up of NN complex numbers (2*NN pairs)
!	and NN must be a power of 2. Spectral components
!	are complex, and ordered from frequency zero to
!	+-NN/2 to -1 in the standard fashion.
!
  use rick_module
  implicit none
  integer, intent(in) :: nn,isign
  real(kind=cp), intent(inout), dimension(2*nn+2) :: rdata
! local
  real(kind=cp) :: tempr,tempi
! this should be double precision locally, regardless
  double precision ::  wr,wi,wpr,wpi,wtemp,theta

  integer :: n,m,i,j,mmax,istep

  n=2*nn 
  j=1    
  do i=1,n,2
     if(j.gt.i)then
        tempr=rdata(j)
        tempi=rdata(j+1)
        rdata(j)=rdata(i)
        rdata(j+1)=rdata(i+1)
        rdata(i)=tempr
        rdata(i+1)=tempi
     endif
     m=n/2
1    if ((m.ge.2).and.(j.gt.m)) then
        j=j-m
        m=m/2
        go to 1
     endif
     j=j+m
  enddo
  mmax=2 
2 if (n.gt.mmax) then
     istep=2*mmax
     theta = 6.28318530717959d0/(isign*mmax)
     wpr=-2.0d0 * sin(0.5d0 * theta)**2
     wpi=sin(theta)
     wr = 1.0d0
     wi = 0.0d0
     
     do m=1,mmax,2
        do i=m,n,istep
           j=i+mmax
           tempr=wr*rdata(j)-wi*rdata(j+1)
           tempi=wr*rdata(j+1)+wi*rdata(j)
           rdata(j)=rdata(i)-tempr
           rdata(j+1)=rdata(i+1)-tempi
           rdata(i)=rdata(i)+tempr
           rdata(i+1)=rdata(i+1)+tempi
        enddo
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
     enddo
     mmax=istep
     go to 2
  endif
  return 
end subroutine rick_f90_four1

