INTEGER, PARAMETER :: dprec = selected_real_kind(12,200)
INTEGER, PARAMETER :: sprec = selected_real_kind(5,30)
#ifdef RICK_DOUBLE_PRECISION

! double
  integer, parameter :: cp = dprec
  real(kind=cp), parameter,public :: eps = 2.0e-15

#else

! single

  integer, parameter :: cp = sprec
  
  real(kind=cp), parameter,public :: eps = 1e-7
  


#endif
