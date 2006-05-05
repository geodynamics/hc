! 
!
! these subroutines deal with the kinematic solution of a Hager &
! O'Connell flow code. they are based on Brad's original code, and later
! Bernhard Steinberger's modifitions
!
! not all the comments reflect all the changes, beware
!
!
!
! $Id: kinsub3.F,v 1.15 2002/03/15 03:42:23 becker Exp $
!
!
!     ****************************************************************
!     * THIS IS THE MAIN PROGRAM FOR THE COMPONENT OF FLOW WITHOUT   *
!     * DENSITY CONTRASTS.  IT USES SEVERAL INPUT/OUTPUT SUBROUTINES *
!     * AND FUNCTIONS TO OBTAIN, CORRECT AND VERIFY A MODEL FROM THE *
!     * USER.  THE FINAL VERSION OF EACH MODEL IS STORED IN A FILE   *
!     * BEFORE THE PROGRAM EXECUTES POLSOL AND TORSOL TO OBTAIN THE  *
!     * POLOIDAL AND TOROIDAL COMPONENTS, RESPECTIVELY, OF THE       *
!     * EQUATIONS OF MOTION.                                         *
!     ****************************************************************
!     Modified such that only toroidal component is calculated
!     Poloidal component is included in densub.f
!
!     ARRAYS:  R: OUTPUT RADII,
!        VISC,RVISC: normalized viscosities and their radii,
!        VISC_LOCAL: local normalized viscosities
!     OTHER VAR:
!     DOUBLE PRECISION:  vis_norm: NORMALIZING VISCOSITY.
!     INTEGERS:  LMAX: MAXIMUM DEGREES,
!           NRADP2: NUMBER OF OUTPUT RADII,
!           n_vislayer: NUMBER OF VISCOSITIES.
!     SUBROUTINES:
!        TORSOL:  CALCULATES THE TOROIDAL VECTOR SOLUTION.

subroutine torsol(nrad,nv,lmax,r,rv,visc,tvec)
  !    
  !     ****************************************************************
  !     * EVALUATES AND PROPAGATES THE TWO TOROIDAL COMPONENTS IN THE  *
  !     * EQUATIONS OF MOTION, AND NORMALIZES THESE SUCH THAT THE      *
  !     * FIRST ELEMENT AT THE SURFACE IS 1.0.                         *
  !     ****************************************************************
  !
  ! in 
  real(kind=cp), intent(in), dimension(nrad+2) :: r
  real(kind=cp), intent(in), dimension(nv+1) :: rv,visc
  integer, intent(in) :: nrad,nv,lmax
  ! out
  real(kind=cp), intent(out), dimension(0:lmax,nrad+2,2) :: tvec
  
  ! local
  real(kind=cp) :: coef,diflog,el,elp2,elm1,exp1,exp2,vecnor,hold, &
       p11,p12,p21,p22,rlast,rnext,tvec1,tvec2
  integer :: l,jvisp1,jvis,i,j,k,nvp1,nradp2
  !
  logical :: qvis,ifp_solvec
  !
  !     PASSED PARAMETERS:  NRADP2: NUMBER OF OUTPUT RADII,
  !        NV: NUMBER OF VISCOSITIES, nvp1 = nv+1
  !        LMAX: MAXIMUM DEGREES.
  !     ARRAYS:  R: OUTPUT RADII,
  !        RV: VISCOSITY RADII,
  !        TVEC: TOROIDAL VECTORS,
  !        VISC: normalized VISCOSITIES.
  !     OTHER VAR:  EXP1,EXP2: EXPONENTIAL FACTORS IN PROPAGATOR,
  !        COEF,ELP2,ELM1: PARAMETERS IN PROPAGATOR,
  !        DIFLOG: DIFFERENCE IN LOGS OF RADII,
  !        EL,L: DEGREE,
  !        VECNOR: NORMALIZES TVEC TO TVEC(N,1),
  !        HOLD: TEMPORARY VAR.,
  !        P11,P12,P21,P22: ELEMENTS OF THE PROPAGATOR MATRIX CORRES-
  !        PONDING TO P(1,1),P(1,2),P(2,1),P(2,2) RESPECTIVELY,
  !        RLAST,RNEXT: RADII FOR PROPAGATOR,
  !        TVEC1,TVEC2: VECTOR COMPONENTS.
  !
  nvp1 = nv + 1
  nradp2 = nrad + 2

  rv(nvp1) = 1.1d0              !why is this 1.1?
  visc(nvp1) = visc(nvp1-1)
!
!     (PREVENTS THE REQUESTING OF NON-EXISTANT VALUES)
!     
!     FOR EACH DEGREE (L) CALCULATE, NORMALIZE AND OUTPUT SOLUTION
!
  tvec(0,:,:) = 0.0_cp          
  
  do l = 1, lmax
     el = dfloat(l)
     !     
     !     SET THE PARAMETERS
     !     
     elp2 = el + 2.0_cp
     elm1 = el - 1.0_cp
     coef = 1.0_cp / (2.0_cp * el + 1.0_cp)
     !
     !     INITIALIZE THE PROPAGATION AT THE CORE
     !
     jvisp1 = 2
     jvis = 1
     rlast = r(1)
     tvec1 = 1.0_cp
     tvec2 = 0.0_cp
     tvec(l,1,1) = tvec1
     tvec(l,1,2) = tvec2
!
!     FIND THE TWO TOROIDAL COMPONENTS AT EACH RADIUS
!     start radius loop
     !
     do i = 2, nradp2
        qvis = .FALSE.
        !
        !     TEST FOR CHANGE IN VISCOSITY IN NEXT LAYER
        !
40      if (rv(jvisp1).gt.r(i)) qvis = .TRUE.
        rnext = rv(jvisp1)
        !
        !     IF NO VISC. CHANGE BEFORE NEXT OUTPUT RADIUS, PROPAGATE DIRECTLY
        !
        if (qvis) rnext = r(i)
        diflog = log(rnext / rlast)
        rlast = rnext
        exp1 = exp(el * diflog)
        exp2 = exp(-(el + 1.0_cp) * diflog)
        !
        !     PROPAGATOR SET UP LINEARLY TO AVOID EXCESS MULTIPLICATIONS
        !
        p11 = (elp2 * exp1 + elm1 * exp2)
        p12 = (exp1 - exp2)   / visc(jvis)
        p21 = elp2 * elm1 * visc(jvis) * (exp1 - exp2)
        p22 = (elm1 * exp1 + elp2 * exp2)
        hold = tvec1
        !
        !     PROPAGATE LAST VECTOR TO GET NEW VECTOR
        !
        tvec1 = coef * (p11 * hold + p12 * tvec2)
        tvec2 = coef * (p21 * hold + p22 * tvec2)
        IF (QVIS) GO TO 50
        jvis = jvisp1
        jvisp1 = jvis + 1
        GO TO 40
!
!     IF AT REQUIRED OUTPUT RADIUS, STORE VECTOR AND CONTINUE
        !
50      tvec(l,i,1) = tvec1
        tvec(l,i,2) = tvec2
            
        !     end of radius loop
     enddo
     ! 
     !     SET TVEC(l,NRADP2,1) = 1.0 AND NORMALIZE ALL VECTORS TO THIS
     !
     vecnor = 1.0_cp / tvec(l,nradp2,1)
     do k = 1, 2
        do j = 1, nradp2
           tvec(l,j,k) = vecnor * tvec(j,k)
        enddo
     enddo
     !     end of l loop         
  enddo
END subroutine torsol

