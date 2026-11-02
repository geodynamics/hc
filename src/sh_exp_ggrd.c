/* 


lower level routines to handle spherical harmonics expansions that
require ggrd compatibility

the higher level routines that operate on single (or triple) 
sets of expansions are in sh_model.c

$Id: sh_exp.c,v 1.15 2006/03/20 05:32:48 becker Exp becker $

*/
#include "hc.h"
#ifndef NO_GMT
#include "hc_ggrd.h"

/* see sh_exp function */
void sh_read_spatial_data_from_grd(struct sh_lms *exp, struct ggrd_gt *ggrd,
				   my_boolean use_3d,int shps, HC_PREC *data, 
				   HC_PREC *z)
{
  int j,k;
  HC_PREC xp[3];
  double dvalue;
  for(j=0;j < exp->npoints;j++){
    /* 
       get expected coordinates to check if the input is OK
    */
    switch(exp->type){
#ifdef HC_USE_HEALPIX
    case SH_HEALPIX:
      switch(exp->heal.ordering){
      case SH_HEALPIX_RING:  
	pix2ang_ring((long)exp->heal.nside,(long)j,
		     (xp+HC_THETA),(xp+HC_PHI));
	break;
      case SH_HEALPIX_NEST:  
	pix2ang_nest((long)exp->heal.nside,(long)j,
		     (xp+HC_THETA),(xp+HC_PHI));
	break;
      default:
	fprintf(stderr,"sh_read_spatial_data_from_grd: error: ordering %i undefined\n",
		exp->heal.ordering);
	exit(-1);
	break;
      }
      break;			/* end Healpix branch */
#endif
    case SH_RICK:
      /* for Rick's type routine */
#ifdef NO_RICK_FORTRAN
      rick_pix2ang(j,exp->lmax,(xp+HC_THETA),(xp+HC_PHI),
		   &exp->rick);
#else
      rick_f90_pix2ang(&j,&exp->lmax,(SH_RICK_PREC)(xp+HC_THETA),
		       (SH_RICK_PREC)(xp+HC_PHI));
#endif
      break;
#ifdef HC_USE_SPHEREPACK
  case SH_SPHEREPACK_GAUSS:
  case SH_SPHEREPACK_EVEN:
    break;
#endif
    default:
      sh_exp_type_error("sh_read_model_spatial_data_from_grd",exp);
      break;
    }	/* end type branch */
    /* 
       read coordinates 
    */
    if(use_3d){
      fprintf(stderr,"sh_read_spatial_data_from_grd: error: 3D not implemented for grd yet\n");
      exit(-1);
    }
    /* 
       interpolate from grd 
    */
    for(k=0;k < shps;k++){
      if(!ggrd_grdtrack_interpolate_tp((double)xp[HC_THETA],(double)xp[HC_PHI],(ggrd+k),
				       &dvalue,FALSE,FALSE)){
	fprintf(stderr,"sh_read_spatial_data_from_grd: interpolation error grd %i, lon %g lat %g\n",
		  k+1,(double)PHI2LON(xp[HC_PHI]),(double)THETA2LAT(xp[HC_THETA]));
	exit(-1);
      }
      data[k*exp[0].npoints+j] = (HC_PREC) dvalue;
    }
  }	/* end points in layer loop */
}


  
#endif
