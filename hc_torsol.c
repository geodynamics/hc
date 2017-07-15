#include "hc.h"
// 
//
// these subroutines deal with the toroidal part of the 
// kinematic solution of a Hager &
// O'Connell flow code. they are based on Brad's original code, and later
// Bernhard Steinberger's modifications
//
// will incorporate the poloidal part of the plate velocities
//
//
// Thorsten Becker, twb@ig.utexas.edu
//
// $Id: hc_torsol.c,v 1.8 2006/01/22 01:11:34 becker Exp becker $
//
//     ****************************************************************
//     * THIS IS THE MAIN PROGRAM FOR THE COMPONENT OF FLOW WITHOUT   *
//     * DENSITY CONTRASTS.  IT USES SEVERAL INPUT/OUTPUT SUBROUTINES *
//     * AND FUNCTIONS TO OBTAIN, CORRECT AND VERIFY A MODEL FROM THE *
//     * USER.  THE FINAL VERSION OF EACH MODEL IS STORED IN A FILE   *
//     * BEFORE THE PROGRAM EXECUTES POLSOL AND TORSOL TO OBTAIN THE  *
//     * POLOIDAL AND TOROIDAL COMPONENTS, RESPECTIVELY, OF THE       *
//     * EQUATIONS OF MOTION.                                         *
//     ****************************************************************
//     Modified such that only toroidal component is calculated
//     Poloidal component is included in densub.f
//
//     input:  r: radii on which output is defined (nrad+2)
//
//        visc,rvisc: normalized viscosities and their radii (nvis)
//
//     lmax: MAXIMUM DEGREES,
//     nrad: NUMBER OF OUTPUT RADII 
//           (without top and bottom layers)
//     nvis: NUMBER OF VISCOSITIES.
//
//     pvel_tor: toroidal part of the plate velocities
//     pkernel: print the two solution vectors to file    
//
//    input/output:
//
//    tvec[nradp2 * lmaxp1 * 2 ]: solution kernel
//
//    output: 
//
//    tor_sol[nradp2 * 2] SHOULD BE PASSED INITIALIZED AS ZEROES
// 
//
//
//
void hc_torsol(struct hcs *hc,
	       int nrad,int nvis,int lmax,HC_PREC *r,
	       HC_PREC **rv,HC_PREC **visc, struct sh_lms *pvel_tor,
	       struct sh_lms *tor_sol,HC_HIGH_PREC *tvec,
	       hc_boolean verbose)
{
  //    
  //     ****************************************************************
  //     * evaluates AND PROPAGATES THE TWO TOROIDAL COMPONENTS IN THE  *
  //     * EQUATIONS OF MOTION, AND NORMALIZES THESE SUCH THAT THE      *
  //     * FIRST ELEMENT AT THE SURFACE IS 1.0.                         *
  //     ****************************************************************
  //
  HC_HIGH_PREC coef,*vecnor,hold,rlast,rnext,tloc[2],*tvec1,*tvec2;
  HC_HIGH_PREC exp_fac[2],p[2][2],diflog,el,elp2,elm1,efdiff;
  int l,jvisp1,jvis,i,j,nvisp1,lmaxp1,os;
  hc_boolean qvis;
  //
  //     PASSED PARAMETERS:  NRADP2: NUMBER OF OUTPUT RADII,
  //        NVIS: NUMBER OF VISCOSITIES, nvisp1 = nvis+1
  //        LMAX: MAXIMUM DEGREES.
  //     ARRAYS:  R: OUTPUT RADII,
  //        RV: VISCOSITY RADII,
  //        TVEC: TOROIDAL VECTORS,
  //        VISC: normalized VISCOSITIES.
  //     OTHER VAR:  EXP_FAC[0],EXP_FAC[1]: EXPONENTIAL FACTORS IN PROPAGATOR,
  //        COEF,ELP2,ELM1: PARAMETERS IN PROPAGATOR,
  //        DIFLOG: DIFFERENCE IN LOGS OF RADII,
  //        EL,L: DEGREE,
  //        VECNOR: NORMALIZES TVEC_LOC TO TVEC(N,1),
  //        HOLD: TEMPORARY VAR.,
  //        P[0][0],P[0][1],P[1][0],P[1][1]: ELEMENTS OF THE PROPAGATOR MATRIX CORRES-
  //        PONDING TO P(1,1),P(1,2),P(2,1),P(2,2) RESPECTIVELY,
  //        RLAST,RNEXT: RADII FOR PROPAGATOR,
  //        TVEC_LOC1,TVEC_LOC2: VECTOR COMPONENTS.
  //
  /* 
     set up some pointers (without those the TVECSOL macro won't
     work!)
  */
  nvisp1 = nvis + 1;		/* length of rv and visc */
  //nradp2 = nrad + 2;		/* radius array */
  lmaxp1 = lmax+1;		/* length of 0:lmax array */

  /* 
     add one item at end of rv and visc arrays 
  */
  hc_dvecrealloc(rv,nvisp1,"hc_torsol: rv");
  hc_dvecrealloc(visc,nvisp1,"hc_torsol: visc");
  /* local reference to viscosity and radii of viscosity */
#define HC_TVISC(i) (*(*visc+(i)))
#define HC_TVR(i) (*(*rv+(i)))

  HC_TVR(nvis) = 1.1;              /* last entry in radius array, why is
				     this 1.1? probably because it has to 
				     be > 1 
				  */
  HC_TVISC(nvis) = HC_TVISC(nvis-1);	/* last entry in viscosity array */
#ifdef DEBUG
  if(hc->nradp2 != nrad + 2){
    fprintf(stderr,"hc_torsol: radius number mismatch\n");
    exit(-1);
  }
  /* 
     test size of expansions 
  */
  j = hc->nradp2 * 2;
  for(i=0;i < j;i++){
    if(tor_sol[i].lmax < pvel_tor->lmax){
      fprintf(stderr,"hc_torsol: error: toroidal expansion %i has lmax %i, plates have %i\n",
	      i+1,tor_sol[i].lmax, pvel_tor->lmax);
      exit(-1);
    }
    if(tor_sol[i].type != pvel_tor->type)
      HC_ERROR("hc_torsol","torsol type error");
  }
#endif
  if(verbose)
    fprintf(stderr,"hc_torsol: toroidal velocities lmax %i and type %i\n",
	    pvel_tor->lmax,pvel_tor->type);
  /* 

  make room for toroidal scaling vectors f(l) and initialize as zeroes

  */
  /* solution factors as f(l,r) */
  /* set local pointes */
  tvec1 = tvec;
  tvec2 = (tvec + hc->nradp2 * lmaxp1);
  //
  //     (PREVENTS THE REQUESTING OF NON-EXISTANT VALUES)
  //     
  //     FOR EACH DEGREE (L) CALCULATE, NORMALIZE AND OUTPUT SOLUTION
  //
  for(l=1;l < lmaxp1;l++){      
    /* 
       loop through all l > 0 
    */
    el = (HC_PREC)l;
    //     
    //     SET THE PARAMETERS
    //     
    elp2 = el + 2.0;
    elm1 = el - 1.0;
    coef = 1.0 / (2.0 * el + 1.0);
    //
    //     INITIALIZE THE PROPAGATION AT THE CORE
    //
    jvisp1 = 1;			/* viscosity layer counters  */
    jvis = 0;

    rlast = r[0];		/* radius of core */
    /* 
       initialize 
    */

    tloc[0] = 1.0; /* there seems to be no best ordering for 
		      addressing this array, later we need l to 
		      be the fastest increasing index */
    tloc[1] = 0.0;
    //
    //     FIND THE TWO TOROIDAL COMPONENTS AT EACH RADIUS
    //     start radius loop
    //
    /* 
       lowest level 
    */
    os = l;
    tvec1[os] = tloc[0];
    tvec2[os] = tloc[1];

    for(i=1;i < hc->nradp2;i++){	/* loop through radii */
      os += lmaxp1;
      //
      //     TEST FOR CHANGE IN VISCOSITY IN NEXT LAYER
      //
      qvis = FALSE;
      do{
	if(HC_TVR(jvisp1) > r[i])
	  qvis = TRUE;
	rnext = HC_TVR(jvisp1);	/*  */
	//
	//     IF NO VISC. CHANGE BEFORE NEXT OUTPUT RADIUS, PROPAGATE DIRECTLY
	//
	if(qvis) 
	  rnext = r[i];
	diflog = log(rnext / rlast);
	exp_fac[0] = exp(         el * diflog);
	exp_fac[1] = exp(-(el + 1.0) * diflog);
	//
	//     PROPAGATOR SET UP LINEARLY TO AVOID EXCESS MULTIPLICATIONS
	//
	efdiff = exp_fac[0] - exp_fac[1];
	p[0][0] = elp2 * exp_fac[0] + elm1 * exp_fac[1];
	p[0][1] = efdiff / HC_TVISC(jvis);
	p[1][0] = elp2 * elm1 * HC_TVISC(jvis) * efdiff;
	p[1][1] = elm1 * exp_fac[0] + elp2 * exp_fac[1];
	//
	//     PROPAGATE LAST VECTOR TO GET NEW VECTOR
	//
	rlast = rnext;
	hold = 	tloc[0];
 	tloc[0] = (p[0][0] * hold + p[0][1] * tloc[1]);
	tloc[1] = (p[1][0] * hold + p[1][1] * tloc[1]);
	tloc[0] *= coef;
	tloc[1] *= coef;
	if(!qvis){
	  jvis = jvisp1;
	  jvisp1++;
	}
      }while(!qvis);
      tvec1[os] = tloc[0];
      tvec2[os] = tloc[1];
    } /* end layer loop */
  } /* end l loop */
  
  // 
  //     set tvec(l,nradp2-1,0) = 1.0 and normalize all vectors to
  //     this
  //
  hc_hvecalloc(&vecnor,lmaxp1,"hc_torsol: vecnor");
  os = (hc->nradp2-1) * lmaxp1;
  vecnor[0] = 1.0;
  for(l=1;l < lmaxp1;l++)
    vecnor[l] = 1.0 / tvec1[os+l];
  /* normalize */
  for(i=os=0;i < hc->nradp2;i++,os+=lmaxp1)
    for(l=0;l < lmaxp1;l++){
      tvec1[os+l] *= vecnor[l];
      tvec2[os+l] *= vecnor[l];
    }
  free(vecnor);
  /* 
     
  the toroidal solution corresponds to the toroidal part of the plate
  motions scaled by the toroidal solution vectors which are functions
  of l and depth

  */
  for(os=i=j=0;i < hc->nradp2;i++,os+=lmaxp1,j+=2){
    /* 
       assign toroidal plate motion fields to solution expansion
    */
    sh_aexp_equals_bexp_coeff((tor_sol+j+0),pvel_tor);
    sh_aexp_equals_bexp_coeff((tor_sol+j+1),pvel_tor);
    /* 
       scale with the toroidal solution at this depth 
    */
    sh_scale_expansion_l_factor((tor_sol+j+0),(tvec1+os));
    sh_scale_expansion_l_factor((tor_sol+j+1),(tvec2+os));
  }
  if(verbose)
    fprintf(stderr,"hc_torsol: done\n");
}
#undef HC_TVISC
#undef HC_TVR



