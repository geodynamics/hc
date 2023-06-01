#include "hc.h"
/* 

compute the solution for the poloidal part of a Hager & O'Connell flow
computation.

the poloidal part has two contributions: density driven flow and plate motions

this routine computes the y_i i=1,...,6 solutions for each layer and
incorporates the poloidal part of the plate motions, which are passed
as pvel_pol



*/
void hc_polsol(struct hcs *hc, 	/* 
				   general Hager & O'Connell solution
				   structure, holds constants and such
				*/
	       /* general output radii  */
	       int nrad, 		/* number of output radii */
	       HC_PREC *rad,		/* output radii, normalized by
					   R */
	       /* 
		  density contribution
	       */
	       int inho,		/* number of density layers */
	       HC_PREC *dfact, 	/* density factors from layer thickness */
	       hc_boolean viscosity_or_layer_changed, /* if TRUE, will
						       (re)compute the arrays 
						       that depend on the spacing
						       of the density anomalies 
						       or the viscosity structure
						      */
	       struct sh_lms *dens_anom, /* 
					expansions of density
					anomalies has to be [inho] 

				     */
	       hc_boolean compressible, /* 
					   if TRUE, will use PREM
					   densities, else, average
					   mantle density
					*/
	       int npb,		/* number of phase boundaries */
	       HC_PREC *rpb,	/* radius and F factor for phase
				   boundaries */
	       HC_PREC *fpb,
	       hc_boolean free_slip,  /* 
					 include plate velocities?
					 possible, if free_slip is
					 FALSE
				      */
	       struct sh_lms *pvel_pol, /* 
					
				     poloidal part of plate motions
				     (only one expansion), only gets accessed
				     if free_slip is false
				     
				     */
	       struct sh_lms *pol_sol,	  /* 
					     poloidal solution
					     expansions 
					     [nout * 6]
					     nout <= nrad
					     
					     SHOULD BE PASSED INITIALIZED AS ZEROES

					  */
	       hc_boolean compute_geoid, /* additionally compute the geoid? */
	       struct sh_lms *geoid, 	/* geoid solution */
	       hc_boolean save_prop_mats, /* 
					     memory intensive speedup
					     by saving all propagator
					     matrices. this makes
					     sense if the density
					     anomalies changes between
					     call, but nothing else
				     */
	       hc_boolean verbose, /* output options */
	       hc_boolean calc_kernel_only /* only compute the
					      kernels */
	       )
{

  //    ****************************************************************
  //    * THIS PROGRAM IS TO USED CALCULATE AND OUTPUT THE POLOIDAL    *
  //    * COMPONENTS OF THE FLOW WITH PLATE MOTIONS                    *
  //    * AND DENSITY CONTRASTS.  THIS PROGRAM REQUIRES THREE          *
  //    * INPUT FILES CONTAINING:  (1) THE MODEL, (2) THE EXPANDED     *
  //    * DENSITY CONTRASTS, AND (3) DENSITY FACTORS (*DR) AND THEIR   *
  //    * RADII.  GENERALLY TO BE USED AFTER PROGRAM NODENC, BUT CAN BE*
  //    * USED BEFORE IT.                                              *
  //    * USES U(A) = P(A,C)*U(C)+SUM OVER I OF (P(A,R(I))*B(I)*DR(I)) *
  //    * PROBLEM IS SEPARATED INTO TWO PARTS:  U(1) THROUGH U(4), THE *
  //    * VELOCITIES AND STRESSES, HAVE BEEN SEPARATED FROM THE        *
  //    * POTENTIAL AND DERIVATIVE USING U(3)NEW = U(3)OLD + RHO * U(5)*
  //    * WHERE RHO IS THE AVERAGE DENSITY OF THE MANTLE.  THE U(3)    *
  //    * IN OUTPUT IS U(3)NEW.  RHO * U(5) MUST BE SUBTRACTED.        *
  //    ****************************************************************
  //    This substraction is not being done in this program. This is,
  //    because vertical stresses are used to calculate stresses 
  //    in the lithosphere, but in this case, not the actual elevation
  //    but the elevation above the equipotential surface, therefore
  //    u(3)new and not u(3)old matters. For the surface jump in gravity
  //    however, u(3)old matters, this is incorporated below.
  //    See e.g. Panasyuk and Hager (1996)
  //
  //
  // this routine has been modified from the orginial version, for which 
  // you find the comments above. not all comments reflect these 
  // changes, so beware. original code by B Hager, then modified by 
  // RJO and Bernhard Steinberger
  //
  //
  // Thorsten Becker twb@ig.utexas.edu
  //
  //
  // $Id: hc_polsol.c,v 1.12 2006/03/20 05:32:48 becker Exp becker $
  //
  //
  //       ARRAYS:
  //       B: SPH. HARM. EXPANSION OF DENSITY CONTRASTS FOR EACH
  //          INHOMOGENEOUS LAYER AND READ IN AT EACH NEW IND1 AND IND2,
  //       pvel_pol: POLOIDAL PART OF PLATE VELOCITIES,        
  //       D: FROM SUBROUTINE GETDEN FOR USE IN U(A) = P(A,C)*U(C)+D,
  //       DEN: THE FACTOR FACT*RDEN*RDEN*ALPHA WHICH IS MULTIPLIED
  //          BY B, THE PRODUCT BEING ADDED TO U(3) AT EACH RADIUS DURING
  //          PROPAGATION WITHIN THE M LOOP, (AT RADII WITH NO DENSITY
  //          CONTRASTS, DEN IS SET TO ZERO), U3(R+)=U3(R-)+DEN*B,
  //       DPOT: FROM GETDEN FOR POTEN(A) = PPOT(A,C)*POTEN(C)+DPOT
  //          WHERE POTEN(A) = [GA,-(L+1)*GA]T, POTEN//= [GC,L*GC]T,
  //       FACT: DENSITY FACTORS (ALLOW VAR. DENS. CONTRAST WITH DEPTH),
  //          THESE ALSO INCLUDE DR, THE RATIO OF RADII MIDPOINTS
  //          BETWEEN RDEN RADII,
  //       PROP: PROPAGATOR,
  //       POTDEN: PRODUCED BY PROPIH FOR USE IN GETDEN (POTENTIALS),
  //       POTEN,POTNEW: THE POTENTIAL AND ITS DERIVATIVE,
  //       PPOT: A POTENTIALS PROPAGATOR (FROM EVPPOT),
  //       PPOTS: THE ARRAY OF ALL POTENTIALS PROPAGATORS TO OBTAIN
  //          POTEN AT EACH DESIRED RADIUS FOR ALL M AT A GIVEN L,
  //       PROPS: THE ARRAY OF ALL PROPAGATORS NECESSARY TO OBTAIN THE
  //          U VECTOR AT EACH DESIRED RADIUS FOR ALL M AT A GIVEN L,
  //       PRPDEN: FACTORS FROM SUBROUTINE PROPIH FOR SUBROUTINE GETDEN,
  //       PVISC: THE VISCOSITY FOR EACH LAYER USED IN EVALUATING PROPS,
  //       QWRITE: A LOGICAL ARRAY USED TO DECIDE WHETHER A GIVEN U
  //          VECTOR IS ONE REQUIRED FOR OUTPUT, (AT AN OUTPUT RADIUS),
  //       RAD: DESIRED OUTPUT RADII,
  //       RDEN: RADII OF INHOMOGENEOUS DENSITY CONTRASTS,
  //       RVISC: RADII OF VISCOSITIES,
  //       U,UNEW: POLOIDAL COMPONENTS OF FLOW,
  //       VISC_LOCAL: VISCOSITIES.
  //    OTHER VAR:
  //       DOUBLE PRECISION:
  //          ALPHA: RE*GACC*180*SECYR*TIMESC*1/(VISNOR*PI)
  //             CONVERSION FACTORS AS USED IN PROPIH,
  //          BETA: -4*PI*G*RE/GACC, CONV. FACTORS AS USED IN PROPIH,
  //          EL: DEGREE (L),
  //          ELIM: PARAMETER USED IN SIMPLE ELIMINATION,
  //          G: UNIVERSAL GRAVITY CONSTANT (SI),
  //          GC: NON-EQUILIBRIUM GRAV. POTENTIAL AT CORE,
  //          GACC: GRAVITATIONAL ACCELERATION AT SURFACE (SI),
  //          RE: RADIUS OF THE EARTH = R_DEF*1e3 (SI),
  //          RNEXT: NEXT RADIUS IN PROPAGATION,
  //          SC: STRESS AT CORE,
  //          SECYR: SECONDS PER YEAR,
  //          TIMESC: TIMESCALE OF MOTION USUALLY 1 M.Y.,
  //          VC: VELOCITY AT CORE,
  //          VISNOR: NORMALIZING VISCOSITY,
  //       INTEGER: 
  //          INDEX: DETERMINES THE ARRAY INDICES FOR DEN,PVISC,QWRITE
  //             AND RPROPS DURING ORDERING AND INITIALIZATION,
  //
  //          INHO: NUMBER OF INHOMOGENEOUS RADII,
  //
  //          IVIS: PRESENT VISCOSITY OR DENSITY LAYER,
  //          L: DEGREE
  //          M: ORDER
  //          NEWPRP: DETERMINES INDEX OF PROPEQ IN STORING PROPAGATORS,
  //          NIH: PRESENT INHOMOGENEOUS LAYER IN EVALUATING DEN,
  //          NINHO: PRESENT INHOMOGENEOUS LAYER IN EVALUATING U,
  //          NPROPS: THE TOTAL NUMBER OF PROPAGATORS TO GET U,
  //          NRADP2: NUMBER OF OUTPUT RADII,
  //          NVIS: NUMBER OF VISCOSITIES,
  //       LOGICAL:  
  //          QINHO: DETERMINES IF NEXT RADIUS IS A DENSITY CONTRAST,
  //          QVIS:  DETERMINES IF NEXT RADIUS IS A VISCOSITY CHANGE,
  //       integer:
  //          a_or_b: ALLOWS CALCULATION OF S(LM) AFTER CALCULATION OF
  //             C(LM) EXCEPT AT M=0.(a_or_b == 0: A, a_or_b == 1: B)
  //
  //    SUBROUTINES AND FUNCTIONS:
  //       SUBROUTINE EVPPOT (L,RATIO,PPOT):  OBTAINS PROPAGATOR FOR
  //          NON-EQUILIBRIUM POTENTIAL AND DERIVATIVE (RATIO IS R(I)/
  //          R(I+1), FOR PROPAGATION FROM R(I) TO R(I+1) AT L),
  int i,i2,i3,i6,j,l,m,nih,nxtv,ivis,os,pos1,pos2,gi,g1,g2,gic,
    prop_s1,prop_s2,nvisp1,nzero,n6,ninho,nl=0,ip1;
  int newprp,newpot,jpb,inho2,ibv,indx[3],a_or_b,ilayer,lmax,
    nprops_max,jsol,mmax;
  int klayer = 1;
  double *xprem;
  HC_HIGH_PREC *b,du1,du2,el,rnext,drho,dadd;
  HC_PREC rbound_kludge;
  HC_HIGH_PREC amat[3][3],bvec[3],u[4],poten[2],
    unew[4],potnew[2],clm[2];
  /* 
     structures which hold u[6][4] type arrays 
  */
  struct hc_sm cmb, *u3;
  hc_boolean qvis,qinho,hit,kludge_warned;
  /*  
      define a few offset and size pointers
  */

#ifdef HC_DEBUG
  if(hc->nradp2 != nrad + 2){
    fprintf(stderr,"hc_polsol: radius number mismatch\n");
    exit(-1);
  }
#endif


  inho2 = inho + 2;
  nvisp1 = hc->nvis+1;
  lmax = pol_sol[0].lmax ;
  /* 
     max number of propagator levels, choose this generously
  */
  nprops_max = hc->nradp2 * 3;
  /* 
     for prop and ppot: one set of propagators for all layers, there
     lmax of those
  */
  prop_s1 = nprops_max * 16;
  prop_s2 = nprops_max * 4;

  /* 
     check if still same general number of layers 
  */
  if((hc->psp.prop_params_init)&&((inho2 != hc->inho2)||
			  (nvisp1 != hc->nvisp1))){
    HC_ERROR("hc_polsol","layer structure changed from last call");
  }
  /*
    
  allocate space for local arrays
  
  */ 				 
  /* inho + 2 */
  u3 = (struct hc_sm *)calloc(inho2,sizeof(struct hc_sm));
  if(!u3)
    HC_MEMERROR("hc_polsol: u3");
  hc_vecalloc(&b,inho2,"hc_polsol");
  if(save_prop_mats){
    /* 
       propagators saved
    */
    if(!hc->psp.prop_mats_init){
      /* 
	 
      we will be saving all propagator matrices. this makes sense if
      the density structure is the only thing that changes
      
      this needs quite a bit more room (array goes from l=1 (not l=0)
      .... lmax)
      
      
      */
      hc_hvecalloc(&hc->props,prop_s1 * lmax,"hc_polsol");
      hc_hvecalloc(&hc->ppots,prop_s2 * lmax,"hc_polsol");
    }
  }else{
    /* 
       propagator recomputed and reallocated each time
    */
    hc_hvecalloc(&hc->props,prop_s1,"hc_polsol");
    hc_hvecalloc(&hc->ppots,prop_s2,"hc_polsol");
  }
  if(!hc->psp.abg_init){
    //
    //    SET alpha, beta and geoid factors
    //
    hc->psp.alpha  =  hc->psp.rho_scale * (hc->re*10.) * hc->gacc / hc->visnor;	/*  */
    hc->psp.alpha *= ONEEIGHTYOVERPI * hc->secyr * hc->timesc; /*  */
    //
    hc->psp.beta   = -4.0 * HC_PI * (hc->g*1e3) * (hc->re*1e2) / hc->gacc;
    if(verbose)
      fprintf(stderr,"hc_polsol: alpha: %.8f beta: %.8f\n",
	      (double)hc->psp.alpha,(double)hc->psp.beta);
    
    /* 
       geoid scaling factor hc->gacc shoud be grav[nprops] for
       compressibility
    */
    hc->psp.geoid_factor = HC_PI * 10.0* hc->visnor/180./hc->secyr/hc->gacc/1.e8;
    hc->psp.abg_init = TRUE;
  }
  if((!hc->psp.prop_params_init) || (viscosity_or_layer_changed)){
    /* 
       
       intialize arrays that depend on viscosity and density layer spacing 
       
    */
    //    
    //    CREATE DEN,PVISC,QWRITE,RPROPS AS FOLLOWS:
    //    1)  INITIALIZE PVISC=VISC(IVIS), DEN=ZERO, QWRITE=FALSE
    //    2)  FIND WHICH RADIUS (RAD,RDEN,RVISC) IS NEXT IN SEQUENCE
    //    TO SURFACE, NOTING THAT ANY TWO OR ALL THREE MAY BE EQUAL
    //    3)  INCREMENT INDEX AND STORE RNEXT IN RPROPS
    //    4)  IF AT RVISC(IVIS) INCREMENT IVIS
    //    5)  IF AT RDEN(NIH) EVALUATE DEN, INCREMENT NIH
    //    6)  IF AT RAD(I) QWRITE = TRUE, INCREMENT I
    //    
    if(!hc->psp.prop_params_init){
      if(verbose)
	fprintf(stderr,"hc_polsol: initializing for %i v layers and %i dens layers\n",
		nrad,inho);
      /* 
	 this is really the first call, allocate arrays

	 arrays that go with nprops
	 
      */
      hc_hvecalloc(&hc->rprops,nprops_max,"hc_polsol: rprop");
      hc_hvecalloc(&hc->pvisc,nprops_max,"hc_polsol");
      hc_hvecalloc(&hc->den,nprops_max,"hc_polsol");
      /* initialize qwrite with zeroes! */
      hc->qwrite = (hc_boolean *)calloc(nprops_max,sizeof(hc_boolean));
      if(!hc->qwrite)
	HC_MEMERROR("hc_polsol: qwrite");
      /* those that go with (inho=nrad)+2 */
      hc_vecrealloc(&hc->rden,inho2,"hc_polsol");
      /* and those for nvis+1 */
      hc_vecrealloc(&hc->rvisc,nvisp1,"hc_polsol");
      hc_vecrealloc(&hc->visc,nvisp1,"hc_polsol");
      /* 
	 save in case we want to check if parameters changed later
      */
      hc->inho2 = inho2;hc->nvisp1=nvisp1;
    }
    //
    //    SET RDEN(INHO+1) = 1.1 TO PREVENT TESTING OF THAT VALUE
    //    
    hc->rden[inho] = 1.1;
    //
    //    APPEND A FINAL RVISC_LOCAL = 1.0 TO PREVENT OUT OF BOUNDS
    //    
    hc->rvisc[hc->nvis] = 1.0;

    //    
    //    INITIALIZE INDEX,IVIS,NIH
    //
    hc->nprops = ivis = nih = 0;
    hc->rprops[0] = rad[0];
    hit = FALSE;
    for(i=1;(i < hc->nradp2)&&(!hit);i++){
      //    
      //    INITIALIZE
      //
      do{
	qinho = TRUE;		/* is next radius a density contrast? */
	qvis = TRUE;
	//    new check, when two radii happen to be the same, exit the
	//    loop
	if((hc->nprops > 0) && 
	   (fabs(hc->rprops[hc->nprops] - hc->rprops[hc->nprops-1])
	    <HC_EPS_PREC)){
	  hit = TRUE;		/* bailout here */
	}else{
	  /* 
	     normal operation 
	  */
	  hc->pvisc[hc->nprops]  = hc->visc[ivis];
	  hc->den[hc->nprops]    = 0.0;
	  hc->qwrite[hc->nprops] = FALSE;	
	  //    
	  //    FIND NEXT RADIUS
	  //
	  nxtv = ivis + 1;
	  if((hc->rden[nih] <= rad[i])&&
	     (hc->rden[nih] <= hc->rvisc[nxtv]))
	    qinho = FALSE;
	  if ((hc->rvisc[nxtv] <= hc->rden[nih])&&
	      (hc->rvisc[nxtv] <= rad[i]) && 
	      (ivis < hc->nvis))
	    qvis = FALSE;
	  rnext = hc->rden[nih];
	  if (!qvis)
	    rnext = hc->rvisc[nxtv];
	  if(qinho && qvis) 
	    rnext = rad[i];
	  //    
	  //    INCREMENT NPROPS, STORE RPROPS
	  //
	  hc->nprops++;
	  if(hc->nprops > nprops_max){ /* check, if we have enough room */
	    fprintf(stderr,"hc_polsol: error: nprops: %i nprops_max: %i\n",
		    hc->nprops,nprops_max);
	    exit(-1);
	  }
	  hc->rprops[hc->nprops] = rnext;
	  //    
	  //    IF RVISC, INCREMENT IVIS
	  //
	  if (!qvis) 
	    ivis = nxtv;
	  if (!qinho) {
	    //
	    //    IF RDEN, EVALUATE DEN, INCREMENT NIH
	    //    
	    hc->den[hc->nprops-1] = dfact[nih] * hc->rden[nih] * hc->rden[nih] * hc->psp.alpha;
	    nih++;
	  }
	}
	//    
	//    IF NOT RAD, DO NOT INCREMENT I
	//    
      }while((!hit) && (fabs(rnext-rad[i])>HC_EPS_PREC));
      if(!hit){
	//    
	//    IF RAD, QWRITE = TRUE
	//
	hc->qwrite[hc->nprops-1] = TRUE;	
      }
    } /* end of nrad loop */
    hc->den[hc->nprops] = 0.0;
    hc->pvisc[hc->nprops]  = hc->pvisc[hc->nprops-1]; /* to look nicer */

    /* 

       number of propagators is now nprops+1

    */
    if(verbose >= 3){
      if(hc->psp.prop_params_init)
	fprintf(stderr,"hc_polsol: using old parameters: %i v layers and %i dens layers\n",
		nrad,inho);
      for(i=i2=0;i < hc->nprops+1;i++){
	if(fabs(hc->den[i]) > HC_EPS_PREC)
	  i2++;
	fprintf(stderr,"hc_polsol: prop: i: %3i(%3i) r: %8.5f v: %8.3f den: %12g ninho: %3i/%3i\n",
		i+1,hc->nprops,(double)hc->rprops[i],
		(double)hc->pvisc[i],(double)hc->den[i],i2,inho);
      }
    }
    if(!hc->psp.rho_init){
      /* 
	 initialize the density factors, for incompressible, those
	 are all constant, else from PREM
      */
      hc_vecalloc(&hc->rho,nprops_max+2,"hc_polsol: rho");
      /* this way, rho_zero can go from -1...nnprops_max */
      hc->rho_zero = (hc->rho+1);
      if(compressible){
	/* 
	   
	for compressible computation, assign densities from PREM, but
	only use the first 10 layers (below crust and ocean, I think)

	densities are in kg/m^3

	*/
	if(!hc->prem_init)
	  HC_ERROR("hc_polsol","PREM wasn't initialized for compressible");
	hc_dvecalloc(&xprem,hc->prem->np,"hc_polsol: rho");
	for(i=0;i < hc->nprops+1;i++){
	  ilayer = prem_find_layer_x((double)hc->rprops[i],1.0,
				     hc->prem->r,
				     10,hc->prem->np, 
				     xprem);
	  hc->rho_zero[i] = prem_compute_pval(xprem,
					      (hc->prem->crho+ilayer*hc->prem->np),
					      hc->prem->np,1.0);
	}
	free(xprem);
      }else{
	/* 

	for the incompressible computation, use average values of
	density for the mantle
	
	densities in kg/m^3

	*/
	hc->rho_zero[-1] = hc->avg_den_core;
	for(i=0;i < hc->nprops+1;i++)
	  hc->rho_zero[i] = hc->avg_den_mantle;
      }
      hc->rho_zero[hc->nprops+1] = 0.0;
      hc->psp.rho_init = TRUE;  
    } /* end rho init */
    hc->psp.prop_params_init = TRUE;
    /* 
       
    end of the propagator factor section, this will only get executed
    once unless the density factors or viscosities change
    
    */
  }
  hc->rprops[hc->nprops+1] = 1.0;


  if(verbose >= 3)
    for(i=0;i < hc->nprops+2;i++)
      fprintf(stderr,"i: %3i nprops: %3i r(i): %11g rho: %11g\n",
	      i,hc->nprops,(double)hc->rprops[i],(double)hc->rho_zero[i]);
  
  //
  //    begin l loop
  //
  if(verbose)
    fprintf(stderr,"hc_polsol: ncalled: %5i for lmax: %i dens lmax: %i, visc or layer %s changed\n",
	    hc->psp.ncalled,pol_sol[0].lmax,dens_anom->lmax,
	    ((viscosity_or_layer_changed)?(""):("not")));
  if(free_slip)			/* select which components of pol solvec to 
				   use */
    nzero = 3;
  else
    nzero = 1;
  pos1 = pos2 = 0;		/*
				  offset pointers for propagators,
				  non-zero only if the propagators are
				  stored
				*/
  for(l = 1;l <= pol_sol[0].lmax;l++){
    /* 
       
    MAIN L LOOP, start at l = 1 (only anomalies)
    
    */
    el = (HC_PREC)l;

    /* 
       this will normally be a very small number so that all
       propagators will be computed above the regular CMB 

       if solver_kludge_l is set to within the [0;L] domain, the depth
       of the bottom will depend on l
       
       (rprops(i).ge.(1.-(1.-0.5448)*50./l))
    */
    rbound_kludge = (1. - (1.-hc->r_cmb)*(HC_PREC)hc->psp.solver_kludge_l/el);
    kludge_warned = FALSE;

    if((!save_prop_mats) || (!hc->psp.prop_mats_init)|| (viscosity_or_layer_changed)){
      //    
      // get all propagators now, as they only depend on l
      //    
      for(newprp = pos1, newpot = pos2, 
	    i = 0;i < hc->nprops;i++, 
	    newprp += 16, newpot += 4){
	/*  
	    obtain and save propagators 
	*/
	hc_evalpa(l,hc->rprops[i],hc->rprops[i+1],
		  hc->pvisc[i],(hc->props+newprp));
	hc_evppot(l,(hc->rprops[i]/hc->rprops[i+1]),
		  (hc->ppots+newpot));
      }	/* i checked the propagator matrices again, those are as in
	   Bernhard's code TWB */
    }
    /* 

    begin m loop

    */
    mmax = (calc_kernel_only)?(0):(l);
    for(m=0;m <= mmax;m++){
      /* 


      START M LOOP 


      */
      //    
      //    CALCULATE C(LM) FOR ALL M, S(LM) FOR M>0
      //
      a_or_b = 0;		/* start with A coefficient */
      do{			/* do loop for A/B evaluation */
	if((!calc_kernel_only) && (l <= dens_anom[0].lmax)){
	  /*
	    
	    obtain the coefficients from the density field expansions

	  */    
	  for(i=0;i < inho;i++)/* 
				  A or B coeff, use the internal
				  convention here, as stored before
			       */
	    sh_get_coeff((dens_anom+i),l,m,a_or_b,FALSE,(b+i));
	  //hc_print_vector(b,inho,stderr);
	}else{
	  /* 
	     density is not expanded to that high an l
	  */
	  for(i=0;i < inho;i++) 
	    b[i] = 0.0;
	}
	b[inho] = 0.0;
	if(calc_kernel_only)
	  b[klayer] = 1.0;
	
	//    
	//    U(C) = [0,VC,SC,0], U(A) = [0,0,SA,SX]
	//    POT(A) = [U5(A),-(L+1)*U5(A)]T, POT(C) = [U5(C),L*U5(C)]T
	//    Find three linear independent solutions of homogeneous eqns and 
	//    one solution of inhomogeneous eqn., all satisfying boundary
	//    conditions at core by integrating from core up to the surface.
	//    Find linear combination that satisfies surface boundary conditions.
	//
	for(i6=0;i6 < 6;i6++)	/* initialize cmb with zeroes */
	  for(ibv=0;ibv < 4;ibv++)
	    cmb.u[i6][ibv] = 0.0;

	if(l > hc->psp.solver_kludge_l){
	  /* 
	     solver trick to ensure stabilty following Steinberger &
	     Torsvik, doi:10.1029/2011GC003808

	     ucmb(4,1)=1.d0
	  */
	  cmb.u[3][1] = 1.0;	/* make core fixed */
	}else{
	  /* regular operation */
	  /*  ucmb(2,1)=1.d0 */
	  cmb.u[1][1] = 1.0;	/* set this to zero for no-slip,
				   in general CMB is free slip */
	}

	cmb.u[2][2] = 1.0;	/* ucmb(3,2)=1.d0 */
	cmb.u[4][3] = 1.0;	/* ucmb(5,3)=1.d0 */
	cmb.u[5][3] = el;	/* ucmb(6,3)=float(l) */

	for(ibv=0;ibv < 4;ibv++){
	  /* 

	  IBV LOOP

	  */
	  for(i=0;i < 4;i++)
	    u[i] = cmb.u[i][ibv];
	  poten[0] = cmb.u[4][ibv];
	  poten[1] = cmb.u[5][ibv];
	  //
	  //    Propagate gravity across CMB. Inside the core, surfaces of
	  //    constant pressure coincide with surfaces of constant potential.
	  //
	  /* 
	     if we were allowing for compressibility, would multi with
	     hc->grav[i]/hc->grav here (beta incorporates 1/grav0)
	  */
	  poten[1] += hc->psp.beta * hc->rprops[0] *
	    (u[2] - (hc->rho_zero[0] - hc->rho_zero[-1]) * 
	     poten[0]);

	  ilayer = 0;
	  

	  u3[ilayer].u[0][ibv] = u[0]; /* flow/stress */
	  u3[ilayer].u[1][ibv] = u[1];
	  u3[ilayer].u[2][ibv] = u[2];
	  u3[ilayer].u[3][ibv] = u[3];
	  u3[ilayer].u[4][ibv] = poten[0]; /* potential solution */
	  u3[ilayer].u[5][ibv] = poten[1];

	  ninho = jpb = 0;
	  for(i=0,ip1=1;i < hc->nprops;i++,ip1++){
	    /* 

	    I NPROPS LOOP

	    */
	    


	    if(hc->rprops[ip1] >= rbound_kludge){
	      //
	      //    PROPAGATE U TO NEXT RADIUS IN RPROPS
	      //    
	      for(os=pos1 + i*16,i2=0;i2 < 4;i2++,os += 4){
		unew[i2] = 0.0;
		for(i3=0;i3 < 4;i3++){
		  unew[i2] += hc->props[os + i3] * u[i3];
		}
	      }
	      hc_a_equals_b_vector(u,unew,4);
	      //    
	      //    PROPAGATE POTEN TO NEXT RADIUS
	      //    
	      os = pos2 + i * 4;
	      potnew[0] = poten[0] * hc->ppots[os+0] + 
		poten[1] * hc->ppots[os+1];
	      poten[1]  = poten[0] * hc->ppots[os+2] + 
		poten[1] * hc->ppots[os+3];
	      poten[0] = potnew[0];
	      if(ibv == 0){
		//    
		//    ADD DEN * B, WHERE DEN = 0 FOR NO DENSITY CONTRAST
		//    
		dadd = hc->den[i] * b[ninho];
		u[2] += dadd;	/* this would have a factor 
				   grav(i)/hc->grav
				*/
		//    
		//    ADD DEN * BETA * B * RDEN
		//    
		//	      fprintf(stderr,"%15.5e %15.5e %15.5e %15.5e\n",
		//      beta, hc->den[i], b[ninho],hc->rden[ninho]);
		poten[1] += hc->psp.beta * dadd * hc->rden[ninho];
	      }
	      //    
	      //    Changes due to radial density variations
	      //
	      drho = hc->rho_zero[i] - hc->rho_zero[ip1];
	      
	      du1 = u[0] * drho/hc->rho_zero[ip1];
	      du2 = du1 * (hc->pvisc[i]+hc->pvisc[ip1]);
	      u[0] +=  du1;
	      u[2] -=  2.0 * du2 + drho * poten[0];
	      u[3] +=        du2;
	      //hc_print_vector(poten,2,stderr);
	      //hc_print_vector(u,4,stderr);
	      //    
	      //    effects of phase boundary deflections
	      //    
	      if ((jpb < npb)&&(hc->rprops[ip1] > rpb[jpb] - 0.0001)){
		if (ibv == 0) {
		  u[2] -= fpb[jpb] * b[ninho];
		  poten[1] -= hc->psp.beta * rpb[jpb] * fpb[jpb] * b[ninho] * hc->psp.rho_scale;
		}
		jpb++;
	      }
	      /* end of l-dependent solver kludge branch */
	    }else{
	      if((verbose)&&(!kludge_warned)){
		fprintf(stderr,"hc_polsol: applying CMB fixed kludge above %3i and shifting CMB to %6.1f km depth for l %3i\n",
			hc->psp.solver_kludge_l,
			(double)HC_Z_DEPTH(rbound_kludge),l);	      
		kludge_warned = TRUE;
	      }
	    }

	    //    
	    //    IF AT A DENSITY CONTRAST, INCREMENT NINHO FOR NEXT ONE
	    
	    if(fabs(hc->den[i]) > HC_EPS_PREC)
	      ninho++;
	    //    
	    //    IF AT OUTPUT RADIUS, assign u, poten
	    //
	    if(hc->qwrite[i]){
	      ilayer++;
	      //fprintf(stderr,"%4i %4i %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e\n",
	      //l,m,u[0],u[1],u[2],u[3],poten[0],poten[1]);
	      u3[ilayer].u[0][ibv] = u[0];
	      u3[ilayer].u[1][ibv] = u[1];
	      u3[ilayer].u[2][ibv] = u[2];
	      u3[ilayer].u[3][ibv] = u[3];
	      u3[ilayer].u[4][ibv] = poten[0];
	      u3[ilayer].u[5][ibv] = poten[1];
	    }
	  } /* 
	       end i,ip1 < nrprops loop  
	    */
	  
	  //    Propagate gravity across surface. Above surface, normal stress
	  //    is zero, which determines the surface elevation. Jump in gravity
	  //    is proportional to total surface elevation (not minus equipotential
	  //    surface)
	  //
	  poten[1] -= hc->psp.beta * hc->rprops[hc->nprops] *
	    (u[2] - hc->rho_zero[hc->nprops] * poten[0]);
	  nl = ilayer;
	  u3[nl].u[5][ibv] = poten[1];

	  //fprintf(stderr,"%3i u %12.4e %12.4e %12.4e %12.4e p %12.4e %12.4e\n",
	  //ilayer, u[0],u[1],u[2],u[3],poten[0],poten[1]);

	  //    end ibv loop
	}
	nl = ilayer+1;
	//    
	//    Here plate motions are incorporated 
	//    Distinguish between free-slip (nzero=4) and no-slip with
	//    optional plate motions (nzero=2)
	//
	//    
	//    AP_l,m = cpol(l+1,m+1), AT_l,m = ctor(l+1,m+1)
	//    and
	//    BP_l,m = cpol(m,l+1),   BT_l,m = ctor(m,l+1)
	//
	//    u_2 = y_2 solution part
	//
	/* 
	   get one coefficient from the poloidal plate motion part
	*/
	if(!free_slip)
	  sh_get_coeff(pvel_pol,l,m,a_or_b,FALSE,clm); /* use internal convention */
	else
	  clm[0] = 0.0;
	/* 
	   B vector 
	*/
	bvec[0]=         u3[ilayer].u[    0][0];
	bvec[1]=         u3[ilayer].u[nzero][0] - clm[0];
	bvec[2]=(el+1.0)*u3[ilayer].u[    4][0] + u3[ilayer].u[5][0];
	/* 
	   A matrix
	*/
	for(i=0,i2=1;i < 3;i++,i2++){
	  amat[0][i] = u3[ilayer].u[    0][i2];
	  amat[1][i] = u3[ilayer].u[nzero][i2];
	  amat[2][i] = (el + 1.0) * u3[ilayer].u[4][i2] + u3[ilayer].u[5][i2];
	}
	/* 
	   solve A x = b, where b will be modified  
	*/
	if(l == 1){
	  jsol = 2;		/* 2x2 solution */
	}else{
	  jsol = 3;		/* 3x3 solution */
	}
	hc_ludcmp_3x3(amat,jsol,indx);
	hc_lubksb_3x3(amat,jsol,indx,bvec);
	/* 
	   assign solution 
	*/
	for(os=ilayer=0;ilayer < nl;ilayer++,os+=6){
	  for(i6=0;i6 < 6;i6++){
	    /* sum up contributions from vector solution */
	    for(i2=1,j=0;j < jsol;j++,i2++){
	      u3[ilayer].u[i6][0] -= bvec[j]*u3[ilayer].u[i6][i2];
	    }
	    //fprintf(stderr,"%i %i %i %i %g\n",l,m,ilayer,i6, u3[ilayer].u[i6][0]);
	    /* 
	       adding vector components to spherical harmonic solution 
	    */
	    /* A or B coefficients */
	    sh_write_coeff((pol_sol+os+i6),l,m,a_or_b,FALSE, /* use internal convention */
			   &u3[ilayer].u[i6][0]);
	  }
 	} /* end layer loop */

	/* 

	   chemical layering (e.g. phase boundaries) go here

	*/
	

	if((!a_or_b) && (m != 0))
	  //    IF S(LM) IS REQUIRED, GO BACK AND CALCULATE IT
	  a_or_b = 1;
	else
	  a_or_b = 0;
     }while(a_or_b);
    } /* end m loop */
    if(save_prop_mats){
      /* 
	 we want save the propagator matrices 
      */
      pos1 += prop_s1;
      pos2 += prop_s2;
    }
  } /* end l loop */
  if(save_prop_mats)
    /* only now can we set the propagator matrix storage scheme to TRUE */
    hc->psp.prop_mats_init = TRUE;
  if(verbose)
    fprintf(stderr,"hc_polsol: assigned nl: %i nprop: %i nrad: %i layers\n",
	    nl,hc->nprops,nrad);
  if(nl != hc->nradp2){
    HC_ERROR("hc_polsol","nl not equal to nrad+2 at end of solution loop");
  }

  if(compute_geoid){
    //
    //    Calculating geoid coefficients. The factor gf comes from
    //    * u(5) is in units of r0 * pi/180 / Ma (about 111 km/Ma) 
    //    * normalizing density is presumably 1 g/cm**3 = 1000 kg / m**3
    //    * geoid is in units of meters
    //    
    if(verbose > 1)
      fprintf(stderr,"hc_polsol: evaluating geoid%s\n",
	      (compute_geoid == 1)?(" at surface"):(", all layers"));
    /* 
       select geoid solution 
    */
    n6 = 4;
    //n6 = -iformat-1;

    /* first coefficients are zero  */
    clm[0] = clm[1] = 0.0;
    switch(compute_geoid){
    case 1:
      g1 = hc->nrad+1;g2=hc->nradp2;	/* only surface */
      break;
    case 2:
      g1 = 0;g2=hc->nradp2;		/* all layers */
      break;
    default:
      fprintf(stderr,"hc_polsol: error, geoid = %i undefined\n",compute_geoid);
      exit(-1);
    }

    for(gic=0,gi=g1;gi < g2;gi++,gic++){			  /* depth loop */
      /* 
	 first coefficients 
      */
      sh_write_coeff((geoid+gic),0,0,0,FALSE,clm); /* 0,0 */
      sh_write_coeff((geoid+gic),1,0,0,FALSE,clm); /* 1,0 */
      sh_write_coeff((geoid+gic),1,1,2,FALSE,clm); /* 1,1 */

      os = gi * 6 + n6;	/* select component */
      for(l=2;l <= pol_sol[0].lmax;l++){
	mmax = (calc_kernel_only)?(0):(l);
	for(m=0;m <= mmax;m++){	/* will typically be <= l, but only 0
				   for kernel computation */
	  if (m != 0){
	    sh_get_coeff((pol_sol+os),l,m,2,FALSE,clm); /* internal convention */
	    clm[0] *= hc->psp.geoid_factor;
	    clm[1] *= hc->psp.geoid_factor;
	    sh_write_coeff((geoid+gic),l,m,2,FALSE,clm);
	  }else{			/* m == 0 */
	    sh_get_coeff((pol_sol+os),l,m,0,FALSE,clm);
	    clm[0] *= hc->psp.geoid_factor;
	    sh_write_coeff((geoid+gic),l,m,0,FALSE,clm);
	  }
	}
      }
    }
    if(verbose > 1)
      fprintf(stderr,"hc_polsol: assigned geoid\n");
  } /* end geoid */
  
  /* 
     free the local arrays 
  */
  free(b);free(u3);
  if(!save_prop_mats){		
    /* 
       destroy individual propagator matrices, if we don't want to
       keep them
    */
    free(hc->props);free(hc->ppots);
  }
  /* all others should be saved */
  hc->psp.ncalled++;
  if(verbose)
    fprintf(stderr,"hc_polsol: done\n");
}

