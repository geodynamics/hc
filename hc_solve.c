#include "hc.h"
/* 

solve the poloidal and toroidal part of a Hager & O'Connell type flow
computation



free_slip: TRUE/FALSE. if false, will either use the plate motions or no-slip,
           depending on how pvel was initialized

solve_mode: solution mode, used for summing the solutions

compute_geoid: compute the geoid (this is done half-way anyway, and
                    should not add much computation time)

geoid: geoid expansion, needs to be initialized

dens_fac_changed: has the density anomaly expansion changed since the last call to 
                  hc_solve?

plate_vel_changed: have the plate motion expansions changed since the last call to
                   hc_solve?

viscosity_or_layer_changed: has the viscosity structure or the layer spacing 
                            of density anomalies changed since the last call to
                            hc_solve?

input/output:
           sol: [3*nradp2] expansions holding radial poloidal and toroidal 
	            components. lmax of expansion has to be >= lmax(plates)
	        has to be initialized before calling this routine


*/

void hc_solve(struct hcs *hc, hc_boolean free_slip, 
	      int solve_mode,
	      struct sh_lms *sol, 
	      hc_boolean dens_anom_changed,
	      hc_boolean plate_vel_changed,
	      hc_boolean viscosity_or_layer_changed,
	      hc_boolean print_pt_sol,
	      hc_boolean compute_geoid,
	      struct sh_lms *pvel, /* plate velocity expansion */
	      struct sh_lms *dens_anom,
	      struct sh_lms *geoid, /* geoid solution, needs to be init */
	      hc_boolean verbose)
{
  int nsh_pol,nsh_tor=0;
  static hc_boolean convert_to_dt = TRUE; /* convert the poloidal and
					     toroidal solution vectors
					     to physical SH convention
					     (if set to FALSE, can
					     compare with Benhard's
					     densub densol output */
  HC_PREC *tvec;

  if(!hc->initialized)
    HC_ERROR("hc_solve","hc structure not initialized");
  if((!free_slip) && (pvel[0].lmax < dens_anom[0].lmax)){
    fprintf(stderr,"hc_solve: error: plate expansion lmax (%i) has to be >= density lmax (%i)\n",
	    pvel[0].lmax,dens_anom[0].lmax);
    exit(-1);
  }
  if(sol[0].lmax < pvel[0].lmax){
    fprintf(stderr,"hc_solve: error: solution lmax (%i) has to be >= plate velocitiy lmax (%i)\n",
	    sol[0].lmax,pvel[0].lmax);
    exit(-1);
  }
  /* 
     POLOIDAL PART  
  */
  /* 
     initialize a bunch of expansions for the poloidal solution 
  */
  nsh_pol = 6 * (hc->nrad+2);	/* u[4] plus poten[2] */
  if((!hc->psp.pol_init)||(!hc->save_solution)){
    /* room for pol solution */
    sh_allocate_and_init(&hc->pol_sol,nsh_pol,
			 dens_anom[0].lmax,hc->sh_type,
			 0,verbose,FALSE); /* irregular grid */
  }
  if((!hc->save_solution) || (!hc->psp.pol_init) || viscosity_or_layer_changed ||
     dens_anom_changed || ((!free_slip) && (plate_vel_changed))){  
    /* 
       
    FIND POLOIDAL SOLUTION 

    if the density anomalies, the viscosity structure, or the plate
    velocities changed
    
    */
    hc_polsol(hc,hc->nrad,hc->r,hc->inho,hc->dfact,
	      viscosity_or_layer_changed,
	      dens_anom,hc->compressible,
	      hc->npb,hc->rpb,hc->fpb,free_slip,
	      (pvel+0),hc->pol_sol,
	      compute_geoid,geoid,hc->save_solution,
	      verbose);
    if(print_pt_sol)		/* print poloidal solution without the
				   scaling factors */
      hc_print_poloidal_solution(hc->pol_sol,hc,31, /* print only up
						       to lmax = 31 or
						       below */
				 HC_POLSOL_FILE,convert_to_dt,verbose);
  }
  if(!free_slip){
    /* 
       
       solve toroidal part only for no-slip surface boundary condition

    */
    if((!hc->psp.tor_init)||(!hc->save_solution)){
      nsh_tor = 2 * (hc->nrad+2);
      sh_allocate_and_init(&hc->tor_sol,nsh_tor,pvel[1].lmax,
			   hc->sh_type,0,verbose,FALSE); /* irregular grid */
    }
    if((!hc->psp.tor_init) || viscosity_or_layer_changed || plate_vel_changed || 
       (!hc->save_solution)){
      /* 
	 if we are not saving solutions, or the velocities or viscosities
	 have changed, we need to (re)compute the toroidal solution
      */
      /* make room for solution kernel */
      hc_vecalloc(&tvec,(hc->nrad+2)*(pvel[1].lmax+1)*2,
		  "hc_solve");
      /* compute kernels, and assign kernel*pvel to tor_sol */
      hc_torsol(hc,hc->nrad,hc->nvis,pvel[1].lmax,hc->r,
		&hc->rvisc,&hc->visc,(pvel+1),hc->tor_sol,tvec,
		verbose);
      if(print_pt_sol)
	hc_print_toroidal_solution(tvec,pvel[1].lmax,
				   hc,pvel[1].lmax,HC_TORSOL_FILE,
				   verbose);
      free(tvec);
    }
  }else{
    nsh_tor = 0;
  }
  switch(solve_mode){
  case HC_VEL:
    if(verbose)
      fprintf(stderr,"hc_solve: computing solution for velocities\n");
    break;
  case HC_RTRACTIONS:
    if(verbose)
      fprintf(stderr,"hc_solve: computing solution for radial tractions\n");
    break;
  case HC_HTRACTIONS:
    if(verbose)
      fprintf(stderr,"hc_solve: computing solution for horizontal tractions\n");
    break;
  default:
    fprintf(stderr,"hc_solve: error: solution mode %i undefined\n",
	    solve_mode);
    exit(-1);
    break;
  }
  /* 

     sum up the poloidal and torodial solutions and set the spectral
     init flag to true for solution expansion
     
  */
  hc_sum(hc,hc->nrad,hc->pol_sol,hc->tor_sol,solve_mode,free_slip,sol,
	 verbose);
  /* 
     free temporary arrays
  */
  if(!hc->save_solution){
    /* 
       POLOIDAL SOLUTION related expansions, those are not saved as they
       change with density anomalies and plate motions
    */
    sh_free_expansion(hc->pol_sol,nsh_pol);
    /* 
       toroidal, maybe save those, since they only depend on plate velocities
       and viscosities 
    */
    if(!free_slip)
      sh_free_expansion(hc->tor_sol,nsh_tor);
  }
  hc->psp.pol_init = TRUE;
  hc->psp.tor_init = TRUE;
  hc->spectral_solution_computed = TRUE;
}
/* 
   
computes the radial, poloidal, and toroidal solution expansions 
as sol[3*nradp2] for each layer

input:

pol_sol[6*nradp2]: y1...y6    (six) poloidal solutions for each layer
tor_sol[2*nradp2]: y9 and y10 (two) toroidal solutions for each layer


THESE SOLUTIONS WILL NEED TO BE SCALED WITH CONSTANTS AS GIVEN IN 

hc_compute_solution_scaling_factors

*/
void hc_sum(struct hcs *hc,
	    int nrad,struct sh_lms *pol_sol, struct sh_lms *tor_sol, 
	    int solve_mode,hc_boolean free_slip, 
	    struct sh_lms *sol, 
	    hc_boolean verbose)
{
  int itchoose,irchoose,ipchoose; /* indices for which solutions to use */
  int i,j,i3,i6;
  if(sol[0].lmax > hc->lfac_init)
    hc_init_l_factors(hc,sol[0].lmax);

  /* 

  pick the right components for the radial, poloidal, and toroidal
  solution
     
  */
  switch(solve_mode){
  case HC_VEL:
    //
    //    velocity output requested 
    //
    irchoose = 0; // y1 for radial
    ipchoose = 1; // y2 for poloidal
    itchoose = 0; // y9 for toroidal
    break;
  case HC_RTRACTIONS:
    //
    //    srr srt srp stress output requested 
    //
    irchoose = 2;// y3  for radial
    ipchoose = 3;// y4  for poloidal
    itchoose = 1;// y10 for toroidal 
    break;
  case HC_HTRACTIONS:
    fprintf(stderr,"hc_sum: horizontal tractions not implemented yet\n");
    exit(-1);
    break;
  default:
    HC_ERROR("hc_sum","solve mode undefined");
    break;
  }
  /* 


  for velocities, this summation is OK. for some other properties,
  might have to rescale by layer radius and such


  */
  for(i=i3=i6=0;i < hc->nradp2;i++,i3+=3,i6+=6){
    /* 
       radial part 
    */
    sh_aexp_equals_bexp_coeff((sol+i3+HC_RAD),(pol_sol+i6+irchoose));
    /* 
       poloidal part, need to scale with sqrt(l(l+1))
    */
    sh_aexp_equals_bexp_coeff((sol+i3+HC_POL),(pol_sol+i6+ipchoose));
    sh_scale_expansion_l_factor((sol+i3+HC_POL),hc->lfac);
    for(j=0;j<3;j++)
      sol[i3+j].spectral_init = TRUE;
    if(!free_slip){
      /* 
	 toroidal part, need to scale with sqrt(l(l+1))
      */
      sh_aexp_equals_bexp_coeff((sol+i3+HC_TOR),(tor_sol+i*2+itchoose));
      sh_scale_expansion_l_factor((sol+i3+HC_TOR),hc->lfac);
    }else{
      /* no toroidal part for free-slip */
      sh_clear_alm((sol+i3+2));
    }
  } /* end layer loop */
}


/* 

given a spherical harmonic solution, compute the spatial 
corresponding solution

sol[nradp2 * 3 ]

data has to be initialized, eg. as NULL
*/
void hc_compute_sol_spatial(struct hcs *hc, struct sh_lms *sol_w,
			    HC_PREC **sol_x, hc_boolean verbose)
{
  int i,i3,np,np2,np3,os;
  static int ntype = 3;
  np = sol_w[0].npoints;
  np2 = np * 2;
  np3 = np2 + np;	/* 
			   number of points per spatial 
			   expansions for r, pol, tor
			*/
  /* allocate space for spatial solution*/
  hc_vecrealloc(sol_x,np3*hc->nradp2,"sol_x");
  /* 
     compute the plm factors 
  */
  sh_compute_plm(sol_w,1,&hc->plm,verbose);
  for(i=i3=0;i < hc->nradp2;i++,i3 += ntype){
    os = i*np3;
    /* radial component */
    sh_compute_spatial((sol_w+i3+HC_RAD),0,TRUE,&hc->plm,
		       (*sol_x+os),verbose);
    os += np;
    /* poloidal/toroidal component */
    sh_compute_spatial((sol_w+i3+HC_POL),1,TRUE,&hc->plm,
		       (*sol_x+os),verbose);
  }
  hc->spatial_solution_computed = TRUE;
}
