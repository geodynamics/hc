#include "hc.h"
/* 
   general routines dealing with the hager & Connell implementation

  
   $Id: hc_init.c,v 1.14 2006/03/21 08:07:18 becker Exp becker $

*/

/* 

initialize the basic operational modes
 */
void hc_init_parameters(struct hc_parameters *p)
{
  /* 
     
  operational modes and parameters

  */
  p->compressible = FALSE;		/* compressibility following Panasyuk
				   & Steinberger */
  /* surface mechanical boundary condition */
  p->free_slip = TRUE;		/* free slip? */
  p->no_slip = FALSE;		/* no slip boundary condition? */
  p->platebc = FALSE;		/* plate velocities? */
  p->compute_geoid = TRUE;	/* compute the geoid?  */
  p->compute_geoid_correlations = FALSE;	/* compute the geoid
						   correlation with
						   refernece only   */
  p->dens_anom_scale = HC_D_LOG_V_D_LOG_D ;	/* default density anomaly scaling to
						   go from PREM percent traveltime
						   anomalies to density anomalies */
  p->scale_dens_anom_with_prem = TRUE; /* scale the input file
					  relative density anomalies
					  with the absolute rho value
					  of PREM at that layer. if
					  set to FALSE, will use the
					  average rho value  */
  p->read_short_dens_sh = FALSE; /* read the density anomaly file in
				    the short format of Becker &
				    Boschi (2002)?  */

  
  p->verbose = 0;		/* debugging output? (0,1,2,3,4...) */
  p->sol_binary_out = TRUE;	/* binary or ASCII output of SH expansion */
  p->solution_mode = HC_VEL;	/* default: velocity */
  p->solver_mode = HC_SOLVER_MODE_DEFAULT ;
  
  p->print_pt_sol = FALSE;
  p->print_spatial = FALSE;	/* by default, only print the spectral solution */
  /* for four layer approaches */
  p->rlayer[0] = HC_ND_RADIUS(660);
  p->rlayer[1] = HC_ND_RADIUS(410);
  p->rlayer[2] = HC_ND_RADIUS(100);
  /* default viscosities for the four layers, in units of reference viscosity*/
  p->elayer[0] = 50.; p->elayer[1] = 1.; p->elayer[2] = 0.1; p->elayer[3] = 50.;
  p->visc_init_mode = HC_INIT_E_FROM_FILE; /* by default, read viscosity from file */
  

  /* 
     depth dependent scaling of density files?
  */
  p->dd_dens_scale = FALSE; /* read a depth-dependent scaling from file? */
  p->ndf = 0;
  p->rdf = p->sdf = NULL;
  /* 

  filenames
  
  */
  strncpy(p->visc_filename,HC_VISC_FILE,HC_CHAR_LENGTH);
  strncpy(p->dens_filename,HC_DENS_SH_FILE,HC_CHAR_LENGTH);
  strncpy(p->prem_model_filename,PREM_MODEL_FILE,HC_CHAR_LENGTH);
}

/* 

   first, call this routine with a blank **hc 
   
*/
void hc_struc_init(struct hcs **hc)
{
  /* this will take care of all flags and such */
  *hc = (struct hcs *)calloc(1,sizeof(struct hcs ));
  if(!(*hc))
    HC_MEMERROR("hc_struc_init: hc");
  /* just to be sure */
  (*hc)->initialized = (*hc)->const_init = (*hc)->visc_init = 
    (*hc)->dens_init = (*hc)->pvel_init = (*hc)->orig_danom_saved = FALSE;
  hc_init_polsol_struct(&((*hc)->psp));
  /* 
     assign NULL pointers to allow reallocating 
  */
  (*hc)->r = (*hc)->visc = (*hc)->rvisc = 
    (*hc)->dfact = (*hc)->rden = (*hc)->dvisc = NULL;
  (*hc)->rpb = (*hc)->fpb= NULL;
  (*hc)->dens_anom = NULL; /* expansions */
  (*hc)->plm = NULL;
  (*hc)->prem_init = FALSE;
}

void hc_init_polsol_struct(struct hc_ps *psp)
{
  psp->ncalled = 0;
   /* scaling factors will only be computed once */
  psp->rho_scale = 1.0;
  psp->rho_init = FALSE;
  psp->prop_params_init = FALSE; 	/* parameters for propagator computation */
  psp->abg_init = FALSE;		/* alpha, beta factors */
  psp->prop_mats_init = FALSE;	/* will be true only if save_prop_mats is  */
  psp->tor_init = psp->pol_init = FALSE;
}
/* 

initialize all variables, after initializing the parameters 


INPUT: sh_type: type of expansion storage/spherical haronics scheme to use

INPUT: hc_parameters: holds all the settings


OUTPUT: hc structure, gets modified
*/
void hc_init_main(struct hcs *hc,int sh_type,
		  struct hc_parameters *p)
{
  int dummy=0;
  HC_PREC dd_dummy[4]={1,1,1,1};
  /* mechanical boundary condition */
  if(p->free_slip){
    if(p->no_slip || p->platebc)
      HC_ERROR("hc_init_main","free slip and no slip does not make sense");
    hc->free_slip = TRUE;
  }
  /* 
     set the default expansion type, input expansions will be 
     converted 
  */
  hc->sh_type = sh_type;
  /* 
     start by reading in physical constants and PREM model
  */
  hc_init_constants(hc,p->dens_anom_scale,p->prem_model_filename,p->verbose);
  /* 
     initialize viscosity structure from file
  */
  hc_assign_viscosity(hc,p->visc_init_mode,p->elayer,p);
  /* 

     initialize possible depth dependent scaling of density model
  */
  hc_assign_dd_scaling(HC_INIT_DD_FROM_FILE,dd_dummy,p,hc->r_cmb);


  if(!p->dd_dens_scale)
    if(p->verbose)
      fprintf(stderr,"hc_init_main: using constant dln\\rho/dln input density scaling of %g\n",
	      hc->dens_scale);
  
  if(p->no_slip && (!p->platebc)){
    /* 

       no slip (zero velocity) surface boundary conditions 

    */
    if(p->free_slip)
      HC_ERROR("hc_init","no slip and free_slip doesn't make sense");
    if(p->verbose)
      fprintf(stderr,"hc_init: initializing for no slip\n");

    /* 
       read in the densities first to determine L from the density expansion
    */
    hc_assign_density(hc,p->compressible,HC_INIT_D_FROM_FILE,
		      p->dens_filename,-1,FALSE,FALSE,p->scale_dens_anom_with_prem,
		      p->verbose,p->read_short_dens_sh,
		      p->dd_dens_scale,p->ndf,p->rdf,p->sdf,
		      (p->solver_mode == HC_SOLVER_MODE_VISC_SCAN)?(TRUE):(FALSE));
    /* 
       assign all zeroes up to the lmax of the density expansion 
    */
    hc_assign_plate_velocities(hc,HC_INIT_P_FROM_FILE,p->pvel_filename,TRUE,hc->dens_anom[0].lmax,
			       FALSE,p->verbose);
  }else if(p->platebc){
    /* 

       surface velocities 

    */
    if(p->free_slip)
      HC_ERROR("hc_init","plate boundary condition and free_slip doesn't make sense");
    if(p->verbose)
      fprintf(stderr,"hc_init: initializing for surface velocities\n");
    /* 
       read in velocities, which will determine the solution lmax 
    */
    hc_assign_plate_velocities(hc,HC_INIT_P_FROM_FILE,p->pvel_filename,FALSE,dummy,FALSE,p->verbose);
    /* then read in the density anomalies */
    hc_assign_density(hc,p->compressible,HC_INIT_D_FROM_FILE,p->dens_filename,hc->pvel[0].lmax,
		      FALSE,FALSE,p->scale_dens_anom_with_prem,
		      p->verbose,p->read_short_dens_sh, p->dd_dens_scale,p->ndf,p->rdf,p->sdf,
		      (p->solver_mode == HC_SOLVER_MODE_VISC_SCAN)?(TRUE):(FALSE));
  }else if(p->free_slip){
    /* 
       
       free slip

    */
    if(p->no_slip)
      HC_ERROR("hc_init","no slip and free slip does not make sense");
    if(p->verbose)
      fprintf(stderr,"hc_init: initializing for free-slip\n");
    /* read in the density fields */
    hc_assign_density(hc,p->compressible,HC_INIT_D_FROM_FILE,p->dens_filename,-1,FALSE,FALSE,
		      p->scale_dens_anom_with_prem,
		      p->verbose,p->read_short_dens_sh, p->dd_dens_scale,p->ndf,p->rdf,p->sdf,
		      (p->solver_mode == HC_SOLVER_MODE_VISC_SCAN)?(TRUE):(FALSE));
  }else{
    HC_ERROR("hc_init","boundary condition logic error");
  }

  /* 
     phase boundaries, if any 
  */
  hc_init_phase_boundaries(hc,0,p->verbose);
  /*  */
  hc->save_solution = TRUE;	/* (don')t save the propagator
				   matrices in hc_polsol and the
				   poloidal/toroidal solutions
				*/
  hc->initialized = TRUE;
}
/* 

   some of those numbers might be a bit funny, but leave them like
   this for now for backward compatibility.

*/
void hc_init_constants(struct hcs *hc, HC_PREC dens_anom_scale,
		       char *prem_filename,hc_boolean verbose)
{
  int ec;
  if(hc->const_init)
    HC_ERROR("hc_init_constants","why would you call this routine twice?")
  if(!hc->prem_init){
    /* PREM constants */
    if((ec=prem_read_model(prem_filename,hc->prem,verbose))){
      fprintf(stderr,"hc_init_constants: error: could not init PREM, error code: %i\n",
	      ec);
      exit(-1);
    }
    hc->prem_init = TRUE;
  }
  /* 
     density scale 
  */
  hc->dens_scale = dens_anom_scale;
  /* 
     constants
  */
  hc->timesc = HC_TIMESCALE_YR;		/* timescale [yr], like 1e6 yrs */
  hc->visnor = HC_VISNOR;		/* normalizing viscosity [Pas]*/
  hc->gacc =  HC_GACC; 		/* gravitational acceleration [cm/s2]*/
  hc->g =  HC_CAPITAL_G;		/* gravitational constant [Nm2/kg2]*/

  /*  

  radius of Earth in [m]

  */
  hc->re = hc->prem->r0;
  if(fabs(hc->re - (HC_RE_KM * 1e3)) > 1e-7)
    HC_ERROR("hc_init_constants","Earth radius mismatch")

  hc->secyr = HC_SECYR;	/* seconds/year  */

  /* 
     those are in g/cm^3
  */
  hc->avg_den_mantle =  HC_AVG_DEN_MANTLE;
  hc->avg_den_core = HC_AVG_DEN_CORE;

  /* 
     take the CMB radius from the Earth model 
  */
  hc->r_cmb = hc->prem->r[1];
  if(fabs(hc->r_cmb - 0.55) > 0.02)
    HC_ERROR("hc_init_constants","Earth model CMB radius appears off");

  /* 

  derived quantities
  
  */
  /* 
  
  velocity scale if input is in [cm/yr], works out to be ~0.11 

  */
  hc->vel_scale = hc->re*PIOVERONEEIGHTY/hc->timesc/HC_VEL_IO_SCALE;
  /* 
  
  stress scaling, will later be divided by non-dim radius, to go 
  to MPa
  
  */
  hc->stress_scale = (PIOVERONEEIGHTY * hc->visnor / hc->secyr)/
    (hc->timesc * HC_TIMESCALE_YR);
  

  hc->const_init = TRUE;
}

/* 
   
     handle command line  parameters
     
     visc_filename[] needs to be [HC_CHAR_LENGTH]

 */
void hc_handle_command_line(int argc, char **argv,
			    struct hc_parameters *p)
{
  int i;
  for(i=1;i < argc;i++){
    if(strcmp(argv[i],"-help")==0 || strcmp(argv[i],"-h")==0 || strcmp(argv[i],"-?")==0){// help
      /* 
	 help page
      */
      fprintf(stderr,"%s - perform Hager & O'Connell flow computation\n\n",
	      argv[0]);
      fprintf(stderr,"This code can compute velocities, tractions, and geoid for simple density distributions\n");
      fprintf(stderr,"and plate velocities using the semi-analytical approach of Hager & O'Connell (1981).\n");
      fprintf(stderr,"This particular implementation illustrates one possible way to combine the HC solver routines.\n");
      fprintf(stderr,"Based on code by Brad Hager, Richard O'Connell, and Bernhard Steinberger.\n");
      fprintf(stderr,"This version by Thorsten Becker and Craig O'Neill\n\n");
      
      fprintf(stderr,"density anomaly options:\n");
      fprintf(stderr,"-dens\tname\tuse name as a SH density anomaly model (%s)\n",
	      p->dens_filename);
      fprintf(stderr,"\t\tAll density anomalies are in units of %g%% of PREM, all SH coefficients\n\t\tin Dahlen & Tromp convention.\n",
	      HC_DENSITY_SCALING*100);
      
      fprintf(stderr,"-dshs\t\tuse the short, Becker & Boschi (2002) format for the SH density model (%s)\n",
	      hc_name_boolean(p->read_short_dens_sh));

      fprintf(stderr,"-ds\t\tdensity scaling factor (%g)\n",
	      p->dens_anom_scale);
      fprintf(stderr,"-dnp\t\tdo not scale density anomalies with PREM but rather mean density (%s)\n",
	      hc_name_boolean(!p->scale_dens_anom_with_prem));
      fprintf(stderr,"-dsf\tfile\tread depth dependent density scaling from file\n");
      fprintf(stderr,"\t\t(overrides -ds, %s), use pdens.py to edit\n\n",
	      hc_name_boolean(p->dd_dens_scale));
      fprintf(stderr,"Earth model options:\n");
      fprintf(stderr,"-prem\tname\tset Earth model to name (%s)\n",
	      p->prem_model_filename);
      fprintf(stderr,"-vf\tname\tviscosity structure filename (%s), use pvisc.py to edit\n",
	      p->visc_filename);
      fprintf(stderr,"\t\tThis file is in non_dim_radius viscosity[Pas] format\n");
      fprintf(stderr,"\t\tif name is \"visc_scan\", will loop through a four layer viscosity scan\n\n");
      fprintf(stderr,"boundary condition options:\n");
      fprintf(stderr,"-fs\t\tperform free slip computation (%s)\n",hc_name_boolean(p->free_slip));
      fprintf(stderr,"-ns\t\tperform no slip computation (%s)\n",hc_name_boolean(p->no_slip));
      fprintf(stderr,"-pvel\tname\tset prescribed surface velocities from file name (%s)\n",
	      hc_name_boolean(p->platebc));
      fprintf(stderr,"\t\tThe file (e.g. %s) is based on a DT expansion of cm/yr velocity fields.\n\n",HC_PVEL_FILE);

      fprintf(stderr,"solution procedure and I/O options:\n");
      fprintf(stderr,"-ng\t\tdo not compute and print the geoid (%s)\n",
	      hc_name_boolean(!p->compute_geoid));
      fprintf(stderr,"-rg\tname\tcompute correlation of computed geoid with that in file \"name\",\n");
      fprintf(stderr,"\t\tthis will not print out the geoid file, but only correlations (%s)\n",
	      hc_name_boolean(p->compute_geoid_correlations));
      fprintf(stderr,"-pptsol\t\tprint pol[6] and tor[2] solution vectors (%s)\n",
	      hc_name_boolean(p->print_pt_sol));
      fprintf(stderr,"-px\t\tprint the spatial solution to file (%s)\n",
	      hc_name_boolean(p->print_spatial));
      fprintf(stderr,"-rtrac\t\tcompute srr,srt,srp tractions [MPa] instead of velocities [cm/yr] (default: vel)\n");
      fprintf(stderr,"-htrac\t\tcompute stt,stp,spp tractions [MPa] instead of velocities [cm/yr] (default: vel)\n");
      fprintf(stderr,"-v\t-vv\t-vvv: verbosity levels (%i)\n",
	      (int)(p->verbose));
      fprintf(stderr,"\n\n");
      exit(-1);
    }else if(strcmp(argv[i],"-ds")==0){	/* density anomaly scaling factor */
      hc_advance_argument(&i,argc,argv);
      sscanf(argv[i],HC_FLT_FORMAT,&p->dens_anom_scale);
    }else if(strcmp(argv[i],"-fs")==0){	/* free slip flag */
      p->free_slip = TRUE;p->no_slip = FALSE;
    }else if(strcmp(argv[i],"-ns")==0){	/* no slip flag */
      p->no_slip = TRUE;p->free_slip = FALSE;
    }else if(strcmp(argv[i],"-px")==0){	/* print spatial solution? */
      hc_toggle_boolean(&p->print_spatial);
    }else if(strcmp(argv[i],"-dshs")==0){ /* use short format? */
      hc_toggle_boolean(&p->read_short_dens_sh);
    }else if(strcmp(argv[i],"-pptsol")==0){	/* print
						   poloidal/toroidal
						   solution
						   parameters */
      hc_toggle_boolean(&p->print_pt_sol);
    }else if(strcmp(argv[i],"-ng")==0){	/* do not compute geoid */
      hc_toggle_boolean(&p->compute_geoid);
    }else if(strcmp(argv[i],"-rg")==0){	/* compute geoid correlations */
      hc_toggle_boolean(&p->compute_geoid_correlations);
      p->compute_geoid = TRUE;
      hc_advance_argument(&i,argc,argv); /* filename */
      strncpy(p->ref_geoid_file,argv[i],HC_CHAR_LENGTH);
      hc_read_geoid(p);
    }else if(strcmp(argv[i],"-vf")==0){ /* viscosity filename */
      hc_advance_argument(&i,argc,argv);
      strncpy(p->visc_filename,argv[i],HC_CHAR_LENGTH);
      if(strcmp(p->visc_filename,"visc_scan")==0)
	{			/* inversion mode */
	  p->solver_mode = HC_SOLVER_MODE_VISC_SCAN;
	  p->visc_init_mode = HC_INIT_E_FOUR_LAYERS;
	}
    }else if(strcmp(argv[i],"-prem")==0){ /* PREM filename */
      hc_advance_argument(&i,argc,argv);
      strncpy(p->prem_model_filename,argv[i],HC_CHAR_LENGTH);
    }else if(strcmp(argv[i],"-dnp")==0){
      hc_toggle_boolean(&p->scale_dens_anom_with_prem);
    }else if(strcmp(argv[i],"-pvel")==0){ /* velocity filename, this will switch off free slip */
      hc_advance_argument(&i,argc,argv);
      strncpy(p->pvel_filename,argv[i],HC_CHAR_LENGTH);
      p->platebc = TRUE;p->no_slip = TRUE;p->free_slip = FALSE;
    }else if(strcmp(argv[i],"-dens")==0){ /* density filename */
      hc_advance_argument(&i,argc,argv);
      strncpy(p->dens_filename,argv[i],HC_CHAR_LENGTH);
    }else if(strcmp(argv[i],"-dsf")==0){ /* density scaling filename */
      hc_toggle_boolean(&p->dd_dens_scale);
      hc_advance_argument(&i,argc,argv);
      strncpy(p->dens_scaling_filename,argv[i],HC_CHAR_LENGTH);
    }else if(strcmp(argv[i],"-visc")==0){ /* viscosity filename */
      hc_advance_argument(&i,argc,argv);
      strncpy(p->visc_filename,argv[i],HC_CHAR_LENGTH);
    }else if(strcmp(argv[i],"-rtrac")==0){	/* compute radial
						   tractions */
      p->solution_mode = HC_RTRACTIONS;
    }else if(strcmp(argv[i],"-htrac")==0){	/* compute horizontal
						   tractions */
      p->solution_mode = HC_HTRACTIONS;
    }else if(strcmp(argv[i],"-v")==0){	/* verbosities */
      p->verbose = 1;
    }else if(strcmp(argv[i],"-vv")==0){	/* verbosities */
      p->verbose = 2;
    }else if(strcmp(argv[i],"-vvv")==0){	
      p->verbose = 3;
    }else{
      fprintf(stderr,"%s: can not use parameter %s, use -h for help page\n",
	      argv[0],argv[i]);
      exit(-1);
    }
  }
}
/* 

assign viscosity structure

mode == 0

read in a viscosity structure with layers of constant viscosity in
format

r visc 

where r is non-dim radius and visc non-dim viscosity, r has to be
ascending

*/
void hc_assign_viscosity(struct hcs *hc,int mode,
			 HC_PREC elayer[4],struct hc_parameters *p)
{
  FILE *in;
  int i;
  char fstring[100];
  HC_PREC mean,mweight,rold,mws;
  switch(mode){
  case HC_INIT_E_FOUR_LAYERS:
    /* initialize a four layer viscosity structure, viscosity values
       for 2871-660, 660-410, 410-100,100-0 should be given in units
       of visnor [1e21] as elayer[4]
     */
    hc_vecrealloc(&hc->rvisc,4,"hc_assign_viscosity");
    hc_vecrealloc(&hc->visc,4,"hc_assign_viscosity");
    /* number of layers */
    hc->nvis = 4;
    /* radii */
    hc->rvisc[0] = hc->r_cmb;
    hc->rvisc[1] = p->rlayer[0];
    hc->rvisc[2] = p->rlayer[1];
    hc->rvisc[3] = p->rlayer[2];

    for(i=0;i < hc->nvis;i++){
      hc->visc[i] = elayer[i];
      //fprintf(stderr,"%11g %11g\n",hc->rvisc[i],hc->visc[i]);
    }
    if(p->verbose)
      fprintf(stderr,"hc_assign_viscosity: assigned four layer viscosity: %.2e %.2e %.2e %.2e\n",
	      hc->visc[0],hc->visc[1],hc->visc[2],hc->visc[3]);
  break;
  case HC_INIT_E_FROM_FILE:		
    /* 
       
       init from file part 
    
    */
    if(hc->visc_init)
      HC_ERROR("hc_assign_viscosity","viscosity already read from file, really read again?");
    /* 
       read viscosity structure from file 

       format:

       r[non-dim] visc[non-dim]

       from bottom to top
    */
    in = ggrd_open(p->visc_filename,"r","hc_assign_viscosity");
    hc_vecrealloc(&hc->rvisc,1,"hc_assign_viscosity");
    hc_vecrealloc(&hc->visc,1,"hc_assign_viscosity");
    hc->nvis = 0;mean = 0.0;mws = 0.0;
    /* read sscanf string */
    hc_get_flt_frmt_string(fstring,2,FALSE);
    rold = hc->r_cmb;
    /* start read loop  */
    while(fscanf(in,"%lf %lf",(hc->rvisc+hc->nvis),(hc->visc+hc->nvis))==2){
      if(hc->visc[hc->nvis] < 1e15)
	fprintf(stderr,"hc_assign_viscosity: WARNING: expecting viscosities in Pas, read %g at layer %i\n",
		hc->visc[hc->nvis],hc->nvis);
      /* normalize viscosity here */
      hc->visc[hc->nvis] /= hc->visnor;
      if(hc->nvis == 0)
	if( hc->rvisc[hc->nvis] < hc->r_cmb-0.01){
	  fprintf(stderr,"hc_assign_viscosity: error: first radius %g is below CMB, %g\n",
		  hc->rvisc[hc->nvis], hc->r_cmb);
	  exit(-1);
	}
      if(p->verbose){
	/* weighted mean, should use volume, really, but this is just
	   for information  */
	mweight = ( hc->rvisc[hc->nvis] - rold); 
	mws += mweight;
	rold = hc->rvisc[hc->nvis];
	mean += log(hc->visc[hc->nvis]) * mweight;
      }
      if(hc->nvis){
	if(hc->rvisc[hc->nvis] < hc->rvisc[hc->nvis-1]){
	  fprintf(stderr,"hc_assign_viscosity: error: radius has to be ascing, entry %i (%g) smaller than last (%g)\n",
		  hc->nvis+1,hc->rvisc[hc->nvis],hc->rvisc[hc->nvis-1]);
	  exit(-1);
	}
      }
      hc->nvis++;
      hc_vecrealloc(&hc->rvisc,hc->nvis+1,"hc_assign_viscosity");
      hc_vecrealloc(&hc->visc,hc->nvis+1,"hc_assign_viscosity");
    }
    fclose(in);
    if(hc->rvisc[hc->nvis-1] > 1.0){
      fprintf(stderr,"hc_assign_viscosity: error: first last %g is above surface, 1.0\n",
	      hc->rvisc[hc->nvis-1]);
      exit(-1);
    }
    if(p->verbose){
      /* last entry */
      mweight = ( 1.0 - hc->rvisc[hc->nvis-1]); 
      mws += mweight;
      rold = hc->rvisc[hc->nvis-1];
      mean += log(hc->visc[hc->nvis-1]) * mweight;
      
      mean = exp(mean/mws);
      fprintf(stderr,"hc_assign_viscosity: read %i layered viscosity[Pas] from %s\n",
	      hc->nvis,p->visc_filename);
      fprintf(stderr,"hc_assign_viscosity: rough estimate of mean viscosity: %g x %g = %g Pas\n",
	      mean, hc->visnor, mean*hc->visnor);
    }
    break;
  default:
    HC_ERROR("hc_assign_viscosity","mode undefined");
    break;
  }
  hc->visc_init = TRUE;
}
/* 

assign/initialize the density anomalies and density factors

if mode==0: expects spherical harmonics of density anomalies [%] with
            respect to the 1-D reference model (PREM) given in SH
            format on decreasing depth levels in [km]
	    

	    spherical harmonics are real, fully normalized as in 
	    Dahlen & Tromp p. 859


this routine assigns the inho density radii, and the total (nrad=inho)+2
radii 

furthermore, the dfact factors are assigned as well

set  density_in_binary to TRUE, if expansion given in binary

nominal_lmax: -1: the max order of the density expansion will either
                  determine the lmax of the solution (free-slip, or vel_bc_zero) or 
		  will have to be the same as the plate expansion lmax (!free_slip && !vel_bc_zero)
              else: will zero out all entries > nominal_lmax

*/
void hc_assign_density(struct hcs *hc,
		       hc_boolean compressible,int mode, /* compressible computation? assignment mode? */
		       char *filename,int nominal_lmax,	 /* input density file name */
		       hc_boolean layer_structure_changed,
		       hc_boolean density_in_binary,
		       hc_boolean scale_dens_anom_with_prem,
		       hc_boolean verbose,
		       hc_boolean use_short_format,
		       hc_boolean dd_dens_scale, /* depth dependent scaling ? */
		       int ndf,HC_PREC *rdf,HC_PREC *sdf, /* depth dependent scaling factors */
		       hc_boolean save_orig_danom) /* save the original density anomalies for later rescaling */
{
  FILE *in;
  int type,lmax,shps,ilayer,nset,ivec,i,j;
  HC_PREC *dtop,*dbot,zlabel,local_scale,dens_scale[1],rho0;
  hc_boolean reported = FALSE,read_on;
  double dtmp[3];
  hc->compressible = compressible;
  hc->inho = 0;
  if(hc->dens_init)			/* clear old expansions, if 
					   already initialized */
    sh_free_expansion(hc->dens_anom,hc->inho);
  /* get PREM model, if not initialized */
  if(!hc->prem_init)
    HC_ERROR("hc_assign_density","assign 1-D reference model (PREM) first");
  switch(mode){
  case HC_RESCALE_D:
    /* resuse old density model, apply new scaling */
    if(!hc->orig_danom_saved)
      HC_ERROR("hc_assign_density","trying to rescale original density anomaly model, but it was not saved");
    for(i=0;i<hc->inho;i++){
      sh_aexp_equals_bexp_coeff((hc->dens_anom+i),(hc->dens_anom_orig+i));
    }
    break;
  case HC_INIT_D_FROM_FILE:
    if(hc->dens_init)
      HC_ERROR("hc_assign_density","really read dens anomalies again from file?");
    /* 
       
    read in density anomalies in spherical harmonics format for
    different layers from file. 

    this assumes that we are reading in anomalies in percent
    
    */

    in = ggrd_open(filename,"r","hc_assign_density");
    if(verbose)
      fprintf(stderr,"hc_assign_density: reading density anomalies in [%g%%] from %s\n",
	      100*HC_DENSITY_SCALING,filename);
    hc->inho = 0;		/* counter for density layers */
    /* get one density expansion */
    hc_get_blank_expansions(&hc->dens_anom,1,0,
			    "hc_assign_density");
    /* 
       read all layers as spherical harmonics assuming real Dahlen &
       Tromp (physical) normalization, short format

    */
    if(use_short_format){
      if(verbose)
	fprintf(stderr,"hc_assign_density: using short format for density SH\n");
      fscanf(in,"%i",&nset);
      ilayer = -1;
    }else{
      if(verbose)
	fprintf(stderr,"hc_assign_density: using default SH format for density\n");
    }
    
    read_on = TRUE;
    while(read_on){
      if(use_short_format){
	/* short format I/O */
	i  = fscanf(in,"%lf",dtmp);zlabel = (HC_PREC)dtmp[0];
	i += fscanf(in,"%i",&lmax);
	read_on = (i == 2)?(TRUE):(FALSE);
	ivec = 0;shps = 1;type = HC_DEFAULT_INTERNAL_FORMAT;
	ilayer++;
      }else{
	read_on = sh_read_parameters_from_file(&type,&lmax,&shps,&ilayer, &nset,
					       &zlabel,&ivec,in,FALSE,density_in_binary,
					       verbose);
      }
      if(read_on){
	if((verbose)&&(!reported)){
	  if(nominal_lmax > lmax)
	    fprintf(stderr,"hc_assign_density: density lmax: %3i filling up to nominal lmax: %3i with zeroes\n",
		    lmax,nominal_lmax);
	  if(nominal_lmax != -1){
	    fprintf(stderr,"hc_assign_density: density lmax: %3i limiting to lmax: %3i\n",
		    lmax,nominal_lmax);
	  }else{
	    fprintf(stderr,"hc_assign_density: density lmax: %3i determines solution lmax\n",
		    lmax);
	  }
	  reported = TRUE;
	  if(verbose >= 2)
	    fprintf(stderr,"hc_assign_density: non_dim radius                 %% factor    PREM \\rho/mean_rho          layer #            depth[km]\n");
	}

	/* 
	   do tests 
	*/
	if((shps != 1)||(ivec))
	  HC_ERROR("hc_assign_density","vector field read in but only scalar expansion expected");
	/* test and assign depth levels */
	hc->rden=(HC_PREC *)
	  realloc(hc->rden,(1+hc->inho)*sizeof(HC_PREC));
	if(!hc->rden)
	  HC_MEMERROR("hc_assign_density: rden");
	/* 
	   assign depth, this assumes that we are reading in depths [km]
	*/
	hc->rden[hc->inho] = HC_ND_RADIUS((double)zlabel);
	if(scale_dens_anom_with_prem){
	  /* 
	     
	     get reference density at this level
	     
	  */
	  prem_get_rho(&rho0,hc->rden[hc->inho],hc->prem);
	  rho0 /= 1000.0;
	  if(rho0 < 3)
	    fprintf(stderr,"\nhc_assign_density: WARNING: using small (%g) density from PREM for layer at depth %g\n\n",
		    rho0*1000,HC_Z_DEPTH(hc->rden[hc->inho]));
	}else{
	  /* mean value */
	  rho0 =  hc->avg_den_mantle;
	}
	/* 
	   density anomaly
	*/
	/* scaling factor without depth dependence */
	dens_scale[0] = ((HC_PREC)HC_DENSITY_SCALING ) * rho0;
	if(verbose >= 2){
	  
	  fprintf(stderr,"hc_assign_density: r: %11g anom scales: %11g x %11g = %11g\t%5i out of %i, z: %11g\n",
		  hc->rden[hc->inho],
		  HC_DENSITY_SCALING,rho0/ hc->avg_den_mantle,dens_scale[0],hc->inho+1,nset,zlabel);
	}
	if(hc->inho){	
	  /* 
	     check by comparison with previous expansion 
	  */
	  if(nominal_lmax == -1)
	    if(lmax != hc->dens_anom[0].lmax)
	      HC_ERROR("hc_assign_density","lmax changed in file");
	  if(hc->rden[hc->inho] <= hc->rden[hc->inho-1]){
	    fprintf(stderr,"hc_assign_density: %i %g %g\n",hc->inho,
		    hc->rden[hc->inho], hc->rden[hc->inho-1]);
	    HC_ERROR("hc_assign_density","depth should decrease, radius increase (give z[km])");
	  }
	}
	/* 
	   make room for new expansion 
	*/
	hc_get_blank_expansions(&hc->dens_anom,hc->inho+1,hc->inho,
				"hc_assign_density");
	/* 
	   initialize expansion on irregular grid
	*/
	sh_init_expansion((hc->dens_anom+hc->inho),
			  (nominal_lmax > lmax) ? (nominal_lmax):(lmax),
			  hc->sh_type,0,verbose,FALSE);
	/* 
	   
	read parameters and scale (put possible depth dependence of
	scaling here)
	
	will assume input is in physical convention
	
	*/
	sh_read_coefficients_from_file((hc->dens_anom+hc->inho),1,lmax,in,density_in_binary,
				       dens_scale,verbose);
	hc->inho++;
      }	/* end actualy read on */
    } /* end while */
    
    if(hc->inho != nset)
      HC_ERROR("hc_assign_density","file mode: mismatch in number of layers");
    fclose(in);
    break;
  default:
    HC_ERROR("hc_assign_density","mode undefined");
    break;
  }
  if(save_orig_danom && (!hc->dens_init)){
    /* make a copy of the original density anomaly before applying
       depth dependent scaling, only done once per run */
    hc_get_blank_expansions(&hc->dens_anom_orig,hc->inho+1,hc->inho,
			    "hc_assign_density");
    for(i=0;i<hc->inho;i++){
      sh_init_expansion((hc->dens_anom_orig+i),hc->dens_anom[0].lmax,
			hc->sh_type,0,FALSE,FALSE);
      sh_aexp_equals_bexp_coeff((hc->dens_anom_orig+i),(hc->dens_anom+i));
    }
    hc->orig_danom_saved=TRUE;
  }
  /* 

     scale with possibly depth dependent factor
     
  */
  for(i=0;i < hc->inho;i++){
    /* depth dependent factor? */
    local_scale = hc_find_dens_scale(hc->rden[i],hc->dens_scale,dd_dens_scale,rdf,sdf,ndf);
    sh_scale_expansion((hc->dens_anom+i),local_scale);
    if(verbose >= 2){
      fprintf(stderr,"hc_assign_density: r: %11g additional %s d\\rho/dinput: %11g \tlayer %5i out of %i\n",
	      hc->rden[i],(dd_dens_scale)?("depth-dependent"):("constant"),local_scale,i,hc->inho);
    }
  }

  if((!hc->dens_init)||(layer_structure_changed)){
    /* 
       
    assign the other radii, nrad + 2
    
    */
    hc->nrad = hc->inho;
    hc->nradp2 = hc->nrad + 2;
    hc_vecrealloc(&hc->r,hc->nradp2,"hc_assign_density");
    /* 
       viscosity at density layers for horizontal stress
       computations */
    hc_vecrealloc(&hc->dvisc,hc->nradp2,"hc_assign_density");

    hc->r[0] = hc->r_cmb;	/* CMB  */
    if(hc->rden[0] <= hc->r[0])
      HC_ERROR("hc_assign_density","first density layer has to be above internal CMD limit");
    for(i=0;i<hc->nrad;i++)	/* density layers */
      hc->r[i+1] = hc->rden[i];
    if(hc->rden[hc->nrad-1] >= 1.0)
      HC_ERROR("hc_assign_density","uppermost density layer has to be below surface");
    hc->r[hc->nrad+1] = 1.0;	/* surface */
    /* 

    assign viscosity at density layers
    
    */
    hc->dvisc[0] = hc->visc[0];
    for(i=1;i < hc->nradp2;i++){
      for(j=hc->nvis-1;j>=0;j--)
	if(hc->rvisc[j] < hc->r[i-1])
	  break;
      hc->dvisc[i] = hc->visc[j];
    }
    /* 

    assign the density jump factors

    */
    /* 
       since we have spherical harmonics at several layers, we assign 
       the layer thickness by picking the intermediate depths
    */
    hc_vecalloc(&dbot,hc->nrad,"hc_assign_density");
    hc_vecalloc(&dtop,hc->nrad,"hc_assign_density");
    //    top boundaries
    j = hc->nrad-1;
    for(i=0;i < j;i++)
      dtop[i] = 1.0 - (hc->rden[i+1] + hc->rden[i])/2.0;
    dtop[j] = 0.0; // top boundary
    //    bottom boundaries
    dbot[0] = 1.0 - hc->r_cmb;  // bottom boundary, ie. CMB 
    for(i=1;i < hc->nrad;i++)
      dbot[i] = dtop[i-1];
    /* 
       density layer thickness factors
    */
    hc_dvecrealloc(&hc->dfact,hc->nrad,"hc_assign_density");
    for(i=0;i<hc->nrad;i++){
      hc->dfact[i] = 1.0/hc->rden[i] *(dbot[i] - dtop[i]);
    }
    if(verbose)
      for(i=0;i < hc->nrad;i++)
	fprintf(stderr,"hc_assign_density: dens %3i: r: %8.6f df: %8.6f |rho|: %8.4f\n",
		i+1,hc->rden[i],hc->dfact[i],
		sqrt(sh_total_power((hc->dens_anom+i))));
    free(dbot);free(dtop);
  } /* end layer structure part */

  

  hc->dens_init = TRUE;
}
/* 

find depth dependent scaling

*/
HC_PREC hc_find_dens_scale(HC_PREC r, HC_PREC s0,hc_boolean depth_dependent, HC_PREC *rd,HC_PREC *sd,int n)
{
  int i;
  if(depth_dependent){
    i=0;
    while((i<n) && (rd[i] < r))
      i++;
    i--;
    return sd[i];
  }else{
    return s0;
  }
}

/* 

assign phase boundary jumps
input:

npb: number of phase boundaries

....



*/
void hc_init_phase_boundaries(struct hcs *hc, int npb,
			      hc_boolean verbose)
{

  hc->npb = npb;		/* no phase boundaries for now */
  if(hc->npb){
    HC_ERROR("hc_init_phase_boundaries","phase boundaries not implemented yet");
    hc_vecrealloc(&hc->rpb,hc->npb,"hc_init_phase_boundaries");
    hc_vecrealloc(&hc->fpb,hc->npb,"hc_init_phase_boundaries");
  }

}

/* 

read in plate velocities, 

vel_bc_zero: if true, will set all surface velocities to zero

lmax will only be referenced if all velocities are supposed to be set to zero

we expect velocities to be in cm/yr, convert to m/yr

*/

void hc_assign_plate_velocities(struct hcs *hc,int mode, char *filename,hc_boolean vel_bc_zero,
				int lmax,hc_boolean pvel_in_binary,hc_boolean verbose)
{
  int type,shps,ilayer,nset,ivec;
  HC_PREC zlabel,vfac[2],t10[2],t11[2];
  FILE *in;
  /* scale to go from cm/yr to internal scale */
  vfac[0] = vfac[1] = 1.0/hc->vel_scale;
  if(hc->pvel_init)
    HC_ERROR("hc_assign_plate_velocities","what to do if called twice?");
  if(!vel_bc_zero){
    /* 

    velocities are NOT all zero


    */
    switch(mode){
    case HC_INIT_P_FROM_FILE:
      /* 
	 read velocities in pol/tor expansion format from file in
	 units of HC_VELOCITY_FILE_FACTOR per year, short format
      */
      if(verbose)
	fprintf(stderr,"hc_assign_plate_velocities: expecting [cm/yr] pol/tor from %s\n",
		filename);
      in = ggrd_open(filename,"r","hc_assign_plate_velocities");
      if(!sh_read_parameters_from_file(&type,&lmax,&shps,&ilayer, &nset,
				       &zlabel,&ivec,in,FALSE,
				       pvel_in_binary,verbose)){
	fprintf(stderr,"hc_assign_plate_velocities: read error file %s\n",
		filename);
	exit(-1);
      } /* check if we read in two sets of expansions */
      if(shps != 2){
	fprintf(stderr,"hc_assign_plate_velocities: two sets expected but found shps: %i in file %s\n",
		shps,filename);
	exit(-1);
      }
      if((nset > 1)||(fabs(zlabel) > 0.01)){
	fprintf(stderr,"hc_assign_plate_velocities: error: expected one layer at surface, but nset: %i z: %g\n",
		nset, zlabel);
	exit(-1);
      }
      /* 
	 initialize expansion using irregular grid
      */
      sh_init_expansion((hc->pvel+0),lmax,hc->sh_type,1,verbose, FALSE);
      sh_init_expansion((hc->pvel+1),lmax,hc->sh_type,1,verbose, FALSE);
      /* 
	 read in expansions, convert to internal format from 
	 physical 
      */
      sh_read_coefficients_from_file(hc->pvel,shps,-1,in,pvel_in_binary,
				     vfac,verbose);
      fclose(in);
      /* 
	 scale by 1/sqrt(l(l+1))
      */
      if(hc->pvel[0].lmax > hc->lfac_init)
	hc_init_l_factors(hc,hc->pvel[0].lmax);
      sh_scale_expansion_l_factor((hc->pvel+0),hc->ilfac);
      sh_scale_expansion_l_factor((hc->pvel+1),hc->ilfac);
      /*  
	  check for net rotation
      */
      sh_get_coeff((hc->pvel+1),1,0,0,TRUE,t10);
      sh_get_coeff((hc->pvel+1),1,0,2,TRUE,t11);
      if(fabs(t10[0])+fabs(t11[0])+fabs(t11[1]) > 1.0e-7)
	fprintf(stderr,"\nhc_assign_plate_velocities: WARNING: toroidal A(1,0): %g A(1,1): %g B(1,1): %g\n\n",
		t10[0],t11[0],t11[1]);
      if(verbose)
	fprintf(stderr,"hc_assign_plate_velocities: read velocities, lmax %i: |pol|: %11g |tor|: %11g\n",
		lmax,sqrt(sh_total_power((hc->pvel+0))),sqrt(sh_total_power((hc->pvel+1))));
      break;
    default:
      HC_ERROR("hc_assign_plate_velocities","op mode undefined");
    }
  }else{
    /* 
       initialize with zeroes
    */
    if(hc->pvel_init){
      sh_clear_alm(hc->pvel);
      sh_clear_alm((hc->pvel+1));
    }else{
      /* use irregular grid */
      sh_init_expansion(hc->pvel,lmax,hc->sh_type, 
			1,verbose,FALSE);
      sh_init_expansion((hc->pvel+1),lmax,hc->sh_type, 
			1,verbose,FALSE);
    }
    if(verbose)
      fprintf(stderr,"hc_assign_plate_velocities: using no-slip surface BC, lmax %i\n",
	      lmax);
  }
  hc->pvel_init = TRUE;
}
 
/* 

initialize an array with sqrt(l(l+1)) factors
from l=0 .. lmax+1

pass lfac initialized (say, as NULL)

*/
void hc_init_l_factors(struct hcs *hc, int lmax)
{
  int lmaxp1,l;
  lmaxp1 = lmax + 1;
  hc_vecrealloc(&hc->lfac,lmaxp1,"hc_init_l_factors");
  hc_vecrealloc(&hc->ilfac,lmaxp1,"hc_init_l_factors");
  /* maybe optimize later */
  hc->lfac[0] = 0.0;
  hc->ilfac[0] = 1.0;		/* shouldn't matter */
  for(l=1;l < lmaxp1;l++){
    hc->lfac[l] = sqrt((HC_PREC)l * ((HC_PREC)l + 1.0));
    hc->ilfac[l] = 1.0/hc->lfac[l];
  }
  hc->lfac_init = lmax;
}

/* 
   
reallocate new spherical harmonic expansion structures 
and initialize them with zeroes

input: 
expansion: expansion **
nold: number of old expansions
nnew: number of new expansions

*/
void hc_get_blank_expansions(struct sh_lms **expansion,
			     int nnew,int nold,
			     char *calling_sub)
{
  struct sh_lms *tmpzero;
  int ngrow;
  ngrow= nnew - nold;
  if(ngrow <= 0){
    fprintf(stderr,"hc_get_blank_expansions: error: ngrow needs to be > 0, ngrow: %i\n",
	    ngrow);
    fprintf(stderr,"hc_get_blank_expansions: was called from %s\n",
	    calling_sub);
    exit(-1);
  }
  /* 
     reallocate space 
  */
  *expansion  = (struct sh_lms *)
    realloc(*expansion,sizeof(struct sh_lms)*nnew);
  if(!(*expansion)){
    fprintf(stderr,"hc_get_blank_expansions: memory error: ngrow: %i\n",
	    nnew-nold);
    fprintf(stderr,"hc_get_blank_expansions: was called from %s\n",
	    calling_sub);
    exit(-1);
  }
  /* zero out new space */
  tmpzero  = (struct sh_lms *)calloc(ngrow,sizeof(struct sh_lms));
  if(!tmpzero)
    HC_MEMERROR("hc_get_blank_expansions: tmpzero");
  /* copy zeroes over */
  memcpy((*expansion+nold),tmpzero,ngrow*sizeof(struct sh_lms));
  free(tmpzero);
}

/* 


be more careful with freeing


 */
void hc_struc_free(struct hcs **hc)
{
  free((*hc)->visc);
  free((*hc)->rvisc);
  free((*hc)->qwrite);

  sh_free_expansion((*hc)->dens_anom,1);
  
  free(*hc);
}

/*  

assign a depth dependent density scale

*/
void hc_assign_dd_scaling(int mode, HC_PREC dlayer[4],struct hc_parameters *p,
			  HC_PREC rcmb)
{
  HC_PREC smean;
  double dtmp[2];
  int i;
  FILE *in;
  if(p->dd_dens_scale){
    /* 
       depth depending scaling
    */
    switch(mode){
      case HC_INIT_DD_FROM_FILE:
	/* 

	   read from file, format same as for viscosity
	*/
	if(p->verbose)
	  fprintf(stderr,"hc_assign_dd_scaling: reading depth dependent  dln\\rho/dln density scaling from %s\n",
		  p->dens_scaling_filename);
	p->ndf=0;smean = 0.0;
	in = ggrd_open(p->dens_scaling_filename,"r","hc_assign_dd_scaling");
	while(fscanf(in,"%lf %lf",dtmp,(dtmp+1)) == 2){
	  hc_vecrealloc(&p->rdf,(1+p->ndf),"hc_assign_dd_scaling");
	  hc_vecrealloc(&p->sdf,(1+p->ndf),"hc_assign_dd_scaling");
	  p->rdf[p->ndf] = dtmp[0];p->sdf[p->ndf] = dtmp[1];
	  smean+=p->sdf[p->ndf];
	  p->ndf++;
	}
	fclose(in);
	if(!p->ndf){
	  fprintf(stderr,"hc_assign_dd_scaling: error: did not read any density scaling factors from %s\n",
		  p->dens_scaling_filename);
	  exit(-1);
	}
	smean /= (HC_PREC)p->ndf;
	if(p->verbose)
	  fprintf(stderr,"hc_assign_dd_density: read scaling on %i layers, rough mean: %g\n",p->ndf,smean);
	break;
	
    case HC_INIT_DD_FOUR_LAYERS:
      p->ndf = 4;
      hc_vecrealloc(&p->rdf,(p->ndf),"hc_assign_dd_scaling");
      hc_vecrealloc(&p->sdf,(p->ndf),"hc_assign_dd_scaling");
      p->rdf[0] = rcmb;p->rdf[1] = p->rlayer[0];p->rdf[2] = p->rlayer[1];p->rdf[3] = p->rlayer[2];
      for(i=0;i<4;i++)
	p->sdf[i] = dlayer[i];
      if(p->verbose)
	fprintf(stderr,"hc_assign_dd_density: assigned four layer density scaling factors %g %g %g %g\n",
		p->sdf[0],	p->sdf[1],	p->sdf[2],	p->sdf[3]);
      break;
    }
    /* end init */
  }else{
    p->ndf = 0;
    hc_vecrealloc(&p->rdf,(1+p->ndf),"hc_assign_dd_scaling");
    hc_vecrealloc(&p->sdf,(1+p->ndf),"hc_assign_dd_scaling");
  }
}

/* read in and assign reference geoid */
void hc_read_geoid(struct hc_parameters *p)
{
  int type,lmax,shps,ilayer,nset,ivec;
  double zlabel;
  FILE *in;
  HC_PREC fac[3]={1,1,1};
  
  in = fopen(p->ref_geoid_file,"r");
  if(!in){
    fprintf(stderr,"hc_read_geoid: error, could not open ref geoid \"%s\", expecting long scalar SH format\n",
	    p->ref_geoid_file);
    exit(-1);
  }
  
    
  /* read in the file */
  sh_read_parameters_from_file(&type,&lmax,&shps,&ilayer,&nset,
			       &zlabel,&ivec,in,FALSE,
			       FALSE,p->verbose);
  if((ivec != 0)||(shps!=1)){
    fprintf(stderr,"hc_read_geoid: error, expecting scalar, long format SH in %s\n",
	    p->ref_geoid_file);
    exit(-1);
  }
  sh_allocate_and_init(&p->ref_geoid,shps,lmax,type,ivec,p->verbose,0);
  sh_read_coefficients_from_file(p->ref_geoid,shps,-1,in,FALSE,fac,p->verbose);
  fclose(in);
  if(p->verbose)
    fprintf(stderr,"hc_read_geoid: read reference geoid from %s, L=%i\n",
	    p->ref_geoid_file,p->ref_geoid->lmax);
  
}

