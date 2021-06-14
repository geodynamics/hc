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
  p->compute_geoid = 1;	/* compute the geoid? 1: surface 2: all layers */
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

  p->read_short_pvel_sh = FALSE; /* read the plate velocities in short format */

  p->verbose = 0;		/* debugging output? (0,1,2,3,4...) */
  p->sol_binary_out = TRUE;	/* binary or ASCII output of SH expansion */
  p->solution_mode = HC_VEL;	/* default: velocity output */

  p->print_density_field = TRUE; /* print the scaled density field (useful for debugging) */


  p->print_kernel_only = FALSE;
  
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
  
  p->solver_kludge_l = INT_MAX;	/* default: no solver tricks */

  p->remove_nr = FALSE;		/* by default, leave NR in for plate motion models */
  /* 
     depth dependent scaling of density files?
  */
  p->dd_dens_scale = HC_DD_CONSTANT; /* can be HC_DD_CONSTANT, HC_DD_READ_FROM_FILE,  HC_DD_POLYNOMIAL  */
  p->ndf = 0;
  p->rdf = p->sdf = NULL;

  /* plate velocities */
  p->pvel_mode = HC_INIT_P_FROM_FILE; /* default is single plate velocity */
  p->pvel_time = -1;

  /* 
     viscosity scan stuff
  */
  p->vscan_n  = 2; /* between 2 and HC_VSCAN_NLAYER_MAX */
  p->vscan_dv =  HC_VSCAN_DV0;
  p->vscan_rlv = FALSE;
  /* 

  filenames
  
  */
  strncpy(p->visc_filename,HC_VISC_FILE,HC_CHAR_LENGTH);
  strncpy(p->dens_filename,HC_DENS_SH_FILE,HC_CHAR_LENGTH);
  strncpy(p->prem_model_filename,PREM_MODEL_FILE,HC_CHAR_LENGTH);
  strncpy(p->ref_geoid_file, HC_GEOID_REF_FILE,HC_CHAR_LENGTH);
  strncpy(p->ref_dtopo_file, HC_DTOPO_REF_FILE,HC_CHAR_LENGTH);
}

/* 

   first, call this routine with a blank **hc 
   
*/
void 
hc_struc_init (hc)
struct hcs **hc;
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

void 
hc_init_polsol_struct (psp)
struct hc_ps *psp;
{
  psp->ncalled = 0;
   /* scaling factors will only be computed once */
  psp->rho_scale = 1.0;
  psp->rho_init = FALSE;
  psp->prop_params_init = FALSE; 	/* parameters for propagator computation */
  psp->abg_init = FALSE;		/* alpha, beta factors */
  psp->prop_mats_init = FALSE;	/* will be true only if save_prop_mats is  */
  psp->tor_init = psp->pol_init = FALSE;
  psp->solver_kludge_l = INT_MAX;  /* every l > psp->solver_kludge_l will
				      have modified core boundary
				      conditions */
}
/* 

initialize all variables, after initializing the parameters 


INPUT: sh_type: type of expansion storage/spherical haronics scheme to use

INPUT: hc_parameters: holds all the settings


OUTPUT: hc structure, gets modified
*/
void 
hc_init_main (hc, sh_type, p)
struct hcs *hc;
int sh_type;
struct hc_parameters *p;
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

  if(!p->print_kernel_only){
    /* 
       
       initialize possible depth dependent scaling of density model
    */
    hc_assign_dd_scaling(HC_INIT_DD_FROM_FILE,dd_dummy,p,hc->r_cmb);
    
    
    if(p->verbose)
      switch(p->dd_dens_scale){
      case HC_DD_CONSTANT:
	fprintf(stderr,"hc_init_main: using constant dln\\rho/dln input density scaling of %g\n",
		(double)hc->dens_scale);
	break;
      case HC_DD_READ_FROM_FILE:
	fprintf(stderr,"hc_init_main: reading density scaling from file\n");
	break;
      case HC_DD_POLYNOMIAL:
	fprintf(stderr,"hc_init_main: using polynomial density scaling (NOT IMPLEMENTED YET)\n");
	exit(-1);
	break;
      default:
	fprintf(stderr,"hc_init_main: error, dd mode %i undefined\n",
		p->dd_dens_scale);
	exit(-1);
      }
  }else{
    fprintf(stderr,"hc_init_main: no density properties read, printing kernels only\n");
  }
  
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
		      (p->solver_mode == HC_SOLVER_MODE_VISC_SCAN)?(TRUE):(FALSE),
		      p->print_kernel_only);
    /* 
       assign all zeroes up to the lmax of the density expansion 
    */
    hc_assign_plate_velocities(hc,p->pvel_mode,p->pvel_filename,TRUE,hc->dens_anom[0].lmax,
			       FALSE,p->read_short_pvel_sh,p->remove_nr,p->verbose);
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
    hc_assign_plate_velocities(hc,HC_INIT_P_FROM_FILE,p->pvel_filename,FALSE,dummy,FALSE,
			       p->read_short_pvel_sh,p->remove_nr,p->verbose);
    /* then read in the density anomalies */
    hc_assign_density(hc,p->compressible,HC_INIT_D_FROM_FILE,p->dens_filename,hc->pvel.p[0].lmax,
		      FALSE,FALSE,p->scale_dens_anom_with_prem,
		      p->verbose,p->read_short_dens_sh, p->dd_dens_scale,p->ndf,p->rdf,p->sdf,
		      (p->solver_mode == HC_SOLVER_MODE_VISC_SCAN)?(TRUE):(FALSE),
		      p->print_kernel_only);
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
		      p->verbose,p->read_short_dens_sh,p->dd_dens_scale,p->ndf,p->rdf,p->sdf,
		      (p->solver_mode == HC_SOLVER_MODE_VISC_SCAN)?(TRUE):(FALSE),
		      p->print_kernel_only);
  }else{
    HC_ERROR("hc_init","boundary condition logic error");
  }
  /* solver tricks */
  hc->psp.solver_kludge_l = p->solver_kludge_l;
  if(hc->psp.solver_kludge_l != INT_MAX)
    if(p->verbose)
      fprintf(stderr,"hc_init_main: WARNING: applying solver CMB kludge for l > %i\n",hc->psp.solver_kludge_l);
  
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
void 
hc_init_constants (hc, dens_anom_scale, prem_filename, verbose)
struct hcs *hc;
HC_PREC dens_anom_scale;
char *prem_filename;
hc_boolean verbose;
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
  if(fabs(hc->re - (HC_RE_KM * 1e3)) > 1e-7){
    fprintf(stderr,"%.7e %.7e\n",(double)(HC_RE_KM * 1e3),(double)hc->re);
    HC_ERROR("hc_init_constants","Earth radius mismatch");
  }

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
void 
hc_handle_command_line (argc, argv, start_from_i, p)
int argc;
char **argv;
int start_from_i;
struct hc_parameters *p;
{
  int i;
  hc_boolean used_parameter;
  HC_PREC tmp;
  
  for(i=start_from_i;i < argc;i++){
    used_parameter = FALSE;			/*  */
    if((p->solver_mode == HC_SOLVER_MODE_VISC_SCAN && argc < 2) || /* need
								      one
								      or
								      two
								      arguments
								      for
								      some
								      of
								      the
								      other
								      modes */
       (p->solver_mode == HC_SOLVER_MODE_DYNTOPO_INVERT && argc < 3) || 
       (strcmp(argv[i],"-help")==0) || 
       (strcmp(argv[i],"--help")==0) || 
       (strcmp(argv[i],"-h")==0) || 
       (strcmp(argv[i],"-?")==0)){// help
      /* 
	 help page
      */
      fprintf(stderr,"%s - perform Hager & O'Connell flow computation\n\n",
	      argv[0]);
      fprintf(stderr,"This code can compute velocities, tractions, and geoid for simple density distributions\n");
      fprintf(stderr,"and plate velocities using the semi-analytical approach of Hager & O'Connell (1981).\n");
      fprintf(stderr,"This particular implementation illustrates one possible way to combine the HC solver routines.\n");
      fprintf(stderr,"Based on code by Brad Hager, Richard O'Connell, and Bernhard Steinberger.\n");
      fprintf(stderr,"This version by Thorsten Becker, with contributions by Craig O'Neill\n");
      fprintf(stderr,"compiled with %s precision ((c) 2017, see README.TXT)\n",(HC_PRECISION==16)?("double"):((HC_PRECISION==8)?"single":"quad"));
      switch(p->solver_mode){
      case HC_SOLVER_MODE_VISC_SCAN:
	fprintf(stderr,"usage example:\n\n");
	fprintf(stderr,"bin/hc_visc_scan geoid.ab\n\n");
	fprintf(stderr,"Scan through viscosity values and compute correlation with the geoid in geoid.ab\n\n");
	break;
      case HC_SOLVER_MODE_DYNTOPO_INVERT:
	fprintf(stderr,"usage example:\n\n");
	fprintf(stderr,"bin/hc_invert_dtopo geoid.ab dtopo.ab\n\n");
	fprintf(stderr,"Invert for density anomalies based on geoid in geoid.ab and res topo wrt to air in dtopo.ab\n\n");
	break;
      default:
	fprintf(stderr,"usage example:\n\n");
	fprintf(stderr,"bin/hc -vvv\n\n");
	fprintf(stderr,"Compute mantle flow solution using the default input files:\n");
	fprintf(stderr,"  viscosity profile visc.dat\n");
	fprintf(stderr,"  density profile   dens.sh.dat\n");
	fprintf(stderr,"  earth model       prem/prem.dat\n");
	fprintf(stderr,"and provide lots of output. Default setting is quiet operation.\n\n");
	break;
      }
      fprintf(stderr,"See README.TXT in the installation directory for example for how to plot output, and\n");
      fprintf(stderr,"http://http://www-udc.ig.utexas.edu/external/becker/seatree/ for a graphical user interface.\n");
      fprintf(stderr,"https://goo.gl/D4E8oy for a VirtualBox install.\n\n");


      fprintf(stderr,"density anomaly options:\n");
      fprintf(stderr,"-dens\tname\tuse name as a SH density anomaly model (%s)\n",
	      p->dens_filename);
      fprintf(stderr,"\t\tAll density anomalies are in units of %g%% of PREM, all SH coefficients\n\t\tin Dahlen & Tromp convention.\n",
	      HC_DENSITY_SCALING*100);
      
      fprintf(stderr,"-dshs\t\tuse the short, Becker & Boschi (2002) format for the SH density model (%s)\n",
	      hc_name_boolean(p->read_short_dens_sh));

      fprintf(stderr,"-ds\tval\tdensity scaling factor (%g)\n",
	      (double)p->dens_anom_scale);
      fprintf(stderr,"-dnp\t\tdo not scale density anomalies with PREM but rather mean density (%s)\n",
	      hc_name_boolean(!p->scale_dens_anom_with_prem));
      fprintf(stderr,"-dsf\tfile\tread depth dependent density scaling from file\n");
      fprintf(stderr,"\t\t(overrides -ds, %s), use pdens.py to edit\n\n",
	      hc_name_boolean((p->dd_dens_scale ==  HC_DD_READ_FROM_FILE)));
      //fprintf(stderr,"-dsp\t\tuse polynomial density scaling (overrides -ds, clashes with -dsf, %s)\n\n", 
      //hc_name_boolean((p->dd_dens_scale ==  HC_DD_POLYNOMIAL)));
    
      fprintf(stderr,"Earth model options:\n");
      fprintf(stderr,"-prem\tname\tset Earth model to name (%s)\n",
	      p->prem_model_filename);
      if((p->solver_mode == HC_SOLVER_MODE_VISC_SCAN)||
	 (p->solver_mode == HC_SOLVER_MODE_DYNTOPO_INVERT)){
	
	fprintf(stderr,"\t\tWARNING: Will loop through a four layer viscosity scan\nscan options:\n");
	
	fprintf(stderr,"-gref\tname\tuse filename for reference geoid (%s)\n",
		p->ref_geoid_file);
	fprintf(stderr,"-vs_n\tn\tuse n layers out of %i for viscosity scane (%i)\n",
		HC_VSCAN_NLAYER_MAX,p->vscan_n);
	fprintf(stderr,"-vs_r\t\trestrict v[%i] to expected relative layer strength (%i)\n",
		HC_VSCAN_NLAYER_MAX,p->vscan_rlv);
	fprintf(stderr,"-vs_dv\tval\tuse val spacing in log space for viscosity scan (%g)\n",
		(double)p->vscan_dv);
	fprintf(stderr,"-vs_zlm\tdepth\tuse depth[km] for the upper/lower mantle boundary (%g)\n",
		(double)HC_Z_DEPTH(p->rlayer[0]));
	fprintf(stderr,"-vs_zau\tdepth\tuse depth[km] for the asthenosphere/upper mantle boundary (%g)\n",
		(double)HC_Z_DEPTH(p->rlayer[1]));

	if(p->solver_mode == HC_SOLVER_MODE_DYNTOPO_INVERT)
	  fprintf(stderr,"-dtref\tname\tuse filename for reference dynamic topography (%s)\n",
		  p->ref_dtopo_file);
      }else{
	fprintf(stderr,"-vf\tname\tviscosity structure filename (%s), use pvisc.py to edit\n",
		p->visc_filename);
	fprintf(stderr,"\t\tThis file is in non_dim_radius viscosity[Pas] format\n");
      }
      
      fprintf(stderr,"\nboundary condition options:\n");
      fprintf(stderr,"-fs\t\tperform free slip computation (%s)\n",hc_name_boolean(p->free_slip));
      fprintf(stderr,"-ns\t\tperform no slip computation (%s)\n",hc_name_boolean(p->no_slip));
      fprintf(stderr,"-pvel\tname\tset prescribed surface velocities from file name (%s)\n",
	      hc_name_boolean(p->platebc));
      fprintf(stderr,"\t\tThe file (e.g. %s) is based on a DT expansion of cm/yr velocity fields.\n",HC_PVEL_FILE);
      fprintf(stderr,"-vshs\t\tuse the short format (only lmax in header) for the plate velocities (%s)\n",
	      hc_name_boolean(p->read_short_pvel_sh));
      fprintf(stderr,"-rnr\t\tremove any net rotation component from the plate velocities (%s)\n",
	      hc_name_boolean(p->remove_nr));
      fprintf(stderr,"-vdir\t\tvelocities are given in files name/vel.1.ab to vel.%i.ab for different times,\n\t\t-%g to -1 Ma before present, where name is from -pvel\n",
	      HC_PVEL_TSTEPS,(double)HC_PVEL_TSTEPS);
      fprintf(stderr,"-vtime\ttime\tuse this particular time step of the plate velocities (%g)\n\n",
	      (double)p->pvel_time);

      fprintf(stderr,"solution procedure and I/O options:\n");
      fprintf(stderr,"-cbckl\tval\twill modify CMB boundary condition for all l > val with solver kludge (%i)\n",
	      p->solver_kludge_l);
      if(p->solver_mode == HC_SOLVER_MODE_DEFAULT){
	/* these only apply to the regular mode */
	fprintf(stderr,"-ng\t\tdo not compute and print the geoid (%i)\n",
		p->compute_geoid);
	fprintf(stderr,"-ag\t\tcompute geoid at all layer depths, as opposed to the surface only\n");
	
	fprintf(stderr,"-pptsol\t\tprint pol[6] and tor[2] solution vectors (%s)\n",
		hc_name_boolean(p->print_pt_sol));
	fprintf(stderr,"-px\t\tprint the spatial solution to file (%s)\n",
		hc_name_boolean(p->print_spatial));
	fprintf(stderr,"-rtrac\t\tcompute srr,srt,srp tractions [MPa] instead of velocities [cm/yr] (default: vel)\n");
	fprintf(stderr,"-htrac\t\tcompute stt,stp,spp tractions [MPa] instead of velocities [cm/yr] (default: vel)\n");
	fprintf(stderr,"-pk\t\tprint kernels only, up to default L = %i\n",HC_LMAX_DEFAULT);
	
      }
      fprintf(stderr,"-v\t-vv\t-vvv: verbosity levels (%i)\n",
	      (int)(p->verbose));
      fprintf(stderr,"\n\n");
      exit(-1);
    }else if(strcmp(argv[i],"-ds")==0){	/* density anomaly scaling factor */
      hc_advance_argument(&i,argc,argv);
      sscanf(argv[i],HC_FLT_FORMAT,&p->dens_anom_scale);
      used_parameter = TRUE;
    }else if(strcmp(argv[i],"-vs_n")==0){	
      hc_advance_argument(&i,argc,argv);
      sscanf(argv[i],"%i",&p->vscan_n);
      if((p->vscan_n < 2) || (p->vscan_n >  HC_VSCAN_NLAYER_MAX)){
	fprintf(stderr,"hc_init: error, vscan layer %i out of bounds, max is %i\n",
		p->vscan_n,HC_VSCAN_NLAYER_MAX);
	exit(-1);
      }
      used_parameter = TRUE;
    }else if(strcmp(argv[i],"-vs_dv")==0){	
      hc_advance_argument(&i,argc,argv);
      sscanf(argv[i],HC_FLT_FORMAT,&p->vscan_dv);
      used_parameter = TRUE;
    }else if(strcmp(argv[i],"-vs_zlm")==0){	
      hc_advance_argument(&i,argc,argv);
      sscanf(argv[i],HC_FLT_FORMAT,&tmp);
      p->rlayer[0] = HC_ND_RADIUS(tmp);
      used_parameter = TRUE;
    }else if(strcmp(argv[i],"-vs_r")==0){	
      hc_toggle_boolean(&p->vscan_rlv);
      used_parameter = TRUE;
      if(p->verbose)
	fprintf(stderr,"hc_init: WARNING: restricting viscosity scan\n");
    }else if(strcmp(argv[i],"-vs_zau")==0){	
      hc_advance_argument(&i,argc,argv);
      sscanf(argv[i],HC_FLT_FORMAT,&tmp);
      p->rlayer[1] = HC_ND_RADIUS(tmp);
      used_parameter = TRUE;
    }else if(strcmp(argv[i],"-cbckl")==0){	/* solver kludge */
      hc_advance_argument(&i,argc,argv);
      sscanf(argv[i],"%i",&p->solver_kludge_l);
      used_parameter = TRUE;
    }else if(strcmp(argv[i],"-vtime")==0){	/* */
      hc_advance_argument(&i,argc,argv);
      sscanf(argv[i],HC_FLT_FORMAT,&p->pvel_time);
      used_parameter = TRUE;
    }else if(strcmp(argv[i],"-fs")==0){	/* free slip flag */
      p->free_slip = TRUE;p->no_slip = FALSE;
      used_parameter = TRUE;
    }else if(strcmp(argv[i],"-ns")==0){	/* no slip flag */
      p->no_slip = TRUE;p->free_slip = FALSE;
      used_parameter = TRUE;
    }else if(strcmp(argv[i],"-pk")==0){	
      p->print_kernel_only = TRUE;
      used_parameter = TRUE;
    }else if(strcmp(argv[i],"-dshs")==0){ /* use short format for densities ? */
      hc_toggle_boolean(&p->read_short_dens_sh);
      used_parameter = TRUE;
    }else if(strcmp(argv[i],"-vshs")==0){ /* use short format for velocities? */
      hc_toggle_boolean(&p->read_short_pvel_sh);
      used_parameter = TRUE;
    }else if(strcmp(argv[i],"-prem")==0){ /* PREM filename */
      hc_advance_argument(&i,argc,argv);
      strncpy(p->prem_model_filename,argv[i],HC_CHAR_LENGTH);
      used_parameter = TRUE;
    }else if(strcmp(argv[i],"-dnp")==0){
      hc_toggle_boolean(&p->scale_dens_anom_with_prem);
      used_parameter = TRUE;
    }else if(strcmp(argv[i],"-rnr")==0){
      hc_toggle_boolean(&p->remove_nr);
      used_parameter = TRUE;
    }else if(strcmp(argv[i],"-pvel")==0){ /* velocity filename, this will switch off free slip */
      hc_advance_argument(&i,argc,argv);
      strncpy(p->pvel_filename,argv[i],HC_CHAR_LENGTH);
      p->platebc = TRUE;p->no_slip = TRUE;p->free_slip = FALSE;
      used_parameter = TRUE;
    }else if(strcmp(argv[i],"-dens")==0){ /* density filename */
      hc_advance_argument(&i,argc,argv);
      strncpy(p->dens_filename,argv[i],HC_CHAR_LENGTH);
      used_parameter = TRUE;
    }else if(strcmp(argv[i],"-dsf")==0){ /* density scaling filename */
      p->dd_dens_scale = HC_DD_READ_FROM_FILE;
      hc_advance_argument(&i,argc,argv);
      strncpy(p->dens_scaling_filename,argv[i],HC_CHAR_LENGTH);
      used_parameter = TRUE;
    }else if(strcmp(argv[i],"-dsp")==0){
      p->dd_dens_scale = HC_DD_POLYNOMIAL;
      hc_advance_argument(&i,argc,argv);
      used_parameter = TRUE;
    }else if(strcmp(argv[i],"-gref")==0){ /* geoid reference */
      hc_advance_argument(&i,argc,argv);
      strncpy(p->ref_geoid_file,argv[i],HC_CHAR_LENGTH);
      used_parameter = TRUE;
    }else if(strcmp(argv[i],"-dtref")==0){ /* dyn topo reference */
      hc_advance_argument(&i,argc,argv);
      strncpy(p->ref_dtopo_file,argv[i],HC_CHAR_LENGTH);
      used_parameter = TRUE;
    }else if(strcmp(argv[i],"-vf")==0){ /* viscosity filename */
      hc_advance_argument(&i,argc,argv);
      strncpy(p->visc_filename,argv[i],HC_CHAR_LENGTH);
      used_parameter = TRUE;
    }else if(strcmp(argv[i],"-vdir")==0){	/* velocities in dir */
      p->pvel_mode = HC_INIT_P_FROM_DIRECTORY;
      used_parameter = TRUE;
    }else if(strcmp(argv[i],"-v")==0){	/* verbosities */
      p->verbose = 1;
      used_parameter = TRUE;
    }else if(strcmp(argv[i],"-vv")==0){	/* verbosities */
      p->verbose = 2;
      used_parameter = TRUE;
    }else if(strcmp(argv[i],"-vvv")==0){	
      p->verbose = 3;
      used_parameter = TRUE;
    }			

    if(p->solver_mode == HC_SOLVER_MODE_DEFAULT){
      /* 
	 those subsequent options only apply for default solver
	 mode  
      */
      if(strcmp(argv[i],"-px")==0){	/* print spatial solution? */
	hc_toggle_boolean(&p->print_spatial);
	used_parameter = TRUE;
      }else if(strcmp(argv[i],"-pptsol")==0){	/* print
						   poloidal/toroidal
						   solution
						   parameters */
	hc_toggle_boolean(&p->print_pt_sol);
	used_parameter = TRUE;
      }else if(strcmp(argv[i],"-ng")==0){	/* do not compute geoid */
	p->compute_geoid = 0;
	used_parameter = TRUE;
      }else if(strcmp(argv[i],"-ag")==0){	/* compute geoid at all layers */
	p->compute_geoid = 2;
	used_parameter = TRUE;		
      }else if(strcmp(argv[i],"-rtrac")==0){	/* compute radial
						   tractions */
	p->solution_mode = HC_RTRACTIONS;
	used_parameter = TRUE;
      }else if(strcmp(argv[i],"-htrac")==0){	/* compute horizontal
						   tractions */
	p->solution_mode = HC_HTRACTIONS;
	used_parameter = TRUE;
      }
    } /* end default operation mode branch  */
    if(!used_parameter){
      fprintf(stderr,"%s: can not use parameter %s, use -h for help page\n",
	      argv[0],argv[i]);
      exit(-1);
    }
  }
  if(p->print_kernel_only){	/* make room to store kernels */
    fprintf(stderr,"%s: print kernel overriding geoid settings\n",argv[0]);
    p->compute_geoid = 2;
  }
  
  if((p->solver_mode == HC_SOLVER_MODE_VISC_SCAN) ||
     (p->solver_mode == HC_SOLVER_MODE_DYNTOPO_INVERT)){
    /* we need a reference  geoid */
    hc_read_scalar_shexp(p->ref_geoid_file,&(p->ref_geoid),"reference geoid",p);
  }
  if(p->solver_mode == HC_SOLVER_MODE_DYNTOPO_INVERT){
    /* we need a reference dynamic topography */
    hc_read_scalar_shexp(p->ref_dtopo_file,&(p->ref_dtopo),
			 "reference dynamic topography",p);
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
void 
hc_assign_viscosity (hc, mode, elayer, p)
struct hcs *hc;
int mode;
HC_PREC elayer[4];
struct hc_parameters *p;
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
	      (double)hc->visc[0],(double)hc->visc[1],
	      (double)hc->visc[2],(double)hc->visc[3]);
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
    in = hc_fopen(p->visc_filename,"r","hc_assign_viscosity");
    hc_vecrealloc(&hc->rvisc,1,"hc_assign_viscosity");
    hc_vecrealloc(&hc->visc,1,"hc_assign_viscosity");
    hc->nvis = 0;mean = 0.0;mws = 0.0;
    /* read sscanf string */
    hc_get_flt_frmt_string(fstring,2,FALSE);
    rold = hc->r_cmb;
    /* start read loop  */
    while(fscanf(in,HC_TWO_FLT_FORMAT,
		 (hc->rvisc+hc->nvis),(hc->visc+hc->nvis))==2){
      if(hc->visc[hc->nvis] < 1e15)
	fprintf(stderr,"hc_assign_viscosity: WARNING: expecting viscosities in Pas, read %g at layer %i\n",
		(double)hc->visc[hc->nvis],hc->nvis);
      /* normalize viscosity here */
      hc->visc[hc->nvis] /= hc->visnor;
      if(hc->nvis == 0)
	if( hc->rvisc[hc->nvis] < hc->r_cmb-0.01){
	  fprintf(stderr,"hc_assign_viscosity: error: first radius %g is below CMB, %g\n",
		  (double)hc->rvisc[hc->nvis], (double)hc->r_cmb);
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
		  hc->nvis+1,(double)hc->rvisc[hc->nvis],
		  (double)hc->rvisc[hc->nvis-1]);
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
	      (double)hc->rvisc[hc->nvis-1]);
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
	      (double)mean, (double)hc->visnor, (double)(mean*hc->visnor));
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
void 
hc_assign_density (hc, compressible, mode, filename, nominal_lmax, layer_structure_changed, density_in_binary, scale_dens_anom_with_prem, verbose, use_short_format, dd_dens_scale, ndf, rdf, sdf, save_orig_danom, print_kernerl_only)
struct hcs *hc;
hc_boolean compressible;
int mode;
char *filename;
int nominal_lmax;
hc_boolean layer_structure_changed;
hc_boolean density_in_binary;
hc_boolean scale_dens_anom_with_prem;
hc_boolean verbose;
hc_boolean use_short_format;
hc_boolean dd_dens_scale;
int ndf;
HC_PREC *rdf;
HC_PREC *sdf;
hc_boolean save_orig_danom;
hc_boolean print_kernerl_only;
{
  FILE *in;
  int type,lmax,shps,ilayer,nset,ivec,i,j;
  HC_PREC *dtop,*dbot,zlabel,local_scale,dens_scale[1];
  double rho0;
  hc_boolean reported = FALSE,read_on;
  HC_PREC dtmp[3];
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
      sh_aexp_equals_bexp_coeff((hc->dens_anom+i),
				(hc->dens_anom_orig+i));
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

    in = hc_fopen(filename,"r","hc_assign_density");
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
      if(fscanf(in,"%i",&nset) != 1)
	HC_ERROR("hc_assign_density","code read error");
      ilayer = -1;
    }else{
      if(verbose)
	fprintf(stderr,"hc_assign_density: using default SH format for density\n");
    }
    
    read_on = TRUE;
    while(read_on){
      if(use_short_format){
	/* short format I/O */
	i  = fscanf(in,HC_FLT_FORMAT,dtmp);zlabel = (HC_PREC)dtmp[0];
	i += fscanf(in,"%i",&lmax);
	read_on = (i == 2)?(TRUE):(FALSE);
	ivec = 0;shps = 1;type = HC_DEFAULT_INTERNAL_FORMAT;
	ilayer++;
      }else{
	read_on = sh_read_parameters_from_stream(&type,&lmax,&shps,&ilayer, &nset,
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
	    fprintf(stderr,"hc_assign_density: non_dim radius                 %% factor    PREM \\rho/mean_rho          layer #             depth[km]  rho[kg/m^3]\n");
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
	hc->rden[hc->inho] = HC_ND_RADIUS((HC_PREC)zlabel);
	if(scale_dens_anom_with_prem){
	  /* 
	     
	     get reference density at this level
	     
	  */


	  prem_get_rho(&rho0,(double)(hc->rden[hc->inho]),
		       hc->prem);

	  rho0 /= 1000.0;
	  if(rho0 < 3)
	    fprintf(stderr,"\nhc_assign_density: WARNING: using small (%g) density from PREM for layer at depth %g\n\n",
		    (double)rho0*1000,(double)HC_Z_DEPTH(hc->rden[hc->inho]));
	}else{
	  /* mean value */
	  rho0 =  (double)hc->avg_den_mantle;
	}
	/* 
	   density anomaly
	*/
	/* scaling factor without depth dependence */
	dens_scale[0] = HC_DENSITY_SCALING  * (HC_PREC)rho0;
	if(verbose >= 2){
	  fprintf(stderr,"hc_assign_density: r: %11g anom scales: %11g x %11g = %11g\t%5i out of %i, z: %11g  %6.1f\n",
		  (double)hc->rden[hc->inho],
		  HC_DENSITY_SCALING,rho0/ (double)hc->avg_den_mantle,(double)dens_scale[0],hc->inho+1,nset,(double)zlabel,rho0*1000);
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
		    (double)hc->rden[hc->inho], 
		    (double)hc->rden[hc->inho-1]);
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
	sh_read_coefficients_from_stream((hc->dens_anom+hc->inho),1,lmax,in,density_in_binary,
					 dens_scale,verbose);
	hc->inho++;
      }	/* end actualy read on */
    } /* end while */
    /* 
       assign actual top layer density 
    */
    hc->rho_top_kg = rho0 * 1000;
    if(hc->inho != nset)
      HC_ERROR("hc_assign_density","file mode: mismatch in number of layers");
    fclose(in);
    break;
  default:
    HC_ERROR("hc_assign_density","mode undefined");
    break;
  }
  if(save_orig_danom && (!hc->dens_init)){
    /* 
       make a copy of the original density anomaly before applying
       depth dependent scaling, only done once per run 
    */
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
    /* 
       depth dependent factor? 
    */
    local_scale = hc_find_dens_scale(hc->rden[i],hc->dens_scale,dd_dens_scale,rdf,sdf,ndf);
    sh_scale_expansion((hc->dens_anom+i),local_scale);
    if(verbose >= 2){
      fprintf(stderr,"hc_assign_density: r: %11g additional %s d\\rho/dinput: %11g \tlayer %5i out of %i\n",
	      (double)hc->rden[i],
	      (dd_dens_scale == HC_DD_READ_FROM_FILE)?("depth-dependent"):((dd_dens_scale==HC_DD_CONSTANT)?("constant"):("polynomial")),(double)local_scale,i,hc->inho);
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
    if(hc->rden[0] <= hc->r[0]){
      fprintf(stderr,"hc_assign_density: rden[0]: %g r[0]: %g\n",(double)hc->rden[0],(double) hc->r[0]);
      HC_ERROR("hc_assign_density","first density layer has to be above internal CMB limit");
    }
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
	fprintf(stderr,"hc_assign_density: dens %3i: r: %8.6f df: %8.6f |rho|: %8.4f dtop: %8.3f\n",
		i+1,(double)hc->rden[i],(double)hc->dfact[i],
		(double)sqrt(sh_total_power((hc->dens_anom+i))),
		(double)dtop[i]);
    free(dbot);free(dtop);
  } /* end layer structure part */

  hc->dens_init = TRUE;
}
/* 

find depth dependent scaling

*/
HC_PREC 
hc_find_dens_scale (r, s0, depth_dependent, rd, sd, n)
HC_PREC r;
HC_PREC s0;
hc_boolean depth_dependent;
HC_PREC *rd;
HC_PREC *sd;
int n;
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
void 
hc_init_phase_boundaries (hc, npb, verbose)
struct hcs *hc;
int npb;
hc_boolean verbose;
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

void hc_assign_plate_velocities (hc, mode, filename, vel_bc_zero, lmax, pvel_in_binary,
				 read_short_pvel_sh, remove_nr, verbose)
     struct hcs *hc;
     int mode;
     char *filename;
     hc_boolean vel_bc_zero;
     int lmax;
     hc_boolean pvel_in_binary;
     hc_boolean read_short_pvel_sh;
     hc_boolean remove_nr;
     hc_boolean verbose;
{
  int i;
  char nfilename[HC_CHAR_LENGTH + 10];
  if(hc->pvel_init)		/* if called twice, trouble with pvels structure */
    HC_ERROR("hc_assign_plate_velocities","what to do if called twice?");
  /* init the plate velo structure  */
  hc->pvel.n = 1;
  hc->pvel.p = (struct sh_lms *)calloc(2,sizeof(struct sh_lms)*2);if(!hc->pvel.p)HC_MEMERROR("pvel");
  hc->pvel.t = (HC_PREC *)calloc(1,sizeof(HC_PREC));if(!hc->pvel.t)HC_MEMERROR("pvel");
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
      hc_init_single_plate_exp(filename,hc,pvel_in_binary,hc->pvel.p,TRUE,read_short_pvel_sh,
			       remove_nr,verbose);
      if(verbose)
	fprintf(stderr,"hc_assign_plate_velocities: read single set of surface velocities, lmax %i: |pol|: %11g |tor|: %11g\n",
		lmax,sqrt(sh_total_power((hc->pvel.p))),sqrt(sh_total_power((hc->pvel.p+1))));
      break;
    case HC_INIT_P_FROM_DIRECTORY:
      /*  */
      hc->pvel.n = HC_PVEL_TSTEPS;
      hc->pvel.p = (struct sh_lms *)realloc(hc->pvel.p,sizeof(struct sh_lms)*2*hc->pvel.n);if(!hc->pvel.p)HC_MEMERROR("pvel");
      hc->pvel.t = (HC_PREC *)realloc(hc->pvel.t,sizeof(HC_PREC)*hc->pvel.n);if(!hc->pvel.t)HC_MEMERROR("pvel");
      for(i=0;i < hc->pvel.n;i++){
	sprintf(nfilename,"%s.%i.ab",filename,HC_PVEL_TSTEPS-i);
	hc_init_single_plate_exp(nfilename,hc,pvel_in_binary,hc->pvel.p+i*2,TRUE,read_short_pvel_sh,
				 remove_nr,verbose);
	hc->pvel.t[i] = (HC_PREC)-(HC_PVEL_TSTEPS-i); /* time, count
							past negative,
							and start with
							earliest
							time */
      }
      break;
    default:
      HC_ERROR("hc_assign_plate_velocities","op mode undefined");
    }
  }else{
    /* 
       initialize with zeroes
    */
    if(hc->pvel_init){
      sh_clear_alm(hc->pvel.p);
      sh_clear_alm((hc->pvel.p+1));
    }else{
      /* use irregular grid */
      sh_init_expansion(hc->pvel.p,    lmax,hc->sh_type, 1,verbose,FALSE);
      sh_init_expansion((hc->pvel.p+1),lmax,hc->sh_type, 1,verbose,FALSE);
    }
    if(verbose)
      fprintf(stderr,"hc_assign_plate_velocities: using no-slip surface BC, lmax %i\n",
	      lmax);
  }
  hc->pvel_init = TRUE;
}
 

/* 

   read in single plate velocity (i.e. poloidal/toroidal) field 

 */
void  hc_init_single_plate_exp (filename, hc, pvel_in_binary, pvel, check_for_nr,
				read_short_pvel_sh, remove_nr,verbose)
     char *filename;
     struct hcs *hc;
     hc_boolean pvel_in_binary;
     struct sh_lms *pvel;
     hc_boolean check_for_nr;
     hc_boolean read_short_pvel_sh;
     hc_boolean verbose;
     hc_boolean remove_nr;
{
  FILE *in;
  int type,shps,ilayer,nset,ivec,lmax;
  HC_PREC zlabel,vfac[2],t10[2],t11[2],nr_amp;
  /* scale to go from cm/yr to internal scale */
  vfac[0] = vfac[1] = 1.0/hc->vel_scale;
  
  in = hc_fopen(filename,"r","hc_init_single_plate_exp");
  if(read_short_pvel_sh){
    ivec = 1;shps = 2;type = HC_DEFAULT_INTERNAL_FORMAT;ilayer=0;zlabel=0;nset=1;
    if(fscanf(in,"%i",&lmax) != 1)
      HC_ERROR("hc_init_single_plate_exp","lmax read error");
  }else{
    if(!sh_read_parameters_from_stream(&type,&lmax,&shps,&ilayer, &nset,
				       &zlabel,&ivec,in,FALSE,
				       pvel_in_binary,verbose)){
      fprintf(stderr,"hc_init_single_plate_exp: read error file %s, expecting long header format\n",
	    filename);
      exit(-1);
    } /* check if we read in two sets of expansions */
  }
  if(shps != 2){
    fprintf(stderr,"hc_init_single_plate_exp: two sets expected but found shps: %i in file %s\n",
	    shps,filename);
    exit(-1);
  }
  if((nset > 1)||(fabs(zlabel) > 0.01)){
    fprintf(stderr,"hc_init_single_plate_exp: error: expected one layer at surface, but nset: %i z: %g\n",
	    nset, (double)zlabel);
    exit(-1);
  }
  /* 
     initialize expansion using irregular grid
  */
  sh_init_expansion((pvel+0),lmax,hc->sh_type,1,verbose, FALSE);
  sh_init_expansion((pvel+1),lmax,hc->sh_type,1,verbose, FALSE);
  /* 
     read in expansions, convert to internal format from 
     physical 
  */
  sh_read_coefficients_from_stream(pvel,shps,-1,in,pvel_in_binary,vfac,verbose);
  fclose(in);
  /* 
     scale by 1/sqrt(l(l+1))
  */
  if(pvel[0].lmax > hc->lfac_init)
    hc_init_l_factors(hc,pvel[0].lmax);
  sh_scale_expansion_l_factor((pvel+0),hc->ilfac);
  sh_scale_expansion_l_factor((pvel+1),hc->ilfac);
  if(check_for_nr || remove_nr){
    /*  
	check for net rotation
    */
    sh_get_coeff((pvel+1),1,0,0,TRUE,t10);
    sh_get_coeff((pvel+1),1,0,2,TRUE,t11);
    nr_amp = fabs(t10[0])+fabs(t11[0])+fabs(t11[1]);
    if(nr_amp > 1.0e-7){
      fprintf(stderr,"\nhc_init_single_plate_exp: WARNING: toroidal NR A(1,0): %g A(1,1): %g B(1,1): %g\n\n",
	      (double)t10[0],(double)t11[0],(double)t11[1]);
      if(remove_nr){		/* set NR coefficients to zero */
	t10[0]=t10[1]=t11[0]=t11[1]=0.0;
	sh_write_coeff((pvel+1),1,0,0,TRUE,t10);
	sh_write_coeff((pvel+1),1,0,2,TRUE,t11);
	fprintf(stderr,"hc_init_single_plate_exp: WARNING: setting those coefficients to zero\n");
      }
    }
  }
}

/* 

initialize an array with sqrt(l(l+1)) factors
from l=0 .. lmax+1

pass lfac initialized (say, as NULL)

*/
void  hc_init_l_factors (hc, lmax)
     struct hcs *hc;
     int lmax;
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
void  hc_get_blank_expansions (expansion, nnew, nold, calling_sub)
     struct sh_lms **expansion;
     int nnew;
     int nold;
     char *calling_sub;
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
void 
hc_struc_free (hc)
struct hcs **hc;
{
  free((*hc)->visc);
  free((*hc)->rvisc);
  free((*hc)->qwrite);

  sh_free_expansion((*hc)->dens_anom,1);
  
  free(*hc);
}

/*  

assign a depth dependent density scale, or assign four layers of scaling

*/
void 
hc_assign_dd_scaling (mode, dlayer, p, rcmb)
int mode;
HC_PREC dlayer[4];
struct hc_parameters *p;
HC_PREC rcmb;
{
  HC_PREC smean;
  HC_PREC dtmp[2];
  int i;
  FILE *in;
  if(p->dd_dens_scale == HC_DD_READ_FROM_FILE){
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
	in = hc_fopen(p->dens_scaling_filename,"r","hc_assign_dd_scaling");
	while(fscanf(in,HC_TWO_FLT_FORMAT,dtmp,(dtmp+1)) == 2){
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
	  fprintf(stderr,"hc_assign_dd_density: read scaling on %i layers, rough mean: %g\n",p->ndf,(double)smean);
	break;
    case HC_INIT_DD_FOUR_LAYERS:
      p->ndf = 4;
      hc_vecrealloc(&p->rdf,(p->ndf),"hc_assign_dd_scaling");
      hc_vecrealloc(&p->sdf,(p->ndf),"hc_assign_dd_scaling");
      p->rdf[0] = rcmb;
      p->rdf[1] = p->rlayer[0];
      p->rdf[2] = p->rlayer[1];
      p->rdf[3] = p->rlayer[2];
      for(i=0;i<4;i++)
	p->sdf[i] = dlayer[i];
      if(p->verbose)
	fprintf(stderr,"hc_assign_dd_density: assigned four layer density scaling factors %g %g %g %g\n",
		(double)p->sdf[0],	
		(double)p->sdf[1],
		(double)p->sdf[2],	
		(double)p->sdf[3]);
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
void  hc_read_scalar_shexp (filename, sh_exp, name, p)
     char *filename;
     struct sh_lms **sh_exp;
     char *name;
     struct hc_parameters *p;
{
  int type,lmax,shps,ilayer,nset,ivec;
  HC_PREC zlabel;
  FILE *in;
  HC_PREC fac[3]={1,1,1};
  
  in = fopen(filename,"r");
  if(!in){
    fprintf(stderr,"hc_read_scalar_shexp: error, could not open %s \"%s\", expecting long scalar SH format\n",
	    name,filename);
    exit(-1);
  }
  
    
  /* read in the file */
  sh_read_parameters_from_stream(&type,&lmax,&shps,&ilayer,&nset,
				 &zlabel,&ivec,in,FALSE,
				 FALSE,p->verbose);
  if((ivec != 0)||(shps!=1)){
    fprintf(stderr,"hc_read_scalar_shexp: error, expecting scalar, long format SH in %s\n",
	    filename);
    exit(-1);
  }
  sh_allocate_and_init(sh_exp,shps,lmax,type,ivec,p->verbose,0);
  sh_read_coefficients_from_stream(*sh_exp,shps,-1,in,FALSE,fac,p->verbose);
  fclose(in);
  if(p->verbose)
    fprintf(stderr,"hc_read_scalar_shexp: read %s from %s, L=%i\n",
	    name,filename,(*sh_exp)->lmax);
  
}

void  hc_select_pvel (time, pvel, p, verbose)
     HC_PREC time;
     struct pvels *pvel;
     struct sh_lms *p;
hc_boolean verbose;
{
  int i;
  hc_boolean hit;
  if(pvel->n == 0){
    /* do nothing, if not initialized proper */
    return;
  }else if(pvel->n == 1){
    if(verbose){
      fprintf(stderr,"hc_select_pvel: only one plate velocity loadeed, disregarding time argument\n");
    }
    /* 
       only one plate velocity 
       
    */
    sh_copy_lms(pvel->p,p);
    sh_copy_lms((pvel->p+1),(p+1));
    
  }else{
    /* this wasn't implemented fully, a bit odd. i just put in some
       simple stuff for now */
    for(hit = FALSE,i=0;i < pvel->n;i++){
      if(fabs(pvel->t[i]-time) < HC_EPS_PREC){
	sh_copy_lms(pvel->p+i*2,p);
	sh_copy_lms((pvel->p+i*2+1),(p+1));
	hit=TRUE;
      }
    }
    if(!hit){
      fprintf(stderr,"hc_select_pvel: was searching for time %g amongst ",
	      (double)time);
      for(i=0;i < pvel->n;i++)fprintf(stderr,"%g ",(double)pvel->t[i]);
      fprintf(stderr,"\n");
      HC_ERROR("hc_select_pvel","interpolation not implemented yet");
    }
  }
}

/* 
   open a file safely and give an error message if there was
   a problem
*/
FILE *
hc_fopen (name, mode, program)
char *name;
char *mode;
char *program;
{
  FILE *in;
  if((in=fopen(name,mode)) == NULL){
    fprintf(stderr,"%s: error: can not open file %s for mode %s access\n",
	    program,name,mode);
    exit(-1);
  }
  return in;
}
