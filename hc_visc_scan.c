#include "hc.h"
/* 

   implementation of Hager & O'Connell (1981) method of solving mantle
   circulation given internal density anomalies, only radially varying
   viscosity, and either free-slip or plate velocity boundary
   condition at surface. based on Hager & O'Connell (1981), Hager &
   Clayton (1989), and Steinberger (2000). the original code is due to
   Brad Hager, Rick O'Connell, and was largely modified by Bernhard
   Steinberger. this version by Thorsten Becker (twb@ig.utexas.edu) for
   additional comments, see hc.c

   scan through viscosities and compute correlation with the geoid


   output viscosities are log10(eta/1e21), top to bottom
   
*/

/* indices for arrays */
#define IVMIN 0
#define IVMAX 1
#define IDV 2



int main(int argc, char **argv)
{
  struct hcs *model;		/* main structure, make sure to initialize with 
				   zeroes */
  struct sh_lms *sol_spectral=NULL, *geoid = NULL;		/* solution expansions */
  struct sh_lms *pvel=NULL;					/* local plate velocity expansion */
  int nsol,lmax;
  struct hc_parameters p[1]; /* parameters */
  hc_boolean solved = FALSE; /* init with FALSE! */
  hc_boolean vary_umlm = FALSE;
  HC_PREC dv_use,vl[HC_VSCAN_NLAYER_MAX][3],v[HC_VSCAN_NLAYER_MAX];			/*  for viscosity scans */
  /* 
     
  
     (1)
  
     initialize the model structure, this is needed to initialize some
     of the default values before callign the parameter handling
     routine this call also involves initializing the hc parameter
     structure
     
  */
  hc_struc_init(&model);
  /* 
  
     (2)
     init parameters to default values

  */
  hc_init_parameters(p);
  /* 

     special options for this computation

  */
  p->solver_mode = HC_SOLVER_MODE_VISC_SCAN;
  p->visc_init_mode = HC_INIT_E_FOUR_LAYERS;
  p->compute_geoid = 1;		/* compute geoid at surface */
  /* 
     handle command line arguments
  */
  hc_handle_command_line(argc,argv,1,p);
  
  fprintf(stderr,"%s: starting scan using reference %s, dv_ref: %g, nlayer: %i/%i, z_ulm: %g z_asth: %g\n",
	  argv[0],p->ref_geoid_file,(double)p->vscan_dv,p->vscan_n,
	  HC_VSCAN_NLAYER_MAX,
	  (double)HC_Z_DEPTH(p->rlayer[0]),
	  (double)HC_Z_DEPTH(p->rlayer[1]));
  if(p->vscan_n < 0){
    p->vscan_n = - p->vscan_n;
    fprintf(stderr,"%s: i.e. %i layers, but additionally varying upper/lower mantle boundary\n",
	    argv[0],p->vscan_n);
    vary_umlm = TRUE;
  }

  /* 

     begin main program part

  */
#ifdef __TIMESTAMP__
  if(p->verbose)
    fprintf(stderr,"%s: starting version compiled on %s\n",
	    argv[0],__TIMESTAMP__);
#else
  if(p->verbose)
    fprintf(stderr,"%s: starting main program\n",argv[0]);
#endif
  /* 

     (3)
  
     initialize all variables
     
     - choose the internal spherical harmonics convention
     - assign constants
     - assign phase boundaries, if any
     - read in viscosity structure
     - assign density anomalies
     - read in plate velocities

  */
  hc_init_main(model,SH_RICK,p);
  nsol = (model->nradp2) * 3;	/* 
				   number of solutions (r,pol,tor) * (nlayer+2) 

				   total number of layers is nlayer +2, 

				   because CMB and surface are added
				   to intermediate layers which are
				   determined by the spacing of the
				   density model

				*/
  if(p->free_slip)		/* maximum degree is determined by the
				   density expansion  */
    lmax = model->dens_anom[0].lmax;
  else				/* max degree is determined by the
				   plate velocities  */
    lmax = model->pvel.p[0].lmax;	/*  shouldn't be larger than that*/
  /* 
     make sure we have room for the plate velocities 
  */
  sh_allocate_and_init(&pvel,2,lmax,model->sh_type,1,p->verbose,FALSE);
  
  /* init done */
  /* 
     SOLUTION PART
  */
  /* 
     make room for the spectral solution on irregular grid
  */
  sh_allocate_and_init(&sol_spectral,nsol,lmax,model->sh_type,HC_VECTOR,
		       p->verbose,FALSE);
  /* make room for geoid solution at surface */
  sh_allocate_and_init(&geoid,1,model->dens_anom[0].lmax,
		       model->sh_type,HC_SCALAR,p->verbose,FALSE);

    
  /* select plate velocity if not free slip */
  if(!p->free_slip)
    hc_select_pvel(p->pvel_time,&model->pvel,pvel,p->verbose);
  
  /* parameter space log bounds */
  dv_use = p->vscan_dv;	/* don't modify */
  
  switch(p->vscan_n){
  case 4:
    /* 
       
       for layer case

     */
    /* uniform "priors" */
    vl[0][IVMIN]=  -HC_VSCAN_VMAX;vl[0][IVMAX]=HC_VSCAN_VMAX+1e-5;vl[0][IDV]=dv_use; /*   0..100
											  layer
											  log
											  bounds
											  and
											  spacing */
    vl[1][IVMIN]=  -HC_VSCAN_VMAX;vl[1][IVMAX]=HC_VSCAN_VMAX+1e-5;vl[1][IDV]=dv_use; /* 100..410 */
    if(p->free_slip){
      vl[2][IVMIN]=  0;vl[2][IVMAX]=0+1e-5;vl[2][IDV]=dv_use; /* for free
								      slip, only relative
								      viscosisites matter
								      for correlation */
      fprintf(stderr,"%s: for free slip, we set upper mantle (layer 2) to unity (only relative viscosities matter)\n",argv[0]);
    }else{
      vl[2][IVMIN]=  -HC_VSCAN_VMAX;vl[2][IVMAX]=HC_VSCAN_VMAX+1e-5;vl[2][IDV]=dv_use; /* need to actually
											       loop 410 .660 */
    }
    vl[3][IVMIN]=  -HC_VSCAN_VMAX;vl[3][IVMAX]=HC_VSCAN_VMAX+1e-5;vl[3][IDV]=dv_use; /* 660 ... 2871 */

    /* loop */
    for(v[0]=vl[0][IVMIN];v[0] <= vl[0][IVMAX];v[0] += vl[0][IDV])
      for(v[1]=vl[1][IVMIN];v[1] <= vl[1][IVMAX];v[1] += vl[1][IDV])
	for(v[2]=vl[2][IVMIN];v[2] <= vl[2][IVMAX];v[2] += vl[2][IDV])
	  for(v[3]=vl[3][IVMIN];v[3] <= vl[3][IVMAX];v[3] += vl[3][IDV])
	    visc_scan_out(v,geoid,sol_spectral,pvel,p,model,&solved,vary_umlm);
    break;
  case 3:
    dv_use /= 3;		/* refine */
    /* 

       three layer case

    */
    vl[0][IVMIN]=  -HC_VSCAN_VMAX;vl[0][IVMAX]=HC_VSCAN_VMAX+1e-5;vl[0][IDV]=dv_use; /*   0..100
											       layer
											       log
											       bounds
											       and
											       spacing */
    /* vl[1] will be same as vl[2] */
    if(p->free_slip){
      vl[2][IVMIN]=  0;vl[2][IVMAX]=0+1e-5;vl[2][IDV]=dv_use; /* for free
								      slip, only relative
								      viscosisites matter
								      for correlation */
      fprintf(stderr,"%s: for free slip, we set upper mantle (layer 2) to unity (only relative viscosities matter)\n",argv[0]);
    }else{
      vl[2][IVMIN]=  -HC_VSCAN_VMAX;vl[2][IVMAX]=HC_VSCAN_VMAX+1e-5;vl[2][IDV]=dv_use; /* need to actually
											       loop 410 .660 */
    }
    vl[3][IVMIN]=  -HC_VSCAN_VMAX;vl[3][IVMAX]=HC_VSCAN_VMAX+1e-5;vl[3][IDV]=dv_use; /* 660 ... 2871 */
    /* loop */
    for(v[0]=vl[0][IVMIN];v[0] <= vl[0][IVMAX];v[0] += vl[0][IDV])
      for(v[2]=vl[2][IVMIN];v[2] <= vl[2][IVMAX];v[2] += vl[2][IDV]){
	v[1] = v[2];
	for(v[3]=vl[3][IVMIN];v[3] <= vl[3][IVMAX];v[3] += vl[3][IDV])
	  visc_scan_out(v,geoid,sol_spectral,pvel,p,model,&solved,vary_umlm);
      }
    break;
  case 2:
    /* 
       two layer case

       vl[0] and vl[1] will be same as vl[2] 

    */
    dv_use /= 10;		/* finer */
    if(p->free_slip){
      vl[2][IVMIN]=  0;vl[2][IVMAX]=0+1e-5;vl[2][IDV]=dv_use; /* for
								      free
								      slip,
								      only
								      relative
								      viscosisites
								      matter
								      for
								      correlation */
      fprintf(stderr,"%s: for free slip, we set upper mantle (layer 2) to unity (only relative viscosities matter)\n",argv[0]);
    }else{
      vl[2][IVMIN]=  -HC_VSCAN_VMAX;vl[2][IVMAX]=HC_VSCAN_VMAX+1e-5;vl[2][IDV]=dv_use; /* need to actually
											       loop 410 .660 */
    }
    vl[3][IVMIN]=  -HC_VSCAN_VMAX;vl[3][IVMAX]=HC_VSCAN_VMAX+1e-5;vl[3][IDV]=dv_use; /* 660 ... 2871 */
    /* loop */
    for(v[2]=vl[2][IVMIN];v[2] <= vl[2][IVMAX];v[2] += vl[2][IDV]){
      v[0] = v[1] = v[2];
      for(v[3]=vl[3][IVMIN];v[3] <= vl[3][IVMAX];v[3] += vl[3][IDV])
	visc_scan_out(v,geoid,sol_spectral,pvel,p,model,&solved,vary_umlm);
    }
    break;
  default:
    fprintf(stderr,"%s: not set up for %i layers\n",argv[0],p->vscan_n);
    exit(-1);
    break;
  }
  
  /*
    
    free memory
    
  */
  sh_free_expansion(sol_spectral,nsol);
  /* local copies of plate velocities */
  sh_free_expansion(pvel,2);
  /*  */
  sh_free_expansion(geoid,1);
  if(p->verbose)
    fprintf(stderr,"%s: done\n",argv[0]);
  hc_struc_free(&model);
  return 0;
}


/* 
   print out a four layer viscosity structure geoid correlation suite,
   or additionally scan through the upper/lower mantle depths
 */
void 
visc_scan_out (v, geoid, sol_spectral, pvel, p, model, solved, vary_umlm)
HC_PREC *v;
struct sh_lms *geoid;
struct sh_lms *sol_spectral;
struct sh_lms *pvel;
struct hc_parameters *p;
struct hcs *model;
hc_boolean *solved;
hc_boolean vary_umlm;
{
  HC_PREC corr[3],r660=660;
  const HC_PREC rtop = 300.1, rbot = 1800+1e-5, dr = 25;
  if(p->vscan_rlv){
    if((v[0] < v[1])||(v[0] < v[2]))		/* lithosphere should be > asth or upper mantle */
      return;
    if(v[1] > v[2])		/* asthenosphere should be < upper mantle */
      return;
    if(v[2] > v[3])
      return;
  }
  if(vary_umlm){
    for(r660=rtop;r660 <= rbot;r660 += dr){
      /* overwrite 660 as first non CMB boundary from the bottom */
      p->rlayer[0] = HC_ND_RADIUS(r660);
      /* print viscosities of 0...100, 100...410, 410 ... 660 and
	 660...2871 layer in log space */ 
      fprintf(stdout,"%14.7e %14.7e %14.7e %14.7e\t",
	      (double)v[0],(double)v[1],(double)v[2],(double)v[3]);
      hc_calc_geoid_corr_four_layer(v,geoid,sol_spectral,pvel,p,model,solved,corr);
      fprintf(stdout,"%10.7f %10.7f %10.7f\t%8.3f\n",(double)corr[0],(double)corr[1],
	      (double)corr[2],(double)r660);
    }
  }else{
    /* no radius scan */
    fprintf(stdout,"%14.7e %14.7e %14.7e %14.7e\t",
	    (double)v[0],(double)v[1],(double)v[2],(double)v[3]);
    hc_calc_geoid_corr_four_layer(v,geoid,sol_spectral,pvel,p,model,solved,corr);
    fprintf(stdout,"%10.7f %10.7f %10.7f\n",(double)corr[0],(double)corr[1],(double)corr[2]);
  }

}
