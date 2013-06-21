#include "hc.h"
/* 

   implementation of Hager & O'Connell (1981) method of solving mantle
   circulation given internal density anomalies, only radially varying
   viscosity, and either free-slip or plate velocity boundary
   condition at surface. based on Hager & O'Connell (1981), Hager &
   Clayton (1989), and Steinberger (2000). the original code is due to
   Brad Hager, Rick O'Connell, and was largely modified by Bernhard
   Steinberger. this version by Thorsten Becker (twb@usc.edu) for
   additional comments, see hc.c

   scan through viscosities and compute correlation with the geoid

*/

int main(int argc, char **argv)
{
  struct hcs *model;		/* main structure, make sure to initialize with 
				   zeroes */
  struct sh_lms *sol_spectral=NULL, *geoid = NULL;		/* solution expansions */
  struct sh_lms *pvel=NULL;					/* local plate velocity expansion */
  int nsol,lmax,solved;
  struct hc_parameters p[1]; /* parameters */
  HC_PREC corr[2];			/* correlations */
  HC_PREC vl[4][3],v[4],dv;			/*  for viscosity scans */
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
  p->compute_geoid = 1;
  p->compute_geoid_correlations = TRUE;
  
  if(argc > 1){
    /* read in the reference geoid */
    strcpy(p->ref_geoid_file,argv[1]);
    hc_read_scalar_shexp(p->ref_geoid_file,&(p->ref_geoid),"reference geoid",p);
  }else{
    fprintf(stderr,"%s: ERROR: need geoid.ab file as an argument\n",argv[0]);
    fprintf(stderr,"%s: usage:\n\n%s geoid.ab\n\n",argv[0],argv[0]);
    fprintf(stderr,"%s: for help, use:\n\n%s geoid.ab -h\n\n",argv[0],argv[0]);
    exit(-1);
  }
  
  /* 
     handle other command line arguments
  */
  hc_handle_command_line(argc,argv,2,p);
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
    
  /* parameter space log bounds */
  
  dv = .1;			/* spacing */
  
  vl[0][0]=  -3;vl[0][1]=3+1e-5;vl[0][2]=dv; /*   0..100 layer log bounds and spacing */
  vl[1][0]=  -3;vl[1][1]=3+1e-5;vl[1][2]=dv; /* 100..410 */
  if(p->free_slip){
    vl[2][0]=  0;vl[2][1]=0+1e-5;vl[2][2]=dv; /* for free slip,
						 only relative
						 viscosisites
						 matter for
						 correlation */
  }else{
    vl[2][0]=  -3;vl[2][1]=3+1e-5;vl[2][2]=dv; /* need to actually
						  loop 410 .660 */
  }
  vl[3][0]=  -3;vl[3][1]=3+1e-5;vl[3][2]=dv; /* 660 ... 2871 */
  
  
  /*  */
  /* select plate velocity */
  if(!p->free_slip)
    hc_select_pvel(p->pvel_time,&model->pvel,pvel,p->verbose);
  
  solved=0;
  for(v[0]=vl[0][0];v[0] <= vl[0][1];v[0] += vl[0][2])
    for(v[1]=vl[1][0];v[1] <= vl[1][1];v[1] += vl[1][2])
      for(v[2]=vl[2][0];v[2] <= vl[2][1];v[2] += vl[2][2])
	for(v[3]=vl[3][0];v[3] <= vl[3][1];v[3] += vl[3][2]){
	  /* layer viscosity structure */
	  p->elayer[0] = pow(10,v[3]); /* bottom */
	  p->elayer[1] = pow(10,v[2]); /* 660..410 */
	  p->elayer[2] = pow(10,v[1]); /* 410..100 */
	  p->elayer[3] = pow(10,v[0]); /* 100..0  */
	  hc_assign_viscosity(model,HC_INIT_E_FOUR_LAYERS,p->elayer,p);
	  /* print viscosities of 0...100, 100...410, 410 ... 660 and 660...2871  layer */
	  fprintf(stdout,"%14.7e %14.7e %14.7e %14.7e\t",
		  (double)p->elayer[3],(double)p->elayer[2],
		  (double)p->elayer[1],(double)p->elayer[0]);
	  /* compute solution */
	  hc_solve(model,p->free_slip,p->solution_mode,sol_spectral,
		   (solved)?(FALSE):(TRUE), /* density changed? */
		   (solved)?(FALSE):(TRUE), /* plate velocity changed? */
		   TRUE,			/* viscosity changed */
		   FALSE,p->compute_geoid,
		   pvel,model->dens_anom,geoid,
		   p->verbose);
	  /* only output are the geoid correlations, for now */
	  hc_compute_correlation(geoid,p->ref_geoid,corr,1,p->verbose);
	  fprintf(stdout,"%10.7f %10.7f ",
		  (double)corr[0],(double)corr[1]);
	  fprintf(stdout,"\n");
	  solved++;
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

