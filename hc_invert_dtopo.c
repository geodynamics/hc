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

   invert for compositional anomalies given geoid anomalies [m] and
   residual topography wrt. to air. [m]

*/

int main(int argc, char **argv)
{
  struct hcs *model;		/* main structure, make sure to initialize with 
				   zeroes */
  struct sh_lms *sol_spectral=NULL, *geoid = NULL, *dtopo = NULL;	/* solution expansions */
  struct sh_lms *pvel=NULL;					/* local plate velocity expansion */
  int nsol,lmax,solved;
  struct hc_parameters p[1]; /* parameters */
  HC_PREC gcorr[3],dcorr[3];			/* correlations */
  hc_struc_init(&model);
  hc_init_parameters(p);

  /* 

     special options for this computation

  */
  p->solver_mode = HC_SOLVER_MODE_DYNTOPO_INVERT;
  p->compute_geoid = 1;
  p->solution_mode = HC_RTRACTIONS; /* make sure to compute tractions */
  /*  */
  p->verbose = 1;


  /* 
     handle other command line arguments
  */
  hc_handle_command_line(argc,argv,3,p);

  fprintf(stderr,"%s: using %s for dyn topo and %s for geoid\n",
	  argv[0],p->ref_dtopo_file,p->ref_geoid_file);
  /* 

     begin main program part

  */
  hc_init_main(model,SH_RICK,p);
  nsol = (model->nradp2) * 3;	
  if(p->free_slip)		/* maximum degree is determined by the
				   density expansion  */
    lmax = model->dens_anom[0].lmax;
  else				/* max degree is determined by the
				   plate velocities  */
    lmax = model->pvel.p[0].lmax;	/*  shouldn't be larger than that*/

  sh_allocate_and_init(&pvel,2,lmax,model->sh_type,1,p->verbose,FALSE);
  sh_allocate_and_init(&sol_spectral,nsol,lmax,model->sh_type,HC_VECTOR,
		       p->verbose,FALSE);
  sh_allocate_and_init(&geoid,1,model->dens_anom[0].lmax,
		       model->sh_type,HC_SCALAR,p->verbose,FALSE);
  if(!p->free_slip)
    hc_select_pvel(p->pvel_time,&model->pvel,pvel,p->verbose);

  /* 


   */
  solved=0;
  {
    /* compute solution */
    hc_solve(model,p->free_slip,p->solution_mode,sol_spectral,
	     TRUE, /* density changed? */
	     (solved)?(FALSE):(TRUE), /* plate velocity changed? */
	     TRUE,			/* viscosity changed */
	     FALSE,p->compute_geoid,pvel,model->dens_anom,geoid,
	     p->verbose,FALSE);
    /* extract the top tractions */
    hc_compute_dynamic_topography(model,sol_spectral,&dtopo,TRUE,p->verbose);
    //sh_single_par_and_exp_to_file(dtopo,"dtopo.ab",TRUE,p->verbose);

    /* geoid correlation */
    hc_compute_correlation(geoid,p->ref_geoid,(gcorr),0,p->verbose); /* full correlation */
    hc_compute_correlation(geoid,p->ref_geoid,(gcorr+1),1,p->verbose); /* up to 20 and 4...9 */
    fprintf(stdout,"geoid full: %10.7f L=20: %10.7f \n",(double)gcorr[0],(double)gcorr[1]);

    hc_compute_correlation(dtopo,p->ref_dtopo,(dcorr),0,p->verbose); /* full correlation */
    hc_compute_correlation(dtopo,p->ref_dtopo,(dcorr+1),1,p->verbose); /* up to 20 and 4..9 */
    fprintf(stdout,"dtopo full: %10.7f L=20: %10.7f \n",(double)dcorr[0],(double)dcorr[1]);

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

