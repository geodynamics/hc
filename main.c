#include "hc.h"
/* 


implementation of Hager & O'Connell (1981) method of solving mantle
circulation given internal density anomalies, only radially varying
viscosity, and either free-slip or plate velocity boundary condition
at surface

based on Hager & O'Connell (1981), Hager & Clayton (1989), and
Steinberger (2000)

the original code is due to Brad Hager, Rick O'Connell, and was
largely modified by Bernhard Steinberger

this version by Thorsten Becker (twb@usc.edu)

$Id: main.c,v 1.13 2006/01/22 01:11:34 becker Exp becker $


*/

int main(int argc, char **argv)
{
  struct hcs *model;		/* main structure, make sure to initialize with 
				   zeroes */
  struct sh_lms *sol_spectral=NULL;		/* solution expansions */
  int nsol,lmax;
  FILE *out;
  struct hc_parameters p[1]; /* parameters */
  char filename[HC_CHAR_LENGTH];
  float *sol_spatial = NULL;	/* spatial solution,
				   e.g. velocities */
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
     handle command line arguments
  */
  hc_handle_command_line(argc,argv,p);
  /* 

  begin main program part

  */
#ifdef __TIMESTAMP__
  if(p->verbose)
    fprintf(stderr,"%s: starting version compiled on %s\n",argv[0],__TIMESTAMP__);
#else
  if(p->verbose)
    fprintf(stderr,"%s: starting\n",argv[0]);
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
  nsol = (model->nrad+2) * 3;	/* number of solution (r,pol,tor)*(nlayer+2) */
  if(p->free_slip)
    lmax = model->dens_anom[0].lmax;
  else
    lmax = model->pvel[0].lmax;	/*  shouldn't be larger than that*/


  /* init done */
  /* 



  SOLUTION PART


  */
  /* 
     make room for the spectral solution
  */
  sh_allocate_and_init(&sol_spectral,nsol,lmax,model->sh_type,1,
		       p->verbose);
  /* 
     solve poloidal and toroidal part and sum
  */
  hc_solve(model,p->free_slip,p->solution_mode,sol_spectral,
	   TRUE,TRUE,TRUE,p->print_pt_sol,p->verbose);
  /* 
     expand velocities to spatial base, compute spatial representation
  */
  hc_compute_sol_spatial(model,sol_spectral,&sol_spatial,
			 p->verbose);
  /* 

  OUTPUT PART
  
  */
  /* 
     output of spherical harmonics solution
  */
  if(p->sol_binary_out)
    sprintf(filename,"%s",HC_SOLOUT_FILE_BINARY);
  else
    sprintf(filename,"%s",HC_SOLOUT_FILE_ASCII);
  if(p->verbose)
    fprintf(stderr,"%s: writing spherical harmonics solution to %s\n",
	    argv[0],filename);
  out = hc_open(filename,"w","main");
  hc_print_spectral_solution(model,sol_spectral,out,
			     p->solution_mode,
			     p->sol_binary_out,p->verbose);
  fclose(out);
  if(0)
    /* 
       output of spatial solution
    */
    /* print lon lat z v_r v_theta v_phi */
    hc_print_spatial_solution(model,sol_spectral,sol_spatial,
			      HC_SPATIAL_SOLOUT_FILE,
			      HC_LAYER_OUT_FILE,
			      p->solution_mode,p->sol_binary_out,
			      p->verbose);
  /* 
     
  free memory

  */
  sh_free_expansion(sol_spectral,nsol);
  free(sol_spatial);
  if(p->verbose)
    fprintf(stderr,"%s: done\n",argv[0]);
  
  return 0;
}

