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

this version by Thorsten Becker (twb@ig.utexas.edu)


>>> SOME COMMENTS FROM THE ORIGINAL CODE <<<

C    * It uses the following Numerical Recipes (Press et al.) routines:
C      four1, realft, gauleg, rk4, rkdumb, ludcmp, lubksb;
C      !!!!!!!!!!!!!!!!!!! rkqc, odeint !!!!!!!!!! take out !!!!!!
C      and the following routines by R.J. O'Connell: 
C      shc2d, shd2c, ab2cs, cs2ab, plmbar, plvel2sh, pltgrid, pltvel,
C      vshd2c, plmbar1.
C      Further subroutines are: kinsub, evalpa, 
C      torsol (all based on "kinflow" by Hager & O'Connell),
C      densub and evppot (based on "denflow" by  Hager & O'Connell),
C      sumsub (based on "sumkd" by  Hager & O'Connell, but
C      substantially speeded up),
C      convert, derivs and shc2dd (based on R.J. O'Connell's shc2d).
C
C      bugs found:
C   * The combination of (1) high viscosity lithosphere 
C     (2) compressible flow (3) kinematic (plate driven) flow
C     doesn't work properly. The problem presumably only occurs
C     at degree 1 (I didn't make this sure) but this is sufficient
C     to screw up everything. It will usually work to reduce the
C     contrast between lithospheric and asthenospheric viscosity.
C     Then make sure that (1) two viscosity structures give similar
C     results for incompressible and (2) incompressible and compressible
C     reduced viscosity look similar (e.g. anomalous mass flux vs. depth)

<<< END OLD COMMENTS 

*/

int main(int argc, char **argv)
{
  struct hcs *model;		/* main structure, make sure to initialize with 
				   zeroes */
  struct sh_lms *sol_spectral=NULL, *geoid = NULL;		/* solution expansions */
  struct sh_lms *pvel=NULL;					/* local plate velocity expansion */
  int nsol,lmax,i;
  FILE *out;
  struct hc_parameters p[1]; /* parameters */
  char filename[HC_CHAR_LENGTH],file_prefix[10];
  HC_PREC *sol_spatial = NULL;	/* spatial solution,
				   e.g. velocities */
  static hc_boolean geoid_binary = FALSE;	/* type of geoid output */
  static HC_CPREC unitya[1] = {1.0};
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
  hc_handle_command_line(argc,argv,1,p);
  /* 

  begin main program part

  */
#ifdef __TIMESTAMP__
  if(p->verbose)
    fprintf(stderr,"%s: starting version compiled on %s\n",argv[0],__TIMESTAMP__);
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
  if(p->compute_geoid == 1)
    /* make room for geoid solution at surface */
    sh_allocate_and_init(&geoid,1,model->dens_anom[0].lmax,
			 model->sh_type,HC_SCALAR,p->verbose,FALSE);
  else if(p->compute_geoid == 2) /* all layers (or kernels) */
    sh_allocate_and_init(&geoid,model->nradp2,model->dens_anom[0].lmax,
			 model->sh_type,HC_SCALAR,p->verbose,FALSE);
  
  /* 
     solve poloidal and toroidal part and sum
  */
  if(!p->free_slip)
    hc_select_pvel(p->pvel_time,&model->pvel,pvel,p->verbose);
  hc_solve(model,p->free_slip,p->solution_mode,sol_spectral,
	   TRUE,TRUE,TRUE,p->print_pt_sol,p->compute_geoid,
	   pvel,model->dens_anom,geoid,
	   p->verbose,p->print_kernel_only);
  /* 
     
     OUTPUT PART
     
  */
  if(!p->print_kernel_only){
    /* 
       
       output of spherical harmonics solution
       
    */
    switch(p->solution_mode){
    case HC_VEL:
      sprintf(file_prefix,"vel");break;
    case HC_RTRACTIONS:
      sprintf(file_prefix,"rtrac");break;
    case HC_HTRACTIONS:
      sprintf(file_prefix,"htrac");break;
    default:
      HC_ERROR(argv[0],"solution mode undefined");break;
    }
    if(p->sol_binary_out)
      sprintf(filename,"%s.%s",file_prefix,HC_SOLOUT_FILE_BINARY);
    else
      sprintf(filename,"%s.%s",file_prefix,HC_SOLOUT_FILE_ASCII);
    if(p->verbose)
      fprintf(stderr,"%s: writing spherical harmonics solution to %s\n",
	      argv[0],filename);
    out = ggrd_open(filename,"w","main");
    hc_print_spectral_solution(model,sol_spectral,out,
			       p->solution_mode,
			       p->sol_binary_out,p->verbose);
    fclose(out);
    /*  */
    if(p->print_density_field){
      /* 
	 print the density field 
      */
      sprintf(file_prefix,"dscaled");
      if(p->sol_binary_out)
	sprintf(filename,"%s.%s",file_prefix,HC_SOLOUT_FILE_BINARY);
      else
	sprintf(filename,"%s.%s",file_prefix,HC_SOLOUT_FILE_ASCII);
      if(p->verbose)
	fprintf(stderr,"%s: writing scaled density anomaly field to %s\n",
		argv[0],filename);
      
      out = ggrd_open(filename,"w","main");
      hc_print_dens_anom(model,out,p->sol_binary_out,p->verbose);
      fclose(out);
    }
    if(p->print_spatial){
      /* 
	 we wish to use the spatial solution
	 
	 expand velocities to spatial base, compute spatial
	 representation
	 
      */
      hc_compute_sol_spatial(model,sol_spectral,&sol_spatial,
			     p->verbose);
      /* 
	 
	 output of spatial solution
	 
      */
      sprintf(filename,"%s.%s",file_prefix,HC_SPATIAL_SOLOUT_FILE);
      /* print lon lat z v_r v_theta v_phi */
      hc_print_spatial_solution(model,sol_spectral,sol_spatial,
				filename,HC_LAYER_OUT_FILE,
				p->solution_mode,p->sol_binary_out,
				p->verbose);
    }
    
  } /* end non-kernel branch */
  
  
  /* compute the geoid? */
  if(p->compute_geoid){
    /* 
       print geoid solution 
    */
    if(p->print_kernel_only)
      sprintf(filename,"%s",HC_GEOID_KERNEL_FILE);
    else
      sprintf(filename,"%s",HC_GEOID_FILE);
    if(p->verbose)
      fprintf(stderr,"%s: writing geoid %sto %s, %s\n",argv[0],
	      ((p->print_kernel_only)?("kernels "):("")),filename,
	      (p->compute_geoid == 1)?("at surface"):("all layers"));
    out = ggrd_open(filename,"w","main");	
    if(p->print_kernel_only){
      hc_print_geoid_kernel(geoid,model->r,model->nradp2,out,p->verbose);
    }else{			/* geoid solutions */
      if(p->compute_geoid == 1) /* surface layer */
	hc_print_sh_scalar_field(geoid,out,FALSE,geoid_binary,p->verbose);
      else{			/* all layers */
	for(i=0;i < model->nradp2;i++){
	  sh_print_parameters_to_stream((geoid+i),1,i,model->nradp2,
					HC_Z_DEPTH(model->r[i]),out,
					FALSE,geoid_binary,p->verbose); 
	  sh_print_coefficients_to_stream((geoid+i),1,out,
					  unitya,geoid_binary,p->verbose); 
	}
      }
    }
    fclose(out);
  }
  /* 
     
  free memory

  */
  sh_free_expansion(sol_spectral,nsol);
  /* local copies of plate velocities */
  sh_free_expansion(pvel,2);
  /*  */
  if(p->compute_geoid == 1)	/* only one layer */
    sh_free_expansion(geoid,1);
  else if(p->compute_geoid == 2) /* all layers */
    sh_free_expansion(geoid,model->nradp2);
  free(sol_spatial);
  if(p->verbose)
    fprintf(stderr,"%s: done\n",argv[0]);
  hc_struc_free(&model);
  return 0;
}

