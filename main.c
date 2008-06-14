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
  int nsol,lmax,solved;
  FILE *out;
  struct hc_parameters p[1]; /* parameters */
  char filename[HC_CHAR_LENGTH],file_prefix[10];
  float *sol_spatial = NULL;	/* spatial solution,
				   e.g. velocities */
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
  nsol = (model->nrad+2) * 3;	/* number of solution (r,pol,tor)*(nlayer+2) */
  if(p->free_slip)		/* maximum degree is determined by the
				   density expansion  */
    lmax = model->dens_anom[0].lmax;
  else				/* max degree is determined by the
				   plate velocities  */
    lmax = model->pvel[0].lmax;	/*  shouldn't be larger than that*/


  /* init done */
  /* 



  SOLUTION PART


  */
  /* 
     make room for the spectral solution on irregular grid
  */
  sh_allocate_and_init(&sol_spectral,nsol,lmax,model->sh_type,1,
		       p->verbose,FALSE);
  if(p->compute_geoid)	
    /* make room for geoid solution */
    sh_allocate_and_init(&geoid,1,model->dens_anom[0].lmax,
			 model->sh_type,0,p->verbose,FALSE);
  switch(p->solver_mode){
  case HC_SOLVER_MODE_DEFAULT:
    /* 
       solve poloidal and toroidal part and sum
    */
    hc_solve(model,p->free_slip,p->solution_mode,sol_spectral,
	     TRUE,TRUE,TRUE,p->print_pt_sol,p->compute_geoid,geoid,
	     p->verbose);
    /* 
       
       OUTPUT PART
       
    */
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
    if(p->compute_geoid){
      if(p->compute_geoid_correlations){
	if(p->verbose)
	  fprintf(stderr,"%s: correlation for geoid with %s\n",argv[0],
		  p->ref_geoid_file);
	hc_compute_correlation(geoid,p->ref_geoid,corr,1,p->verbose);
	fprintf(stdout,"%10.7f %10.7f\n",corr[0],corr[1]);
      }else{
	/* 
	   print geoid solution 
	*/
	sprintf(filename,"%s",HC_GEOID_FILE);
	if(p->verbose)
	  fprintf(stderr,"%s: writing geoid to %s\n",argv[0],filename);
	out = ggrd_open(filename,"w","main");
	hc_print_sh_scalar_field(geoid,out,FALSE,FALSE,p->verbose);
	fclose(out);
      }
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
    break;			/* end default mode */
  case HC_SOLVER_MODE_VISC_SCAN: /* scan through viscosities for a
				    four layer problem */
 
    
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
	    fprintf(stdout,"%14.7e %14.7e %14.7e %14.7e\t",p->elayer[3],p->elayer[2],p->elayer[1],p->elayer[0]);
	    /* compute solution */
	    hc_solve(model,p->free_slip,p->solution_mode,sol_spectral,
		     (solved)?(FALSE):(TRUE), /* density changed? */
		     (solved)?(FALSE):(TRUE), /* plate velocity changed? */
		     TRUE,			/* viscosity changed */
		     FALSE,p->compute_geoid,geoid,
		     p->verbose);
	    /* only output are the geoid correlations, for now */
	    hc_compute_correlation(geoid,p->ref_geoid,corr,1,p->verbose);
	    fprintf(stdout,"%10.7f %10.7f ",corr[0],corr[1]);
	    fprintf(stdout,"\n");
	    solved++;
	  }
    break;
  default:
    HC_ERROR("hc","solver mode undefined");
  }
   
  /* 
     
  free memory

  */
  sh_free_expansion(sol_spectral,nsol);
  if(p->compute_geoid)
    sh_free_expansion(geoid,1);

  free(sol_spatial);
  if(p->verbose)
    fprintf(stderr,"%s: done\n",argv[0]);
  hc_struc_free(&model);
  return 0;
}

