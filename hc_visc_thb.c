#include "hc.h"
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

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

//helper functions
int randInt(gsl_rng *rng,int n)//random integer in the range 0-n-1, inclusive
{
  int r = (int) gsl_rng_uniform_int(rng, (unsigned long int) n);
  return r;
}

double randDouble(gsl_rng *rng)//random double in the range 0-1
{
  double r = gsl_rng_uniform(rng);// gsl routine to generate random double in [0-1)
  return r;
}

double randn(gsl_rng *rng)//normally distributed random number
{
  double r = gsl_ran_gaussian(rng, 1.0);
  return r;
}

void print_solution(FILE *fh,struct thb_solution *sol){
  /* print a user-readable form of the solution */
  fprintf(fh,"Solution:\nTotal residual, Likelihood %le,%le\nVar, mdist= %le,%le\nNLayers=%d\n",(double) sol->total_residual,(double) sol->likeprob,sol->var,sol->mdist,sol->nlayer);
  fprintf(fh,"r-values:\t");
  for(int i=0;i<sol->nlayer;i++){
    fprintf(fh,"%.6f\t",sol->r[i]);
  }
  fprintf(fh,"\nvisc:\t");
  for(int i=0;i<sol->nlayer;i++){
    fprintf(fh,"%.6f\t",sol->visc[i]);
  }
  fprintf(fh,"\n");
}

void interpolate_viscosity(struct thb_solution *solution,HC_PREC *rvisc,HC_PREC *visc, struct hc_parameters *p){
  const int nlayer = HC_INTERP_LAYERS;
  int i=0;
  int j=0;
  for(i=0;i<nlayer;i++){
    double this_r = rvisc[i];
    // find j such that this_r <= layer_r[j]
    while( solution->r[j+1] < this_r && (j+1) < solution->nlayer )
      j++;
    visc[i] = solution->visc[j]+(solution->visc[j+1]-solution->visc[j])/(solution->r[j+1]-solution->r[j])*(this_r-solution->r[j]);
    visc[i] = pow(10.0,visc[i]);
  }
  if( p->verbose >= 3 ){
    fprintf(stderr,"interpolate_viscosity input:\n");
    for(i=0;i<solution->nlayer;i++){
      fprintf(stderr,"%.3e %.3e\n",(double) solution->r[i],(double) solution->visc[i]);
    }
    fprintf(stderr,"interpolate viscosity output:\n");
    for(i=0;i<HC_INTERP_LAYERS;i++){
      fprintf(stderr,"%.3e %.3e\n",(double) rvisc[i],(double) visc[i]);
    }

  }
}

int thb_max_layers(int iter){
  // Calculate the maximum number of control points/layers allowed according to a burn-in schedule
  // and the maximum number of allowed layers.
  int burnin_steps[MAX_NUM_VOR];

  int i;
  for(i=0;i<MAX_NUM_VOR;i++){
    if(i<3){
      burnin_steps[i] = 0;
    }else{
      burnin_steps[i] = burnin_steps[i-1] + 10000*(i-2);
    }
  }
  int max_layers = 3;
  while(max_layers < MAX_NUM_VOR && burnin_steps[max_layers] < iter){
    max_layers++;
  }
  return max_layers;
}

void propose_solution(struct hcs *model, struct thb_solution *old_solution, struct thb_solution *new_solution,gsl_rng *rng, int iter, struct hc_parameters *p){
  const double visc_min = 19.0;
  const double visc_max = 25.0;
  const double visc_range = visc_max - visc_min;
  const double visc_change = 0.1;
  
  const double rad_min = model->r_cmb;
  const double rad_max = 1.0;
  const double rad_change = 0.05;      // shape parameter for proposal distribution
  const double rad_range = rad_max-rad_min;
  const double drmin = 0.0157;           // minimum layer thickness
  
  const double var_min = 1e-3;
  const double var_max = 1e6;
  const double var_change = 0.05;

  const int max_vor = thb_max_layers( iter );
  const int min_vor = (iter > 10000) ? 3 : 2;//forbid N=2 solutions

  // choose one of five options at random
  int random_choice = 0;

  /* Choose the random option and use rejection sampling to restrict forbidden cases */
  do{
    if( old_solution->nlayer <= min_vor ){ /* change probabilities for the case where N=2 to remove biases */
      double tmp = randDouble(rng);
      if( 0 && p->thb_no_hierarchical ){
	if( tmp < 0.25 ){
	  random_choice = 0;
	}else{
	  random_choice = 3;
	}
      }else{/* hierarchical case */
	if( tmp < 0.2 ){
	  random_choice = 0;
	}else if(tmp<0.6){
	  random_choice = 3;
	}else{
	  random_choice = 4;
	}
      }    
    }else if( old_solution->nlayer >= max_vor ){/* special case for N=Nmax */
      double tmp = randDouble(rng);
      if( 0 && p->thb_no_hierarchical ){
	if( tmp < 0.25 ){
	  random_choice = 1;
	}else if(tmp < 0.25 + 0.75/2.0){
	  random_choice = 2;
	}else{
	  random_choice = 3;
	}
      }else{
	/* hierarchical case */	
	if( tmp < 0.2 ){
	  random_choice = 1;
	}else if(tmp < 0.2+0.8/3.0){
	  random_choice = 2;
	}else if(tmp < 0.2+2.0*0.8/3.0){
	  random_choice = 3;
	}else{
	  random_choice = 4;
	}
      }
    }else{
      /* Normal case - all options have equal probability */
      if( 0 && p->thb_no_hierarchical ){
	random_choice = randInt(rng,4);
      }else{
	random_choice = randInt(rng,5);
      }
    }
  } while( !p->thb_no_hierarchical || (p->thb_no_hierarchical && random_choice == 4 ) );
    
  int success = 0;
  int failcount = 0;
  while(!success){
    failcount++;
    success=1;
    new_solution[0] = old_solution[0]; // copy old to new

    if(random_choice == 0){
      // Add a control point at random
      double new_rad = rad_min + rad_range*randDouble(rng);
      double new_visc = visc_min + visc_range*randDouble(rng);
      int i=0;
      while(new_solution->r[i] < new_rad && i < new_solution->nlayer)
	i++;
      int j;
      for(j=new_solution->nlayer;j>i;j--){
	new_solution->r[j] = new_solution->r[j-1];
	new_solution->visc[j] = new_solution->visc[j-1];
      }
      new_solution->r[i] = new_rad;
      new_solution->visc[i] = new_visc;
      new_solution->nlayer++;      
    }else if(random_choice == 1){
      // kill a layer (at random)
      if( new_solution->nlayer == 3){
	int kill_layer = 1;
	new_solution->r[1] = new_solution->r[2];
	new_solution->visc[1] = new_solution->visc[2];
	new_solution->nlayer--;
      }else{
	int kill_layer = randInt(rng,new_solution->nlayer-2) + 1;
	int j;
	for(j=kill_layer;j<new_solution->nlayer-1;j++){
	  new_solution->r[j] = new_solution->r[j+1];
	  new_solution->visc[j] = new_solution->visc[j+1];
	}
	new_solution->nlayer--;
      }
    }else if(random_choice == 2){
      // perturb the location of a control point
      if(new_solution->nlayer == 2){
	success = 0;
      }else{
	int perturb_layer = randInt(rng,new_solution->nlayer-2)+1;
	new_solution->r[ perturb_layer ] += rad_change*randn(rng);
	if(new_solution->r[ perturb_layer] <= model-> r_cmb || new_solution->r[ perturb_layer ] >= 1.0){
	  success = 0;
	}
      }
    }else if(random_choice == 3){
      //perturb viscosity
      int perturb_layer = randInt(rng,new_solution->nlayer);
      new_solution->visc[perturb_layer] += visc_change*randn(rng);
      if( new_solution->visc[perturb_layer] > visc_max || new_solution->visc[perturb_layer] < visc_min){
	success = 0;
      }
    }else if(random_choice == 4){
      //change variance
      new_solution->var = pow(10.0, log10(new_solution->var) + var_change*randn(rng));
      if( new_solution->var > var_max || new_solution->var < var_min){
	success = 0;
      }
    }else{
      fprintf(stderr,"Warning: this random choice isn't implemented\n");
      success = 0;
    }
    // any additional sanity checks on the proposed solution:
    int j;
    for(j=1;j<new_solution->nlayer;j++){
      double dr = new_solution->r[j]-new_solution->r[j-1];
      if( dr < drmin ){
	success = 0; // This ensures that the nodes remain in increasing order
      }
    }
    if( failcount > 10000){
      fprintf(stderr,"Error: Fail count %d, random option %d\n",failcount,random_choice);
      exit(-1);
    }
  }// End while not successful
  if( p-> verbose){
    fprintf(stderr,"Proposed solution:\n");
    int i;
    for(i=0;i<new_solution->nlayer;i++){
      fprintf(stderr,"%.3e %.3e\n",new_solution->r[i],new_solution->visc[i]);
    }
  }
}

/* COVARIANCE MATRIX STUFF */
gsl_matrix * load_inv_covariance_matrix(struct hc_parameters *p){
  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  gsl_matrix *matrix;
  FILE *fh = NULL;
  int i1=0;
  if(p->verbose){
    fprintf(stdout,"Reading covariance matrix from %s\n",p->thb_covmat_filename);
  }
  fh = fopen(p->thb_covmat_filename,"r");
  double det;
  i1 += fscanf(fh,"%le\n",&det);
  int nshdeg;
  i1 += fscanf(fh,"%d\n",&nshdeg);
  if( nshdeg != p->thb_nl || nshdeg>10 ){
    fprintf(stderr,"Error: mismatch in number of spherical harmonic degrees (%d vs %d)!\n",nshdeg,p->thb_nl);
    exit(-1);
  }

  int sh_degs[10];
  for(int i=0;i<nshdeg;i++){
    i1 += fscanf(fh,"%d\t",&sh_degs[i]);
    if(sh_degs[i] != p->thb_ll[i]){
      i1 += fprintf(stderr,"Error: mismatch in spherical harmonic degrees!\n");
      fclose(fh);
      exit(-1);
    }
  }
  int M,N;
  i1 += fscanf(fh,"%d %d\n",&M,&N);
  matrix = gsl_matrix_alloc(M,N);
  /* load the data from the file */
  for(int j=0;j<N;j++){
    for(int i=0;i<M;i++){
      i1 += fscanf(fh,"%le\t",&(matrix->data[i*matrix->tda+j]));     
    }
    i1 += fscanf(fh,"\n");
  }
  fclose(fh);
  if(p->verbose){
    for(int irank=0;irank<size;irank++){
      MPI_Barrier(MPI_COMM_WORLD);
      if(irank == rank){
	fprintf(stdout,"[%d] Cdinv:\n",rank);
	//gsl_matrix_fprintf(stdout,matrix,"%e");
	for(int i=0;i<N;i++){
	  for(int j=0;j<N;j++){
	    fprintf(stdout,"%le,",matrix->data[i*matrix->tda+j]);
	  }
	  fprintf(stdout,"\n");
	}
      }
    }
  }
  return matrix;
}


void thb_checkpoint(struct thb_solution *sol, double chain_temperature, int iter, struct hc_parameters *p){
  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  struct thb_solution *all_solutions;
  double *all_temperatures;
  if(!rank){   
    char checkpoint_file[HC_CHAR_LENGTH+7];
    sprintf(checkpoint_file,"%s.chkpt",p->thb_ensemble_filename);
    
    FILE *fh = NULL;
    fh = fopen(checkpoint_file,"w");
    /* gather solutions from all ranks */
    all_solutions = (struct thb_solution *) malloc( sizeof(struct thb_solution)*size );
    all_temperatures = (double *) malloc( sizeof(double)*size );
    MPI_Gather(sol,sizeof(struct thb_solution),MPI_BYTE,all_solutions,sizeof(struct thb_solution),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Gather(&chain_temperature,1,MPI_DOUBLE,all_temperatures,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    fprintf(fh,"%d\n",iter);
    for(int i=0;i<size;i++){
      struct thb_solution *sol = all_solutions+i;
      fprintf(fh,"%d\t%le\t",i,all_temperatures[i]);
      fprintf(fh,"%le\t%le\t",sol->var,sol->mdist);
      fprintf(fh,"%d\t",sol->nlayer);
      fprintf(fh,"%le\t%le\t",(double) sol->total_residual,(double) sol->likeprob);
      for(int j=0;j<sol->nlayer;j++){
	/* r-values */
	fprintf(fh,"%le\t",sol->r[j]);
      }
      for(int j=0;j<sol->nlayer;j++){
	/* visc values */
	fprintf(fh,"%le\t",sol->visc[j]);       
      }
      fprintf(fh,"\n");
    }
    fclose(fh);
    free(all_solutions);
    free(all_temperatures);
  }else{
    /* send solution to rank 0 */
    MPI_Gather(sol,sizeof(struct thb_solution),MPI_BYTE,all_solutions,sizeof(struct thb_solution),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Gather(&chain_temperature,1,MPI_DOUBLE,all_temperatures,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  }
  if( p->verbose ){
    fprintf(stderr,"Checkpoint written\n");
  }
}

void thb_checkpoint_resume(struct thb_solution *sol, double *chain_temperature, int *iter, struct hc_parameters *p){
  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  struct thb_solution *all_solutions;
  double *all_temperatures;
  if(!rank){
    /* open checkpoint file*/
    char checkpoint_file[HC_CHAR_LENGTH+7];
    sprintf(checkpoint_file,"%s.chkpt",p->thb_ensemble_filename);
    FILE *fh = NULL;
    fh = fopen(checkpoint_file,"r");
    if(fh == NULL){
      perror("Error opening checkpoint file");
      exit(-1);
    }
    all_solutions = (struct thb_solution *) malloc( sizeof(struct thb_solution)*size );
    all_temperatures = (double *) malloc( sizeof(double)*size );
    int i1=0;
    i1 += fscanf(fh,"%d\n",iter);
    for(int i=0;i<size;i++){
      struct thb_solution *sol = all_solutions+i;
      int rank_tmp;
      i1+=fscanf(fh,"%d\t%le\t",&rank_tmp,all_temperatures+i);
      i1+=fscanf(fh,"%le\t%le\t",&(sol->var),&(sol->mdist));
      i1+=fscanf(fh,"%d\t",&(sol->nlayer));
      double residual_tmp, likeprob_tmp;
      i1+=fscanf(fh,"%le\t%le\t", &residual_tmp, &likeprob_tmp);
      sol->total_residual = (HC_PREC) residual_tmp;
      sol->likeprob = (HC_PREC) likeprob_tmp;
      for(int j=0;j<sol->nlayer;j++){
	/* r-values */
	i1+=fscanf(fh,"%le\t",&(sol->r[j]));
      }
      for(int j=0;j<sol->nlayer;j++){
	/* visc values */
	i1+=fscanf(fh,"%le\t",&(sol->visc[j]));       
      }
      i1+=fscanf(fh,"\n");
    }
    chain_temperature[0] = all_temperatures[rank];
    sol[0] = all_solutions[rank];
    fclose(fh);
  }
  /* Collective scatter of solutions */
  MPI_Bcast(iter,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Scatter(all_solutions,sizeof(struct thb_solution),MPI_BYTE,sol,sizeof(struct thb_solution),MPI_BYTE,0,MPI_COMM_WORLD);
  MPI_Scatter(all_temperatures,1,MPI_DOUBLE,chain_temperature,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  if(!rank){
    free(all_temperatures);
    free(all_solutions);
  }
  if(p->verbose){
    fprintf(stdout,"Checkpoint restored!\n");
    for(int i=0;i<size;i++){
      MPI_Barrier(MPI_COMM_WORLD);
      if( i== rank ){
	fprintf(stdout,"Solution for rank %d:\n",rank);
	fprintf(stdout,"iter=%d\n",*iter);
	print_solution(stdout,sol);
      }
    }
  }
}



int main(int argc, char **argv)
{
  MPI_Init(&argc,&argv);
  int rank,size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if(size % 2){
    fprintf(stderr,"Error: MPI size must be even!\n");
    exit(-1);
  }

  
  struct hcs *model;		/* main structure, make sure to initialize with 
				   zeroes */
  struct sh_lms *sol_spectral=NULL, *geoid = NULL;		/* solution expansions */
  struct sh_lms *pvel=NULL;					/* local plate velocity expansion */
  int nsol,lmax,solved;
  struct hc_parameters p[1]; /* parameters */
  
  /* Initialize random number generation using GNU Scientific Library */
  
  const gsl_rng_type *rng_type;
  gsl_rng *rng;
  gsl_rng_default_seed = (unsigned long int) rank; /* ensure that each mpi rank has a unique random seed */
  gsl_rng_env_setup();
  rng_type = gsl_rng_default;
  rng = gsl_rng_alloc (rng_type);
  for(int i=0;i<rank;i++) gsl_rng_get(rng); /* burn rank-many random numbers to ensure that each MPI process takes a different path */
  printf ("[%d]: first random value = %lu\n",rank, gsl_rng_get (rng));

  
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
  p->visc_init_mode = HC_INIT_E_INTERP;
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

  /* BEGIN THB VISCOSITY LOOP */
  
  struct thb_solution sol1,sol2; /* accepted and proposed solutions */
  /* THB residual vectors */
  int thb_nlm=0;
  for(int i=0;i<p->thb_nl;i++)
    thb_nlm += 2*(p->thb_ll[i])+1;

  gsl_vector *residual1;
  gsl_vector *residual2;
  gsl_vector *tmp_vector;
  residual1 = gsl_vector_calloc(thb_nlm);
  residual2 = gsl_vector_calloc(thb_nlm);
  tmp_vector = gsl_vector_calloc(thb_nlm);
  
  //  double *residual1 = (double *) malloc(sizeof(double)*thb_nlm);
  //double *residual2 = (double *) malloc(sizeof(double)*thb_nlm);
  fprintf(stdout,"Allocated residual [%dx1]\n",thb_nlm);
  gsl_matrix *inverse_covariance_matrix;
  if( p->thb_use_covmat)
    inverse_covariance_matrix = load_inv_covariance_matrix(p);

  
  double chain_temperature;
  int iter;
  // initialize solution
  if( !p->thb_checkpoint_resume){
    iter = 0;
    sol1.nlayer = 2;
    sol1.r[0] = model->r_cmb;
    sol1.r[1] = 1.0;
    sol1.visc[0] = 22.0;
    sol1.visc[1] = 22.0;
    sol1.var = 1.0;
    sol1.iter_accepted = -1;
    if(p->thb_parallel_tempering){
      /* space chain temperatures log-uniformly from 1 to 10 */
      if( !rank ){
	chain_temperature = 1.0;
      }else{      
	chain_temperature = pow(10.0, ((double) rank)/((double) (size-1))*(log10(10.0)-log10(1.0)) );
      }
    }else{
      chain_temperature = 1.0;
    }
  }else{
    thb_checkpoint_resume(&sol1,&chain_temperature,&iter,p);
    if(p->verbose){
      print_solution(stderr,&sol1);
    }
  }  
  interpolate_viscosity(&sol1,model->rvisc,model->visc,p);

  fprintf(stdout,"[%d]: made it here, temperature is %le\n",rank,chain_temperature); fflush(stdout); MPI_Barrier(MPI_COMM_WORLD);
  p->thb_save_skip = (p->thb_iter - p->thb_save_start)/p->thb_sample_target;
  FILE *thb_ensemble_file = NULL;
  
  if(!rank && !p->thb_checkpoint_resume){
    thb_ensemble_file = fopen(p->thb_ensemble_filename,"w");
    {
      fprintf(thb_ensemble_file,"#HC version %s\n",HC_VERSION);
      fprintf(thb_ensemble_file,"#HC invoked with the following options\n#");
      for(int i=0;i<argc;i++)
	fprintf(thb_ensemble_file,"%s ",argv[i]);
      fprintf(thb_ensemble_file,"\n");
      fprintf(thb_ensemble_file,"#THB Viscosity parameters:\n");
      fprintf(thb_ensemble_file,"#Maximum iterations: %d\n",p->thb_iter);
      fprintf(thb_ensemble_file,"#Save start %d\n#Save skip %d\n#Sample target %d\n",p->thb_save_start,p->thb_save_skip,p->thb_sample_target);
      fprintf(thb_ensemble_file,"#Number of chains: %d\n",size);
      fprintf(thb_ensemble_file,"#Parallel tmpering: %s\n",p->thb_parallel_tempering ? "True" : "False");
      fprintf(thb_ensemble_file,"#Covariance Matrix: %s\n",p->thb_use_covmat ? p->thb_covmat_filename : "False");
      
      fprintf(thb_ensemble_file,"#Parallel tempering start %d\n",p->thb_swap_start);
      fprintf(thb_ensemble_file,"#Parallel tempering steps for swap: %d\n",p->thb_steps_for_swap);
      fprintf(thb_ensemble_file,"#NL=%d, L=[%d",p->thb_nl,p->thb_ll[0]);

      int i;
      for(i=1;i<p->thb_nl;i++){
	fprintf(thb_ensemble_file,",%d",p->thb_ll[i]);
      }
      fprintf(thb_ensemble_file,"]\n");
      
      fprintf(thb_ensemble_file,"#Burnin schedule follows:\n#iter max layers\n");
      int max_layers=0;
      i=0;
      while(max_layers < MAX_NUM_VOR){
	if( thb_max_layers(i) > max_layers ){
	  max_layers = thb_max_layers(i);
	  fprintf(thb_ensemble_file,"#%d, %d\n",i,max_layers);
	}
	i++;
      }
      fflush(thb_ensemble_file);
    }   
    /* Write viscosity file header */
    fprintf(thb_ensemble_file,"#rank\t,iter\t,iter_accepted\t,total_residual\t,variance_reduction\t,likeprob\t,var\t,nlayer\t,rr\t,visc\n");
  }else if(!rank && p->thb_checkpoint_resume){
    thb_ensemble_file = fopen(p->thb_ensemble_filename,"a");
    fprintf(thb_ensemble_file,"# Resumed from checkpoint\n");
  }

  /* Calculate total variance at the included spherical harmonic degrees */
  double total_variance = sh_residual_vector(geoid,p->ref_geoid,p->thb_ll,p->thb_nl,residual1->data,0);
  if(!rank) fprintf(stdout,"Total Variance %le\n",total_variance);
  
  /* select plate velocity */
  if(!p->free_slip)
    hc_select_pvel(p->pvel_time,&model->pvel,pvel,p->verbose);
  // do the initial solution
  solved=0;
  if( !p->thb_sample_prior)
    hc_solve(model,p->free_slip,p->solution_mode,sol_spectral,
	     (solved)?(FALSE):(TRUE),   /* density changed? */
	     (solved)?(FALSE):(TRUE),   /* plate velocity changed? */
	     TRUE,			/* viscosity changed */
	     FALSE,p->compute_geoid,
	     pvel,model->dens_anom,geoid,
	     p->verbose);
  solved=1;
  //hc_compute_correlation(geoid,p->ref_geoid,corr,1,p->verbose);
  sol1.total_residual = sh_residual_vector(geoid,p->ref_geoid,p->thb_ll,p->thb_nl,residual1->data,0);
  
  /* Calculate the Mahalanobis Distance Phi */
  {
    double mdist = 0.0;    
    int i;
    if( p->thb_sample_prior ){
      sol1.mdist = 0.0;
    }else{
      if( p->thb_use_covmat ){
	/* form r'*Cd^-1*r */
	double rtr;
	gsl_blas_dgemv(CblasNoTrans,1.0,inverse_covariance_matrix,residual1,0.0,tmp_vector);
	gsl_blas_ddot(residual1,tmp_vector,&rtr);
	sol1.mdist = rtr/sol1.var;
      }else{
	double rtr;
	gsl_blas_ddot(residual1,residual1,&rtr);
	sol1.mdist = rtr/sol1.var;
      }
    }
    sol1.likeprob = -0.5 * sol1.mdist/chain_temperature;
  }

  /* BEGIN THB LOOP */
  while(iter < p->thb_iter){
    /* ------------------- */
    /* 1. Propose Solution */
    /* ------------------- */
    propose_solution( model, &sol1, &sol2, rng, iter, p);
    interpolate_viscosity(&sol2, model->rvisc, model->visc,p);

    /* --------------------------- */
    /* 2A. Calculate forward model */
    /* --------------------------- */
    if(!p->thb_sample_prior)
      hc_solve(model,p->free_slip,p->solution_mode,sol_spectral,
	       (solved)?(FALSE):(TRUE), /* density changed? */
	       (solved)?(FALSE):(TRUE), /* plate velocity changed? */
	       TRUE,			/* viscosity changed */
	       FALSE,p->compute_geoid,
	       pvel,model->dens_anom,geoid,
	       p->verbose); 

    sol2.total_residual = sh_residual_vector(geoid,p->ref_geoid,p->thb_ll,p->thb_nl,residual2->data,0);
    /* ------------------------------------------ */
    /* 2B. Calculate the Mahalanobis Distance Phi */
    /* ------------------------------------------ */
    {
      double mdist = 0.0;    
      int i;
      if( p->thb_sample_prior ){
	sol2.mdist = 0.0;
      }else{
	if( p->thb_use_covmat ){
	  /* form r'*Cd^-1*r */
	  double rtr;
	  gsl_blas_dgemv(CblasNoTrans,1.0,inverse_covariance_matrix,residual2,0.0,tmp_vector);
	  gsl_blas_ddot(residual2,tmp_vector,&rtr);
	  sol2.mdist = rtr/sol2.var;
	}else{
	  double rtr;
	  gsl_blas_ddot(residual2,residual2,&rtr);
	  sol2.mdist = rtr/sol2.var;
	}
      }
      sol2.likeprob = -0.5 * sol2.mdist/chain_temperature;
    }
    /* ------------------------------------------ */
    /* 3. Calculate the probablity of acceptance: */
    /* ------------------------------------------ */
    int k2 = sol2.nlayer;
    int k1 = sol1.nlayer;
    double varfakt = sol2.var / sol1.var;
    double prefactor = p->thb_no_hierarchical ? 0.0 : -0.5 * ((double) thb_nlm)*log(varfakt);
    double probAccept = prefactor + sol2.likeprob - sol1.likeprob + log((double) (k1)) - log((double) (k2));
    if( probAccept > 0 || probAccept > log(randDouble(rng))){
      // Accept the proposed solution
      sol1 = sol2;
      if( p->verbose ){
	fprintf(stdout,"Accepted new solution:\n");
	int i;
	for(i=0;i<sol2.nlayer;i++){
	  fprintf(stdout,"%le %le\n",sol2.r[i],sol2.visc[i]);
	}
	fprintf(stdout,"Residual: %le\n",(double) sqrt(sol2.total_residual));
      }
      sol1.iter_accepted = iter;
    }

    /* --------------------------------------- */
    /* 4. Consider swap for parallel tempering */
    /* --------------------------------------- */
    
    if( p->thb_parallel_tempering && iter > p->thb_swap_start && !(iter % p->thb_steps_for_swap)){
      int *ranks;
      ranks = (int *) malloc(sizeof(int)*size);
      if(!rank){
	for(int i=0;i<size;i++) ranks[i] = i;
	gsl_ran_shuffle(rng,ranks,size,sizeof(int));/* random permutation of the ranks */
      }
      int *partners;
      partners = (int *) malloc(sizeof(int)*size);
      if(!rank){
	for(int i=0;i<size;i++) partners[i] = i;
	int i=0;
	for(int i=0;i<size-1;i+=2){
	  int swapi = ranks[i];
	  int swapj = ranks[i+1];	    
	  partners[swapi] = swapj;
	  partners[swapj] = swapi;	    
	  if( p->verbose >= 3 )
	    fprintf(stderr,"Propose swap %d <--> %d\n",swapi,swapj);
	}	
      }
      MPI_Bcast(partners,size,MPI_INT,0,MPI_COMM_WORLD);
      int swapi = rank;
      int swapj = partners[rank];
      double data[4];
	 
      if( swapi != swapj ){
	if( swapi < swapj ){
	  /* chain with the lower rank calculates the probability of swap */
	  /* get Varj, temperaturej, mdistj */
	  MPI_Recv(data,4,MPI_DOUBLE,swapj, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  double varj = data[0];
	  double vari = sol1.var;
	  double mdistj = data[1];
	  double mdisti = sol1.mdist;
	  double Tj = data[2];
	  double Ti = chain_temperature;
	  double kj = data[3];
	  double ki = (double) sol1.nlayer;
	  double VarfaktS = vari/varj;
	  double prefactor = (p->thb_no_hierarchical) ? 0.0 : (1.0/Ti-1.0/Tj)*thb_nlm/2.0*log(VarfaktS);
	  double probAcceptSwap = prefactor + (1.0/Ti-1.0/Tj)*(-0.5*(mdistj-mdisti) + log(ki) - log(kj));
	  if( probAcceptSwap > 0 || probAcceptSwap > log(randDouble(rng)) ){
	    /* accept the swap */
	    if( p->verbose >= 3) fprintf(stderr,"Swapping %d <--> %d. mdists (%le,%le) prob=%le\n",swapi,swapj,sol1.mdist,mdistj,probAcceptSwap);
	    chain_temperature = Tj;
	    Tj = Ti;
	  }
	  MPI_Send(&Tj,1,MPI_DOUBLE,swapj,0,MPI_COMM_WORLD);	 
	}else{
	  /* chain with the higher rank sends stuff for swap to be calculated */
	  data[0] = sol1.var;
	  data[1] = sol1.mdist;
	  data[2] = chain_temperature;
	  data[3] = (double) sol1.nlayer;
	  MPI_Send(data,4,MPI_DOUBLE,swapj, 0, MPI_COMM_WORLD);             // Send data required to evaluate the swap probability
	  MPI_Recv(&chain_temperature,1,MPI_DOUBLE,swapj,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE); // Receive the temperature from partner
	}
	free(partners);
	free(ranks);
      }
      /* update the likeprob for the accepted solution using the new chain temperature */
      sol1.likeprob = -0.5*sol1.mdist/chain_temperature;
    }
    /* ------------------------------------ */
    /* 5. Output Ensemble and Checkpointing */
    /* ------------------------------------ */
    
    /* update progress */
    if(!rank && (!(iter%1000) || p->verbose>=2))
      fprintf(stderr,"[%d]: Finished iteration %d\n",rank,iter);
    
    //save ensemble   
    if( iter >= p->thb_save_start && !(iter%p->thb_save_skip)){
      struct thb_solution *all_solutions;
      double *all_temperatures;
      if(!rank){
	all_solutions = (struct thb_solution *) malloc( sizeof(struct thb_solution)*size );
	all_temperatures = (double *) malloc( sizeof(double)*size );
      }
      MPI_Gather(&sol1,sizeof(struct thb_solution),MPI_BYTE,all_solutions,sizeof(struct thb_solution),MPI_BYTE,0,MPI_COMM_WORLD);
      MPI_Gather(&chain_temperature,1,MPI_DOUBLE,all_temperatures,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      if(!rank){
	int write_success = 0;
	for(int irank=0;irank<size;irank++){
	  if( all_temperatures[irank] == 1.0 ){
	    double variance_reduction = 1.0 - pow((sqrt(all_solutions[irank].total_residual)/sqrt(total_variance)),2.0);
	    fprintf(thb_ensemble_file,"%02d,%08d,%08d,%.6le,%.6le,%le,%le,%02d",irank,iter,all_solutions[irank].iter_accepted,sqrt(all_solutions[irank].total_residual),variance_reduction,(double) all_solutions[irank].likeprob,all_solutions[irank].var,all_solutions[irank].nlayer);
	    for(int i=0;i<all_solutions[irank].nlayer;i++) fprintf(thb_ensemble_file,",%le",all_solutions[irank].r[i]);
	    for(int i=0;i<all_solutions[irank].nlayer;i++) fprintf(thb_ensemble_file,",%le",all_solutions[irank].visc[i]);
	    fprintf(thb_ensemble_file,"\n");
	    write_success++;
	  }
	}
	if( !write_success ){
	  fprintf(stderr,"Warning: Nothing written to ensemble file\n");
	  fprintf(stderr,"Chain temperatures:");
	  for(int i=0;i<size;i++) fprintf(stderr,"%.8le (%d) ",all_temperatures[i],(int) all_temperatures[i]==1.0);
	  fprintf(stderr,"\n");
	}	
	free(all_solutions);
	free(all_temperatures);
      }
    }
    if( !(iter % THB_CHECKPOINT_INTERVAL) )
      thb_checkpoint(&sol1,chain_temperature,iter,p);
    
    iter++;
  }
  if( !rank ){
    fflush(thb_ensemble_file);
    fclose(thb_ensemble_file);
  }
  /* END THB VISCOSITY LOOP */  
  /*
    
    free memory
    
  */
  gsl_vector_free(residual1);
  gsl_vector_free(residual2);
  gsl_vector_free(tmp_vector);
  if( p->thb_use_covmat )
    gsl_matrix_free(inverse_covariance_matrix);
  
  sh_free_expansion(sol_spectral,nsol);
  /* local copies of plate velocities */
  sh_free_expansion(pvel,2);
  /*  */
  sh_free_expansion(geoid,1);
  if(p->verbose)
    fprintf(stderr,"%s: done\n",argv[0]);
  hc_struc_free(&model);

  MPI_Finalize();
  return 0;
}

