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

/* This version performs a grid-search to demonstrate the non-linearity of the inverse problem.
   The code is setup to run a search over viscosity structures defined by four control points.
   
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
  {
    int i1=0;
    while(i1<solution->nlayer-1){
      if(solution->r[i1] >= solution->r[i1+1]){
	fprintf(stderr,"Layer %d cannot have r>= layer %d\n",i1,i1+1);
	exit(-10);
      }
      i1++;
    }

  }
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
  const double drmin = 0.00785;           // minimum layer thickness
  
  const double var_min = 1e-3;
  const double var_max = 1e6;
  const double var_change = 0.05;

  const int max_vor = thb_max_layers( iter );
  const int min_vor = (iter > 10000) ? 3 : 2;//forbid N=2 solutions

  // choose one of five options at random
  int random_choice = 0;

  /* Choose the random option and use rejection sampling to restrict forbidden cases */
  
  if( old_solution->nlayer <= min_vor ){ /* change probabilities for the case where N=2 to remove biases */
    double tmp = randDouble(rng);
    if( p->thb_no_hierarchical ){
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
    if( p->thb_no_hierarchical ){
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
    if( p->thb_no_hierarchical ){
      random_choice = randInt(rng,4);
    }else{
      random_choice = randInt(rng,5);
    }
  }

    
  int success = 0;
  int failcount = 0;
  while(!success){
    failcount++;
    success=1;
    new_solution[0] = old_solution[0]; // copy old to new

    if(random_choice == 0){
      // Add a control point at random
      double new_rad = rad_min + rad_range*randDouble(rng);
      //double new_visc = visc_min + visc_range*randDouble(rng);
      int i=0;
      while(new_solution->r[i] < new_rad && i < new_solution->nlayer)
	i++;
      int j;
      for(j=new_solution->nlayer;j>i;j--){
	new_solution->r[j] = new_solution->r[j-1];
	new_solution->visc[j] = new_solution->visc[j-1];
      }
      if( i<1 || i>new_solution->nlayer-1 ){
	fprintf(stderr,"Error - inserted layer index doesn't make sense.\n");
	exit(-1);
      }
      double dr = new_solution->r[i+1]-new_solution->r[i-1];
      double dvisc = new_solution->visc[i+1]-new_solution->visc[i-1];
      double new_visc = new_solution->visc[i-1] + (new_rad - new_solution->r[i-1])*dvisc/dr + visc_change*randn(rng);
      
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
    if(!success && failcount > 10000){
      fprintf(stderr,"Warning: Fail count %d, random option %d\n",failcount,random_choice);
      new_solution[0] = old_solution[0];
      return;
    }
  }// End while not successful
  if( p-> verbose>=3){
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


int main(int argc, char **argv)
{
  MPI_Init(&argc,&argv);
  int rank,size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);


  
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
  
  struct thb_solution sol1; /* accepted and proposed solutions */
  /* THB residual vectors */
  int thb_nlm=0;
  for(int i=0;i<p->thb_nl;i++)
    thb_nlm += 2*(p->thb_ll[i])+1;

  gsl_vector *residual1;
  gsl_vector *tmp_vector;
  residual1 = gsl_vector_calloc(thb_nlm);
  tmp_vector = gsl_vector_calloc(thb_nlm);
  
  //  double *residual1 = (double *) malloc(sizeof(double)*thb_nlm);
  //double *residual2 = (double *) malloc(sizeof(double)*thb_nlm);
  fprintf(stdout,"Allocated residual [%dx1]\n",thb_nlm);
  gsl_matrix *inverse_covariance_matrix;
  if( p->thb_use_covmat)
    inverse_covariance_matrix = load_inv_covariance_matrix(p);

  
  double chain_temperature;
  // initialize solution
  
  sol1.nlayer = 2;
  sol1.r[0] = model->r_cmb;
  sol1.r[1] = 1.0;
  sol1.visc[0] = 22.0;
  sol1.visc[1] = 22.0;
  sol1.var = 1.0;
  sol1.iter_accepted = -1;
  
  interpolate_viscosity(&sol1,model->rvisc,model->visc,p);

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
  
  /* Calculate the Mahalanobis Distance Phi for the starting model */
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

  /* Sweep through viscosity structures */
  const int visc_layers = 4;
  const int nvisc = 16;
  const int ndepth = 20;
  double *visc_list;
  double *depth_list;
  visc_list = (double *) malloc(nvisc*sizeof(double));
  depth_list = (double *) malloc(ndepth*sizeof(double));
  for(int i=0;i<nvisc;i++){
    visc_list[i] = 19.0 + (25.0-19.0)/((double) nvisc-1)*((double) i);
  }
  for(int i=0;i<ndepth;i++){
    depth_list[i] = model->r_cmb + (1.0-model->r_cmb)/((double) ndepth+1)*((double) i+1);
  }
  
  /* loop over layers */
  /* loop over viscosities for first layer */
  sol1.nlayer=visc_layers;
  sol1.r[0] = model->r_cmb;
  sol1.r[3] = 1.0;
  sol1.var = 1.0;
  FILE *scan_file;
  if( !rank ){
    scan_file = fopen("scan_results.txt","w");
    fprintf(scan_file,"#Results of visc scan: %d layers\n",visc_layers);
    fprintf(scan_file,"#loop order visc0,visc1,r1,visc2,r2,visc3\n");
    fprintf(scan_file,"# %d %d %d %d %d %d\n",nvisc,nvisc,ndepth,nvisc,ndepth,nvisc);
    fflush(scan_file);
  }
  unsigned long int ntotal = nvisc*nvisc*nvisc*nvisc*ndepth*ndepth;
  unsigned long int iter=0;
  for(int v0=0;v0<nvisc;v0++){
    sol1.visc[0] = visc_list[v0];
    for(int v1=0;v1<nvisc;v1++){
      sol1.visc[1] = visc_list[v1];
      for(int r1=0;r1<ndepth;r1++){
	sol1.r[1] = depth_list[r1];
	for(int v2=0;v2<nvisc;v2++){
	  sol1.visc[2] = visc_list[v2];
	  for(int r2=0;r2<ndepth;r2++){
	    sol1.r[2] = depth_list[r2];	   
	    for(int v3=0;v3<nvisc;v3++){
	      sol1.visc[3] = visc_list[v3];
	      
	      if( iter % size == rank ){// obtain this solution
		struct thb_solution sol2 = sol1;
		if( sol2.r[2] < sol2.r[1] ){
		  double tmp = sol2.r[2];
		  sol2.r[2] = sol2.r[1];
		  sol2.r[1] = tmp;
		  tmp = sol2.visc[2];
		  sol2.visc[2] = sol2.visc[1];
		  sol2.visc[1] = tmp;
		}else if(sol2.r[1] == sol2.r[2]){
		  sol2.nlayer = 3;
		  sol2.r[2] = sol2.r[3];
		  sol2.visc[1] = 0.5*(sol2.visc[1]+sol2.visc[2]);
		  sol2.visc[2] = sol2.visc[3];
		}
		//if( sol1.r[2] <= sol1.r[1] ){
		//  sol1.mdist = NAN;
		//  sol1.total_residual = NAN;
		//}else{
		{
		  /* Interpolate solution */
		  interpolate_viscosity(&sol2, model->rvisc, model->visc,p);
		  /* Solve */
		  hc_solve(model,p->free_slip,p->solution_mode,sol_spectral,
			   (solved)?(FALSE):(TRUE), /* density changed? */
			   (solved)?(FALSE):(TRUE), /* plate velocity changed? */
			   TRUE,		    /* viscosity changed */
			   FALSE,p->compute_geoid,
			   pvel,model->dens_anom,geoid,
			   p->verbose); 
		  
		  sol1.total_residual = sh_residual_vector(geoid,p->ref_geoid,p->thb_ll,p->thb_nl,residual1->data,0);
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
		}
		struct thb_solution *solutions = NULL;
		if( !rank ){
		  /* allocate space for the solutions */
		  solutions = (struct thb_solution *) malloc(size*sizeof(struct thb_solution));
		  MPI_Gather(&sol1,sizeof(struct thb_solution),MPI_BYTE,solutions,sizeof(struct thb_solution),MPI_BYTE,0,MPI_COMM_WORLD);
		  for(int irank=0;irank<size;irank++){		
		    /* output */
		    for(int i=0;i<solutions[irank].nlayer;i++)
		      fprintf(scan_file,"%le\t",solutions[irank].r[i]);
		    for(int i=0;i<solutions[irank].nlayer;i++)
		      fprintf(scan_file,"%le\t",solutions[irank].visc[i]);
		    fprintf(scan_file,"%le\t%le\n",(double) solutions[irank].total_residual,solutions[irank].mdist);
		  }
		  free(solutions);
		}else{/* send solution to rank 0*/
		  MPI_Gather(&sol1,sizeof(struct thb_solution),MPI_BYTE,solutions,sizeof(struct thb_solution),MPI_BYTE,0,MPI_COMM_WORLD);
		}
	      }
	      
	      if(!rank && !(iter%1000) )
		fprintf(stderr,"Finished %ld/%ld\n",iter,ntotal);
	      
	      iter++;
	    }
	  }
	}
      }
    }
  }
  if( !rank )
    fclose(scan_file);
  free(depth_list);
  free(visc_list);
  
  /* END THB VISCOSITY LOOP */  
  /*
    
    free memory
    
  */
  gsl_vector_free(residual1);

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

