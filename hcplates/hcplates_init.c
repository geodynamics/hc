#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "hcplates.h"


/* 
   
   input file for Lithgow-Bertelloni's & Richards Plates extensions to Hager & O'Connell experimental code
   O'Neill (coneill@els.qm.edu.au) & Becker (twb@usc.edu).

   This has routine to: set initial defaults
			handle the command line args (ie. just the input file name)
			read the input file, and assign parameters
  
*/
void hc_init_params(struct hc_plates_params *plates) 
{
	int NPLT, LDIM, dd, m, i;
   	
	HC_PREC r2[] = {0.0, 3480.0e5, 5700.0e5, 6140.0e5, 6240.0e5,6370.0e5};  /* Note we adopt the fotran indices scheme so the 0 entry in the array is just a decoy! */
 	HC_PREC visc2[] = {0, 50, 1, 1, 10};
	HC_PREC pparea[] = {0, 0.77691E+14, 0.59987E+14, 0.49031E+13, 0.36125E+13, 0.31173E+13, 0.68647E+14, 0.61026E+14, 0.58440E+14, 0.16462E+14, 0.10806E+15, 0.57208E+13, 0.42244E+14};
	

	NPLT = 11;
	LDIM = 25;
	printf("Assigning default parameters first\n");
	plates->RCore = 3480.0e5;	/* Why the e5? Its in cm's? */
	plates->erad = 6370.0e5;	/*Earth radius */
	plates->visc0 = 1.0e22;
	plates->NPLT = NPLT;
	plates->LDIM = LDIM; 	/* Absolute max */
	plates->NPDIM = 3 * NPLT;
	plates->LDIM1 = LDIM + 1;
	plates->KDIM = (LDIM + 2)*(LDIM+1)/2;
	
	plates->dd=1;
	dd = plates->dd;
	plates->NLAT = 180/dd;
	plates->NLONG = 360/dd;

	plates->stoy = 365.0*24.0*60.0*60.0;

	/* EArth information - to be imported. Defaults are from CRLB present day file */
	plates->Lmax = 20;
	plates->n = 4;
	plates->iba = 0; /*Always want no-slip for the plates calculation */

	plates->Lload = plates->Lmax;
	plates->Lplt = plates->Lmax;

	m=plates->n + 1;
	plates->r[0] = 0.0;
	plates->visc[0] = 0.0;
	for (i=1;i<m;i++) {
		plates->r[i] = r2[i];
		plates->visc[i]= visc2[i];
	}
	plates->r[m] = r2[m];
	
	for (i=1; i<=plates->NPLT; i++) {
		plates->parea[i] = pparea[i];
	}
	
	/* Init filename defaults */
	strncpy(plates->loadfile,"point.j",HC_CHAR_LENGTH);
	strncpy(plates->unitrotfile,"unitrot.coeffs",HC_CHAR_LENGTH);
	strncpy(plates->platemapfile,"map.plate",HC_CHAR_LENGTH);
	strncpy(plates->outputfile,"hcplates.output",HC_CHAR_LENGTH);
	strncpy(plates->polesfile,"hcpoles.output",HC_CHAR_LENGTH);
	strncpy(plates->velgridfile,"hc_velgrid.output",HC_CHAR_LENGTH);
	
	plates->ratio = 1.0;  /* This is rigged from CRLB's final answer which depends on observed plate rotations*/
								/* We want our answers independent of observations, but will keep this scaling factor for consistency with old fortran */
	printf(" NLAT = %i \n",plates->NLAT);
}

double *c_alloc1d(int n)
{
	double *dp = (double *)malloc(n * sizeof(double));
	if(dp == NULL) {
		fprintf(stderr, "out of memory\n");
		exit(1);
	}
	return dp;
}

double **c_alloc2d(int m, int n)
{
	int i;
	double **dpp = (double **)malloc(m * sizeof(double *));
	if(dpp == NULL) {
		fprintf(stderr, "out of memory\n");
		exit(1);
	}
	for(i = 0; i < m; i++) {
		if((dpp[i] = (double *)malloc(n * sizeof(double))) == NULL) {
			fprintf(stderr, "out of memory\n");
			exit(1);
		}
	}
	return dpp;
}

int **i_alloc2d(int m, int n)
{
	int i;
	int **dpp = (int **)malloc(m * sizeof(int *));
	if(dpp == NULL) {
		fprintf(stderr, "out of memory\n");
		exit(1);
	}
	for(i = 0; i < m; i++) {
		if((dpp[i] = (int *)malloc(n * sizeof(int))) == NULL) {
			fprintf(stderr, "out of memory\n");
			exit(1);
		}
	}
	return dpp;
}


void Free1D(double *data)
 {
  if(data != NULL)
   {
    free(data);
    data = NULL;
   }
 }

void c_free2d(double **dpp, int m)
{
	int i;
	for(i = 0; i < m; i++)
		free((void *)dpp[i]);

	free((void *)dpp);
}


void Free2D(double **data)
 {
  Free1D(*data);
  Free1D((double *) data);
 }


void Free3D(double ***data)
 {
  Free1D(**data);
  Free2D((double **) data);
 }


/* 
  Initializes the allocated space to 0 .
  Returns NULL if there is not enough memory 
  or dim1 <= 0 or dim2 <= 0 or dim3 <= 0 .
*/

double ***c_alloc3d(int dim1, int dim2, int dim3)
 {
  int i, j;
  int offset1, offset2;
  double ***data;

  if(dim1 <= 0 || dim2 <= 0 || dim3 <= 0)
   return NULL;

  data = (double ***) malloc(dim1 * sizeof(double **));

  if(data == NULL)
   return NULL;

  *data = (double **) malloc(dim1 * dim2 * sizeof(double *));

  if(*data == NULL)
   return NULL;
  
  **data = (double *) calloc(dim1 * dim2 * dim3, sizeof(double));  //typo?

  if(**data == NULL)
   return NULL;

  for(i = 0; i < dim1; i++)
   {
    offset1 = i * dim2;
    offset2 = offset1 * dim3;
    data[i] = (*data) + offset1;

    for(j = 0; j < dim2; j++)
     (*data)[offset1 + j] = (**data) + offset2 + j * dim3;
   }

  return data;
 }


void hc_init_arrays(struct hc_plates_params *plates, struct hc_plates_arrays *A) 
{
  /* Allocate shear stress coefficient arrays */
  A->y4in = c_alloc2d(plates->KDIM+1,5);
  A->y2in = c_alloc2d(plates->KDIM+1,5);
  A->y3in = c_alloc2d(plates->KDIM+1,5);
  A->z3in = c_alloc2d(plates->KDIM+1,5);

  A->y4ex = c_alloc3d(plates->NPDIM+1,plates->KDIM+1,4);
  A->y10ex = c_alloc3d(plates->NPDIM+1,plates->KDIM+1,4);

  /* Allocate spherical harmonic arrays */

  A->p = c_alloc1d(plates->KDIM+1);
  A->dpdt = c_alloc1d(plates->KDIM+1);
  A->pbyst = c_alloc1d(plates->KDIM+1);
  A->torp = c_alloc2d(plates->NLAT+1,plates->KDIM+1);
  A->cm = c_alloc1d(plates->LDIM1+1);
  A->sm = c_alloc1d(plates->LDIM1+1);

  /* Allocate torque arrays */

  A->fin = c_alloc1d(plates->NPDIM+1);
  A->rots = c_alloc1d(plates->NPDIM+1);
  A->finet = c_alloc1d(plates->NPDIM+1);
  A->fex = c_alloc2d(plates->NPDIM+1,plates->NPDIM+1);
  A->srt = c_alloc1d(plates->NPDIM+1);
  A->srp = c_alloc1d(plates->NPDIM+1);
  A->fins = c_alloc1d(plates->NPDIM+1);
  A->fexs = c_alloc2d(plates->NPDIM+2,plates->NPDIM+2);
  A->ww = c_alloc1d(plates->NPDIM+1);
  A->vv = c_alloc2d(plates->NPDIM+1,plates->NPDIM+1);

  A->trta = c_alloc2d(plates->NLAT+1,plates->NLONG+1);
  A->trpa = c_alloc2d(plates->NLAT+1,plates->NLONG+1);
  A->trra = c_alloc2d(plates->NLAT+1,plates->NLONG+1);

  printf("   Now for the big torque arrays: srpa * srta...\n");
  A->srpa = c_alloc3d(plates->NLAT+1,plates->NLONG+1,plates->NPDIM+1);
  A->srta = c_alloc3d(plates->NLAT+1,plates->NLONG+1,plates->NPDIM+1);

  A->idp = i_alloc2d(plates->NLAT+1,plates->NLONG+1);

  A->WW = c_alloc1d(plates->NPDIM+1);
  A->VV = c_alloc2d(plates->NPDIM+1,plates->NPDIM+1);
  

  A->pointx = c_alloc1d(plates->NPDIM+1);
  A->pointy = c_alloc1d(plates->NPDIM+1);
  A->pointz = c_alloc1d(plates->NPDIM+1);
  
  A->pointx2 = c_alloc1d(plates->NPDIM+1);
  A->pointy2 = c_alloc1d(plates->NPDIM+1);
  A->pointz2 = c_alloc1d(plates->NPDIM+1);
  
  A->vtheta = c_alloc2d(plates->NLAT+1,plates->NLONG+1);
  A->vphi = c_alloc2d(plates->NLAT+1,plates->NLONG+1);
  
  //printf(" alloc: NLAT = %i \n",plates->NLAT);

}
   	
	
void hcplates_command_line(int argc, char **argv,
			    struct hc_plates_params *p)
{
  int i;
  //printf("In here?\n");
  for(i=1;i < argc;i++){
   //printf("In here? 2\n");
    if(strcmp(argv[i],"-h")==0 || strcmp(argv[i],"-?")==0){// help
      /* 
	 help page
      */
      fprintf(stderr,"%s - HCPLATES  \n",
	      argv[0]);
      fprintf(stderr,"Options:\n\n");
      fprintf(stderr,"-L load file. Default point.j.\n");
      fprintf(stderr,"-U Unitrot file (output from ptrot). Default unitrot.coeffs. \n");
	   fprintf(stderr,"-T map plate file for torque integration. Default map.plate. \n");
	  fprintf(stderr,"-O output file (default: hcplates.output) containing ... stuff. \n");
	  fprintf(stderr,"-Op POLES output file (default: hcpoles.output) containing ... stuff. \n");
	  fprintf(stderr,"-Ov VEL GRID output file (default: hc_velgrid.output) containing ... stuff. \n");
		fprintf(stderr,"-P parameter file name to override defaults. Be a little careful here, file structure counts. \n");  
      fprintf(stderr,"\n\n");
      exit(-1);
	}
    else if(strcmp(argv[i],"-L")==0){ /* load filename */
	  // printf("In here? 3\n");
      hcplates_advance_argument(&i,argc,argv);
      strncpy(p->loadfile,argv[i],HC_CHAR_LENGTH);
	}
    else if(strcmp(argv[i],"-U")==0){ /* data filename */
	//   printf("In here? 4\n");
      hcplates_advance_argument(&i,argc,argv);
      strncpy(p->unitrotfile,argv[i],HC_CHAR_LENGTH);
    }
	else if(strcmp(argv[i],"-T")==0){ /* torque filename */
	//   printf("In here? 4\n");
      hcplates_advance_argument(&i,argc,argv);
      strncpy(p->platemapfile,argv[i],HC_CHAR_LENGTH);
    }
	else if(strcmp(argv[i],"-O")==0){ /* output filename */
	//   printf("In here? 4\n");
      hcplates_advance_argument(&i,argc,argv);
      strncpy(p->outputfile,argv[i],HC_CHAR_LENGTH);
    }
	else if(strcmp(argv[i],"-Op")==0){ /* output filename */
	//   printf("In here? 4\n");
      hcplates_advance_argument(&i,argc,argv);
      strncpy(p->polesfile,argv[i],HC_CHAR_LENGTH);
    }
	else if(strcmp(argv[i],"-Ov")==0){ /* output filename */
	//   printf("In here? 4\n");
      hcplates_advance_argument(&i,argc,argv);
      strncpy(p->velgridfile,argv[i],HC_CHAR_LENGTH);
    }
	else if(strcmp(argv[i],"-P")==0){ /* parameter filename */
	//   printf("In here? 4\n");
	  printf(" P!: NLAT = %i \n",p->NLAT);
      hcplates_advance_argument(&i,argc,argv);
      strncpy(p->parameterfile,argv[i],HC_CHAR_LENGTH);  /*prob string to integer? */
	  read_parameter_file(p);
	  printf(" P2: NLAT = %i \n",p->NLAT);
    }
	else{
	//   printf("In here? 5\n");
      fprintf(stderr,"%s: can not use parameter %s, use -h for help page\n",
	      argv[0],argv[i]);
      exit(-1);
    }
  }
}

/* 

//
// check, if we can read a value for the option flag in a chain of command line
// arguments
//  */
void hcplates_advance_argument(int *i,int argc, char **argv)
{
  if(argc <= *i + 1){// no arguments left
    fprintf(stderr,"%s: input parameters: error: option \"%s\" needs a value\n",
	    argv[0],argv[*i]);
    exit(-1);
  }
  *i += 1;
}

void read_parameter_file(struct hc_plates_params *p) {
	FILE *filePtr_P;
	char file_P[HC_CHAR_LENGTH];
	char toss1[3000];
	double toss2,toss3;
	int i,nlayer,m;
	HC_PREC r2[30],visc2[30],p_area[30];
	
	fprintf(stderr,"Reading parameter file %s\n",p->parameterfile);	
	
	strncpy(file_P,p->parameterfile,HC_CHAR_LENGTH);
	filePtr_P = fopen(file_P, "r");
   if (filePtr_P == NULL)
   { printf("Error loading parameter file.  Will use defaults. Ciao. \n"); }
   else {

	fscanf(filePtr_P,"%[^\n]%*1[\n]",toss1);	//Skip first 4 lines
	//fprintf(stderr,"%s \n",toss1);
	fscanf(filePtr_P,"%[^\n]%*1[\n]",toss1);
	//fprintf(stderr,"%s \n",toss1);
	fscanf(filePtr_P,"%[^\n]%*1[\n]",toss1);
	//fprintf(stderr,"%s \n",toss1);
	fscanf(filePtr_P,"%[^\n]%*1[\n]",toss1);
	//fprintf(stderr,"%s \n",toss1);
	fscanf(filePtr_P," %[^\n]%*1[\n]",toss1);  //This will read blank line and next text line as one

	// read no. distinct layers
	fscanf(filePtr_P,"%s %i",toss1,&nlayer);
	fprintf(stderr,"	nlayer = %i \n",nlayer);
	fscanf(filePtr_P," %[^\n]%*1[\n]",toss1);  //This will read blank line and next text line as one
	
	//Scan in viscosity for different layers
	fscanf(filePtr_P,"%s",toss1);
	for (i=0; i<=nlayer; i++) {
		fscanf(filePtr_P,"%lf",&r2[i]);
	}
	
	//fscanf(filePtr_P,"%s %lf %lf %lf %lf %lf %lf",&toss1,&p->r[0],&p->r[1],&p->r[2],&p->r[3],&p->r[4],&p->r[5]);
	for (i=0; i<=nlayer; i++) {
		r2[i] *=1e5;
	}

	//fprintf(stderr,"	%s  %lf, %lf, %lf, %lf %lf, %lf\n",toss1,p->r[0],p->r[1],p->r[2],p->r[3],p->r[4],p->r[5]);
	fscanf(filePtr_P," %[^\n]%*1[\n]",toss1);  //This will read blank line and next text line as one
	
	// Read viscosity
	fscanf(filePtr_P,"%s %lf",toss1,&p->visc0);
	fprintf(stderr,"	visc0 = %lf \n",p->visc0);
	fscanf(filePtr_P," %[^\n]%*1[\n]",toss1);  //This will read blank line and next text line as one
	
	// Viscosity layers
	fscanf(filePtr_P,"%s",toss1);
	for (i=0; i<nlayer; i++) {
			fscanf(filePtr_P,"%lf",&visc2[i]);
			fprintf(stderr,"		%i  %lf\n",i,visc2[i]);
	}

	fscanf(filePtr_P," %[^\n]%*1[\n]",toss1);  //This will read blank line and next text line as one
	// core radius
	fscanf(filePtr_P,"%s %lf",toss1,&p->RCore);
	fprintf(stderr,"	Rcore = %lf km \n",p->RCore);
	p->RCore *= 1e5;		//units not kms
	fscanf(filePtr_P," %[^\n]%*1[\n]",toss1);  //This will read blank line and next text line as one
	// earth radius
	fscanf(filePtr_P,"%s %lf",toss1,&p->erad);
	fprintf(stderr,"	Rearth = %lf km \n",p->erad);
	p->erad *= 1e5;			//units not kms
	fscanf(filePtr_P," %[^\n]%*1[\n]",toss1);  //This will read blank line and next text line as one
	//No .plates
	fscanf(filePtr_P,"%s %i",toss1,&p->NPLT);
	fprintf(stderr,"	NPLT = %i \n",p->NPLT);
	fscanf(filePtr_P," %[^\n]%*1[\n]",toss1);  //This will read blank line and next text line as one
	// Plate areas
	// Viscosity layers
	fscanf(filePtr_P,"%s",toss1);
	for (i=0; i<=p->NPLT; i++) {
			fscanf(filePtr_P,"%lf",&p_area[i]);
			fprintf(stderr,"		%i  %lf\n",i,p_area[i]);
	}
	fscanf(filePtr_P," %[^\n]%*1[\n]",toss1);  //This will read blank line and next text line as one
	//LDIM
	fscanf(filePtr_P,"%s %i",toss1,&p->LDIM);
	fprintf(stderr,"	LDIM = %i \n",p->LDIM);
	fscanf(filePtr_P," %[^\n]%*1[\n]",toss1);  //This will read blank line and next text line as one
	
	// Dependents
	p->NPDIM = 3 * p->NPLT;
	p->LDIM1 = p->LDIM + 1;
	p->KDIM = (p->LDIM + 2)*(p->LDIM+1)/2;

	// Geopgraphic divisor
	fscanf(filePtr_P,"%s %i",toss1,&p->dd);
	fprintf(stderr,"	dd = %i \n",p->dd);
	fscanf(filePtr_P," %[^\n]%*1[\n]",toss1);  //This will read blank line and next text line as one
	// time conversion factor

	fscanf(filePtr_P,"%s %lf",toss1,&p->stoy);
	fprintf(stderr,"	stoy = %lf \n",p->stoy);
	fscanf(filePtr_P," %[^\n]%*1[\n]",toss1);  //This will read blank line and next text line as one
	//Lmax

	fscanf(filePtr_P,"%s %i",toss1,&p->Lmax);
	fprintf(stderr,"	Lmax = %i \n",p->Lmax);
	fscanf(filePtr_P," %[^\n]%*1[\n]",toss1);  //This will read blank line and next text line as one
	//n
	fscanf(filePtr_P,"%s %i",toss1,&p->n);
	fprintf(stderr,"	n = %i \n",p->n);
	fscanf(filePtr_P," %[^\n]%*1[\n]",toss1);  //This will read blank line and next text line as one
	//no-slip
	fscanf(filePtr_P,"%s %i",toss1,&p->iba);
	fprintf(stderr,"	iba = %i \n",p->iba);
	// scaling ratio for velocities at end
	//no-slip
	fscanf(filePtr_P," %[^\n]%*1[\n]",toss1);  //This will read blank line and next text line as one
	fscanf(filePtr_P,"%s %lf",toss1,&p->ratio);
	fprintf(stderr,"	ratio = %lf \n",p->ratio);
	//DONE
	
	
	p->Lload = p->Lmax;
	p->Lplt = p->Lmax;

	m=p->n + 1;
	p->r[0] = 0.0;
	p->visc[0] = 0.0;
	for (i=1;i<m;i++) {
		p->r[i] = r2[i];
		p->visc[i]= visc2[i];
	}
	p->r[m] = r2[m];
	
	for (i=1; i<=p->NPLT; i++) {
		p->parea[i] = p_area[i];
	}
	
	printf(" NLAT = %i \n",p->NLAT);
	
	fclose(filePtr_P);
   }

}

