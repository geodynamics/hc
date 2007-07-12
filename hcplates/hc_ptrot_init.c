#include "hc_ptrot.h"

/* Initialization routine */

void hc_ptrot_init (struct hc_ptrot *A) 
{
	int ip,l,m;

	/* Default filenames */
	strncpy(A->boundary_file,"data",HC_CHAR_LENGTH);
	strncpy(A->data_filename,"enes",HC_CHAR_LENGTH);
	strncpy(A->plateids_file,"plates_ids.ixz",HC_CHAR_LENGTH);
	strncpy(A->unitrot_file,"unitrot.coeffs",HC_CHAR_LENGTH);
	
	A->nlm = 20;
	A->ylm1 = 0;
	A->ylm2 = 0;
	A->plm = 0;
	
	A->tvelxm = 0;
	A->tvelym = 0;
	A->tvelzm = 0;
	A->tvelxp = 0;
	A->tvelyp = 0;
	A->tvelzp = 0;

	/* Allocate arrays */
	A->dlmx = c_alloc4d(16,3,50,50);
	A->dlmy = c_alloc4d(16,3,50,50);
	A->dlmz = c_alloc4d(16,3,50,50);
	
	A->vlmx = c_alloc4d(17,4,50,50);
	A->vlmy = c_alloc4d(16,3,50,50);
	A->vlmz = c_alloc4d(16,3,50,50);
	
	A->dlm = c_alloc5d(16,4,3,50,50);
	A->vlm = c_alloc5d(16,4,3,50,50);
	
	A->num = c_alloc1d(16);
	
	for (ip = 1; ip<=15; ip++) {
		for (l=0; l<=49; l++) {
			for (m=0; m<=l; m++) {
				
				
				A->dlmx[ip][1][l][m] = 0.0;
				A->dlmx[ip][2][l][m] = 0.0;
				
				A->dlmy[ip][1][l][m] = 0.0;
				A->dlmy[ip][2][l][m] = 0.0;
				
				A->dlmz[ip][1][l][m] = 0.0;
				A->dlmz[ip][2][l][m] = 0.0;
				
				A->vlmx[ip][1][l][m] = 0.0;
				A->vlmx[ip][2][l][m] = 0.0;
				
				A->vlmy[ip][1][l][m] = 0.0;
				A->vlmy[ip][2][l][m] = 0.0;
				
				A->vlmz[ip][1][l][m] = 0.0;
				A->vlmz[ip][2][l][m] = 0.0;
				
				A->dlm[ip][1][1][l][m] = 0.0;
				A->dlm[ip][1][2][l][m] = 0.0;
				A->dlm[ip][2][1][l][m] = 0.0;
				A->dlm[ip][2][2][l][m] = 0.0;
				A->dlm[ip][3][1][l][m] = 0.0;
				A->dlm[ip][3][2][l][m] = 0.0;
				
				A->vlm[ip][1][1][l][m] = 0.0;
				A->vlm[ip][1][2][l][m] = 0.0;
				A->vlm[ip][2][1][l][m] = 0.0;
				A->vlm[ip][2][2][l][m] = 0.0;
				A->vlm[ip][3][1][l][m] = 0.0;
				A->vlm[ip][3][2][l][m] = 0.0;
			}
		}
	}
				
}

void hc_ptrot_command_line(int argc, char **argv,
			    struct hc_ptrot *p)
{
  int i;
  //printf("In here?\n");
  for(i=1;i < argc;i++){
   //printf("In here? 2\n");
    if(strcmp(argv[i],"-h")==0 || strcmp(argv[i],"-?")==0){// help
      /* 
	 help page
      */
      fprintf(stderr,"%s - Routine to create unitrot.coeffs files\n",
	      argv[0]);
      fprintf(stderr,"Options:\n\n");
      fprintf(stderr,"-B plate boundary file (%s)\n",
	      p->boundary_file);
      fprintf(stderr,"-D data file (enes) containg no. plates (first line), then no. points on each plate (rest). Both files needed.\n");
	   fprintf(stderr,"-P plate id file (default: plate_ids.ixz) containg plates id, lat & long. \n");
	  fprintf(stderr,"-O output file (default: unitrots.coeff) containg coefficients for unit rotations on each plate \n");
	  fprintf(stderr,"-d degree of spherical harmonics (defaults to 20) \n");  
      fprintf(stderr,"\n\n");
      exit(-1);
	}
    else if(strcmp(argv[i],"-B")==0){ /* boundary filename */
	  // printf("In here? 3\n");
      hc_ptrot_advance_argument(&i,argc,argv);
      strncpy(p->boundary_file,argv[i],HC_CHAR_LENGTH);
	}
    else if(strcmp(argv[i],"-D")==0){ /* data filename */
	//   printf("In here? 4\n");
      hc_ptrot_advance_argument(&i,argc,argv);
      strncpy(p->data_filename,argv[i],HC_CHAR_LENGTH);
    }
	else if(strcmp(argv[i],"-P")==0){ /* output filename */
	//   printf("In here? 4\n");
      hc_ptrot_advance_argument(&i,argc,argv);
      strncpy(p->plateids_file,argv[i],HC_CHAR_LENGTH);
    }
	else if(strcmp(argv[i],"-O")==0){ /* output filename */
	//   printf("In here? 4\n");
      hc_ptrot_advance_argument(&i,argc,argv);
      strncpy(p->unitrot_file,argv[i],HC_CHAR_LENGTH);
    }
	else if(strcmp(argv[i],"-d")==0){ /* output filename */
	//   printf("In here? 4\n");
      hc_ptrot_advance_argument(&i,argc,argv);
      strncpy((void *)p->nlm,argv[i],HC_CHAR_LENGTH);  /*prob string to integer? */
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
void hc_ptrot_advance_argument(int *i,int argc, char **argv)
{
  if(argc <= *i + 1){// no arguments left
    fprintf(stderr,"%s: input parameters: error: option \"%s\" needs a value\n",
	    argv[0],argv[*i]);
    exit(-1);
  }
  *i += 1;
}
