#include "hc_findplate.h"


/* 
   
   stand-alone routine to take plate boundary file, and metadata file with no. plates + no. points on each plate
   And then create a file with plate_id lat long for every lat/long on a 180x360 grid.
   
   For Lithgow-Bertelloni's & Richards Plates extensions to Hager & O'Connell experimental code
   O'Neill (coneill@els.qm.edu.au) & Becker (twb@usc.edu).

   This has routine to: set initial defaults
			handle the command line args (ie. just the input file name)
			read the input file, and assign parameters
  
*/

int main(int argc, char **argv)
{
struct hc_findplates p[1];
FILE *filePtr, *filePtr2, *filePtr3;
int nplate,ii,jj,npoint[100],tmp,k,l,ifind,i,ncross,j,mod1,test1,m,mm,id[182][362],done_long;
double x[60][3000],y[60][3000],dlat,dlong, dlat1,dlong1,xp,yp;
double yo,xo,xj,yj,xj1,yj1,dely,a1,b1,a2,b2,xint,yint,plat,plong;
int cross1,cross2,cross3,cross4,lll,idlast;

char file1[HC_CHAR_LENGTH], file2[HC_CHAR_LENGTH];
char file3[HC_CHAR_LENGTH];

/*	   Note: had to shrink x,y & id arrays in size to stop seg faulting
       Here: handle command line  parameters
   */
   //printf("Initializing defaults a \n");
   dlat = 1.0;
   dlong = 1.0;
   
   printf("Initializing defaults \n");
   hc_findplates_init(p);
   printf("Reading arguments\n");
   hc_findplates_command_line(argc, argv, p);
   fprintf(stderr,"Plate boundary file is %s\n",p->boundary_file);
   fprintf(stderr,"Plate data file is %s\n",p->data_filename);
   fprintf(stderr,"Output file is %s\n",p->output_file);
   
	 //Opening files
   strncpy(file3,p->output_file,HC_CHAR_LENGTH);
   filePtr3 = fopen(file3, "w");
   if (filePtr3 == NULL)
   { printf("Error loading output file. \n"); }
   else {		 
	 
	strncpy(file2,p->boundary_file,HC_CHAR_LENGTH);
	filePtr2 = fopen(file2, "r");
	if (filePtr2 == NULL)
	{ printf("Error loading boundary file. \n"); }
	else {	
	
	 strncpy(file1,p->data_filename,HC_CHAR_LENGTH);
	 filePtr = fopen(file1, "r");
	 if (filePtr == NULL)
	   { printf("Error loading Data file. \n"); }
     else
      {
	  /* Get nplate on first line */
	  fscanf(filePtr, "%i", &nplate);
	  printf("nplate: %i \n",nplate);
	
	  for (ii=1; ii<=nplate; ii++) {
		/* Read number of points on each plate */
		fscanf(filePtr, "%i", &npoint[ii]);
		printf("npoint: %i \n", npoint[ii]);
		
		for (jj=1; jj<=npoint[ii]; jj++) {
			fscanf(filePtr2, "%lf %lf", &x[ii][jj], &y[ii][jj]);
			//printf("		x, y: %lf %lf \n", x[ii][jj],y[ii][jj]);
			if (y[ii][jj] < 0.0) {
				y[ii][jj] += 360.00;
			}
		}
		
		tmp = npoint[ii]+1;
		x[ii][tmp]=x[ii][1];
		y[ii][tmp]=y[ii][1];
	  }		//finish nplates loop
	  
	  dlat1 = 180/dlat;
	  dlong1 = 360/dlong;
	  
	  for (k=1; k <= dlat1; k++) {
		xp = 0.0001+dlat*(k-0.5)-90.0;
		for (l=1; l<=dlong1; l++) {
			ifind = 0;
			done_long = 1;
			yp = 0.0001+dlong*(l-0.5);
			for (i=1; i<=nplate; i++) {
				ncross=0;
				for (j=1; j<=npoint[i]; j++) {
					yo=yp;
					xo=90;
					xj=x[i][j];
					yj=y[i][j];
					xj1=x[i][j+1];
					yj1=y[i][j+1];
					dely = fabs(yj-yj1);
					
					if (dely > 100) {
						if (yp < 180) {
							if (yj1 > 180) {
								yj1 -= 360;
							}
							if (yj > 180) {
								yj -= 360;
							}
						}
						else {
							if (yj1 < 180) {
								yj1 += 360;
							}
							if (yj < 180) {
								yj += 360;
							}
						}
					}  /*end ifs */
					
					/* Why the special cases for these plates? */
					if (xp > 0) {
						xo = -90.0;
					}
					if (i == 2) {
						xo = 90.0;
					}
					if (i == 7) {
						xo = -90.0;
					}
					
					a1 = (yp-yo)/(xp-xo);
					b1 = yo - a1*xo;
					
					test1=0;
					if (xj == xj1) {
						cross3 = (((yp>yj)&&(yp<yj1))||((yp>yj1)&&(yp<yj)));
						cross4 = (((xj>xo)&&(xj<xp))||((xj>xp)&&(xj<xo)));
						xint = xj;
						yint = yp;
						if (cross3 && cross4) {
							ncross++;
							test1=1;
						}
					}
					if (test1 != 1) {
						a2 = (yj-yj1)/(xj-xj1);
						b2 = yj - a2*xj;
						xint = (b2-b1)/(a1-a2);
						yint = a1*xint + b1;
						cross1 = (((xint>xj)&&(xint<xj1))||((xint>xj1)&&(xint<xj)));
						cross2 = (((xint>xo)&&(xint<xp))||((xint>xp)&&(xint<xo)));
						if (cross1 && cross2) {
							ncross++;
						}
					}
				}  /* npoint for loops end */
				mod1 = ncross%2;
				
				if (mod1 != 0 && id[k][l]==0 && ifind==0) {  /*This is a big fudge to stop doubling up at same points */
					id[k][l]=i;
					ifind=1;
					done_long=0;
				}
				else {
					id[k][l]=0;
				}
				if (id[k][l] > 0 && done_long<1){
					plat = 0.0001 + dlat*(k-0.5) - 90.0;
					plong = 0.0001 + dlong*(l-0.5);
					fprintf(filePtr3,"%i %lf %lf\n",id[k][l],plat,plong);	//output here
					done_long = 1;
					idlast = id[k][l];
				}
				
			}  /*nplate loop ends */
			if (ifind==0) {
						id[k][l]=idlast;
						fprintf(filePtr3,"%i %lf %lf \n",id[k][l],plat,plong);					
			}
		}	// long loop
			
	   }			/*Last for loop (dlat1) ends */
	   
	  } // files open loops
	}
	}	//output file
	 
	 fclose(filePtr3);
     fclose(filePtr);
	 fclose(filePtr2);

 printf("Exiting ok \n");
	 return 0;
}


void hc_findplates_init(struct hc_findplates *p) {
	strncpy(p->boundary_file,"data",HC_CHAR_LENGTH);
	strncpy(p->data_filename,"enes",HC_CHAR_LENGTH);
	strncpy(p->output_file,"plates_ids.ixz",HC_CHAR_LENGTH);
}

void hc_findplates_command_line(int argc, char **argv,
			    struct hc_findplates *p)
{
  int i;
  //printf("In here?\n");
  for(i=1;i < argc;i++){
   //printf("In here? 2\n");
    if(strcmp(argv[i],"-h")==0 || strcmp(argv[i],"-?")==0){// help
      /* 
	 help page
      */
      fprintf(stderr,"%s - Routine to create lat/lon/plate_id dile\n\n",
	      argv[0]);
      fprintf(stderr,"options:\n\n");
      fprintf(stderr,"-B plate boundary file (%s)\n",
	      p->boundary_file);
      fprintf(stderr,"-D data file (enes) containg no. plates (first line), then no. points on each plate (rest). Both files needed. \n");
	   fprintf(stderr,"-O output file (default: plate_ids.ixz) containg plates id, lat & long\n");	  
      fprintf(stderr,"\n\n");
      exit(-1);
	}
    else if(strcmp(argv[i],"-B")==0){ /* boundary filename */
	  // printf("In here? 3\n");
      hc_findplates_advance_argument(&i,argc,argv);
      strncpy(p->boundary_file,argv[i],HC_CHAR_LENGTH);
	}
    else if(strcmp(argv[i],"-D")==0){ /* data filename */
	//   printf("In here? 4\n");
      hc_findplates_advance_argument(&i,argc,argv);
      strncpy(p->data_filename,argv[i],HC_CHAR_LENGTH);
    }
	else if(strcmp(argv[i],"-O")==0){ /* output filename */
	//   printf("In here? 4\n");
      hc_findplates_advance_argument(&i,argc,argv);
      strncpy(p->output_file,argv[i],HC_CHAR_LENGTH);
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
void hc_findplates_advance_argument(int *i,int argc, char **argv)
{
  if(argc <= *i + 1){// no arguments left
    fprintf(stderr,"%s: input parameters: error: option \"%s\" needs a value\n",
	    argv[0],argv[*i]);
    exit(-1);
  }
  *i += 1;
}

