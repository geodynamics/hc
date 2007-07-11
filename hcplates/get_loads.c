#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "hcplates.h"

/* routine to read a load file and add up the products of the internal loads and the poloidal shear
c**********     stress Green functions for all load radii
c**********
c
c*** read first line of internal load file: number of load radii, max L
      read(22,*) nrad,Lload
c*** loop on load radius, harmonic degree and order: read load radius, read load 
c***     coefficients, calculate Green functions, and accumulate poloidal shear
c***     stress coefficients ('y4in')
c*** note that a factor of visc0/erad has been left out of 'y4in',
c***     since it ultimately multiplies both sides of the final matrix equation
c*** tomography conversion factor=tconv*****

Makes  y4in, y2in, y3in, z3in

Inputs dens_file: load input file (points.j)
	iba:  flag for poload, 0=no slip, 1=free slip

Make structure for plates earth model:
	tconv (tomography converions factor)
	n: no. radial layers
	erad, Earth radius 
	r: radial structure
	visc: viscosity layers
	visc0: ref visc ~1e22

TESTED AGAINST CRLB'S FORTRAN CODE AND WORKING -- SEE COMMENTS FOR POSSIBLE TYPOS IN FORTRAN

 */

void get_loads(struct hc_plates_params *plates, struct hc_plates_arrays *A, char dens_file[])
{
   FILE *filePtr;
   int nrad,Lload,irad,L,m,Ldum,mdum,k;
   double rload,rload1,cload,sload,cload1,sload1,tconv;
   double ua[5],uc[5];

/* Open load file */
   printf("Get_loads: density file is %s \n",dens_file);
   filePtr = fopen(dens_file, "r");
  if (filePtr == NULL)
	{ printf("Error load opening file. \n"); }
  else
	{
	/* Get nrad & Lload on first line */
       	fscanf(filePtr, "%i %i", &nrad, &Lload);
	printf("Get nrad * Lload: %i %i \n",nrad,Lload);
	plates->Lload = Lload;
	fscanf(filePtr, "%lf", &tconv);
	printf("Get tomography conversion factor: %lf \n",tconv);
	for (irad=1; irad<=nrad; irad++) {
		fscanf(filePtr, "%lf",&rload1); /* check this input */
		/*printf("Rload1: %lf 	",rload1); */
		rload = plates->erad - rload1*1.0e5;
		/* printf("Rload: %lf \n",rload); */
		for (L=0; L<=Lload; L++) {
			if (L==0) {
				ua[2]=0.0;
				ua[3]=0.0;
				ua[4]=0.0; /* uc not initialized - why?? */
			}
			else {
				poload(plates,L,plates->n,plates->iba,1,rload,plates->visc0,ua,uc);
				//fprintf(stderr,"LOADS: %i/%i ua2:%le ua3:%le ua4:%le\n",L,m,ua[2],ua[3],ua[4]);
				//fprintf(stderr," LOADS: %i/%i uc2:%le uc3:%le uc4:%le\n",L,m,uc[2],uc[3],uc[4]);
				/*printf("%i ",L); */
			}
			for (m=0; m <= L; m++) {
				fscanf(filePtr, " %d %d %lf %lf",&Ldum,&mdum,&cload1,&sload1); /* check this input */
				if (irad > 0) {
					cload=cload1*98.10*tconv;
					sload=sload1*98.10*tconv;
				}
				else {
					cload=cload1*98.10;
					sload=sload1*98.10;
				}
				k=(L+2)*(L+1)/2-m;
				/* In CRLBs code, y4in & y3in are intially set to zero
				   Additionally, ua[2-4] are zero'ed above.
				   However, y2in and z3in are not intially zero'd here - is this a problem?
				   Also, uc[3] is not initially set to zero if L=0 - how come? */
				A->y4in[k][1] += cload*ua[4];
				A->y4in[k][2] += sload*ua[4];
				A->y2in[k][1] += cload*ua[2];
				A->y2in[k][2] += sload*ua[2];
				A->y3in[k][1] += cload*ua[3];
				A->y3in[k][2] += sload*ua[3];
				A->z3in[k][1] += cload*uc[3];
				A->z3in[k][2] += cload*uc[3];		//CRLB's code has this as cload not sload - is this right? 
													// k,1 & k,2 are now the same  - typo??
			
				//printf("irad=%d L=%d m=%d y2in:%lf	%lf \n",irad,L,m,y2in[k][1],y2in[k][2]);
			//	printf("LOADS: irad=%d L=%d m=%d z3in:%le	%le \n",irad,L,m,A->z3in[k][1],A->z3in[k][2]);
				
				if (irad == nrad) {
					/* Ouput routines here - intstress and intvel */
					
				}
			}
		}
	}
  }
  fclose(filePtr);

}


