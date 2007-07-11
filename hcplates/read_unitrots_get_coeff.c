#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "hcplates.h"

/* This subroutine reads in a file with all the unit rotation velocities for each plate
add up the products of the unit rotation velocities and the shear stress
c**********     Green functions (poloidal and toroidal) for each plate and each
c**********     orthogonal rotation direction

unitrots_filename= name of file containing unitrots SH coefficients.
n=no. layers
r = radii
visc = viscosity structure

L= degree
m = SH index
ua, uc poloidal green's functions
va vc toroidal green's functions

Needs polvel & torvel. And a file with unit rotation coefficients.

TESTED AGAINST FORTRAN AND SEEMS OK

*/

void read_unitrots_get_coeff(struct hc_plates_params *plates, struct hc_plates_arrays *A, char unitrots_filename[])
{
  FILE *filePtr;
  int Lplt,count;
  double gpol[plates->LDIM1],gtor[plates->LDIM1],gfact;
  double ua[5], uc[5], va[3], vc[3];
  int Ldum, mdum;
  double cvel,svel;
  int L,iplt,idir,ii,k,m;

 /* read file unitrot.coeffs - firstline is Lplt=max L. Then format is Ldum, mdum, cvel, svel */
  fprintf(stderr,"Unitrots filename=%s\n",unitrots_filename);
  filePtr = fopen(unitrots_filename, "r");
  if (filePtr == NULL)
	{ printf("Error opening file. \n"); }
  else
	{
	/* Get Lplt on first line */
	fscanf(filePtr, "%d", &Lplt);
	plates->Lplt = Lplt;
	printf("  Lplt: %d %d\n",plates->Lplt,Lplt);
	/* Before we go any further, get toroidal and poloidal Green functions
	  up to degree Lplt */
	gpol[1]=0.0;
	gtor[1]=0.0;
	for (L=1; L<=Lplt; L++) {
		polvel(plates,L,plates->n,ua,uc);
		torvel(plates,L,plates->n,va,vc);
		gpol[L+1]=ua[4];
		gtor[L+1]=va[2];
		//printf("		setting gfact: L=%d gpol=%lf gtor=%lf ua4=%lf va4=%lf\n",L,gpol[L+1],gtor[L+1],ua[4],va[2]); //good till here
	}	
	/* c*** loop on plate number, rotation direction: accumulate shear poloidal and
c***     toroidal shear stress coefficients ('y4ex' and 'y10ex')
c*** note that a factor of visc0/erad has been left out of 'gfact', i.e., 'y4ex' and 'y10ex',
c***   since it ultimately multiplies both sides of the final matrix equation
	*/
    count = 0;
	for (iplt=1; iplt <= plates->NPLT; iplt++) {
		for (idir=1; idir<=3; idir++) {
			ii=3*(iplt-1)+idir;
			/*Poloidal */
			gfact=0.0;
			for (L=0; L<=Lplt; L++) {
				gfact=gpol[L+1];
				for (m=0; m<=L; m++) {
					fscanf(filePtr, " %d %d %lf %lf", &Ldum, &mdum, &cvel, &svel); /*check this input */
					k=(L+2)*(L+1)/2-m; 	/*ambiguous*/
					/*printf("P: iplt %d L %d m %d k %d __ Ldum %d mdum %d __ cvel %lf svel %lf\n",iplt,L,m,k,Ldum,mdum,cvel,svel); */

					A->y4ex[ii][k][1]=cvel*gfact;
					A->y4ex[ii][k][2]=svel*gfact;
					//printf("%d  %d/%d/%dP/%d_%d/%d cvel=%lf svel=%lf gfact=%lf y4ex1=%lf y4ex2=%lf\n",count,iplt,idir,L,m,ii,k,cvel,svel,gfact,A->y4ex[ii][k][1],A->y4ex[ii][k][2]); 
					count++;
				}
			}
			/*Toroidal*/
			for (L=0; L<=Lplt; L++) {
				gfact=gtor[L+1];
				for (m=0; m<=L; m++) {
					fscanf(filePtr, "%d %d %lf %lf",&Ldum, &mdum, &cvel, &svel); /* check */
					/*printf("T: iplt %d L %d m 	%d Ldum %d mdum %d 	cvel %lf svel %lf\n",iplt,L,m,Ldum,mdum,cvel,svel); */
					k=(L+2)*(L+1)/2-m; 	/*ambiguous*/
					A->y10ex[ii][k][1]=cvel*gfact;
					A->y10ex[ii][k][2]=svel*gfact;
				//	printf("%d  %d/%d/%dT/%d_%d/%d y10ex[1,2]= [%le] [%le] cvel[%le] sv[%le] gfact[%le]\n",count,iplt,idir,L,m,ii,k,A->y10ex[ii][k][1],A->y10ex[ii][k][2],cvel,svel,gfact);
					count++;
				}
			}
		}
	}
  }
  printf("\n");
  printf("   Line count is %d \n",count);
  printf("\n");
  fclose(filePtr);

}
		

 

