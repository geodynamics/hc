#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "hcplates.h"


/* 
   These are just a few simple routines to zero all the major arrays we will encounter
   Note: all these arrays are being initialized with plates->KDIM not MAX_KDIM
   First do shear stress coefficient arrays
     
*/

void zero_arrays(struct hc_plates_params *plates,struct hc_plates_arrays *A)
{
  int i,j,k,KDIM,NPDIM,LDIM1;

  KDIM = plates->KDIM;
  NPDIM = plates->NPDIM;
  LDIM1 = plates->LDIM1;

  printf("SS arrays:  1\n");
  for (i=1;i<=KDIM;i++) {
	for (j=1;j<3;j++) {
		A->y4in[i][j] = 0.0;
		A->y3in[i][j] = 0.0;
		A->y2in[i][j] = 0.0;
		A->z3in[i][j] = 0.0;
	}
   }
  printf("            2\n");
  for (i=1;i<=NPDIM;i++) {
	for (j=1;j<=KDIM;j++) {
		for (k=1;k<3;k++) {
			A->y4ex[i][j][k] = 0.0;
			A->y10ex[i][j][k] = 0.0;
		}
	}
  }

 
   printf("SH arrays:  1\n");

  for (i=1;i<=KDIM;i++) {
	A->p[i] = 0.0;
	A->dpdt[i] = 0.0;	
	A->pbyst[i] = 0.0;
	for (j=1;j<plates->NLAT;j++) {
		A->torp[j][i] = 0.0;
	}
   }
   printf("            2\n");

  for (i=1; i<=LDIM1; i++) {
	A->cm[i] = 0.0;
	A->sm[i] = 0.0;
  }

  printf("Torque arrays:   1\n");
  for (i=1;i <= plates->NPDIM; i++) {
	A->srt[i] = 0.0;
	A->srp[i] = 0.0;
	A->fin[i] = 0.0;
	A->finet[i] = 0.0;
	A->fins[i] = 0.0;
	A->rots[i] = 0.0;
	A->ww[i] = 0.0;
	for (j=1; j <= plates->NPDIM; j++)  {
		A->fex[i][j] = 0.0;
		A->fexs[i][j] = 0.0;
		A->vv[i][j] = 0.0;
	/*	printf("		Zeroing: i=%d j=%d fexs=%lf\n",i,j,A->fexs[i][j]);*/
	}
  }
  printf("                 2\n");
  for (i=1;i <= plates->NLAT; i++) {
	for (j=1; j <= plates->NLONG; j++)  {
		A->trta[i][j] = 0.0;
		A->trpa[i][j] = 0.0;
		A->trra[i][j] = 0.0;
		for (k=1; k<= plates->NPDIM;k++) {
			A->srpa[i][j][k] = 0.0;
			A->srta[i][j][k] = 0.0;
		}
	}
  }
}
