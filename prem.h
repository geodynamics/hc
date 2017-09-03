/* 

header file for PREM utilities

$Id: prem.h,v 1.4 2006/04/18 01:08:39 twb Exp twb $


*/

#ifndef __READ_PREM_HEADER__
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PREM_F_STRING "%lf"

/* default number of layers */
#define PREM_N 13
/* number of coefficients  */
#define PREM_NP 4		


#ifndef hc_boolean
#define hc_boolean unsigned short
#endif

#ifndef TRUE 
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

/* 
   PREM model structure 
*/
struct prem_model{
  double crho[PREM_N*PREM_NP];	/* density parameters */
  double cvp[PREM_N*PREM_NP];		/* v_p  */
  double cvs[PREM_N*PREM_NP]; 	/* v_s */
  double cvpv[PREM_N*PREM_NP];	/* anisotropic velocities */
  double cvph[PREM_N*PREM_NP];
  double cvsv[PREM_N*PREM_NP];
  double cvsh[PREM_N*PREM_NP];
  double ceta[PREM_N*PREM_NP];	/* anisotropy parameter */
  double cqmu[PREM_N];		/* q factors */
  double cqkappa[PREM_N];
  
  int n;			/* number of layers */
  int np;			/* number of polynomial coefficients
				   in general (except qmu and qkappa) 
				*/
  double rb[PREM_N];		/* top boundary of layers, in meters */
  double r[PREM_N];		/* non-dimensionalized version  */
  double r0;			/* surface r in meters */
  
  char model_filename[2000];	/* model file name */

  hc_boolean init;
};
/* 
   constants:
   
   radius of earth in km 

*/
#define PREM_RE_KM 6371.0087714

/* 

functions

*/
int prem_find_layer_x(double , double, double *,int ,int, double *);
double prem_compute_pval(double  *, double *, int , double );
double prem_compute_dpval(double *, double *, int , double );
double prem_vs_voigt(double , double , double , double , 
		    double);
void prem_get_rhodrho(double *, double *,  double , struct prem_model *);
void prem_get_rho(double *,double , struct prem_model *);

void prem_get_values(double *, double *, double *, double *, 
		     double *, double *, double *, double *, 
		     double *, double *, double , 
		     struct prem_model *);



int prem_read_model(char *,struct prem_model *, 
		     hc_boolean );

int prem_read_para_set(double *, int , int ,FILE *);
#define __READ_PREM_HEADER__

#endif
