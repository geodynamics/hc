/* 

header files for modified GMT grd interpolation routines dealing
with grd interpolation

*/


#ifndef __GGRD_READ_GMT__
#include "gmt.h"
void GMT_grdio_init (void);
#define __GGRD_READ_GMT__
#endif


#include "ggrd_base.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef USE_GMT3		/* new GMT */
#define GGRD_GMT_GLOBAL_STRING "-fg"
#define GGRD_GMT_GEOGRAPHIC_STRING "-fg"
#define GGRD_GMT_XPERIODIC_STRING "-Lx"
#else  /* old GMT */
#define GGRD_GMT_GLOBAL_STRING "-Lg"
#define GGRD_GMT_GEOGRAPHIC_STRING ""
#define GGRD_GMT_XPERIODIC_STRING ""
#endif

/* 

wrappers

*/
int ggrd_grdtrack_init_general(ggrd_boolean ,char *, char *,char *,
			       struct ggrd_gt *,ggrd_boolean ,
			       ggrd_boolean, ggrd_boolean);
ggrd_boolean ggrd_grdtrack_interpolate_lonlatz(double ,double ,double ,
					       struct ggrd_gt *,double *,
					       ggrd_boolean );
ggrd_boolean ggrd_grdtrack_interpolate_rtp(double ,double ,double ,
					    struct ggrd_gt *,double *,
					   ggrd_boolean,ggrd_boolean,double);
ggrd_boolean ggrd_grdtrack_interpolate_xyz(double ,double ,double ,
					    struct ggrd_gt *,double *,
					    ggrd_boolean);
ggrd_boolean ggrd_grdtrack_interpolate_xy(double ,double ,
					   struct ggrd_gt *,
					   double *,
					   ggrd_boolean );
ggrd_boolean ggrd_grdtrack_interpolate_tp(double ,double ,
					   struct ggrd_gt *,
					   double *,
					   ggrd_boolean ,ggrd_boolean);

void ggrd_grdtrack_free_gstruc(struct ggrd_gt *);

int ggrd_grdtrack_rescale(struct ggrd_gt *,ggrd_boolean , ggrd_boolean , 
			  ggrd_boolean ,double);

void ggrd_init_vstruc(struct ggrd_master *);

/* 

moderately external

*/
void ggrd_init_master(struct ggrd_master *);

int ggrd_init_thist_from_file(struct ggrd_t *,char *,ggrd_boolean ,ggrd_boolean);
int ggrd_read_vel_grids(struct ggrd_master *, double, unsigned short, unsigned short, char *,ggrd_boolean);

#ifndef USE_GMT3
/* GMT >4.1.2 */
ggrd_boolean ggrd_grdtrack_interpolate(double *, ggrd_boolean , struct GRD_HEADER *, float *,
					struct GMT_EDGEINFO *, int, float *, int ,	double *,ggrd_boolean,
					struct GMT_BCR *);
// GMT < 4.5.1
//int ggrd_grdtrack_init(double *, double *, double *, double *, float **, int *, char *, struct GRD_HEADER **, struct GMT_EDGEINFO **, char *, ggrd_boolean *, int *, ggrd_boolean, char *, float **, int *, ggrd_boolean, ggrd_boolean, ggrd_boolean, struct GMT_BCR *);
// GMT >= 4.5.1
int ggrd_grdtrack_init(double *, double *, double *, double *, float **, int *, char *, 
		       struct GRD_HEADER **, struct GMT_EDGEINFO **, char *, 
		       ggrd_boolean *, GMT_LONG *, ggrd_boolean, char *, 
		       float **, int *,GMT_LONG, ggrd_boolean, 
		       ggrd_boolean, struct GMT_BCR *);

/* < 4.5.1 */
//void ggrd_print_layer_avg(float *,float *,int , int ,int, FILE *,int *);
/* >= 4.5.1 */
void ggrd_print_layer_avg(float *,float *,int , int ,int, FILE *,GMT_LONG *);

#else
void ggrd_print_layer_avg(float *,float *,int , int ,int, FILE *,int *);

ggrd_boolean ggrd_grdtrack_interpolate(double *, ggrd_boolean , struct GRD_HEADER *, float *,
				       struct GMT_EDGEINFO *, int, 
				       float *, int ,	
				       double *,ggrd_boolean,
				       struct BCR *);

int ggrd_grdtrack_init(double *, double *,double *, double *, /* geographic bounds,
								 set all to zero to 
								 get the whole range from the
								 input grid files
							      */
			float **,	/* data */
			int *,  /* size of data */
			char *,	/* name, or prefix, of grd file with scalars */
			struct GRD_HEADER **,
			struct GMT_EDGEINFO **,
			char *,ggrd_boolean *,
			int *,	/* [4] array with padding (output) */
			ggrd_boolean _d, char *, 	/* depth file name */
			float **,	/* layers, pass as NULL */
			int *,		/* number of layers */
			ggrd_boolean , /* linear/cubic? */
		       ggrd_boolean ,ggrd_boolean,
		       struct BCR *);

void ggrd_global_bcr_assign(struct BCR *);

void my_GMT_bcr_init (struct GRD_HEADER *, int *, 
		      int ,struct BCR *);

#endif
/* 

local 

 */

void ggrd_gt_interpolate_z(double,float *,int ,
			   int *, int *, double *, double *,ggrd_boolean,
			   ggrd_boolean *); /*  */
float ggrd_gt_rms(float *,int );
float ggrd_gt_mean(float *,int );
FILE *ggrd_open(char *, char *, char *);


void ggrd_find_spherical_vel_from_rigid_cart_rot(double *,
						 double *,
						 double *,
						 double *, 
						 double *);
						 


void ggrd_vecalloc(double **,int,char *);
void ggrd_vecrealloc(double **,int,char *);
void ggrd_indexx(int ,GGRD_CPREC *, int *);
void ggrd_calc_mean_and_stddev(GGRD_CPREC *, GGRD_CPREC *,int ,GGRD_CPREC *,
			       GGRD_CPREC *,GGRD_CPREC *, 
			       ggrd_boolean , ggrd_boolean,GGRD_CPREC *);

