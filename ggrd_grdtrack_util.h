/* 

header files for modified GMT grd interpolation routines dealing
with grd interpolation

*/

#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "gmt.h"

#include "ggrd.h"
/* 

wrappers

*/
int ggrd_grdtrack_init_general(unsigned char ,char *, char *,char *,
			       struct ggrd_gt *,unsigned char ,
			       unsigned char);
unsigned char ggrd_grdtrack_interpolate_rtp(double ,double ,double ,
					    struct ggrd_gt *,double *,
					    unsigned char);
unsigned char ggrd_grdtrack_interpolate_xyz(double ,double ,double ,
					    struct ggrd_gt *,double *,
					    unsigned char);
unsigned char ggrd_grdtrack_interpolate_xy(double ,double ,
					   struct ggrd_gt *,
					   double *,
					   unsigned char );
unsigned char ggrd_grdtrack_interpolate_tp(double ,double ,
					   struct ggrd_gt *,
					   double *,
					   unsigned char );

void ggrd_grdtrack_free_gstruc(struct ggrd_gt *);

int ggrd_grdtrack_rescale(struct ggrd_gt *,unsigned char , unsigned char , 
			  unsigned char ,double);


/* 

moderately external

*/
int ggrd_read_time_intervals(struct ggrd_t *,char *,unsigned char ,unsigned char);
int ggrd_read_vel_grids(struct ggrd_vel *, double, unsigned short, unsigned short, char *);
#ifdef USE_GMT4
unsigned char ggrd_grdtrack_interpolate(double *, unsigned char , struct GRD_HEADER *, float *,
					struct GMT_EDGEINFO *, int, float *, int ,	double *,unsigned char,
					struct GMT_BCR *);
int ggrd_grdtrack_init(double *, double *, double *, double *, float **, int *, char *, struct GRD_HEADER **, struct GMT_EDGEINFO **, char *, unsigned char *, int *, unsigned char, char *, float **, int *, unsigned char, unsigned char, unsigned char, struct GMT_BCR *);

#else
unsigned char ggrd_grdtrack_interpolate(double *, unsigned char , struct GRD_HEADER *, float *,
				  struct GMT_EDGEINFO *, int, float *, int ,	double *,unsigned char);

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
			char *,unsigned char *,
			int *,	/* [4] array with padding (output) */
			unsigned char _d, char *, 	/* depth file name */
			float **,	/* layers, pass as NULL */
			int *,		/* number of layers */
			unsigned char , /* linear/cubic? */
			unsigned char ,unsigned char);
#endif
/* 

local 

 */

void ggrd_gt_bcr_init_loc(void);
void ggrd_gt_interpolate_z(double,float *,int ,
			   int *, int *, double *, double *,unsigned char); /*  */
float ggrd_gt_rms(float *,int );
float ggrd_gt_mean(float *,int );

void ggrd_print_layer_avg(float *,float *,int , int ,FILE *);

void ggrd_find_spherical_vel_from_rigid_cart_rot(double *,
						 double *,
						 double *,
						 double *, 
						 double *);
						 

/* hc related utility functions */

void hc_vecalloc(double **,int,char *);
void hc_vecrealloc(double **,int,char *);
