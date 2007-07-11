#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "hc.h"

#define HC_MAX_LAYERS 100 /* Hope this is enough radial layers */
#define MAX_NO_PLATE 30
#define MAX_LDIM 100
#define MAX_NPDIM 90
#define MAX_KDIM 515	/* Note: having this too high gives seg faults But we need it... weird.*/

/* 
   
   header file for Lithgow-Bertelloni's & Richards Plates extensions to Hager & O'Connell experimental code
   O'Neill (coneill@els.qm.edu.au) & Becker (twb@usc.edu).
  
*/
/*void hc_init_params(struct hc_plates_params *);*/

struct hc_plates_params{
  HC_PREC RCore; /* Radius of core - following format in fotran (why e5) */
  HC_PREC visc0;   /* Normalizing viscosity */
  HC_PREC erad;

  /* Arrays, (harmonic coefficients, number of plate components etc)
     Note: LDIM should be the largest harmonic degree used in ANY input file */ 
  int NPLT;
  int LDIM;
  int NPDIM;
  int LDIM1;
  int KDIM;

  /*Grid discretization - changing these to integers for loops within routine (CO'N). */
  int dd;
  int NLAT;
  int NLONG;

  /* Conversion factors */
  HC_PREC stoy;  /* seconds to years */

  /* Earth information - to be imported */
  int Lmax; 	/*Maximum order for the final torque integration for both internal and external forces*/
		/*Should not be larger than Lplt OR Lload, but may be smaller */
  int n;	/* Number of radial layers */
  int iba;	/* Note - already use free-slip for this in hc, but not importing that structure here */
  HC_PREC r[HC_MAX_LAYERS+1];   /* Radial layer radii */
  HC_PREC visc[HC_MAX_LAYERS];  /* Viscosity layers */

  int Lload;	/* Max L of load file - will init to Lmax, should be reset when load file read */
  int Lplt;	/* Max order of SH expansion of unitrot velocities */
  
   /* File names */
   char loadfile[HC_CHAR_LENGTH];
   char unitrotfile[HC_CHAR_LENGTH];
   char platemapfile[HC_CHAR_LENGTH];
   char outputfile[HC_CHAR_LENGTH];
   char parameterfile[HC_CHAR_LENGTH];
   char polesfile[HC_CHAR_LENGTH];
   char velgridfile[HC_CHAR_LENGTH];
   
   double ratio;
   HC_PREC parea[MAX_NO_PLATE];

  
};

struct hc_plates_arrays{
   /* Stress arrays */
   HC_PREC **y4in;
   HC_PREC **y2in;
   HC_PREC **y3in;
   HC_PREC **z3in;
   HC_PREC ***y4ex;
   HC_PREC ***y10ex;

   /* Spherical harmonic arrays */

   HC_PREC *p;
   HC_PREC *dpdt;
   HC_PREC *pbyst;
   HC_PREC **torp;
   HC_PREC *cm;
   HC_PREC *sm;

   /* Torque arrays */

   HC_PREC *fin;
   HC_PREC *rots;
   HC_PREC *finet;
   HC_PREC **fex;
   HC_PREC *srt;
   HC_PREC *srp;
   HC_PREC *fins;
   HC_PREC **fexs;
   HC_PREC *ww;
   HC_PREC **vv;

   HC_PREC **trta;
   HC_PREC **trpa;
   HC_PREC **trra;
   HC_PREC ***srpa;
   HC_PREC ***srta;

  /* Plate id's */

   int **idp;

  /* For solve rots routines */

   HC_PREC *WW;
   HC_PREC **VV;
   
   /* velgrid */

   HC_PREC *pointx;
   HC_PREC *pointy;
   HC_PREC *pointz;
   HC_PREC *pointx2;
   HC_PREC *pointy2;
   HC_PREC *pointz2;
   HC_PREC **vtheta;
   HC_PREC **vphi;
   
  
};



/* hcplates_init.c */
void hc_init_params(struct hc_plates_params *);
double *c_alloc1d(int);
double **c_alloc2d(int, int);
int **i_alloc2d(int,int);
void Free1D(double *);
/*where's the rest */

void hc_init_arrays(struct hc_plates_params *, struct hc_plates_arrays *);
void hcplates_command_line(int, char **, struct hc_plates_params *);
void hcplates_advance_argument(int *, int, char **);
void read_parameter_file(struct hc_plates_params *);


/* zero_arrays.c */
void zero_arrays(struct hc_plates_params *,struct hc_plates_arrays *);

/*plates_propagator.c */ 
void mamlt2(double [3][3], double [3][3], double [3][3]);
void mamlt4(double [5][5], double [5][5], double [5][5]);
void prop4(int ,double, double, double [5][5]);
void prop2(int ,double, double, double [3][3]);
void torvel(struct hc_plates_params *,int, int, double [3], double [3]);
void polvel(struct hc_plates_params *,int,int, double [5], double [5]);
void poload(struct hc_plates_params *,int,int,int,int,double,double, double [5], double [5]);

/* get_loads.c */
void get_loads(struct hc_plates_params *, struct hc_plates_arrays *, char []);

/* Read unitrots */
void read_unitrots_get_coeff(struct hc_plates_params *, struct hc_plates_arrays *, char []);

/*crlb_vector_plm.c	*/
void crlb_vecplm(struct hc_plates_params *, struct hc_plates_arrays *, double, int);
void crlb_polyfact(int,double []);
void crlb_plm(double,int,struct hc_plates_arrays *);

/* integrate_torques.c */
void integrate_torques(struct hc_plates_params *, struct hc_plates_arrays *, char []);

/* solve_rot.c */
void solve_rot(struct hc_plates_params *, struct hc_plates_arrays * );
double pythag(double,double);
void svdcmp(struct hc_plates_params *,struct hc_plates_arrays *,int,int);
void svbksb(struct hc_plates_params *, struct hc_plates_arrays *, int,int);
/* in crlb_residual.c */
void crlb_residual(struct hc_plates_arrays *, int,int, double *, double *);

/* in hc_plates_write_output.c */
void hc_poles(struct hc_plates_params *, struct hc_plates_arrays *, char[]);
void hcplates_velgrid(struct hc_plates_params *, struct hc_plates_arrays *, char[]);
void vspher(struct hc_plates_params *, struct hc_plates_arrays *);
void vellin(struct hc_plates_params *, struct hc_plates_arrays *, int *, double, double, double *, double *, double *);
void point(struct hc_plates_params *, struct hc_plates_arrays *, int *, double, double, double *, double *, double *);
void point2(struct hc_plates_params *, struct hc_plates_arrays *, int *, double, double, double *, double *, double *);



