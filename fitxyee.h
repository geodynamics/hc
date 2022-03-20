/* fitxyee.c */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define PRECISION double
//#define NR_ACC 1e-7
#define NR_ACC 1e-15
#define NR_EPS 5e-15
#define FOURFORMAT "%lf %lf %lf %lf"
#define TWOFORMAT "%lf %lf"

/*

  fit a line through x y data with uncertainties in both x and y

  $Id: fitxyee.c,v 1.2 2005/07/04 16:47:43 becker Exp becker $

  From Numerical Recipes in C, p.668


  major modifications:

  - using data structure instead of x,y,sigx,sigy
  - removed global variables and put those into fit structure

*/
static PRECISION _tmp_sqrarg;
#define SQUARE(a) ((((_tmp_sqrarg=(a))) == 0.0) ? (0.0) : (_tmp_sqrarg*_tmp_sqrarg))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
static PRECISION _tmp_maxarg1,_tmp_maxarg2;
#define FMAX(a,b) (_tmp_maxarg1=(a),_tmp_maxarg2=(b),(_tmp_maxarg1) > (_tmp_maxarg2) ?\
        (_tmp_maxarg1) : (_tmp_maxarg2))
/* 
   data structure 
*/
struct dp{
  PRECISION x,y,sigx,sigy;
};
/* 
   fit structure 
*/
struct fits{
  int nn;
  PRECISION *xx,*yy,*sx,*sy,*ww,aa,offs;
};



void fit(struct dp *, int, PRECISION *, PRECISION *, PRECISION *, PRECISION *, PRECISION *, PRECISION *);
void fitexy(struct dp *, int, PRECISION *, PRECISION *, PRECISION *, PRECISION *, PRECISION *, PRECISION *);
PRECISION brent(PRECISION, PRECISION, PRECISION, PRECISION (*)(void), struct fits *, PRECISION, PRECISION *);
PRECISION chixy(PRECISION, struct fits *);
PRECISION gammq(PRECISION, PRECISION);
void gcf(PRECISION *, PRECISION, PRECISION, PRECISION *);
void gser(PRECISION *, PRECISION, PRECISION, PRECISION *);
void nrerror(char []);
PRECISION *vector(long, long);
int *ivector(long, long);
unsigned char *cvector(long, long);
unsigned long *lvector(long, long);
PRECISION *dvector(long, long);
PRECISION **matrix(long, long, long, long);
PRECISION **dmatrix(long, long, long, long);
int **imatrix(long, long, long, long);
PRECISION **submatrix(PRECISION **, long, long, long, long, long, long);
PRECISION **convert_matrix(PRECISION *, long, long, long, long);
PRECISION ***f3tensor(long, long, long, long, long, long);
void free_vector(PRECISION *, long, long);
void free_ivector(int *, long, long);
void free_cvector(unsigned char *, long, long);
void free_lvector(unsigned long *, long, long);
void free_dvector(PRECISION *, long, long);
void free_matrix(PRECISION **, long, long, long, long);
void free_dmatrix(PRECISION **, long, long, long, long);
void free_imatrix(int **, long, long, long, long);
void free_submatrix(PRECISION **, long, long, long, long);
void free_convert_matrix(PRECISION **, long, long, long, long);
void free_f3tensor(PRECISION ***, long, long, long, long, long, long);
void avevar(struct dp *, unsigned long, PRECISION *, PRECISION *);
void fitline(PRECISION *, PRECISION *, int, PRECISION *, int, PRECISION *, PRECISION *, PRECISION *, PRECISION *, PRECISION *, PRECISION *);
void mnbrak(PRECISION *, PRECISION *, PRECISION *, PRECISION *, PRECISION *, PRECISION *, PRECISION (*)(void), struct fits *);
PRECISION zbrent(PRECISION (*)(void), struct fits *, PRECISION, PRECISION, PRECISION);
PRECISION gammln(PRECISION);
int comparef(struct dp *, struct dp *);
