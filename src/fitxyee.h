/* fitxyee.c */
#include "hc.h"
#define NR_ACC (HC_EPS_PREC*10)
#define NR_EPS HC_EPS_PREC
/*

  fit a line through x y data with uncertainties in both x and y

  Copyright Numerical Recipes in C, p.668, do not distribute without
  permission


  minor modifications:

  - using data structure instead of x,y,sigx,sigy
  - removed global variables and put those into fit structure

*/
static HC_PREC _tmp_sqrarg;
#define NR_SQUARE(a) ((((_tmp_sqrarg=(a))) == 0.0) ? (0.0) : (_tmp_sqrarg*_tmp_sqrarg))
#define NR_SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
static HC_PREC _tmp_maxarg1,_tmp_maxarg2;
#define NR_FMAX(a,b) (_tmp_maxarg1=(a),_tmp_maxarg2=(b),(_tmp_maxarg1) > (_tmp_maxarg2) ? \
        (_tmp_maxarg1) : (_tmp_maxarg2))
/* 
   data structure 
*/
struct nr_dp{
  HC_PREC x,y,sigx,sigy;
};
/* 
   fit structure 
*/
struct nr_fits{
  int nn;
  HC_PREC *xx,*yy,*sx,*sy,*ww,aa,offs;
};



void nr_fit(struct nr_dp *, int, HC_PREC *, HC_PREC *, HC_PREC *, HC_PREC *, HC_PREC *, HC_PREC *);
void nr_fitexy(struct nr_dp *, int, HC_PREC *, HC_PREC *, HC_PREC *, HC_PREC *, HC_PREC *, HC_PREC *);
HC_PREC nr_brent(HC_PREC, HC_PREC, HC_PREC, HC_PREC (*)(void), struct nr_fits *, HC_PREC, HC_PREC *);
HC_PREC nr_chixy(HC_PREC, struct nr_fits *);
HC_PREC nr_gammq(HC_PREC, HC_PREC);
void nr_gcf(HC_PREC *, HC_PREC, HC_PREC, HC_PREC *);
void nr_gser(HC_PREC *, HC_PREC, HC_PREC, HC_PREC *);
void nr_error(char []);
HC_PREC *nr_vector(long, long);
int *nr_ivector(long, long);
unsigned char *nr_cvector(long, long);
unsigned long *nr_lvector(long, long);
HC_PREC *nr_dvector(long, long);
HC_PREC **nr_matrix(long, long, long, long);
HC_PREC **nr_dmatrix(long, long, long, long);
int **nr_imatrix(long, long, long, long);
HC_PREC **nr_submatrix(HC_PREC **, long, long, long, long, long, long);
HC_PREC **nr_convert_matrix(HC_PREC *, long, long, long, long);
HC_PREC ***nr_f3tensor(long, long, long, long, long, long);
void nr_free_vector(HC_PREC *, long, long);
void nr_free_ivector(int *, long, long);
void nr_free_cvector(unsigned char *, long, long);
void nr_free_lvector(unsigned long *, long, long);
void nr_free_dvector(HC_PREC *, long, long);
void nr_free_matrix(HC_PREC **, long, long, long, long);
void nr_free_dmatrix(HC_PREC **, long, long, long, long);
void nr_free_imatrix(int **, long, long, long, long);
void nr_free_submatrix(HC_PREC **, long, long, long, long);
void nr_free_convert_matrix(HC_PREC **, long, long, long, long);
void nr_free_f3tensor(HC_PREC ***, long, long, long, long, long, long);
void nr_avevar(struct nr_dp *, unsigned long, HC_PREC *, HC_PREC *);
void nr_fitline(HC_PREC *, HC_PREC *, int, HC_PREC *, int, HC_PREC *, HC_PREC *, HC_PREC *, HC_PREC *, HC_PREC *, HC_PREC *);
void nr_mnbrak(HC_PREC *, HC_PREC *, HC_PREC *, HC_PREC *, HC_PREC *, HC_PREC *, HC_PREC (*)(void), struct nr_fits *);
HC_PREC nr_zbrent(HC_PREC (*)(void), struct nr_fits *, HC_PREC, HC_PREC, HC_PREC);
HC_PREC nr_gammln(HC_PREC);
int nr_comparef(struct nr_dp *, struct nr_dp *);
