#include <math.h>
#include <stdio.h>
#include <stdlib.h>
// for numrec
#define NR_END 1
#define FREE_ARG char*

//bunch of correlation and cross correlation routines from Numerical Recipes
//see their copyright 



void compute_correl(double *, double *, double **, int, int *, int);
void correl(double *, double *, int, double *);
void twofft(double *, double *, double *, double *, int);
void four1(double *, int, int);
void realft(double *, int, int);
void spear(double *, double *, unsigned long, double *, double *, double *, double *, double *);
double *vector(int, int);
void free_vector(double *, int, int);
void pearsn(double *, double *, unsigned long, double *, double *, double *);
double betai(double, double, double);
double betacf(double, double, double);
void nrerror(char *);
double gammln(double);
double erfcc(double);
void crank(unsigned long, double *, double *);
void sort2(unsigned long, double *, double *);
int *ivector(long, long);
void free_ivector(int *, long, long);


#define ME {fprintf(stderr,"memory error\n");exit(-1);}
#define INVERSE_LN_TWO  1.44269504088896341
#define log2(x)         (log(x) * INVERSE_LN_TWO)
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
