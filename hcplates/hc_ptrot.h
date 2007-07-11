#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


#define HC_MAX_LAYERS 100 /* Hope this is enough radial layers */
#define MAX_NO_PLATE 30
#define MAX_LDIM 100
#define MAX_NPDIM 90
#define MAX_KDIM 515	/* Note: having this too high gives seg faults But we need it... weird.*/

#define HC_CHAR_LENGTH 300		/* length of char arrays */
#define HC_PREC double
#define HC_PI 3.1415926535897932384626433832795


struct hc_ptrot{
	char boundary_file[HC_CHAR_LENGTH];
	char data_filename[HC_CHAR_LENGTH];
	char plateids_file[HC_CHAR_LENGTH];
	char unitrot_file[HC_CHAR_LENGTH];
	
	int nlm;
	
	double plm;
	
	double tvelxm;
	double tvelym;
	double tvelzm;
	
	double tvelxp;
	double tvelyp;
	double tvelzp;
	
	HC_PREC ylm1;
	HC_PREC ylm2;
	
	/* Define arrays */
	HC_PREC ****dlmx;
	HC_PREC ****dlmy;
	HC_PREC ****dlmz;
	
	HC_PREC ****vlmx;
	HC_PREC ****vlmy;
	HC_PREC ****vlmz;
	
	HC_PREC *****dlm;
	HC_PREC *****vlm;
	
	int *num;
};

/* In ptrot.c */
void polcoeff(struct hc_ptrot *, int);
void torcoeff(struct hc_ptrot *);

void spharm(int, int, double, double, struct hc_ptrot *);
void plgndr(int, int, double, struct hc_ptrot *);
void vel(double, double, struct hc_ptrot *);

/* in hc_ptrot_mem.c */
int *c_alloc1d(int);
double ****c_alloc4d(int, int, int, int);
double *****c_alloc5d(int, int, int, int, int);

/* in hc_ptrot_init.c */
void hc_ptrot_init(struct hc_ptrot *);

void hc_ptrot_command_line(int, char **, struct hc_ptrot *);
void hc_ptrot_advance_argument(int *, int, char **);
