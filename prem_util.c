/* 


routines to deal with PREM earth model by Dziewonki & Anderson (1981)

$Id: prem_util.c,v 1.1 2005/10/19 22:52:57 becker Exp $

twb@ig.utexas.edu



*/
#include "prem.h"
#include <math.h>
#include <stdlib.h>
/* 

given a radius in ruinits and n layers at radii, find the layer index
ilayer (0...n-1) and assigns x[0...3], where

output:

x[0] = 1, x[1] = r/r_unit (nondim radius), x[2]=x[1]^2, x[3]=x[1]^3

pass r and rb in the same units, and r_unit as r0 if units are [m], else as 1

*/
int prem_find_layer_x(double r, double r_unit,double *rb,
		      int n,int np, double *x)
{
  int i,ilayer,nm1;
  //
  // find layer
  //
  ilayer = 0;nm1=n-1;
  while((rb[ilayer] < r) && (ilayer < nm1))
    ilayer ++;
  //
  // get radii and radii to some power
  //
  x[0] = 1.0;
  x[1] = r/r_unit;// non dimensionalize r
  for(i=2;i < np;i++)
    x[i] = x[i-1]*x[1]; //r^i
  return ilayer;
}

/* 
  
compute the polynomial expansion of x at y[0,...np-1] and scale by
fac. y[0] = 1, y[1]=r', y[2]=r'^2

val = 

fac * (par[0]*y[0] + par[1]*y[1] + par[2]*y[2] + .. par[np-1]*y[np-1]) 

=

fac * (par[0]*1+ par[1]*r  + par[2]*r^2  + .. par[np-1]*r^(np-1)) 

*/
double prem_compute_pval(double *y, double *par, int np, 
			 double fac)
{
  int i;
  double val;
  val = par[0];
  for(i=1;i < np;i++)
    val += par[i] * y[i];
  val *= fac;
  return val;
}
/*

compute the derivative with  respect to y (r)

fac * (par[1]+ par[2]*r^1  + .. par[np-1]*r^(np-2)) 


*/
double prem_compute_dpval(double *y, double *par, 
			  int np, double fac)
{
  int i;
  double val;
  val = 0.0;
  for(i=1;i < np;i++)
    val += par[i] * y[i-1];
  val *= fac;
  return val;
}

/*

  calculate the mean S velocity from the Voigt shear modulus, see eg. Dahlen and
  Tromp:

  \mu = 1/15(A+C-2F+5N+6L) with
  \eta = F/(A-2L)
  and N, L, A, C corresponding to vsh, vsv, vph, and vpv respectively 
  according to, e.g, A = vph^2 \rho

  then

  vs=sqrt(1/15((1-2\eta)vph^2 + vpv^2 +5vsh^2+(6+4\eta)vsv^2))

  a simplified version without information about P and eta 
  
  
  vs = sqrt(1/3 * vsh^2 + 2/3 * vsv^2);
*/
double prem_vs_voigt(double vsh, double vsv, double vph, double vpv, 
		     double eta)
{
  double v;
  v = sqrt((1./15.) * ( (1.0-2.0*eta) * vph*vph + vpv*vpv + 
			5.0*vsh*vsh + 
			(6.0+4.0*eta)*vsv*vsv));
  return v;
}

/* 

compute density and drho/dr ([kg/m^3] and [kg/m^4], resp.)

rnd is the non-dimensionalized radius

*/
void prem_get_rhodrho(double *rho, double *drhodr,  
		      double rnd, struct prem_model *prem)
{
  int ilayer,os;
  double *x;
   if(!prem->init){
    fprintf(stderr,"prem_get_rhodrho: error, struc not init \n");
    exit(-1);
  }
   x=(double *)
    malloc(sizeof(double)*prem->np);
  if(!x){fprintf(stderr,"prem_get_rhodrho:mem error\n");exit(-1);}
  ilayer = prem_find_layer_x(rnd,1.0,prem->r,prem->n,prem->np,x);
  os = ilayer * prem->np;
  *rho = prem_compute_pval(x,(prem->crho+os),prem->np,1e3);
  *drhodr = prem_compute_dpval(x,(prem->crho+os),prem->np,1e3);
  free(x);
}
/* 
   get density in SI units [kg/m^3] at non-dim radius rnd
*/
void prem_get_rho(double *rho,double rnd, 
		  struct prem_model *prem)
{
  int ilayer,os;
  double *x;
  if(!prem->init){
    fprintf(stderr,"prem_get_rho: error, struc not init \n");
    exit(-1);
  }
  x=(double *)malloc(sizeof(double)*prem->np);
  if(!x){
    fprintf(stderr,"prem_get_rhodrho: mem error\n");
    exit(-1);
  }
  ilayer = prem_find_layer_x(rnd,1.0,prem->r,prem->n,prem->np,x);
  os = ilayer * prem->np;
  *rho = prem_compute_pval(x,(prem->crho+os),prem->np,1e3);
  free(x);
}



//
//
// compute several values from PREM at radius r [m]
// output of densities in kg/m^3 and velocities m/s
//
void prem_get_values(double *rho, double *vs,  double *vsv, double *vsh,
		     double *vp,  double *vpv, double *vph, double *eta,
		     double *qmu, double *qkappa,
		     double r, struct prem_model *prem)
{
  int ilayer,os;
  double *x;
   if(!prem->init){
    fprintf(stderr,"prem_get_values: error, struc not init \n");
    exit(-1);
  }
   x=(double *)malloc(sizeof(double)*prem->np);
  if(!x){fprintf(stderr,"prem_get_values:mem error\n");exit(-1);}
  ilayer = prem_find_layer_x(r,prem->r0,prem->rb,prem->n,prem->np,x);
  /* velocities and densities */
  os = ilayer * prem->np;
  /* density */
  *rho = prem_compute_pval(x, (prem->crho+os),prem->np,1e3);
  /* velocities */
  *vp =  prem_compute_pval(x,(prem->cvp+os), prem->np,1e3);
  *vs =  prem_compute_pval(x,(prem->cvs+os), prem->np,1e3);
  *vpv = prem_compute_pval(x,(prem->cvpv+os),prem->np,1e3);
  *vph = prem_compute_pval(x,(prem->cvph+os),prem->np,1e3);
  *vsv = prem_compute_pval(x,(prem->cvsv+os),prem->np,1e3);
  *vsh = prem_compute_pval(x,(prem->cvsh+os),prem->np,1e3);
  /* anisotropy ellipticity factor */
  *eta = prem_compute_pval(x,(prem->ceta+os),prem->np,1.0);
  /* q factors */
  *qmu = prem_compute_pval(x,   (prem->cqmu   +ilayer),1,1.0);
  *qkappa = prem_compute_pval(x,(prem->cqkappa+ilayer),1,1.0);
  free(x);
}

/* 

read in a PREM type model of Earth structure for velocities and densities 

for format, see example below and comments. n= number of layers,
most parameters p1, p2 ... have three coefficients, the last two
(q factors) only one

returns error code 

*/

int prem_read_model(char *filename,struct prem_model *prem, 
		    hc_boolean verbose)
{
  int i,rcnt;
  FILE *in;
  rcnt = 0;
  sprintf(prem->model_filename,"%s",filename);

  in=fopen(filename,"r");
  if(!in){
    fprintf(stderr,"prem_read_model: error: can't open model file %s\n",filename);
    return 1;
  }
  rcnt += fscanf(in,"%i",&prem->n);// nr layers
  rcnt += fscanf(in,"%i",&prem->np);// nr of coeff for most parameters
  if(prem->n > PREM_N){
    fprintf(stderr,"prem_read_model: error: read %i layers from file %s but max is %i\n",
	    prem->n,filename,PREM_N);
    fprintf(stderr,"prem_read_model: increase PREM_N\n");
    return 2;
  }
  if(prem->np > PREM_NP){
    fprintf(stderr,"prem_read_model: error: read %i parameters from file %s but max is %i\n",
	    prem->np,filename,PREM_NP);
    fprintf(stderr,"prem_read_model: increase PREM_NP\n");
    return 3;
  }
  rcnt += fscanf(in,PREM_F_STRING,&prem->r0); // read top level
  /* 
     upper layer boundary radii in [m]
  */
  rcnt += prem_read_para_set(prem->rb,prem->n,1,in);
  for(i=0;i < prem->n;i++)	/* non dim version */
    prem->r[i] = prem->rb[i]/prem->r0;
  if(rcnt != 3+prem->n){
    fprintf(stderr,"prem_read_model: read error: after layer input\n");
    return 4;
  }
  /* for each parameter, read in four polynomial coefficients for each
     layer  */
  /* 
     
  we expect the velocities to be in km/s and the densities in 
  g/cm^3
  
  */
  rcnt += prem_read_para_set(prem->crho,prem->n,prem->np,in);
  rcnt += prem_read_para_set(prem->cvp,prem->n,prem->np,in);
  rcnt += prem_read_para_set(prem->cvs,prem->n,prem->np,in);
  if(rcnt != 3+prem->n+3*prem->np*prem->n){
    fprintf(stderr,"prem_read_model: read error: after vp/vs/rho: %i\n",
	    rcnt);
    return 5;
  }
  rcnt += prem_read_para_set(prem->cvpv,prem->n,prem->np,in);
  rcnt += prem_read_para_set(prem->cvph,prem->n,prem->np,in);
  rcnt += prem_read_para_set(prem->cvsv,prem->n,prem->np,in);
  rcnt += prem_read_para_set(prem->cvsh,prem->n,prem->np,in);
  rcnt += prem_read_para_set(prem->ceta,prem->n,prem->np,in);
  if(rcnt != 3+prem->n+8*prem->np*prem->n){
    fprintf(stderr,"prem_read_model: read error: after anisotropic parameters\n");
    exit(-1);
  }
  /* 
     Q factors only are constant, only one parameter
  */
  for(i=0;i< prem->n;i++){// Q_mu
    rcnt += fscanf(in,PREM_F_STRING,&prem->cqmu[i]);
    if(fabs(prem->cqmu[i] - 1e10) < 1e-3){ /* == 1e10 */
#ifdef DBL_MAX
      prem->cqmu[i] = DBL_MAX;
#else
      prem->cqmu[i] = 3.40282347E+38;
#endif
    }
  }
  rcnt += prem_read_para_set(prem->cqkappa,prem->n,1,in);
  fclose(in);
  if(rcnt != 3 + prem->n+ 8*prem->np*prem->n + 2 * prem->n){
    fprintf(stderr,"prem_read_model: read error: after Q factors\n");
    return 6;
  }
  prem->init = TRUE;
  if(verbose)
    fprintf(stderr,"prem_read_model: initialized PREM model from file %s\n",
	    filename);
  return 0;
}
/* 
   read in a set of np parameters for n layers
*/
int prem_read_para_set(double *par, int n, int np,FILE *in)
{
  int i,j,rc,os;
  rc=0;
  for(i=os=0;i<n;i++,os+=np)
    for(j=0;j<np;j++)
      rc += fscanf(in,PREM_F_STRING,(par+os+j));
  return rc;
}
/* 

the following comment block holds the prem.dat file, for clarification

note that radii are given in [m]

this is the transversely isotropic version of PREM with the ocean
layer

*/
/*

        13 4
          6371000.
          1221500.  3480000.  3630000.  5600000.  5701000.  5771000.
            5971000.  6151000.  6291000.  6346600.  6356000.  6368000.
            6371000.
          
	   13.0885  0.0000 -8.8381  0.0000
           12.5815 -1.2638 -3.6426 -5.5281
            7.9565 -6.4761  5.5283 -3.0807
            7.9565 -6.4761  5.5283 -3.0807
            7.9565 -6.4761  5.5283 -3.0807
            5.3197 -1.4836  0.0000  0.0000
           11.2494 -8.0298  0.0000  0.0000
            7.1089 -3.8045  0.0000  0.0000
            2.6910  0.6924  0.0000  0.0000
            2.6910  0.6924  0.0000  0.0000
            2.9000  0.0000  0.0000  0.0000
            2.6000  0.0000  0.0000  0.0000
            1.0200  0.0000  0.0000  0.0000
            
           11.2622   0.0000 -6.3640   0.0000
           11.0487  -4.0362  4.8023 -13.5732
           15.3891  -5.3181  5.5242  -2.5514
           24.9520 -40.4673 51.4832 -26.6419
           29.2766 -23.6027  5.5242  -2.5514
           19.0957  -9.8672  0.0000   0.0000
           39.7027 -32.6166  0.0000   0.0000
           20.3926 -12.2569  0.0000   0.0000
            4.1875   3.9382  0.0000   0.0000
            4.1875   3.9382  0.0000   0.0000
            6.8000   0.0000  0.0000   0.0000
            5.8000   0.0000  0.0000   0.0000
            1.4500   0.0000  0.0000   0.0000

            3.6678   0.0000 -4.4475  0.0000
            0.0000   0.0000  0.0000  0.0000
            6.9254   1.4672 -2.0834  0.9783
           11.1671 -13.7818 17.4575 -9.2777
           22.3459 -17.2473 -2.0834  0.9783
            9.9839  -4.9324  0.0000  0.0000
           22.3512 -18.5856  0.0000  0.0000
            8.9496  -4.4597  0.0000  0.0000
            2.1519   2.3481  0.0000  0.0000
            2.1519   2.3481  0.0000  0.0000
            3.9000   0.0000  0.0000  0.0000
            3.2000   0.0000  0.0000  0.0000
            0.0000   0.0000  0.0000  0.0000

           11.2622   0.0000 -6.3640   0.0000
           11.0487  -4.0362  4.8023 -13.5732
           15.3891  -5.3181  5.5242  -2.5514
           24.9520 -40.4673 51.4832 -26.6419
           29.2766 -23.6027  5.5242  -2.5514
           19.0957  -9.8672  0.0000   0.0000
           39.7027 -32.6166  0.0000   0.0000
           20.3926 -12.2569  0.0000   0.0000
            0.8317   7.2180  0.0000   0.0000
            0.8317   7.2180  0.0000   0.0000
            6.8000   0.0000  0.0000   0.0000
            5.8000   0.0000  0.0000   0.0000
            1.4500   0.0000  0.0000   0.0000

           11.2622   0.0000 -6.3640   0.0000
           11.0487  -4.0362  4.8023 -13.5732
           15.3891  -5.3181  5.5242  -2.5514
           24.9520 -40.4673 51.4832 -26.6419
           29.2766 -23.6027  5.5242  -2.5514
           19.0957  -9.8672  0.0000   0.0000
           39.7027 -32.6166  0.0000   0.0000
           20.3926 -12.2569  0.0000   0.0000
            3.5908   4.6172  0.0000   0.0000
            3.5908   4.6172  0.0000   0.0000
            6.8000   0.0000  0.0000   0.0000
            5.8000   0.0000  0.0000   0.0000
            1.4500   0.0000  0.0000   0.0000

            3.6678   0.0000 -4.4475  0.0000
            0.0000   0.0000  0.0000  0.0000
            6.9254   1.4672 -2.0834  0.9783
           11.1671 -13.7818 17.4575 -9.2777
           22.3459 -17.2473 -2.0834  0.9783
            9.9839  -4.9324  0.0000  0.0000
           22.3512 -18.5856  0.0000  0.0000
            8.9496  -4.4597  0.0000  0.0000
            5.8582  -1.4678  0.0000  0.0000
            5.8582  -1.4678  0.0000  0.0000
            3.9000   0.0000  0.0000  0.0000
            3.2000   0.0000  0.0000  0.0000
            0.0000   0.0000  0.0000  0.0000

            3.6678   0.0000 -4.4475  0.0000
            0.0000   0.0000  0.0000  0.0000
            6.9254   1.4672 -2.0834  0.9783
           11.1671 -13.7818 17.4575 -9.2777
           22.3459 -17.2473 -2.0834  0.9783
            9.9839  -4.9324  0.0000  0.0000
           22.3512 -18.5856  0.0000  0.0000
            8.9496  -4.4597  0.0000  0.0000
           -1.0839   5.7176  0.0000  0.0000
           -1.0839   5.7176  0.0000  0.0000
            3.9000   0.0000  0.0000  0.0000
            3.2000   0.0000  0.0000  0.0000
            0.0000   0.0000  0.0000  0.0000

            1.0000   0.0000  0.0000  0.0000    
            1.0000   0.0000  0.0000  0.0000
            1.0000   0.0000  0.0000  0.0000    
            1.0000   0.0000  0.0000  0.0000
            1.0000   0.0000  0.0000  0.0000    
            1.0000   0.0000  0.0000  0.0000
            1.0000   0.0000  0.0000  0.0000    
            1.0000   0.0000  0.0000  0.0000    
            3.3687  -2.4778  0.0000  0.0000
            3.3687  -2.4778  0.0000  0.0000    
            1.0000   0.0000  0.0000  0.0000
            1.0000   0.0000  0.0000  0.0000          
	    1.0000   0.0000  0.0000  0.0000

	84.6
	1e10
	312
	312
	312
	143
	143
	143
	80
	600
	600
	600
	1e10

	1327.7
	57823
	57823
	57823
	57823
	57823
	57823
	57823
	57823
	57823
	57823
	57823
	57823

*/
