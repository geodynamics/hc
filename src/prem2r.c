#include "hc.h"
/* 

read the PREM model and print density with depth, compute mass,
gravity and pressure

this is a super sloppy routine and should not be used for any real
work

*/
double fm(double ,double *, struct prem_model *);

int main(int argc,char **argv)
{
  int i,n;
  double rho,rho2,r,z,rnd,dr,drh,m,old_r,mfac,old_mfac,*zp,*gp,*mp,*p;
  struct prem_model prem[1];
  char filename[HC_CHAR_LENGTH];
  if(argc>1)
    strncpy(filename,argv[1],HC_CHAR_LENGTH);
  else
    strncpy(filename,PREM_MODEL_FILE,HC_CHAR_LENGTH);

  prem_read_model(filename,prem,TRUE);

  dr=.00001;			/* in km */
  
  drh=dr/2;

  rnd = old_r = old_mfac=0;
  m=0;
  zp = (double *)malloc(sizeof(double));
  mp = (double *)malloc(sizeof(double));
  gp = (double *)malloc(sizeof(double));
  
  n=0;
  /* 
     integrate mass and compute M(r) and g(r) 
  */
  for(r=dr;r <= HC_RE_KM;r += dr){
    z = HC_RE_KM-r;
    rnd = r/HC_RE_KM;
    mfac = fm(rnd,&rho,prem);
    m += drh * (mfac+old_mfac);
    if(r-old_r >= .05-1e-7){
      zp[n] = z;mp[n] = m;
      gp[n] = HC_CAPITAL_G * m/(r*r*1e6); /* compute g */
      n++;
      zp = (double *)realloc(zp,sizeof(double)*(n+1));
      mp = (double *)realloc(mp,sizeof(double)*(n+1));
      gp = (double *)realloc(gp,sizeof(double)*(n+1));
      old_r = r;
    }
    old_mfac = mfac;
  }
  zp[n] = z;mp[n] = m;
  gp[n] = HC_CAPITAL_G * m/(r*r*1e6); 
  n++;

  
  p = (double *)calloc(n,sizeof(double));
  /* 
     compute pressure from overburden, mid point method
  */
  for(i=n-2;i>=0;i--){
    prem_get_rho(&rho, (HC_RE_KM - zp[i]  )/HC_RE_KM,prem);
    prem_get_rho(&rho2,(HC_RE_KM - zp[i+1])/HC_RE_KM,prem);
    p[i] = p[i+1] + (gp[i]*rho + gp[i+1]*rho2)/2.*(zp[i]-zp[i+1])*1e3;
  }

  
  for(i=0;i<n;i++){
    r = HC_RE_KM-zp[i];
    rnd = r/HC_RE_KM;
    prem_get_rho(&rho,rnd,prem);
    /* output is z[km] r rho mp gp p  */
    if(fabs(zp[i]-(int)zp[i]) < 0.01)
      fprintf(stdout,"%11g %11g %12.8e %12.8e %12.8e %12.8e\n",zp[i],r,rho,mp[i],gp[i],p[i]);
  }
  
  free(mp);free(gp);free(p);free(zp);
  return 0;
}

double fm(double rnd,double *rho, struct prem_model *prem)
{
  /* 4pi/3 r^2 \rho(r)  */
  const double fac = 4*M_PI*HC_RE_KM*HC_RE_KM*1e9;

  prem_get_rho(rho,rnd, prem);
  return fac*rnd*rnd*(*rho);
}
