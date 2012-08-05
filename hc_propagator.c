/* 

propagator generation routines for hc_polsol part of HC code

$Id: hc_propagator.c,v 1.1 2004/07/01 23:52:23 becker Exp $

*/
#include "hc.h"

void hc_evalpa(int l,HC_HIGH_PREC r1,HC_HIGH_PREC r2,HC_HIGH_PREC visc, HC_HIGH_PREC *p)
{
  //
  //    ****************************************************************
  //    * THIS SUBROUTINE CALCULATES THE 4X4 PROPAGATOR MATRIX THROUGH *
  //    * USE OF AN ALGEBRAI//FORM (IE: ALL MATRIX MULTIPLICATIONS     *
  //    * PERFORMED ALGEBRAICLY, THE RESULT BEING CODED HEREIN).  THIS *
  //    * CODE IS A REVISION OF M. RICHARDS' SUBROUTINE, DESIGNED TO   *
  //    * INTERFACE WITH THE EXISTING PROGRAMS.                        *
  //    ****************************************************************
  //
  HC_HIGH_PREC den1,den2,f[4],r,rlm1,rlp1,rmlm2,rml,v2,rs;
  long int np[4][4][4];
  long int lp1,lp2,lp3,lm1,lm2,lpp,lmm,l2p3,l2p1,l2m1,lltp1,lltp2;
  long int i,j,k,os1,os2;
  //fprintf(stderr,"hc_evalpa: %i %g %g %g\n",l,r1,r2,visc);
  //
  //    PASSED PARAMETERS:  L: DEGREE,
  //    R1,R2: PROPAGATE FROM R1 TO R2 (RADII),
  //    VISC: VISCOSITY
  //    P: THE PROPAGATOR [4*4]
  //
  //    OTHER VARIABLES:
  //       NP: INTEGER FACTORS (FUNCTIONS OF L) IN THE EXPRESSIONS
  //          FOR THE ELEMENTS OF P,
  //       F: HC_HIGH_PREC PRECISION FACTORS (FUNCTIONS OF R AND L) FOR P,
  //       R IS THE RADIUS RATIO R/R0 (HC_HIGH_PREC PRECISION).
  //
  r = r2 / r1;

  lp1 = l + 1;
  lp2 = l + 2;
  lp3 = l + 3;
  lm1 = l - 1;
  lm2 = l - 2;

  lpp = l * lp3 - 1;		
  lmm = l * l - lp3;
  
  lltp1 = l * lp1;
  lltp2 = l * lp2;

  l2p3 = lp3 + l;
  l2p1 = lp1 + l;
  l2m1 = lm1 + l;

  v2 = visc * 2.0;
  rs = r * r;

  rlm1 = pow(r,(HC_HIGH_PREC)lm1);
  rlp1 = rlm1 * rs;
  rmlm2 = pow(r,(HC_HIGH_PREC)(-lp2));
  rml = rmlm2 * rs;

  den1 = (HC_HIGH_PREC)(l2p1 * l2p3);
  den2 = (HC_HIGH_PREC)(l2p1 * l2m1);

  f[0] = rlp1 / den1;
  f[1] = rlm1 / den2;
  f[2] = rml  / den2;
  f[3] = rmlm2/ den1;

  // DO INTEGER MULTIPLIES

  // FIRST, SET UP 16 'REFERENCE' ELEMENTS
  
  np[0][2][0] = -lltp1;
  np[0][2][1] = -np[0][2][0];
  np[0][2][2] =  np[0][2][0];
  np[0][2][3] =  np[0][2][1];

  np[1][2][0] = -lp3;
  np[1][2][1] =  lp1;
  np[1][2][2] =  lm2;
  np[1][2][3] = -l;

  np[2][2][0] = - lp1  * lmm;
  np[2][2][1] =  lltp1 * lm1;
  np[2][2][2] =      l * lpp;
  np[2][2][3] = -lltp1 * lp2;

  np[3][2][0] = -lltp2;
  np[3][2][1] =  lm1 * lp1;
  np[3][2][2] = -np[3][2][1];
  np[3][2][3] = -np[3][2][0];
 
  //generate other elements

  for(i=0;i < 4;i++){
    np[i][0][0] =  np[i][2][0] * lp2;
    np[i][1][0] = -np[i][0][0] * l;
    np[i][3][0] = -np[i][2][0] * l;
    
    np[i][0][1] = (np[i][2][1] * lpp) / lp1;
    np[i][1][1] = -np[i][2][1] * lp1 * lm1;
    np[i][3][1] = -np[i][2][1] * lm2;
    
    np[i][0][2] = -np[i][2][2] * lm1;
    np[i][1][2] =  np[i][0][2] * lp1;
    np[i][3][2] =  np[i][2][2] * lp1;

    np[i][0][3] =(-np[i][2][3] * lmm) / l;
    np[i][1][3] = -np[i][2][3] * lltp2;
    np[i][3][3] =  np[i][2][3] * lp3;
    
  }
  //put IN COMMON FACTORS
  
  for(os1=i=0;i < 4;i++,os1+=4){
    for(j=0;j<4;j++){
      os2 = os1 + j;
      p[os2] = 0.0;
      for(k=0;k < 4;k++)
	p[os2] += ((HC_HIGH_PREC)(np[i][j][k])) * f[k];
    }
  }
  
  // PUT IN VISCOSITIES

  p[0*4+2] /= v2;
  p[0*4+3] /= v2;
  p[1*4+2] /= v2;
  p[1*4+3] /= v2;

  p[2*4+0] *= v2;
  p[2*4+1] *= v2;
  p[3*4+0] *= v2;
  p[3*4+1] *= v2;
}

void hc_evppot(int l,HC_HIGH_PREC ratio, HC_HIGH_PREC *ppot)
{
  //    ********************************************
  //    * THIS SUBROUTINE OBTAINS THE POTENTIALS   *
  //    * PROPAGATOR FROM R1 TO R2 AT L, WHERE     *
  //    * RATIO = R1/R2.                           *
  //    ********************************************
  //
  HC_HIGH_PREC  c,expf1,expf2,x,xp1;
  
  //    PASSED PARAMETERS:  L: DEGREE,
  //       RATIO: R1 / R2 (R2>R1),
  //       PPOT: THE POTENTIALS PROPAGATOR [4]
  //
  //    HC_HIGH_PREC PRECISION:  C: COEFFICIENT MULTIPLYING PROPAGATOR,
  //       EXPF1: RATIO**L, EXPONENTIAL FACTOR IN PPOT,
  //       EXPF2: (1/RATIO)**(L+1), EXP. FACTOR IN PPOT,
  //       X: DEGREE (L),
  //       XP1: X+1.
  
  x = (HC_HIGH_PREC)l;
  c = 1.0 / (2.0 * x + 1.0);
  xp1 = x + 1.0;
  expf1 = pow(ratio,x);
  expf2 = 1.0 / (expf1*ratio);
  ppot[0] = c * (x * expf1 + xp1 * expf2);
  ppot[1] = c * (expf2 - expf1);
  ppot[2] = x * xp1 * ppot[1];
  ppot[3] = ppot[0] - ppot[1];
}
