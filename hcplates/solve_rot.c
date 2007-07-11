#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "hcplates.h"

/* Routine to solve for plate rotation components
*fin, *fexs, *srta, *srpa, *trta, *trpa, *trra
 */

void solve_rot( struct hc_plates_params *plates, struct hc_plates_arrays *A )
{
 HC_PREC tol,wmax,wmin,wdum,resx,bmag,ratw;
 int jk, jkk, mkk, NPDIM, i, izero, NLAT, NLONG, ilat, ilon,NPLT;
 HC_PREC sumt,sump,dtheta,theta,ctheta,stheta,dphi,phi,area,rtheta,id2,rphi,cphi,sphi,rrt,rrp;
 int iplt,idir,ii,j,ind,ind1;
 HC_PREC tnt,tnp,tnr,con,e1,e2,e3,e4,e5;

 NPDIM = plates->NPDIM;
 NLAT = plates->NLAT;
 NLONG = plates->NLONG;
 NPLT = plates->NPLT;

/* First condition the ww matrix, after singular value decomp */
printf("	Conditiong the ww matrix \n");
 tol=1.0e-4;
 for (jk=1; jk<=NPDIM; jk++) {
	A->fins[jk]=A->fin[jk];
	A->rots[jk]=A->fins[jk];
	}
 for (jkk=1; jkk<=NPDIM; jkk++) {
	for (mkk=1; mkk<=NPDIM; mkk++) {	
		A->fexs[jkk][mkk] = A->fex[jkk][mkk];
	}
	}


  /* These are gonna take some work...*/
  printf("	Calling svdcmp, fexs11=%lf fexs[%i][%i]=%lf fexs[1][22]=%lf\n",A->fexs[1][1],plates->NPDIM,plates->NPDIM,A->fexs[plates->NPDIM][plates->NPDIM],A->fexs[1][22]);
  //printf(" WW=%le \n",A->WW[30]); 
  svdcmp(plates,A,NPDIM,NPDIM);	/*This had 2 too many NPDIMS again ... */
  //printf(" WW=%le \n",A->WW[30]); 
  printf("	Done svdcmp\n");
  wmax = -1.0e15;
  wmin = 1.0e15;	

  for (i=1; i<=NPDIM; i++) {
	if (wmax < A->WW[i]) {
		wmax = A->WW[i];
	}
	if (wmin < A->WW[i]) {
		wmin = A->WW[i];
	}
   }
  izero = 0;
  for (j=1; j<=NPDIM; j++) {
	wdum = A->WW[j]/wmax;
	if (wdum <= tol ) {
		/* Output to infomat (? wanna ?) */
		//printf(" j:%i WW:%le wmax:%le wdum:%le",j,A->WW[j],wmax,wdum);
		A->WW[j] = 0.0;
		izero++;
	}
  }

  ratw = wmin/wmax;
  /* Write to infomat ?? */

  printf("	Solve for rotation components: svbsb\n");
  /* Now solve for rotation components */
  //printf("  1. rots[3]=%le fexs22=%le WW=%le VV=%le\n",A->rots[3],A->fexs[2][2],A->WW[30],A->VV[2][2]);
  svbksb(plates,A,NPDIM,NPDIM);  /* This had too many args (NPDIMS) - check for consistency... */
//printf("  2. rots[3]=%le fexs22=%le WW=%le VV=%le\n",A->rots[3],A->fexs[2][2],A->WW[30],A->VV[2][2]);
  printf("	Calling residuals\n");
  crlb_residual(A,NPDIM,NPDIM,&resx,&bmag);
  /* Write infomat output ?? */
  
  sumt = 0.0;
  sumt = 0.0;
  dtheta = (plates->dd)*HC_PI/180;
  dphi = dtheta;
  theta = 180.5*dtheta;

  printf("	Looping through lats/longs\n");
  for (ilat=1;ilat<=NLAT;ilat++) {
	theta -= dtheta;
	ctheta = cos(theta);
	stheta = sin(theta);
	phi = -0.5*dphi;
	area = stheta*dtheta*dphi;
	rtheta = ilat - 90.5;
	printf("	  Solving rots - lat=%d\n",ilat);
	for (ilon=1;ilon<=NLONG;ilon++) {
		id2 = A->idp[ilat][ilon];
		phi += dphi;
		rphi = ilon - 0.5;
		cphi = cos(phi);
		sphi = sin(phi);
		rrt = 0.0;
		rrp = 0.0;

		/*Compute coordinate conversion factors */

		e1 = area*ctheta*cphi;
		e2 = area*sphi;
		e3 = area*ctheta*sphi;
		e4 = -area*cphi;
		e5 = -area*stheta;

		for (iplt=1;iplt<=NPLT;iplt++) {
			for (idir=1;idir<=3;idir++) {
				ii = 3*(iplt-1) + idir;
				rrt += A->srta[ilat][ilon][ii] * A->rots[ii];
				rrp += A->srpa[ilat][ilon][ii] * A->rots[ii];
			}
		}
		tnt = A->trta[ilat][ilon] - rrt;
		tnp = A->trpa[ilat][ilon] - rrp;
		tnr = A->trra[ilat][ilon];
		//printf("  %i/%i tnp: %le trp:%le rrp:%le\n",ilat,ilon,tnp,A->trpa[ilat][ilon],rrp);
		con = 1.0e-7;
		ind1 = 3*(id2 - 1);
		A->finet[ind1+1] += tnp*e1 + tnt*e2;
		A->finet[ind1+2] += tnp*e3 + tnt*e4;
		A->finet[ind1+3] += tnp*e5;

		/* Accumulate torques to see if net torque is zero */

		/* Write net shear stress: out.netstress (?) */
		/* Write resisting shear stress : out.resstress (?) */
	}
	//printf("   tnp=%le e1:%le e2:%le e3:%le e4:%le e5:%le\n",tnp,e1,e2,e3,e4,e5); 
  }

  for (i=1;i<=NPLT;i++) {
	ind = 3*(i-1);
	printf("FINET %i %le %le %le \n",i,A->finet[ind+1],A->finet[ind+2],A->finet[ind+3]);
  }
}
