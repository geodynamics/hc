#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "hcplates.h"

/*

Writing output.
FEX -> out.forces
fin*ratio ->out.forces

rot*ratio -> out.predrot

*/

void hc_poles(struct hc_plates_params *p, struct hc_plates_arrays *A, char pole_file[]) {
	FILE *filePtr;
	int i,tmp,NPLT;
	double x,y,z,angvel,rlat,rlong,w,rtheta,rphi,tmp1,tmp2,tmp3;
	
	NPLT = p->NPLT;
	filePtr = fopen(pole_file, "w");
	if (filePtr == NULL)
	{ printf("Error writing pole file. \n"); }
  else
	{
	printf("	Opened pole file file ok\n");

	for (i=1; i<=NPLT; i++) {
	tmp = 3*(i-1);
	tmp1 = A->rots[tmp+1]*p->stoy;
	tmp2 = A->rots[tmp+2]*p->stoy;
	tmp3 = A->rots[tmp+3]*p->stoy;
	
	x = tmp1*p->ratio;
	y = tmp2*p->ratio;
	z = tmp3*p->ratio;
	
	//printf("%i x %le y %le z %le ratio %le stoy %le \n",i,x,y,z,p->ratio,p->stoy);
	//printf(" Rots: %lf %lf %lf \n",tmp1,tmp2,tmp3);
	
	angvel = sqrt(x*x + y*y + z*z);
	rlat = acos(z/angvel);
	rlong = atan2(y,x);
	w = angvel*1e6*180.0/(M_PI*p->erad);
	rtheta = 90.0 - rlat*180.0/M_PI;
	rphi = rlong*180.0/M_PI;
	//printf(" angvel %le rlat %lf rlong %lf w %lf rtheta %lf rphi %lf \n",angvel,rlat,rlong,w,rtheta,rphi);
	
	/* OUTPUT */
	fprintf(filePtr," %i %lf %lf %lf \n",i,w,rtheta,rphi);
	}
	}
	fclose(filePtr);
}


void hcplates_velgrid(struct hc_plates_params *p, struct hc_plates_arrays *A, char velgrid_file[]) 
{
	FILE *filePtr;
	double pxlat,pxlong,ratio;
	int jl,jm;
			
	/* Open plates_map file for reading */
  filePtr = fopen(velgrid_file, "w");
  if (filePtr == NULL)
	{ printf("Error writing velgrid file. \n"); }
  else
	{
	printf("	Opened velgrid file ok\n");
	vspher(p,A);	
	printf("	Done vspher \n");	
	for (jl=1;jl<=p->NLAT;jl++) {
		pxlat = p->dd*(jl - 0.5) - 90.0;
		for (jm=1;jm<p->NLONG;jm++) {
			pxlong = p->dd*(jm - 0.5);
			A->vtheta[jl][jm] *= p->ratio;
			A->vphi[jl][jm] *= p->ratio;
			fprintf(filePtr," %lf %lf %le %le \n",pxlong,pxlat,A->vphi[jl][jm],A->vtheta[jl][jm]);
		}
	 }
	fclose(filePtr);
	}
}

void vspher(struct hc_plates_params *p, struct hc_plates_arrays *A) 
{
	int ik,k,i,idx;
	double npdiv3,plat,plong,ethex,ethey,ethez;
	double ephix,ephiy,ephiz,velx,vely,velz,velx2,vely2,velz2;
	
	//Open output files...
	npdiv3 = p->NPDIM/3;
	
	for (ik=1; ik<=npdiv3; ik++) {
		A->pointx[ik] = 0.0;
		A->pointy[ik] = 0.0;
		A->pointz[ik] = 0.0;
		A->pointx2[ik] = 0.0;
		A->pointy2[ik] = 0.0;
		A->pointz2[ik] = 0.0;
	}
	
	for (k=1; k <=p->NLAT; k++) {
		plat = (90.0 - (p->dd*(k-0.5) - 90.0))*M_PI/180.0;
		for (i=1; i<=p->NLONG; i++) {
			plong = (p->dd*(i-0.5))*M_PI/180.0;
			idx = A->idp[k][i];
			
			ethex = -sin(M_PI/2 - plat)*cos(plong);
			ethey = -sin(M_PI/2 - plat)*sin(plong);
			ethez = cos(M_PI/2 - plat);
			
			ephix = -sin(plong);
			ephiy = cos(plong);
			ephiz = 0.0;
			
			//printf(" ... plong:%lf plat:%lf\n",plong,plat);
			//printf(" calling vellin %lf \n",plat);
			vellin(p,A,&idx,plat,plong,&velx,&vely,&velz);
			//printf(" %i/%i ethex:%lf ethey:%lf ethez:%lf velx:%lf vely:%lf velz:%lf \n",k,i,ethex,ethey,ethez,velx,vely,velz);
			A->vtheta[k][i] = ethex*velx + ethey*vely + ethez*velz;
			//printf("       ephix:%lf ephiy:%lf    vtheta:%lf vphi:%lf \n",ephix,ephiy,A->vtheta[k][i],A->vphi[k][i]);
			A->vphi[k][i] = ephix*velx + ephiy*vely;
			
			/*   NOW DO WEIGHTING OF VELX VELY VELZ */
			//printf(" calling point %lf \n",vely);
			point(p,A,&idx,plat,plong,&velx,&vely,&velz);
			
			velx2 = pow(velx,2);
			vely2 = pow(vely,2);
			velz2 = pow(velz,2);
			//printf(" calling point2 %lf \n",vely2);
			point2(p,A,&idx,plat,plong,&velx2,&vely2,&velz2);
			
			/* DoNE WEIGHTING */
			
			/* OUTPUT LINEAR VEL COMPONENTS HERE */
			
		}	
		
	}
	
	
}

void vellin(struct hc_plates_params *p, struct hc_plates_arrays *A, int *idp, double plat, double plong, double *velx, double *vely, double *velz)
{
	double px,py,pz;
	int jj;
	
	px = sin(plat)*cos(plong);
	py = sin(plat)*sin(plong);
	pz = cos(plat);
	
	jj = 3*(*idp - 1);
	
	*velx = p->stoy*(A->rots[jj+2]*pz - A->rots[jj+3]*py);
	*vely = p->stoy*(A->rots[jj+3]*px - A->rots[jj+1]*pz);
	*velz = p->stoy*(A->rots[jj+1]*py - A->rots[jj+2]*px);
}

void point(struct hc_plates_params *p, struct hc_plates_arrays *A, int *idx, double plat, double plong, double *velx, double *vely, double *velz)
{
	int npdiv3;
	double dtheta,dphi,earea2;
	int iidx;
	
	iidx = (int) (*idx);
	
	dtheta = p->dd*M_PI/180.0;
	dphi = p->dd*M_PI/180.0;
	earea2 = pow(p->erad,2);
	
	A->pointx[iidx] += (*velx)*earea2*sin(plat)*dtheta*dphi;
	A->pointy[iidx] += (*vely)*earea2*sin(plat)*dtheta*dphi;
	A->pointz[iidx] += (*velz)*earea2*sin(plat)*dtheta*dphi;
}

void point2(struct hc_plates_params *p, struct hc_plates_arrays *A, int *idx, double plat, double plong, double *velx, double *vely, double *velz)
{
	int npdiv3;
	double dtheta,dphi,earea2;
	
	dtheta = p->dd*M_PI/180.0;
	dphi = p->dd*M_PI/180.0;
	earea2 = pow(p->erad,2);
	
	A->pointx2[*idx] += (*velx)*earea2*sin(plat)*dtheta*dphi;
	A->pointy2[*idx] += (*vely)*earea2*sin(plat)*dtheta*dphi;
	A->pointz2[*idx] += (*velz)*earea2*sin(plat)*dtheta*dphi;
}