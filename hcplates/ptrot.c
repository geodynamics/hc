#include <stdio.h>
#include "hc_ptrot.h"

/* 
   
   stand-alone routine to take plate boundary file, and metadata file with no. plates + no. points on each plate
   and a file with plate_id lat long for every lat/long on a 180x360 grid, then create a unitrot.coeff file
   (ie coeffs for a unit rotation on each plate - needed for hcplate routine).
   
   For Lithgow-Bertelloni's & Richards Plates extensions to Hager & O'Connell experimental code
   O'Neill (coneill@els.qm.edu.au) & Becker (twb@usc.edu).

   This has routines to: set initial defaults
			handle the command line args (ie. just the input file name)
			read the input file, and assign parameters
  
  Tested against CRLB's fortran version, and working ok
  
*/

int main(int argc, char **argv)
{
struct hc_ptrot p[1];
FILE *filePtr1, *filePtr2;
char file1[HC_CHAR_LENGTH], file2[HC_CHAR_LENGTH];
int nplate,ii,ip,m,l,k,kk,ll,mm;
double factor,factor2,dlm1,dlm2,xx,yy,vlm1,vlm2;

fprintf(stderr," Initializing....\n");
hc_ptrot_init(p);
xx = 0.0;
yy = 0.0;

fprintf(stderr," Reading command line args....\n");
hc_ptrot_command_line(argc,argv,p);

fprintf(stderr,"Plate boundary file is %s\n",p->boundary_file);
fprintf(stderr,"Plate data file is %s\n",p->data_filename);
fprintf(stderr,"Plate IDs file is %s\n",p->plateids_file);
fprintf(stderr,"Output file is %s\n",p->unitrot_file);
fprintf(stderr,"Default degree is %i\n",p->nlm);

/* Opening up plate info file (enes) for a looksies */
	strncpy(file1,p->data_filename,HC_CHAR_LENGTH);
	filePtr1 = fopen(file1, "r");
	if (filePtr1 == NULL)
	{ printf("Error loading data file. \n"); }
	else {
		fscanf(filePtr1, "%i", &nplate);
		printf("nplate: %i \n",nplate);
		for (ii=1; ii<=nplate; ii++) {
			/* Read number of points on each plate */
			fscanf(filePtr1, "%i", &p->num[ii]);
			printf("npoint: %i \n", p->num[ii]);
		}
	}

fclose(filePtr1);
polcoeff(p,nplate);
/*fprintf(stderr,"p->dlmx: %lf\n",p->dlmx[1][1][1][1]);
fprintf(stderr,"p->dlmy: %lf\n",p->dlmx[2][1][4][5]);
fprintf(stderr,"p->vlmz: %lf\n",p->vlmz[8][2][5][7]); 
polcoeef seems to be working ok*/
torcoeff(p);

for (ip=1; ip<= nplate; ip++) {
	for (l=0; l<=p->nlm; l++) {
		for (m=0; m<=l; m++) {
			p->dlm[ip][1][1][l][m] = p->dlmx[ip][1][l][m];
			p->dlm[ip][1][2][l][m] = p->dlmx[ip][2][l][m];
			p->dlm[ip][2][1][l][m] = p->dlmy[ip][1][l][m];
			p->dlm[ip][2][2][l][m] = p->dlmy[ip][2][l][m];
			p->dlm[ip][3][1][l][m] = p->dlmz[ip][1][l][m];
			p->dlm[ip][3][2][l][m] = p->dlmz[ip][2][l][m];
			
			p->vlm[ip][1][1][l][m] = p->vlmx[ip][1][l][m];
			p->vlm[ip][1][2][l][m] = p->vlmx[ip][2][l][m];
			p->vlm[ip][2][1][l][m] = p->vlmy[ip][1][l][m];
			p->vlm[ip][2][2][l][m] = p->vlmy[ip][2][l][m];
			p->vlm[ip][3][1][l][m] = p->vlmz[ip][1][l][m];
			p->vlm[ip][3][2][l][m] = p->vlmz[ip][2][l][m];
		}
	}
}
strncpy(file2,p->unitrot_file,HC_CHAR_LENGTH);
	filePtr2 = fopen(file2, "w");
	if (filePtr2 == NULL)
	{ printf("Error loading data file. \n"); }
	else {
		fprintf(filePtr2,"%i \n", p->nlm);
		for (k=1; k<= nplate; k++) {
			for (kk=1; kk <=3; kk++) {
				for (l=0; l <= p->nlm; l++) {
					for (m=0; m <= l; m++) {
						factor = l*(l+1);
						if (l != 0) {
							dlm1 = p->dlm[k][kk][1][l][m]/factor;
							dlm2 = p->dlm[k][kk][2][l][m]/factor;
							fprintf(filePtr2,"  %i	%i	%le	%le\n",l,m,dlm1,dlm2);
						}
						else {
							fprintf(filePtr2,"  %i	%i	%le	%le \n",l,m,xx,yy);
						}
					}
				}
				for (ll=0; ll <= p->nlm; ll++) {
					for (mm = 0; mm <= ll; mm++) {
						factor2 = ll*(ll+1);
						vlm1 = p->vlm[k][kk][1][ll][mm]/factor2;
						vlm2 = p->vlm[k][kk][2][ll][mm]/factor2;
						if (ll != 0) {
							fprintf(filePtr2,"%i %i %lf %lf\n",ll,mm,vlm1,vlm2);
						}
						else {
							fprintf(filePtr2,"%i %i %lf %lf \n",ll,mm,xx,yy);
						}
					}
				}
			}
		}
	
}
fclose(filePtr2);
fprintf(stderr," Done!\n");
return 0;
}



/*													*/

void polcoeff(struct hc_ptrot *p, int nplate)
{
	FILE *filePtr3;
	char file3[HC_CHAR_LENGTH];
	double pbla[16][901],pblo[16][901],px[901],py[901],pz[901],decoy1,decoy2;
	double ex,ey,ez,rlinx,rliny,rlinz,theta,phi,rlmag,rnnx,rnny,rnnz;
	double rnmag,vxm,vym,vxp,vyp,vzm,vzp,vdotnxp,vdotnxm,vdotnyp,vdotnym,vdotnzp,vdotnzm;
	double rx1,ry1,rz1,vranx,vrany,vranz,x;
	int mj,mmj,np,i,k,m,l;
	
	/*Opening up boundary file */
	strncpy(file3,p->boundary_file,HC_CHAR_LENGTH);
	filePtr3 = fopen(file3, "r");
	if (filePtr3 == NULL)
	{ printf("Error loading data file. \n"); }
	else {
		for (mj=1; mj <=nplate; mj++) {
			for (mmj=1; mmj<= p->num[mj]; mmj++) {
				//fprintf(stderr,"	Scanning file... %i %i \n",mj,mmj);
				fscanf(filePtr3,"%lf %lf",&pbla[mj][mmj],&pblo[mj][mmj]);
				//printf("	Scanning file...Done \n");
				decoy1 = pbla[mj][mmj];
				decoy2 = pblo[mj][mmj];
				pbla[mj][mmj] =  HC_PI*(0.5 - decoy1/180.0);
				pblo[mj][mmj] = HC_PI*decoy2/180.0;
			}
		}
	}
	/* Finished loading boundary data */
	fclose(filePtr3);
	
	for (np=1; np<=nplate; np++) {
			ex=1.0;
			ey=1.0;
			ez=1.0;
			for (i=1; i<= p->num[np]; i++) {
					px[i] = sin(pbla[np][i])*cos(pblo[np][i]);
					py[i] = sin(pbla[np][i])*sin(pblo[np][i]);
					pz[i] = cos(pbla[np][i]);
			}
			for (k=1; k<= p->num[np]; k++) {
				if (k==1) {
					rlinx = (px[2] - px[p->num[np]])/2;
					rliny = (py[2] - py[p->num[np]])/2;
					rlinz = (pz[2] - pz[p->num[np]])/2;
				}
				if (k == p->num[np]) {
					rlinx = (px[1] - px[p->num[np]-1])/2;
					rliny = (py[1] - py[p->num[np]-1])/2;
					rlinz = (pz[1] - pz[p->num[np]-1])/2;
				}
				if ((k != 1)&&(k != p->num[np])) {
					rlinx = (px[k+1] - px[k-1])/2;
					rliny = (py[k+1] - py[k-1])/2;
					rlinz = (pz[k+1] - pz[k-1])/2;
				}
				theta = pbla[np][k];
				phi = pblo[np][k];
				rlmag = sqrt(rlinx*rlinx + rliny*rliny + rlinz*rlinz);
				rnnx = rliny*pz[k] - rlinz*py[k];
				rnny = rlinz*px[k] - rlinx*pz[k];
				rnnz = rlinx*py[k] - rliny*px[k];
				rnmag = sqrt(pow(rnnx,2) + pow(rnny,2) + pow(rnnz,2)); 
				
				vxm = -ey*pz[k];
				vxp = py[k]*ez;
				vym = -ez*px[k];
				vyp = pz[k]*ez;
				vzm = -ex*py[k];
				vzp = px[k]*ey;
				
				if (rnmag == 0.0) {
					vdotnxp = 0.0;
					vdotnxm = 0.0;
					vdotnyp = 0.0;
					vdotnym = 0.0;
					vdotnzp = 0.0;
					vdotnzm = 0.0;
				}
				else {
					rx1 = rnnx/rnmag;
					ry1 = rnny/rnmag;
					rz1 = rnnz/rnmag;
					
					vdotnxp = rlmag*(vxp*rx1);
					vdotnxm = rlmag*(vxm*rx1);
					vdotnyp = rlmag*(vyp*ry1);
					vdotnym = rlmag*(vym*ry1);
					vdotnzp = rlmag*(vzp*rz1);
					vdotnzm = rlmag*(vzm*rz1);
				}
				vranx = (vdotnyp+vdotnzm);
				vrany = (vdotnxm+vdotnzp);
				vranz = (vdotnxp+vdotnym);
				x = cos(theta);
				for (l=0; l<=p->nlm; l++) {
					for (m=0; m<=l; m++) {
						spharm(l,m,x,phi,p);
						p->dlmx[np][1][l][m] += vranx*p->ylm1;
						p->dlmx[np][2][l][m] += vranx*p->ylm2;
						p->dlmy[np][1][l][m] += vrany*p->ylm1;
						p->dlmy[np][2][l][m] += vrany*p->ylm2; 
						p->dlmz[np][1][l][m] += vranz*p->ylm1;
						p->dlmz[np][2][l][m] += vranz*p->ylm2;
					}
				}
			}		
		
	}
	
	fprintf(stderr,"Finished polcoeffs \n");
}

void torcoeff(struct hc_ptrot *p) {
	FILE *filePtr4;
	char file4[HC_CHAR_LENGTH];
	double vthetax[201][401],vphix[201][401],vthetay[201][401],vphiy[201][401];
	double vthetaz[201][401],vphiz[201][401],toss1,toss2;
	int idj,idk,idp,k,i,kj,ij,m,l,cnt0,cnt1,cnt2,id[181][361],count,toss3,toss4;
	double dtheta,dphi,plat,plong,ethex,ethey,ethez,ephix,ephiy,ephiz,eradx,erady,eradz;
	double dvthetaxm,dvthetaym,dvthetazm,dvthetaxp,dvthetayp,dvthetazp;
	double dvphixp,dvphixm,dvphiyp,dvphiym,dvphizp,dvphizm;
	double term1x,term1y,term1z,term2xp,term2yp,term2zp,term2xm,term2ym,term2zm;
	double term3xp,term3yp,term3zp,term3xm,term3ym,term3zm,theta,avephi,x,facx;
	/* added this TWB */
	toss4 = 360;
	/*Opening up boundary file */
	strncpy(file4,p->plateids_file,HC_CHAR_LENGTH);
	filePtr4 = fopen(file4, "r");
	if (filePtr4 == NULL)
	{ printf("Error loading plate ids file. \n"); }
	else {
		count=0;
		for (idj=1; idj<=180; idj++) {
			for (idk=1; idk<=360; idk++) {
				//fscanf(filePtr4,"%i %lf %lf %i %i",&id[idj][idk],&toss1,&toss2,&toss3,&toss4);	
				fscanf(filePtr4,"%i %lf %lf ",&id[idj][idk],&toss1,&toss2);
				count++;
				/*if (idj != toss3 || idk != toss4 || id[idj][idk]==0) {
					fprintf(stderr," %i %i %i %i %i %i %lf %lf \n",count,idj,idk,toss3,toss4,id[idj][idk],toss1,toss2);
					}*/
			}
		}
	}
	if (count != 64800 || toss4 != 360) {
		fprintf(stderr,"Weird line numbers in %s, id array probably off, this sometimes cause seg faults \n",p->plateids_file);
	}
	fprintf(stderr,"In torcoeffs ...working\n");
	//fprintf(stderr,"	id(179,355):%i		%i\n",id[179][355],id[180][355]);
	dtheta = HC_PI/180.0;
	dphi = HC_PI/180.0;
	for (k=1; k<=180; k++) {
		plat = (90.0 - ((k-0.5) - 90.0))*HC_PI/180.0;
		for (i=1; i<=360; i++) {
				plong = (i-0.5)*HC_PI/180.0;
				idp = id[k][i];
				//fprintf(stderr," 1a.		k=%i i=%i  id=%i\n",k,i,id[k][i]); //id ok here
				
				ethex = cos(plat)*cos(plong);
				ethey = cos(plat)*sin(plong);
				ethez = -sin(plat);
				
				ephix = -sin(plong);
				ephiy = cos(plong);
				ephiz = 0.0;
				
				eradx = sin(plat)*cos(plong);
				erady = sin(plat)*sin(plong);
				eradz = cos(plat);
				
				vel(plat,plong,p);
				
				/*fprintf(stderr,"Testing vel: p->tvelzm: %lf\n",p->tvelzm);
				fprintf(stderr,"Testing vel: p->tvelxp: %lf\n",p->tvelxp);
				seems ok to here */
				
				vthetax[k][i] = ethez*p->tvelzm + ethey*p->tvelyp;
				vthetay[k][i] = ethex*p->tvelxm + ethez*p->tvelzp;
				vthetaz[k][i] = ethex*p->tvelxp + ethey*p->tvelym;
				
				vphix[k][i] = ephiy*p->tvelyp;
				vphiy[k][i] = ephix*p->tvelxm;
				vphiz[k][i] = ephix*p->tvelxp + ephiy*p->tvelym;
		}
	}
	//fprintf(stderr,"Testing torcoeff: vphiz: %lf\n",vphiz[1][1]);  -OK
	
	for (kj=1; kj<=180; kj++) {			//3
		plat = (90.0 - ((kj-0.5) - 90.0))*HC_PI/180.0;
		for (ij=1; ij <=360; ij++) {	//2
			plong = (ij-0.5)*HC_PI/180.0;
//			fprintf(stderr,"plat = %lf plong = %lf \n",plat,plong);
/*  kj conditions ###################################################	*/
			if ((kj==1)||(kj==180)) {
				if (kj == 1) {
					dvphixp = vphix[1][ij]/dtheta;
					dvphixm = -vphix[2][ij]/dtheta;
					dvphiyp = vphiy[1][ij]/dtheta;
					dvphiym = -vphiy[2][ij]/dtheta;
					dvphizp = vphiz[1][ij]/dtheta;
					dvphizm = -vphiz[2][ij]/dtheta;
				}
				if (kj == 180) {
					dvphixp = vphix[179][ij]/dtheta;
					dvphixm = -vphix[180][ij]/dtheta;
					dvphiyp = vphiy[179][ij]/dtheta;
					dvphiym = -vphiy[180][ij]/dtheta;
					dvphizp = vphiz[179][ij]/dtheta;
					dvphizm = -vphiz[180][ij]/dtheta;
				}
			}
			else {
					dvphixp = 0.5*vphix[kj-1][ij]/dtheta;
					dvphixm = -0.5*vphix[kj+1][ij]/dtheta;
					dvphiyp = 0.5*vphiy[kj-1][ij]/dtheta;
					dvphiym = -0.5*vphiy[kj+1][ij]/dtheta;
					dvphizp = 0.5*vphiz[kj-1][ij]/dtheta;
					dvphizm = -0.5*vphiz[kj+1][ij]/dtheta;
			}

/* end kj conditions #############################################
ij conditions #################################################### */		
			if ((ij==1)||(ij==360)) {
				if (ij == 1) {
					dvthetaxp = 0.5*(vthetax[kj][2])/dphi;
					dvthetayp = 0.5*(vthetay[kj][2])/dphi;
					dvthetazp = 0.5*(vthetaz[kj][2])/dphi;
					
					dvthetaxm = 0.5*(-vthetax[kj][360])/dphi;
					dvthetaym = 0.5*(vthetay[kj][360])/dphi;
					dvthetazm = 0.5*(vthetaz[kj][360])/dphi;
				}
				if (ij == 360) {
					dvthetaxp = 0.5*(vthetax[kj][1])/dphi;
					dvthetayp = 0.5*(vthetay[kj][1])/dphi;
					dvthetazp = 0.5*(vthetaz[kj][1])/dphi;
					
					dvthetaxm = 0.5*(-vthetax[kj][359])/dphi;
					dvthetaym = 0.5*(vthetay[kj][359])/dphi;
					dvthetazm = 0.5*(vthetaz[kj][359])/dphi;
				}
			}
			else {
					dvthetaxp = 0.5*(vthetax[kj][ij+1])/dphi;
					dvthetayp = 0.5*(vthetay[kj][ij+1])/dphi;
					dvthetazp = 0.5*(vthetaz[kj][ij+1])/dphi;
					
					dvthetaxm = 0.5*(-vthetax[kj][ij-1])/dphi;
					dvthetaym = 0.5*(vthetay[kj][ij-1])/dphi;
					dvthetazm = 0.5*(vthetaz[kj][ij-1])/dphi;
			}
			
/* end ij conditions ##############################################*/
			facx = 1.0;
			if ((kj ==1) || (kj == 180)) {
				facx = 0.75;
				}
			term1x = cos(plat)*vphix[kj][ij]*dtheta*dphi*facx;
			term1y = cos(plat)*vphiy[kj][ij]*dtheta*dphi*facx;
			term1z = cos(plat)*vphiz[kj][ij]*dtheta*dphi*facx;
			term2xp = sin(plat)*dvphixp*dphi*dtheta*facx;
			term2yp = sin(plat)*dvphiyp*dphi*dtheta*facx;
			term2zp = sin(plat)*dvphizp*dphi*dtheta*facx;
			term2xm = sin(plat)*dvphixm*dphi*dtheta*facx;
			term2ym = sin(plat)*dvphiym*dphi*dtheta*facx;
			term2zm = sin(plat)*dvphizm*dphi*dtheta*facx;
			term3xp = -dvthetaxp*dtheta*dphi*facx;
			term3yp = -dvthetayp*dtheta*dphi*facx;
			term3zp = -dvthetazp*dtheta*dphi*facx;
			term3xm = -dvthetaxm*dtheta*dphi*facx;
			term3ym = -dvthetaym*dtheta*dphi*facx;
			term3zm = -dvthetazm*dtheta*dphi*facx;

			theta = plat;
			avephi = plong;
			x = cos(theta);
			//fprintf(stderr,"Testing torcoeff: term3zm: %le\n",term3zm);  
			
			cnt0 = id[kj][ij];
			
			for (l=0; l <= p->nlm; l++) {		//33
				for (m=0; m <= l; m++) {
//					fprintf(stderr,"%i/%i  %i/%i  x:%lf avephi:%lf \n",kj,ij,l,m,x,avephi);  s'ok
					spharm(l,m,x,avephi,p);
					
					p->vlmx[cnt0][1][l][m] += (term1x*p->ylm1);
					p->vlmx[cnt0][2][l][m] += (term1x*p->ylm2);
					p->vlmy[cnt0][1][l][m] += term1y*p->ylm1;
					p->vlmy[cnt0][2][l][m] += term1y*p->ylm2;
					p->vlmz[cnt0][1][l][m] += term1z*p->ylm1;
					p->vlmz[cnt0][2][l][m] += term1z*p->ylm2;
					
					
					
					
/* kj conditions #################################################*/					
					if ((kj==1)||(kj==180)) {
						if (kj ==1) {
							cnt1 = id[1][ij];
							cnt2 = id[2][ij];
							p->vlmx[cnt1][1][l][m] += (term2xp*p->ylm1);
							p->vlmx[cnt2][1][l][m] += (term2xm*p->ylm1);
							p->vlmx[cnt1][2][l][m] += (term2xp*p->ylm2);
							p->vlmx[cnt2][2][l][m] += (term2xm*p->ylm2);
							
							p->vlmy[cnt1][1][l][m] += term2yp*p->ylm1;
							p->vlmy[cnt2][1][l][m] += term2ym*p->ylm1;
							p->vlmy[cnt1][2][l][m] += term2yp*p->ylm2;
							p->vlmy[cnt2][2][l][m] += term2ym*p->ylm2;
							
							p->vlmz[cnt1][1][l][m] += term2zp*p->ylm1;
							p->vlmz[cnt2][1][l][m] += term2zm*p->ylm1;
							p->vlmz[cnt1][2][l][m] += term2zp*p->ylm2;
							p->vlmz[cnt2][2][l][m] += term2zm*p->ylm2;
							
						}
						if (kj ==180) {
							cnt1 = id[179][ij];
							cnt2 = id[180][ij];
							p->vlmx[cnt1][1][l][m] += (term2xp*p->ylm1);
							p->vlmx[cnt2][1][l][m] += (term2xm*p->ylm1);
							p->vlmx[cnt1][2][l][m] += (term2xp*p->ylm2);
							p->vlmx[cnt2][2][l][m] += (term2xm*p->ylm2);
							
							p->vlmy[cnt1][1][l][m] += term2yp*p->ylm1;
							p->vlmy[cnt2][1][l][m] += term2ym*p->ylm1;
							p->vlmy[cnt1][2][l][m] += term2yp*p->ylm2;
							p->vlmy[cnt2][2][l][m] += term2ym*p->ylm2;
							
							p->vlmz[cnt1][1][l][m] += term2zp*p->ylm1;
							p->vlmz[cnt2][1][l][m] += term2zm*p->ylm1;
							p->vlmz[cnt1][2][l][m] += term2zp*p->ylm2;
							p->vlmz[cnt2][2][l][m] += term2zm*p->ylm2;

						}
					}
					else {
						
						cnt1 = id[kj-1][ij];
						cnt2 = id[kj+1][ij];
						
						p->vlmx[cnt1][1][l][m] += (term2xp*p->ylm1);
						p->vlmx[cnt2][1][l][m] += (term2xm*p->ylm1);	// this bad boy seg faults if the input file is not long enough  at kj/ij:179/352 cnt2:7 l:20 m:20 
						p->vlmx[cnt1][2][l][m] += (term2xp*p->ylm2);
						p->vlmx[cnt2][2][l][m] += (term2xm*p->ylm2);
						
						p->vlmy[cnt1][1][l][m] += term2yp*p->ylm1;
						p->vlmy[cnt2][1][l][m] += term2ym*p->ylm1;
						p->vlmy[cnt1][2][l][m] += term2yp*p->ylm2;
						p->vlmy[cnt2][2][l][m] += term2ym*p->ylm2;
						
						p->vlmz[cnt1][1][l][m] += term2zp*p->ylm1;
						p->vlmz[cnt2][1][l][m] += term2zm*p->ylm1;
						p->vlmz[cnt1][2][l][m] += term2zp*p->ylm2;
						p->vlmz[cnt2][2][l][m] += term2zm*p->ylm2;
					}
/* end kj conditions ###########################################*/
/*ij conditions ################################################*/					
					if ((ij==1)||(ij==360)) {
						if (ij==1) {
							cnt1 = id[kj][2];
							cnt2 = id[kj][360]; 
							//fprintf(stderr,"	2. (ij:%i)	cnt2:%i l:%i m:%i vlmx:%le term3xm:%le ylm1:%le \n",ij,cnt2,l,m,p->vlmx[cnt2][1][l][m],term3xm,p->ylm1);
							p->vlmx[cnt1][1][l][m] += (term3xp*p->ylm1);
							p->vlmx[cnt2][1][l][m] += (term3xm*p->ylm1);
							p->vlmx[cnt1][2][l][m] += (term3xp*p->ylm2);
							p->vlmx[cnt2][2][l][m] += (term3xm*p->ylm2);
						
							p->vlmy[cnt1][1][l][m] += term3yp*p->ylm1;
							p->vlmy[cnt2][1][l][m] += term3ym*p->ylm1;
							p->vlmy[cnt1][2][l][m] += term3yp*p->ylm2;
							p->vlmy[cnt2][2][l][m] += term3ym*p->ylm2;
						
							p->vlmz[cnt1][1][l][m] += term3zp*p->ylm1;
							p->vlmz[cnt2][1][l][m] += term3zm*p->ylm1;
							p->vlmz[cnt1][2][l][m] += term3zp*p->ylm2;
							p->vlmz[cnt2][2][l][m] += term3zm*p->ylm2;
						}
						if (ij == 360) {
							cnt1 = id[kj][1];
							cnt2 = id[kj][359];
							p->vlmx[cnt1][1][l][m] += (term3xp*p->ylm1);
							p->vlmx[cnt2][1][l][m] += (term3xm*p->ylm1);
							p->vlmx[cnt1][2][l][m] += (term3xp*p->ylm2);
							p->vlmx[cnt2][2][l][m] += (term3xm*p->ylm2);
						
							p->vlmy[cnt1][1][l][m] += term3yp*p->ylm1;
							p->vlmy[cnt2][1][l][m] += term3ym*p->ylm1;
							p->vlmy[cnt1][2][l][m] += term3yp*p->ylm2;
							p->vlmy[cnt2][2][l][m] += term3ym*p->ylm2;
						
							p->vlmz[cnt1][1][l][m] += term3zp*p->ylm1;
							p->vlmz[cnt2][1][l][m] += term3zm*p->ylm1;
							p->vlmz[cnt1][2][l][m] += term3zp*p->ylm2;
							p->vlmz[cnt2][2][l][m] += term3zm*p->ylm2;
						}
					}
					else {
							cnt1 = id[kj][ij+1];
							cnt2 = id[kj][ij-1];
							p->vlmx[cnt1][1][l][m] += (term3xp*p->ylm1);
							p->vlmx[cnt2][1][l][m] += (term3xm*p->ylm1);
							p->vlmx[cnt1][2][l][m] += (term3xp*p->ylm2);
							p->vlmx[cnt2][2][l][m] += (term3xm*p->ylm2);
						
							p->vlmy[cnt1][1][l][m] += term3yp*p->ylm1;
							p->vlmy[cnt2][1][l][m] += term3ym*p->ylm1;
							p->vlmy[cnt1][2][l][m] += term3yp*p->ylm2;
							p->vlmy[cnt2][2][l][m] += term3ym*p->ylm2;
						
							p->vlmz[cnt1][1][l][m] += term3zp*p->ylm1;
							p->vlmz[cnt2][1][l][m] += term3zm*p->ylm1;
							p->vlmz[cnt1][2][l][m] += term3zp*p->ylm2;
							p->vlmz[cnt2][2][l][m] += term3zm*p->ylm2;
					}
					
					/*if (((cnt0==7)||(cnt1=7))&&((l=15)&&(m==15))){
					fprintf(stderr,"	%i/%i  %i/%i %i ->vlmx=%lf  term1x=%le  ylm1=%le \n",kj,ij,l,m,cnt0,p->vlmx[7][1][l][m],term1x,p->ylm1);	// term1x & ylm1 are good at end but bad at start 
					}*/
/* end ij conditions #######################################################*/
				}
			} //33
		} //2
	} //3
	//fprintf(stderr,"	->vlmx=%le	vlmy=%le\n",p->vlmx[1][1][1][1],p->vlmy[1][1][3][3]);  //Not ok so far, test above
	fprintf(stderr,"Finished torcoeffs \n");

}

void vel(double plat, double plong, struct hc_ptrot *p) 
{
double tex,tey,tez,tpx,tpy,tpz;

tex = 1.0;
tey = 1.0;
tez = 1.0;
tpx = sin(plat)*cos(plong);
tpy = sin(plat)*sin(plong);
tpz = cos(plat);

p->tvelxm = -tey*tpz;
p->tvelym = -tez*tpx;
p->tvelzm = -tex*tpy;

p->tvelxp = tez*tpy;
p->tvelyp = tex*tpz;
p->tvelzp = tey*tpx;

}


void spharm(int l, int m, double x, double phi, struct hc_ptrot *p) 
{
	int lmm, lpm,jj,ii,kk,sign;
	double top,down,fact,arg,R,opi;

	lmm = l-m;
	lpm = l+m; 
	
	top = 1.0;
	down = 1.0;
	jj=0;
	for (ii=1; ii<=lmm; ii++) {
		//fprintf(stderr,"		-> top %lf jj %i \n",top,jj);
		jj++;
		top *= jj;
		//fprintf(stderr,"		---> top %lf jj %i \n",top,jj);
	}
	
	kk=0;
	for (ii=1; ii<=lpm; ii++) {
		kk++;
		down *= kk;
	}
	
	fact = top/down;
	arg = (2*l+1.0)*fact;
	R = sqrt(arg);
	//fprintf(stderr,"	lmm %i lpm %i top: %lf down:  %lf fact: %lf arg: %lf  \n",lmm,lpm,top,down,fact,arg);
	if (m > 0) {
		R = sqrt(2.0*arg);
	}
	sign = pow(-1,m);
	opi = 1/(4.0*HC_PI);
	plgndr(l,m,x,p);
	//fprintf(stderr,"sign: %i opi %lf R %lf plm: %lf m: %i phi %lf \n",sign,opi,R,p->plm,m,phi);
	p->ylm1 = sign*opi*R*p->plm*cos(m*phi);
	p->ylm2 = sign*opi*R*p->plm*sin(m*phi);
	//fprintf(stderr,"	... l: %i m: %i x: %lf plm: %lf ylm1: %lf ylm2: %lf \n",l,m,x,p->plm,p->ylm1,p->ylm2);
	
}

void plgndr(int l, int m, double x, struct hc_ptrot *p) {
	int i,ll,n;
	double somx2,fact,plgnd,pmmp1,pmm,pmm2,pll;
	
	if ((m < 0)||(m>l)||(fabs(x)>1)) {
		printf("Screwing argument in plgndr - ABORT!!\n");
		exit(-1);
	}
	
	pmm = 1.0;
	if (m > 0) {
		somx2 = sqrt((1.0-x)*(1.0+x));
		fact = 1.0;
		for (i=1; i<=m; i++) {
			pmm2 = pmm;
			pmm = (-pmm2)*fact*somx2;
			fact += 2;
		}
	}
		//fprintf(stderr,"	pmm %lf fact %lf somx %lf m %i \n",pmm,fact,somx2,m);
	if (l == m) {
			plgnd = pmm;
	}
	else {
			pmmp1 = x*(2*m+1)*pmm;
			n = m+1;
			if (l == n) {
				plgnd = pmmp1;
			}
			else {
				n = m+2;
				for (ll=n; ll<=l; ll++) {
					pll = (x*(2*ll - 1)*pmmp1 - (ll+m-1)*pmm)/(ll-m);
					pmm = pmmp1;
					pmmp1=pll;
				}
				plgnd = pll;
			}
  }
	p->plm = plgnd;
}
	
