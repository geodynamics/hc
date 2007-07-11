#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "hcplates.h"

/* Routine to integrate the torques on the plates due to 1) internal loads, and
2)  externally imposed plate motions

loop on latitude, longitude coordinates: read plate index number (from the file map.plate - assumes nlat=0:180 & NLONG=0:360, if not, better mkae these loops and map.plate match up...), compute
coordinate angles, Legendre polynomials, and spherical harmonic functions
For the record, at the moment dd=1, nlat=180/dd=180, and nlong=360/dd=360.

Returns arrays fin & fex, as well as srta, srpa, trta, trpa & trra.

 */

void integrate_torques(struct hc_plates_params *plates, struct hc_plates_arrays *A, char platesmap_filename[])
{
  FILE *filePtr;
  double theta,dtheta,dd,dphi,phi,ctheta,stheta,rtheta,pi;
  double cphi,sphi,rphi;
  double e1,e2,e3,e4,e5;
  double trt,trp,trr,trrc,srt[plates->NPDIM],srp[plates->NPDIM];
  double cc,ss,dtylm1,dtylm2,dpylm1,dpylm2,pk1,pk2;
  double area, torp[plates->NLAT][102], cm[101],sm[101];
  double toss1,toss2;
  int ilat,ilong,kkmax,kk,lineNumber,cursorPosition,id,m,L,k,ind1,ii,iplt,idir,ind2,ilon;
  char ch,inputLine[500];

  pi=M_PI;

  /* Open plates_map file for reading */
  filePtr = fopen(platesmap_filename, "r");
  if (filePtr == NULL)
	{ printf("Error opening file. \n"); }
  else
	{
	printf("	Opened platesmap file ok\n");
	/* Initialize these degree partition counters */
	dd=1.0; /* 1x1 degree mesh */
	dtheta=dd*pi/180.0;	
	dphi=dtheta;
	theta=180.5*dtheta;
	printf("	Pre lat loop: NLAT=%d, theta=%lf dtheta=%lf \n",plates->NLAT,theta,dtheta);
	for (ilat=1; ilat <= plates->NLAT; ilat++) {
		theta -= dtheta;	/* decrements on each lat loop */
		printf("	Integrating - looping though lats: ilat=%d theta=%lf %i \n",ilat,theta,plates->NPLT);
		ctheta=cos(theta);
		stheta=sin(theta);  //already in rad
		phi=-0.5*dphi;
		rtheta= ilat - 90.5;

		/*  note that a factor of erad*erad*erad has been left out of 'area', since it ultimately
		c***     multiplies both sides of the final matrix equation */
		area=stheta*dtheta*dphi;

		/* 	c*** compute Legendre functions and derivatives */
		//printf("	Calling crlb_vecplm, theta=%lf Lmax=%d A->p[0]=%lf\n",theta,plates->Lmax,A->p[0]);
		crlb_vecplm(plates, A,theta,plates->Lmax); /*WORKING */
		//printf("torpy: \n");
		
		kkmax=(plates->Lmax+1)*(plates->Lmax+2)/2;
		
		//printf(" Lmax= %i  kkmax = %i ilat=%i A->p[kk]=%lf torp=%lf \n",plates->Lmax,kkmax,ilat,A->p[kkmax]);
		
		for (kk=1; kk<=kkmax; kk++) {
		   // printf("  ilat=%i, kk=%i p[kk]=%lf torp=%lf \n",ilat,kk,A->p[kk],A->torp[ilat][kk],A->torp[ilat][1]);
			A->torp[ilat][kk]=A->p[kk];
		}
		//printf("torpy -done:\n");
		for (ilon=1; ilon <= plates->NLONG; ilon++) {
			/* read plate index number */
			//lineNumber=0;	
			//cursorPosition=0;
			//ch=fgetc(filePtr);
			/*while (ch!='\n' || ch!='^M') {
				if (cursorPosition >= 500) {
					printf("Line in platesmap file unexpectedly long. Ciao.\n");
					exit(0);
				}
				inputLine[cursorPosition] = ch;
				cursorPosition++;
				ch=fgetc(filePtr); 
			}*/
			//printf("Reading line:\n");
			//fscanf(filePtr,"%s",&inputLine);
			//fscanf(filePtr,"%i %lf %lf",&id,&toss1,&toss2);
			//printf("Reading line1: (line): %s",inputLine);
			fgets(inputLine, 100, filePtr);
			id = atoi(inputLine); /*atoi "should" read the first integer, and ignore all else */
			//printf("Reading line: id=%i  \n",id); 
			//printf("Reading line: id=%i, %lf %lf \n",id,toss1,toss2); 
			A->idp[ilat][ilon]=id;
		 	phi += dphi;
			cphi = cos(phi);
			sphi = sin(phi);
			rphi = ilon - 0.5;	
			//fprintf(stderr,"id:%i  phi:%le cphi:%le sphi:%le ctheta:%le stheta:%le area:%le\n",id,phi,cphi,sphi,ctheta,stheta,area); 
			/* Compute coordinate conversion factors (r, theta,phi) --> (x,y,z) */
			e1 = area*ctheta*cphi;
			e2 = area*sphi;
			e3 = area*ctheta*sphi;
			e4 = -area*cphi;
			e5 = -area*stheta;
//			printf("   e1=%lf e2=%le e3=%le e4=%le e5=%le\n",e1,e2,e3,e4,e5); //seem ok
			/*  compute and store cos(m*phi) and sin(m*phii) functions using simple
			c***     trigonometric recurrence relation */
			cm[1]=1.0;
			sm[1]=0.0;
			cm[2]=cphi;
			sm[2]=sphi;
			for (m = 2; m <= plates->Lmax; m++) {
				cm[m+1] = 2.0*cm[m]*cphi - cm[m-1]; 
				sm[m+1] = 2.0*sm[m]*cphi - sm[m-1];
			}
			/* loop on L,m: compute spherical harmonics and accumulate torques */
			trt = 0.0;
			trp = 0.0;
			trr = 0.0;
			trrc = 0.0;
			for (ii=1; ii <= plates->NPDIM; ii++) {
				A->srt[ii]=0.0;
				A->srp[ii]=0.0;
			}
			for (L=0; L <= plates->Lmax; L++) {
				for (m=0; m<=L; m++) {
					k=(L+2)*(L+1)/2-m;
					cc = cm[m+1];
					ss = sm[m+1];
					dtylm1 = A->dpdt[k]*cc;		/* Problem is that at L=1,m=0,k=3, dpdt=inf */
				 	dtylm2 = A->dpdt[k]*ss;
// WORKING					printf("		Testing dty's: L=%d, m=%d, k=%d, cc=%lf ss=%lf dtylm1=%lf A->dpdt=%lf\n",L,m,k,cc,ss,dtylm1,A->dpdt[k]); 
					dpylm1 = -m*(A->pbyst[k])*ss;
					dpylm2 = m*(A->pbyst[k])*cc;
					pk1 = A->p[k]*cc;
					pk2 = A->p[k]*ss;
				/*	printf("		dts: dtylm1=%lf dtylm2=%lf dpylm1=%lf dpylm2=%lf\n");*/
					/* Accumulate shear stresses due to internal load */
					if (L < plates->Lload) {
					/*	printf("		Testing const. of trt: trt=%lf y4in=%lf,dtylm1=%lf y4in2=%lf dtylm2=%lf\n",trt,A->y4in[k][1],dtylm1,A->y4in[k][2],dtylm2); */
							trt += A->y4in[k][1]*dtylm1 + A->y4in[k][2]*dtylm2;
							trp += A->y4in[k][1]*dpylm1 + A->y4in[k][2]*dpylm2;
	      					trr += A->y3in[k][1]*pk1 + A->y3in[k][2]*pk2;
	      					trrc += A->z3in[k][1]*pk1 + A->z3in[k][2]*pk2;
// Seem ok 							printf("  tr_: %le %le %le %le\n",trt,trp,trr,trrc);
					}
					/* loop on imposed plate rotation direction */
					if (L <= plates->Lplt) {
						for (iplt=1; iplt <= plates->NPLT; iplt++) {
							for (idir=1; idir <=3; idir++) {
								ii = 3*(iplt-1) + idir;
								/* Accumulate shear stresses due to external imposed plate motions  - ROUNDING DIFFERENCES BW FORTRAN & C HERE*/
								A->srt[ii] += A->y4ex[ii][k][1]*dtylm1 + A->y4ex[ii][k][2]*dtylm2;
								A->srt[ii] += A->y10ex[ii][k][1]*dpylm1 + A->y10ex[ii][k][2]*dpylm2;
								A->srp[ii] += A->y4ex[ii][k][1]*dpylm1 + A->y4ex[ii][k][2]*dpylm2;
								A->srp[ii] += (-1)*A->y10ex[ii][k][1]*dtylm1 - A->y10ex[ii][k][2]*dtylm2;
		//						printf(" ii=%d k=%d y10ex1=%lf y10ex2=%lf y4ex1=%lf y4ex2=%lf\n",ii,k,A->y10ex[ii][k][1],A->y10ex[ii][k][2],A->y4ex[ii][k][1],A->y4ex[ii][k][2]);
								//printf("E: srt=%le srp=%le\n",A->srt[1],A->srp[1]);
							}
						}
					}
				}
			}
//			printf("E: srt=%le srp=%le\n",A->srt[1],A->srp[1]);
			/* Reassign driving and resisting tractions for net calculation */
			for (iplt=1; iplt <= plates->NPLT; iplt++) {
				for (idir=1; idir <=3; idir++) {
					ii = 3*(iplt-1) + idir;
					A->srta[ilat][ilon][ii] = A->srt[ii];
					A->srpa[ilat][ilon][ii] = A->srp[ii];
				}
			}
			A->trta[ilat][ilon] = trt;
			A->trpa[ilat][ilon] = trp;
			A->trra[ilat][ilon] = trr;
			/* Write output out.normstress here? */
	
			/* Finish reassigning and writing out */
			ind1 = 3*(id-1);
			/* accumulate torques: fin, fex should already be zero'd by here, so += as we go */
			A->fin[ind1+1] += trp*e1 + trt*e2;
			A->fin[ind1+2] += trp*e3 + trt*e4;
			A->fin[ind1+3] += trp*e5;
			//printf("%i  FIN123: [%le] [%le] [%le]\n",ilon,A->fin[ind1+1],A->fin[ind1+2],A->fin[ind1+3]);
		/*	printf("		Testing fin: ilat=%d ilon=%d %lf %lf %lf id=%d ind1=%d \n",ilat,ilon,A->fin[ind1+1],A->fin[ind1+2],A->fin[ind1+3],id,ind1);*/
			for (iplt=1; iplt<=plates->NPLT; iplt++) {
				for (idir=1; idir <= 3; idir++) {
					ind2 = 3*(iplt-1) + idir;
					A->fex[ind2][ind1+1] += A->srp[ind2]*e1 + A->srt[ind2]*e2;
					A->fex[ind2][ind1+2] += A->srp[ind2]*e3 + A->srt[ind2]*e4;
					A->fex[ind2][ind1+3] += A->srp[ind2]*e5;
//					printf("%i	Integrating: fex[%d][%d]=%le fex[%d][%d]=%le fex[%d][%d]=%le     \n",ilat,ind2,ind1+1,A->fex[ind2][ind1+1],ind2,ind1+2,A->fex[ind2][ind1+2],ind2,ind1+3,A->fex[ind2][ind1+3]);
//					printf("	srp[%d]=%lf srt[%d]=%lf \n",ind2,A->srp[ind2],ind2,A->srt[ind2]);
				}
			}
		}
	}

  }
fclose(filePtr);

}

