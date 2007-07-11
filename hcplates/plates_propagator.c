#include <math.h>
#include "hcplates.h"

/* Subroutines to calculate green functions */

/* 2x2 toroidal propagator subroutine */


/* 
poload: subroutine for internal loads (poloidal flow)
Calls prop4, mamlt4, returns ua & uc	tested - working */

void poload(struct hc_plates_params *plates,int L,int n,int iba,int nrb,double rb, double visc0,double *ua,double *uc)
{

int i,j,layr,nrad;
double pba[5][5],pca[5][5],p1[5][5],p2[5][5];
double rnorm,fact1,fact2,bsigma,tmp1,tmp2,tmp3,tmp4;

/* Want all arrays fotran style 1->4, so initialized to size 5 */

//fprintf(stderr,"POLOAD: n:%i iba:%i nrb:%i rb:%lf visc0:%le\n",n,iba,nrb,rb,visc0);

/* zero arrays */
	for(i=1;i<=4;i++) {
		for (j=1;j<=4;j++) {
			pba[i][j]=0.0;
			pca[i][j]=0.0;
			p1[i][j]=0.0;
			p2[i][j]=0.0;
		}
	}
/* Calculate pca
Initialize pca to an identity matrix */
	for (i=1;i<=4;i++) {
		pca[i][i]=1.0;
	}

/* Multiply pca by each matrix layer */

	for (layr=1;layr<=n;layr++) {
		rnorm=plates->r[layr+1]/plates->r[layr];
		prop4(L,rnorm,plates->visc[layr],p1);	
	/*	fprintf(stderr,"	POLOAD, L:%i   p1: %le %le %le %le \n",L,p1[1][1],p1[2][1],p1[3][1],p1[4][1]);  //
		fprintf(stderr,"	POLOAD, L:%i   p2: %le %le %le %le \n",L,p2[1][1],p2[2][1],p2[3][1],p2[4][1]);  //
		fprintf(stderr,"	 POLOAD, L:%i   pca: %le %le %le %le \n",L,pca[1][1],pca[2][1],pca[3][1],pca[4][1]);   */
		mamlt4(p1,pca,p2);
		//fprintf(stderr,"	  POLOAD, L:%i   p2: %le %le %le %le \n",L,p2[1][1],p2[2][1],p2[3][1],p2[4][1]);  //ok till here
		for (i=1;i<=4;i++) {
			for (j=1; j<=4; j++) {
				pca[i][j]=p2[i][j];
				//fprintf(stderr,"		%i/%i pca=%le \n",i,j,pca[i][j]);
			}
		}
	}
/* Find layer number of load */

	for (i=1; i<=4; i++) {
		if ( rb > plates->r[i] && rb <= plates->r[i+1] ) { /* O'Neill: this is only called for nrb=1, so setting rb to just double */
			nrad=i;
		}
	}
/*printf("  Found layer number of load; rb =  %lf, nrad=%d \n",rb,nrad); */
/* Calculate pba
Initialize pba to identity matrix */
	for (i=1;i<=4;i++) {
		pba[i][i]=1.0;
	}
/*Mulitply pba by each layer matrix */ /* for soe reason memory is corrupted here - the layr countered is everywhere... */

	for (layr=nrad; layr<=n; layr++) {
		rnorm=plates->r[layr+1]/plates->r[layr];
		if (layr == nrad) {
			rnorm=plates->r[layr+1]/rb;
		}
		prop4(L,rnorm,plates->visc[layr],p1);
		mamlt4(p1,pba,p2);
		for (i=1;i<=4;i++) {
			for (j=1; j<=4; j++) {
				pba[i][j]=p2[i][j];
				
			}
		}
		//fprintf(stderr,"	  POLOAD, L:%i   p2: %le %le %le %le \n",L,p2[1][1],p2[2][1],p2[3][1],p2[4][1]); //good till here
	}
/* Solve for the unknown boundary vectors
Branch according to surface boundary condition 
iba=1 -> free slip; iba=0 -> no-slip */
	if (iba == 0) {
		/* no slip */
		fact1=pca[1][2]/pca[1][3] - pca[2][2]/pca[2][3];
		fact2=pba[1][3]/pca[1][3] - pba[2][3]/pca[2][3];
		//fprintf(stderr,"	POLOAD: visc0= %le\n",visc0);
		bsigma=rb/visc0;
		/* Get CMB solution */
		uc[1]=0.0;
		uc[2]= -bsigma*fact2/fact1;
		uc[3]=(-uc[2]*pca[1][2] - bsigma*pba[1][3])/pca[1][3];
		uc[4]=0.0;
		/*Get surface solution */
		ua[1]=0.0;
		ua[2]=0.0;
		tmp1 = pca[3][2]*uc[2];
		tmp2 = pca[3][3]*uc[3];
		tmp3 =  pba[3][3]*bsigma;
		tmp4 = tmp1+tmp2+tmp3;
		ua[3] = pca[3][2]*uc[2] + pca[3][3]*uc[3] + pba[3][3]*bsigma;  //Note - small rounding error c/w fortran version for v small ua[3]
	/*	printf("	  POLOAD: pca32=%le pca33=%le upba33=%le  \n",pca[3][2],pca[3][3],pba[3][3]); //- working
		printf("	  POLOAD: uc2=%le uc3=%le bsigma=%le  ua3=%le \n",uc[2],uc[3],bsigma,ua[3]);
		printf("	  POLOAD: tmps = %le %le %le %le \n",tmp1,tmp2,tmp3,tmp4); */
		ua[4]=pca[4][2]*uc[2] + pca[4][3]*uc[3] + pba[4][3]*bsigma;
	}
	else {
		/* Free slip */
		fact1=pca[1][2]/pca[1][3] - pca[4][2]/pca[4][3];
		fact2=pba[1][3]/pca[1][3] - pba[4][3]/pca[4][3];
		bsigma=rb/visc0;
		/* Get CMB solution */
		uc[1]=0.0;
		uc[2]= -bsigma*fact2/fact1;
		uc[3]=(-uc[2]*pca[1][2] - bsigma*pba[1][3])/pca[1][3];
		uc[4]=0.0;
		/* Get surface solution */
		ua[1]=0.0;
      		ua[2]=pca[2][2]*uc[2]+pca[2][3]*uc[3]+pba[2][3]*bsigma;		//different - slightly...
     		ua[3]=pca[3][2]*uc[2]+pca[3][3]*uc[3]+pba[3][3]*bsigma;
      		ua[4]=0.0;
	}
	//printf("	  POLOAD: uc2=%le uc3=%le ua3=%le ua4=%le \n",uc[2],uc[3],ua[3],ua[4]); //- working


}

/* subroutine for imposed surface poloidal velocity field */

void polvel(struct hc_plates_params *plates, int L, int n,double  *ua, double *uc)
{
	int i,j,layr;
	double pca[5][5],p1[5][5],p2[5][5];
	double rnorm;
	
	/* Zero arrays */
	for (i=1;i<=4;i++) {
		for (j=1;j<=4;j++) {
			pca[i][j]=0.0;
			p1[i][j]=0.0;
			p2[i][j]=0.0;
		}
	}
	/* Calculate pca
	Initialize pca to an identity matrix */
	for (i=1;i<=4;i++) {	
		pca[i][i]=1.0;
	}
	/*Multiply pca by each layer matrix */
	for (layr=1;layr<=n;layr++) {
		rnorm=plates->r[layr+1]/plates->r[layr];
	/*	printf("	Prop4 input: layr=%d L=%d rnorm=%lf visc=%lf pca11=%lf\n",layr,L,rnorm,plates->visc[layr],pca[1][1]); */
		prop4(L,rnorm,plates->visc[layr],p1);
		mamlt4(p1,pca,p2);
		for (i=1;i<=4;i++) {
			for (j=1;j<=4;j++) {
			/*	printf("	   i=%d j=%d p1=%lf pca=%lf p2=%lf\n",i,j,p1[i][j],pca[i][j],p2[i][j]);*/
				pca[i][j]=p2[i][j];
			}
		}
	}
	/* Get CMB solution */
/*	printf("		Am I dividing by zero: pca22=%lf pca23=%lf pca12=%lf pca13=%lf\n",pca[2][2],pca[2][3],pca[1][2],pca[1][3]); */
	uc[1]=0.0;
	uc[2]=1.0/(pca[2][2]-pca[2][3]*pca[1][2]/pca[1][3]);
	uc[3]=-uc[2]*pca[1][2]/pca[1][3];
	uc[4]=0.0;
	/*Get surface solution */
	ua[1]=0.0;
	ua[2]=1.0;
	ua[3]=pca[3][2]*uc[2]+pca[3][3]*uc[3];
	ua[4]=pca[4][2]*uc[2]+pca[4][3]*uc[3];		/*something wrong with ua[4]=nan */

/*	printf("		uc2=%lf uc3=%lf ua3=%lf ua4=%lf pca42=%lf pca43=%lf\n",uc[2],uc[3],ua[3],ua[4],pca[4][2],pca[4][3]);*/
}

/* Subroutine for imposed toroidal velocity field */
void torvel(struct hc_plates_params *plates, int L, int n, double *va,double *vc)
{
	double pca[3][3],p1[3][3],p2[3][3];
	int i,j,layr;
	double rnorm;

	/* Zero arrays */
	for (i=1;i<=2;i++) {
		for (j=1;j<=2;j++) {
			pca[i][j]=0.0;
			p1[i][j]=0.0;
			p2[i][j]=0.0;
		}
	}
	/* Calculate pca
	Initialize pca to indentity matrix */
	for (i=1;i<=2;i++) {
		pca[i][i]=1.0;
	}
	/* Mulitply pca by each layer matrix */
	for (layr=1;layr<=n;layr++) {
		rnorm=plates->r[layr+1]/plates->r[layr];
		prop2(L,rnorm,plates->visc[layr],p1);
		mamlt2(p1,pca,p2);
		for (i=1;i<=2;i++) {
			for (j=1;j<=2;j++) {	
				pca[i][j]=p2[i][j];
			}
		}
	}
	/* Get CMB solution */	
	vc[1]=1.0/pca[1][1];
	vc[2]=0.0;
	/* Get surface solution */
	va[1]=1.0;
	va[2]=vc[1]*pca[2][1];
}

/* Matrix multiply routine, where matrices are 4x4, but the indices are from 1:4 
(as opposed to the C convention 0:3) */

void mamlt4(double A[5][5], double B[5][5], double C[5][5])
{
	int i,j,k;
	double sum1,sum2;
	for (i=1;i<=4;i++) {
		for (j=1;j<=4;j++) {
			sum2=0.0;
			for (k=1;k<=4;k++) {
				sum1 = A[i][k]*B[k][j];
				sum2 += sum1; 
				C[i][j]=sum2;
			}
		}
	}
}

/* Matrix multiply routine, where matrices are 2x2, but the indices are from 1:2 
(as opposed to the C convention 0:1) */
	
void mamlt2(double A[3][3], double B[3][3], double C[3][3])
{
	int i,j,k;
	double sum1,sum2;
	for (i=1;i<=2;i++) {
		for (j=1;j<=2;j++) {
			sum2=0;
			for (k=1;k<=2;k++) {
				sum1 = A[i][k]*B[k][j];
				sum2 += sum1;
				C[i][j] = sum2;
			}
		}
	}
}

void prop2(int l, double r,double visc, double p[3][3])
{
	int LL,rLL;
	double alpha,rnu,x1,x2,x3,ss,cc;

	LL=l*(l+1);
	rLL=LL;
	alpha=sqrt(1 + 4*rLL);
	rnu=log(r);
	x2=exp(-0.5*rnu);
	x1=x2/alpha;
	x2=exp(0.5*rnu*alpha);
	x3=exp(-0.5*rnu*alpha);
	ss=(x2-x3)/2.0;
	cc=(x2+x3)/2.0;
	p[1][1]=x1*(3.0*ss + alpha*cc);
	p[1][2]=x1*2.0*ss/visc;
	p[2][1]=x1*2.0*(rLL-2.0)*visc*ss;
	p[2][2]=x1*(-3.0*ss + alpha*cc);
}

/* 4x4 poloidal propagator subroutine. tested - works */
void prop4(int l, double r, double v, double p[5][5]) 
{
 int lp1,lm1,lp2,lm2,lp3,lm3,l2,i,j;
 int rl,rlp1,rlm1,rlp2,rlm2,rlp3,rlm3,rl2;
 double c[5][5],rlambda[5][5],b[5][5],x[5][5];
 
 lp1=l+1;
 lm1=l-1;
 lp2=l+2;
 lm2=l-2;
 lp3=l+3;
 lm3=l-3;
 l2=l*l;
 

 rl=l;
 rlp1=l+1;
 rlm1=l-1;
 rlp2=l+2;
 rlm2=l-2;
 rlp3=l+3;
 rlm3=l-3;
 rl2=l*l;

/* Calculate c matrix components */

 c[1][1]=l*lp1;
 c[1][2]=l*lp1;
 c[1][3]=l*lp1;
 c[1][4]=l*lp1;
 c[2][1]=lp3;
 c[2][2]=lp1;
 c[2][3]=-lm2;
 c[2][4]=-l;
 c[3][1]=2*v*lp1*(l2-lp3);
 c[3][2]=2*v*l*lm1*lp1;
 c[3][3]=-2*v*l*(l2+3*l-1);
 c[3][4]=-2*v*l*lp2*lp1;
 c[4][1]=2*v*l*lp2;
 c[4][2]=2*v*lm1*lp1;
 c[4][3]=2*v*lm1*lp1;
 c[4][4]=2*v*l*(lp2);

/* Calculate lambda a=matrix */
	for (i=1;i<=4;i++) {
		for (j=1;j<=4;j++) {
			rlambda[i][j]=0.0;
		}
	}
	rlambda[1][1]=pow(r,lp1);
	rlambda[2][2]=pow(r,lm1);
	rlambda[3][3]=pow(r,(-l));
	rlambda[4][4]=pow(r,(-lp2));
/* Calculate b(=cinv) matrix */
	b[1][1]=-rlp2;
	b[1][2]=rl*rlp2;
	b[1][3]=-1.0/(2.0*v);
	b[1][4]=rl/(2.0*v);
	b[2][1]=(rl2+3*rl-1.0)/rlp1;
      	b[2][2]=-rlm1*rlp1;
      	b[2][3]=1.0/(2.0*v);
      	b[2][4]=-rlm2/(2.0*v);
      	b[3][1]=rlm1;
      	b[3][2]=rlm1*rlp1;
      	b[3][3]=-1.0/(2.0*v);
      	b[3][4]=-rlp1/(2.0*v);
      	b[4][1]=(3.0+rl-rl2)/rl;
      	b[4][2]=-rl*rlp2;
      	b[4][3]=1.0/(2.0*v);
      	b[4][4]=rlp3/(2.0*v);
	
	for (i=1;i<=4;i++) {
		b[1][i]=b[1][i]/(2.0*rl + 3.0);
		b[4][i]=b[4][i]/(2.0*rl + 3.0);
	}
	for (i=1;i<=4;i++) {
		b[2][i]=b[2][i]/(2.0*rl - 1.0);
		b[3][i]=b[3][i]/(2.0*rl - 1.0);
	}
	for (i=1;i<=4;i++) {
		for (j=1;j<=4;j++) {
			b[i][j]=b[i][j]/(2.0*rl + 1.0);
		}
	}
	/* Get product of c*rlambda */
	mamlt4(c,rlambda, x);
	/* Get product of (c*rlambda)*b=p */
	mamlt4(x, b, p);

}



