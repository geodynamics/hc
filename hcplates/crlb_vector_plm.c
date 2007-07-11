#include <math.h>
#include "hcplates.h"

/* 
This subroutine computes the plm, the theta derivative of plm,
c        and plm/sin(theta). It uses the code plm.f to compute
c        fully normalized Legendre polynomials. It uses the
c        code polyfact.f to fix up the normalization factors.
c        The formula for the derivative of plm w.r.t. theta
c        is given by
c
c            d(plm)/dtheta=plm*L*cos(theta)/sin(theta)-pl(m-1)*(L+m)/sin(theta)
c
c        for associated Legendre functions. Hence, the need to take care of the
c        full normalization.
c
c    Theta is the co-latitude in radians, n is the maximum harmonic degree,
c        and p, dpdt, and pbyst are the spherical harmonic quantities,
c        with an indexing scheme given by
c
c            k=(l+1)*(l+2)/2-m
c
 */

/* Which are pointers? 
Assumes p, dpdt & pbyst are passed in as pointers
Size of pfact? */

void crlb_vecplm(struct hc_plates_params *plates, struct hc_plates_arrays *A, double theta, int n) 
{
	double pfact[2000];
	double stheta, ctheta;
	int l, m, k, k1;
	double pi;

	pi=M_PI;
	
	/*printf("	   Compute normalization factors, n=%d \n",n);
	printf("	   Check plates structure: plates->NPDIM=%d \n",plates->NPDIM);
	printf("	   Check A struct: A->y4in[1][1]=%lf\n",A->y4in[1][1]);*/
/* Compute normalization factors */
	crlb_polyfact(n,pfact);		//WORKING
	/*printf(" 	   Done polyfact \n");
	printf(" 	   Compute the plm, A->y4in[1][1]=%lf\n",A->y4in[1][1]);

	
	printf(" 	   Compute the plm, A->y4in[0][0]=%lf\n",A->y4in[1][1]);*/
/* Compute the plm */
	crlb_plm(theta,n,A);		// WORKING!

/* Do not allow sin(theta) to be zero since it is a divisor */
	if(theta < 0.001) {
		theta=0.001;
	}
	if ((M_PI - theta) < 0.001) {
		theta= M_PI - 0.001;
	}

	stheta=sin(theta);
	ctheta=cos(theta);
//printf(" 	   vecplm 2 \n");
/* Compute the derivative of the plm */
 for (l=0; l<=n; l++) {
	for (m=0; m<=l; m++) {
		k=(l+1)*(l+2)/2-m;
		k1=(l)*(l+1)/2-m;	/*ambiguous */
	//	printf(" 	   vecplm 3 \n");
		A->dpdt[k]=l*(A->p[k])*ctheta;

		/*Make sure the second term in the derivative formula exists before adding it in */
		if (l > m ) {
		//	printf(" 	      if  vecplm 4 \n");
			A->dpdt[k] += (-1)*(l+m)*(A->p[k1])*pfact[k]/pfact[k1];	
		}
		A->dpdt[k] /= stheta;	/*Think this is causing probs */
		
		/* Compute plm.sin(theta) */
		A->pbyst[k] = A->p[k]/stheta;
	//	printf(" 	   vecplm 5:%i/%i   k:%i  \n",l,m,k);

	}
  }
   // printf(" Finished crlb_vecplm \n");
}


/* 
c****This subroutine computes the normalization factors to convert
c        associated Legendre functions into fully normalized functions.
c        Note: These factors are not used by plm.f even though it
c        does compute fully normalized functions.
Assumes pfact is passed in as pointer, n is just by value
*/
	
void crlb_polyfact(int n, double *pfact) 
{	 
	int l, k, m;
/*	pfact[0]=0; */
	pfact[1]=1.0;
/*	printf("	      In polyfact, pfact[1]=%lf\n",pfact[1]); */
	for (l=1; l <= n ; l++) {
		k=(l+1)*(l+2)/2;
		pfact[k]=sqrt(2.0*(2.0*l+1.0));
		for (m=1; m <= l; m++) {
			k=(l+1)*(l+2)/2-m;
			pfact[k]=pfact[k+1]*sqrt(1.0/((l+m)*(l-m+1.0)));
//			printf("  %i/%i  %le\n",l,m,pfact[k]);
		}
	}
/*Fix up the m=0 terms here */
	for (l=1; l<=n; l++) {
		k=(l+1)*(l+2)/2;	
		pfact[k]=sqrt(2.0*l+1.0);
	}
/*	printf("	      In end polyfact, pfact[1]=%lf, k=%d\n",pfact[1],k); */
}

/* c
c****This subroutine efficiently computes fully normalized Legendre funtions
c        according to the normalization given in Stacey's book, Appendix C.
c ASsumes p is pointer, the rest passed in by value
theta is in radians already
*/

void crlb_plm(double theta, int n, struct hc_plates_arrays *A)
{
	int k,l,m,kold;
	double stheta, ctheta,pi,top,den;
	
	kold=1;
	pi = M_PI;	
	if(theta < 0.001) {
		theta=0.001;
	}
	if ((M_PI - theta) < 0.001) {
		theta= M_PI - 0.001;
	}

	stheta=sin(theta);  //*pi/180 ??
	ctheta=cos(theta);
	
	//printf("  stheta %lf ctheta %lf\n",stheta,ctheta);

	A->p[0]=0.0;
	A->p[1]=sqrt(2.0);
	for (l=1; l<=n; l++) {
		k=(l+1)*(l+2)/2-l;
		A->p[k]=stheta*sqrt((2.0*l+1.0)/(2.0*l))*(A->p[kold]);
		kold=k;	
		for (m=l-1;m>=0;m--) {
			k=(l+1)*(l+2)/2-m;
			den=sqrt((l+m+1.0)*(l-m));
			top=sqrt((l-m-1.0)*(l+m+2.0));
			A->p[k]=2.0*(m+1.0)*ctheta/(stheta*den)*(A->p[k-1])-(top/den)*(A->p[k-2]); 
			if (m==0) {
				A->p[k] /= sqrt(2.0);
			}
	//		printf("   %i/%i %le \n",l,m,A->p[k]);
		}
	}
	A->p[1]=1.0;
}
