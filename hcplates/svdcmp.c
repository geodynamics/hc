#include <math.h>
#define NRANSI
#include "nrutil.h"
#include "hcplates.h"

double pythag( double side_a, double side_b )
{
   return( sqrt( (side_a*side_a) + (side_b*side_b ))) ;
}


void svdcmp(struct hc_plates_params *plates, struct hc_plates_arrays *A, int m, int n)
{
	int flag,i,ii,its,j,jj,k,l,nm;
	double anorm,c,f,g,h,s,scale,x,y,z, *rv1;

	rv1 = (double *) c_alloc1d(n+1); /*Complains abour cast */
	printf("	   Allocated rv1 ok\n");
	g=scale=anorm=0.0;
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) {
				scale += fabs(A->fexs[k][i]);
			}
			if (scale != 0.0) {
				for (k=i;k<=m;k++) {
					A->fexs[k][i] /= scale;
					s += A->fexs[k][i]*A->fexs[k][i];
				}
				f=A->fexs[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				A->fexs[i][i]=f-g;
				for (j=l;j<=n;j++) {
					s=0.0;
					for (k=i;k<=m;k++) {
						if (A->fexs[k][i] && A->fexs[k][j]) {
							s += A->fexs[k][i]*A->fexs[k][j];
						}
						else {
							printf("	bad fexs here\n");
						}
					}
					f=s/h;
					for (k=i;k<=m;k++) {
						A->fexs[k][j] += f*A->fexs[k][i];
					}
				}
				for (k=i;k<=m;k++) {
					A->fexs[k][i] *= scale;
				}
			}
		}
		A->WW[i] = scale*g;
		g=s=scale=0.0;
		if ((i <= m) && (i != n)) {
			for (k=l;k<=n;k++) {	
				scale += fabs(A->fexs[i][k]);
			}
			if (scale != 0.0) {
				for (k=l;k<=n;k++) {
					A->fexs[i][k] /= scale;
					s += A->fexs[i][k]*A->fexs[i][k];
				}
				f=A->fexs[i][l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				A->fexs[i][l]=f-g;
				for (k=l;k<=n;k++) {
					rv1[k]=A->fexs[i][k]/h;
				}
				for (j=l;j<=m;j++) {
					s=0.0;
					for (k=l;k<=n;k++) {
						s += A->fexs[j][k]*A->fexs[i][k];
					}
					for (k=l;k<=n;k++) {
						A->fexs[j][k] += s*rv1[k];
					}
				}
				for (k=l;k<=n;k++) {
					A->fexs[i][k] *= scale;
				}
			}
		}
		anorm=FMAX(anorm,(fabs(A->WW[i])+fabs(rv1[i])));
	}
	for (i=n;i>=1;i--) {
		if (i < n) {
			if (g != 0.0) {
				for (j=l;j<=n;j++){
					A->VV[j][i]=(A->fexs[i][j]/A->fexs[i][l])/g;
				}
				for (j=l;j<=n;j++) {
					s=0.0;
					for (k=l;k<=n;k++) {
						s += A->fexs[i][k]*A->VV[k][j];
					}
					for (k=l;k<=n;k++) {
						A->VV[k][j] += s*A->VV[k][i];
					}
				}
			}
			for (j=l;j<=n;j++) {
				A->VV[i][j]=A->VV[j][i]=0.0;
			}
		}
		A->VV[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	ii = IMIN(m,n);
	for (i=ii;i>=1;i--) {
		l=i+1;
		g=A->WW[i];
		for (j=l;j<=n;j++) {
			A->fexs[i][j]=0.0;
		}
		if (g !=0) {
			g=1.0/g;
			for (j=l;j<=n;j++) {
				s=0.0;
				for (k=l;k<=m;k++) {
					s += A->fexs[k][i]*A->fexs[k][j];
				}
				f=(s/A->fexs[i][i])*g;
				for (k=i;k<=m;k++) {
					A->fexs[k][j] += f*A->fexs[k][i];
				}
			}
			for (j=i;j<=m;j++) A->fexs[j][i] *= g;
		} else for (j=i;j<=m;j++) A->fexs[j][i]=0.0;
		++(A->fexs[i][i]);
	}
	for (k=n;k>=1;k--) {
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=1;l--) {
				nm=l-1;
				if ((float)(fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((float)(fabs(A->WW[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((float)(fabs(f)+anorm) == anorm) break;		//is this working ok?
					g=A->WW[i];
					h=pythag(f,g);
					A->WW[i]=h;		//ici?
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=1;j<=m;j++) {
						y=A->fexs[j][nm];
						z=A->fexs[j][i];
						A->fexs[j][nm]=y*c+z*s;
						A->fexs[j][i]=z*c-y*s;
					}
				}
			}
			flag = 0;
			z=A->WW[k];
			if (l == k) {
				if (z < 0.0) {
					A->WW[k] = -z;
					for (j=1;j<=n;j++) A->VV[j][k] = -A->VV[j][k];
				}
				flag = 1;
				break;  /// This should break out completely!
			}
			if (its == 30) {
				printf("no convergence in 30 svdcmp iterations");
				return;
			}	
			if (flag < 1) {
			x=A->WW[l];
			nm=k-1;
			y=A->WW[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=A->WW[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=1;jj<=n;jj++) {
					x=A->VV[jj][j];
					z=A->VV[jj][i];
					A->VV[jj][j]=x*c+z*s;
					A->VV[jj][i]=z*c-x*s;
				}
				z=pythag(f,h);
				A->WW[j]=z;
				if (z != 0.0) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=1;jj<=m;jj++) {
					y=A->fexs[jj][j];
					z=A->fexs[jj][i];
					A->fexs[jj][j]=y*c+z*s;
					A->fexs[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			A->WW[k]=x;
			}
		}
	}
	Free1D(rv1);
}
#undef NRANSI
