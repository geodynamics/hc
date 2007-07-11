#define NRANSI
#include "nrutil.h"
#include "hcplates.h"

void svbksb(struct hc_plates_params *plates, struct hc_plates_arrays *A, int m, int n)
{
	int jj,j,i;
	double s,*tmp;

	tmp= (double *) c_alloc1d(n+1); /*Complaing about cast */
	for (j=1;j<=n;j++) {
		s=0.0;
		if (A->WW[j] > 1e-6) {
			for (i=1;i<=m;i++) {
				s += A->fexs[i][j]*A->fins[i];
			}
			s /= A->WW[j];
		}
		tmp[j]=s;
		printf("%i...tmp=%le ww=%le fx=%le fin=%le\n",j,tmp[j],A->WW[j],A->fexs[j][j],A->fins[j]);
	}
	for (j=1;j<=n;j++) {
		s=0.0;
		for (jj=1;jj<=n;jj++) {
			s += A->VV[j][jj]*tmp[jj];
		}
		A->rots[j]=s;
	}
	Free1D(tmp);
}
#undef NRANSI
