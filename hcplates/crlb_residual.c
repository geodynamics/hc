#include <math.h>
#include "hcplates.h"

/* Calculates the residuals, from CRLB 
Assumes res & bmag are passed in as pointers*/

void crlb_residual(struct hc_plates_arrays *A, int n, int np, double *res, double *bmag)
{
	int res1, bmag1, i, j;
	double sum;

	*res=0.0;
	*bmag=0.0;

	for (i=1; i<=n; i++) {
		*bmag += A->fins[i]*A->fins[i];
		sum=0.0;	
		for (j=1; j<=n; j++) {	 
			sum += A->fex[i][j]*A->rots[j];
		}
		sum += -A->fins[i];
		*res += sum*sum;
	}
	res1=sqrt(*res);
	*res = res1;
	bmag1=sqrt(*bmag);
	*bmag=bmag1;
}
			
