#include "hc_ptrot.h"

int *c_alloc1d(int n)
{
	int *data;
	
	data = (int *) malloc(n * sizeof(int));
	if(data == NULL) {
		fprintf(stderr, "out of memory\n");
		exit(1);
	}
	return data;
}
 
 double ****c_alloc4d(int dim1, int dim2, int dim3, int dim4) 
 {
 int i,k,j;
 double ****data;
 
 data = (double ****) malloc(dim1*sizeof(double ***));
 
 for (i=0; i<dim1; i++) {
	data[i] = (double ***) malloc(dim2*sizeof(double **));
	for (j=0; j<dim2; j++) {
		data[i][j] = (double **) malloc(dim3*sizeof(double *));
		for (k=0;k<dim3;k++) {
			data[i][j][k] = (double *) malloc(dim4*sizeof(double));
		}
	}
  }
  if (data == NULL) {fprintf(stderr, "out of memory\n");
		exit(1);
	}
	
  return data;
 }

double *****c_alloc5d(int dim1, int dim2, int dim3, int dim4, int dim5) 
 {
 int i,k,j,l;
 double *****data;
 
 data = (double *****) malloc(dim1*sizeof(double ****));
 
 for (i=0; i<dim1; i++) {
	data[i] = (double ****) malloc(dim2*sizeof(double ***));
	for (j=0; j<dim2; j++) {
		data[i][j] = (double ***) malloc(dim3*sizeof(double **));
		for (k=0;k<dim3;k++) {
			data[i][j][k] = (double **) malloc(dim4*sizeof(double *));
			for (l=0;l<dim4;l++) {
				data[i][j][k][l] = (double *) malloc(dim5*sizeof(double));
			}
		}
	}
  }
  return data;
 }


