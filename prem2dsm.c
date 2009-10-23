#include "hc.h"
/* 

read the PREM model in our format and convert to DSM format

 */

int main(int argc,char **argv)
{
  int i,j;
  double rlast;
  struct prem_model prem[1];
  char filename[HC_CHAR_LENGTH];
  if(argc>1)
    strncpy(filename,argv[1],HC_CHAR_LENGTH);
  else
    strncpy(filename,PREM_MODEL_FILE,HC_CHAR_LENGTH);

  prem_read_model(filename,prem,TRUE);
  //fprintf(stderr,"nl %i np %i\n",prem->n,prem->np);
  printf("%i\tn_structure_zone\n",prem->n);
  rlast =0;
  for(i=0;i < prem->n;i++){
    printf("%12.1f %12.1f\t",rlast,prem->rb[i]/1000.);
    for(j=0;j < prem->np;j++)	/* rho */
      printf("%9.4f ",prem->crho[i*prem->np+j]);
    printf("\n\t\t\t\t");
    for(j=0;j < prem->np;j++)	/* vpv */
      printf("%9.4f ",prem->cvpv[i*prem->np+j]);
    printf("\n\t\t\t\t");
    for(j=0;j < prem->np;j++)	/* vph */
      printf("%9.4f ",prem->cvph[i*prem->np+j]);
    printf("\n\t\t\t\t");
    for(j=0;j < prem->np;j++)	/* vsv */
      printf("%9.4f ",prem->cvsv[i*prem->np+j]);
    printf("\n\t\t\t\t");
    for(j=0;j < prem->np;j++)	/* vsh */
      printf("%9.4f ",prem->cvsh[i*prem->np+j]);
    printf("\n\t\t\t\t");
    for(j=0;j < prem->np;j++)	/* eta */
      printf("%9.4f ",prem->ceta[i*prem->np+j]);

    printf("%9.1f %9.1f\n",	/* qmu qkappa */
	   (prem->cqmu[i]> 1e8)?(-1.0):(prem->cqmu[i]),prem->cqkappa[i]);

    
    rlast = prem->rb[i]/1000.;
  }
  return 0;
}
