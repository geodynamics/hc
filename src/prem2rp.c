#include "hc.h"
/* 

read the PREM model and print density with depth

 */

int main(int argc,char **argv)
{
  double rho, r,z,rnd;
  struct prem_model prem[1];
  char filename[HC_CHAR_LENGTH];
  if(argc>1)
    strncpy(filename,argv[1],HC_CHAR_LENGTH);
  else
    strncpy(filename,PREM_MODEL_FILE,HC_CHAR_LENGTH);

  prem_read_model(filename,prem,TRUE);

  for(r=0;r< HC_RE_KM;r+=0.5){
    rnd = r/HC_RE_KM;
    z = HC_RE_KM-r;
    prem_get_rho(&rho,rnd, prem);
    fprintf(stdout,"%11g %11g %12.4e %12.4e\n",z,r,rho);
  }
  return 0;
}
