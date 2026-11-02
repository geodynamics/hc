#include "hc.h"

#define N 5
int main(void)
{
  int i,n,m;
  SH_RICK_PREC x[3*N],y[3*N];

  n=N;
  m=2*n+2;

  for(i=0;i<3*N;i++){
    x[i] = y[i] = i;
  }
  i=1;
  rick_realft_nr((x-1),n,i);
  rick_f90_realft(y,&n,&i);
  for(i=0;i<m;i++)
    printf("%11g %11g\n",x[i],y[i]);

  return 0;
}
