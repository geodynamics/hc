/*


  bunch of correlation and cross correlation routines from Numerical Recipes
  see their copyright 


  calculates the cross-correlation 
  of two timeseries
  using Numerical Recipes routines



*/
#include "correl_nr.h"
#define COMP_PRECISION double
/* 

   compute cross correlation for different lags using FFT, and normalize by
   the two RMS values to actually get correlation

   corr nn on output
   

 */
void compute_correl(COMP_PRECISION *x, COMP_PRECISION *y, COMP_PRECISION **corr, int n, int *nn,int verbose)
{

  int lag;
  int i;
  COMP_PRECISION xm,ym,sum1,sum2,rmst,nn2,*xl,*yl;
  // max lag
  lag = n / 2;
  // zero padding
  *nn=pow(2,ceil(log2((COMP_PRECISION)(n+lag))));
  nn2 = *nn * 2;
  
  /* remove means */
  xm = ym = 0;
  sum1 = sum2 = 0.0;
  for(i=0;i<n;i++){
    xm += x[i];sum1 += x[i] * x[i];
    ym += y[i];sum2 += y[i] * y[i];
  }
  xm /= (COMP_PRECISION)n;
  ym /= (COMP_PRECISION)n;
  sum1 /= (COMP_PRECISION)n;
  sum2 /= (COMP_PRECISION)n;

  xl=(COMP_PRECISION *)calloc(*nn,sizeof(COMP_PRECISION));
  yl=(COMP_PRECISION *)calloc(*nn,sizeof(COMP_PRECISION));
  *corr=(COMP_PRECISION *)realloc(*corr,sizeof(COMP_PRECISION)*nn2);
  if(!xl || !yl || !(*corr))ME;
  
  for(i=0;i < n;i++){
    xl[i] = x[i] -= xm;
    yl[i] = y[i] -= ym;
  }
  if(verbose)
    fprintf(stderr,"compute_correl: using %i data pairs (removed means %g %g), next power of two padding for max lag of %i: %i\n",
	    n,xm,ym,lag,*nn);

  // calculate correlation, numrec style call
  correl((xl-1),(yl-1),(*nn),(*corr-1));
  /* normalize to get correlation */
  //fprintf(stderr,"normalizing by %g\n",rmst);
  rmst = sqrt(sum1*sum2)*(n);
  for(i=0;i<nn2;i++)
    *(*corr+i) /= rmst;
  free(xl);free(yl);
}
/*

  from here, numerical recipes routines
  cross correlation
*/
void correl(COMP_PRECISION *data1,COMP_PRECISION *data2,int n,COMP_PRECISION *ans)
{
  unsigned long no2,i;
  COMP_PRECISION dum,*fft;
 
  
  fft=vector(1,n<<1);
  twofft(data1,data2,fft,ans,n);
  no2=n>>1;
  for (i=2;i<=n+2;i+=2) {
    ans[i-1]=(fft[i-1]*(dum=ans[i-1])+fft[i]*ans[i])/no2;
    ans[i]=(fft[i]*dum-fft[i-1]*ans[i])/no2;
  }
  ans[2]=ans[n+1];
  realft(ans,n,-1);
  
  free_vector(fft,1,n<<1);
}

void twofft(COMP_PRECISION *data1,COMP_PRECISION *data2,COMP_PRECISION *fft1,COMP_PRECISION *fft2,
	    int n)
{
  unsigned long nn3,nn2,jj,j;
  COMP_PRECISION rep,rem,aip,aim;
  
  nn3=1+(nn2=2+n+n);
  for (j=1,jj=2;j<=n;j++,jj+=2) {
    fft1[jj-1]=data1[j];
    fft1[jj]=data2[j];
  }
  four1(fft1,n,1);
  fft2[1]=fft1[2];
  fft1[2]=fft2[2]=0.0;
  for (j=3;j<=n+1;j+=2) {
    rep=0.5*(fft1[j]+fft1[nn2-j]);
    rem=0.5*(fft1[j]-fft1[nn2-j]);
    aip=0.5*(fft1[j+1]+fft1[nn3-j]);
    aim=0.5*(fft1[j+1]-fft1[nn3-j]);
    fft1[j]=rep;
    fft1[j+1]=aim;
    fft1[nn2-j]=rep;
    fft1[nn3-j] = -aim;
    fft2[j]=aip;
    fft2[j+1] = -rem;
    fft2[nn2-j]=aip;
    fft2[nn3-j]=rem;
  }
}
void four1(COMP_PRECISION *data,int nn,int isign)
{
  unsigned long n,mmax,m,j,istep,i;
  COMP_PRECISION wtemp,wr,wpr,wpi,wi,theta;
  COMP_PRECISION tempr,tempi;
  
  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2) {
    if (j > i) {
      SWAP(data[j],data[i]);
      SWAP(data[j+1],data[i+1]);
    }
    m=n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax=2;
  while (n > mmax) {
    istep=mmax << 1;
    theta=isign*(6.28318530717959/mmax);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2) {
      for (i=m;i<=n;i+=istep) {
	j=i+mmax;
	tempr=wr*data[j]-wi*data[j+1];
	tempi=wr*data[j+1]+wi*data[j];
	data[j]=data[i]-tempr;
	data[j+1]=data[i+1]-tempi;
	data[i] += tempr;
	data[i+1] += tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}
void realft(COMP_PRECISION *data,int n,int isign)
{
  unsigned long i,i1,i2,i3,i4,np3;
  COMP_PRECISION c1=0.5,c2,h1r,h1i,h2r,h2i;
  COMP_PRECISION wr,wi,wpr,wpi,wtemp,theta;
  
  theta=3.141592653589793/(COMP_PRECISION) (n>>1);
  if (isign == 1) {
    c2 = -0.5;
    four1(data,n>>1,1);
  } else {
    c2=0.5;
    theta = -theta;
  }
  wtemp=sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi=sin(theta);
  wr=1.0+wpr;
  wi=wpi;
  np3=n+3;
  for (i=2;i<=(n>>2);i++) {
    i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
    h1r=c1*(data[i1]+data[i3]);
    h1i=c1*(data[i2]-data[i4]);
    h2r = -c2*(data[i2]+data[i4]);
    h2i=c2*(data[i1]-data[i3]);
    data[i1]=h1r+wr*h2r-wi*h2i;
    data[i2]=h1i+wr*h2i+wi*h2r;
    data[i3]=h1r-wr*h2r+wi*h2i;
    data[i4] = -h1i+wr*h2i+wi*h2r;
    wr=(wtemp=wr)*wpr-wi*wpi+wr;
    wi=wi*wpr+wtemp*wpi+wi;
  }
  if (isign == 1) {
    data[1] = (h1r=data[1])+data[2];
    data[2] = h1r-data[2];
  } else {
    data[1]=c1*((h1r=data[1])+data[2]);
    data[2]=c1*(h1r-data[2]);
    four1(data,n>>1,-1);
  }
}




void  spear(COMP_PRECISION *data1, COMP_PRECISION *data2, unsigned long n,
	    COMP_PRECISION *d, COMP_PRECISION *zd, COMP_PRECISION *probd,
	    COMP_PRECISION *rs, COMP_PRECISION *probrs)
{
  //COMP_PRECISION betai(),erfcc();
  //void crank(),sort2();
	unsigned long j;
	COMP_PRECISION vard,t,sg,sf,fac,en3n,en,df,aved,*wksp1,*wksp2;

	wksp1=vector(1,n);
	wksp2=vector(1,n);
	for (j=1;j<=n;j++) {
		wksp1[j]=data1[j];
		wksp2[j]=data2[j];
	}
	sort2(n,wksp1,wksp2);
	crank(n,wksp1,&sf);
	sort2(n,wksp2,wksp1);
	crank(n,wksp2,&sg);
	*d=0.0;
	for (j=1;j<=n;j++)
		*d += SQR(wksp1[j]-wksp2[j]);
	en=n;
	en3n=en*en*en-en;
	aved=en3n/6.0-(sf+sg)/12.0;
	fac=(1.0-sf/en3n)*(1.0-sg/en3n);
	vard=((en-1.0)*en*en*SQR(en+1.0)/36.0)*fac;
	*zd=(*d-aved)/sqrt(vard);
	*probd=erfcc(fabs(*zd)/1.4142136);
	*rs=(1.0-(6.0/en3n)*(*d+(sf+sg)/12.0))/sqrt(fac);
	fac=(*rs+1.0)*(1.0-(*rs));
	if (fac > 0.0) {
		t=(*rs)*sqrt((en-2.0)/fac);
		df=en-2.0;
		*probrs=betai(0.5*df,0.5,df/(df+t*t));
	} else
		*probrs=0.0;
	free_vector(wksp2,1,n);
	free_vector(wksp1,1,n);
}




/* allocate a COMP_PRECISION vector with subscript range v[nl..nh] */
COMP_PRECISION *vector(int nl,int nh)
{
  COMP_PRECISION *v;
  
  v=(COMP_PRECISION *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(COMP_PRECISION)));
  if (!v) {fprintf(stderr,"allocation failure in vector\n");exit(-1);}
  return v-nl+NR_END;
}

/* free a COMP_PRECISION vector allocated with vector() */
void free_vector(COMP_PRECISION *v,int nl,int nh)
{
  free((FREE_ARG) (v+nl-NR_END));
}


#define NR_TINY 1.0e-20

void pearsn(COMP_PRECISION *x,COMP_PRECISION *y,unsigned long n,
	    COMP_PRECISION *r,COMP_PRECISION *prob,COMP_PRECISION *z)
{
  //COMP_PRECISION betai(),erfcc();
  unsigned long j;
  COMP_PRECISION yt,xt,t,df;
  COMP_PRECISION syy=0.0,sxy=0.0,sxx=0.0,ay=0.0,ax=0.0;
  
  for (j=1;j<=n;j++) {
    ax += x[j];
    ay += y[j];
  }
  ax /= (COMP_PRECISION)n;
  ay /= (COMP_PRECISION)n;
  for (j=1;j<=n;j++) {
    xt=x[j]-ax;
    yt=y[j]-ay;
    sxx += xt*xt;
    syy += yt*yt;
    sxy += xt*yt;
  }
  *r=sxy/sqrt(sxx*syy);
  *z=0.5*log((1.0+(*r)+NR_TINY)/(1.0-(*r)+NR_TINY));
  df=(COMP_PRECISION)n-2.;
  
  /* this is 14.5.5 */
  
  t=(*r)*sqrt(df/((1.0-(*r)+NR_TINY)*(1.0+(*r)+NR_TINY)));
  *prob=betai(0.5*df,0.5,df/(df+t*t));
}
#undef NR_TINY


COMP_PRECISION betai(a,b,x)
COMP_PRECISION a,b,x;
{
	COMP_PRECISION betacf(),gammln();
	void nrerror();
	COMP_PRECISION bt;

	if (x < 0.0 || x > 1.0) nrerror("Bad x in routine betai");
	if (x == 0.0 || x == 1.0) bt=0.0;
	else
		bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
	if (x < (a+1.0)/(a+b+2.0))
		return bt*betacf(a,b,x)/a;
	else
		return 1.0-bt*betacf(b,a,1.0-x)/b;
}

#define MAXIT 100
#define EPS 5.0e-15
#define FPMIN 1.0e-30

COMP_PRECISION betacf(a,b,x)
COMP_PRECISION a,b,x;
{
	void nrerror();
	int m,m2;
	COMP_PRECISION aa,c,d,del,h,qab,qam,qap;

	qab=a+b;
	qap=a+1.0;
	qam=a-1.0;
	c=1.0;
	d=1.0-qab*x/qap;
	if (fabs(d) < FPMIN) d=FPMIN;
	d=1.0/d;
	h=d;
	for (m=1;m<=MAXIT;m++) {
		m2=2*m;
		aa=m*(b-m)*x/((qam+m2)*(a+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (m > MAXIT) nrerror("a or b too big, or MAXIT too small in betacf");
	return h;
}
#undef MAXIT
#undef EPS
#undef FPMIN
void nrerror(error_text)
char error_text[];
/* Numerical Recipes standard error handler */
{


	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}


COMP_PRECISION gammln(xx)
COMP_PRECISION xx;
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}


COMP_PRECISION 
erfcc (double x)
{
	COMP_PRECISION t,z,ans;

	z=fabs(x);
	t=1.0/(1.0+0.5*z);
	ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
		t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
		t*(-0.82215223+t*0.17087277)))))))));
	return x >= 0.0 ? ans : 2.0-ans;
}
void 
crank (unsigned long n, COMP_PRECISION w[], COMP_PRECISION *s)
{
	unsigned long j=1,ji,jt;
	COMP_PRECISION t,rank;

	*s=0.0;
	while (j < n) {
		if (w[j+1] != w[j]) {
			w[j]=j;
			++j;
		} else {
			for (jt=j+1;jt<=n && w[jt]==w[j];jt++);
			rank=0.5*(j+jt-1);
			for (ji=j;ji<=(jt-1);ji++) w[ji]=rank;
			t=jt-j;
			*s += t*t*t-t;
			j=jt;
		}
	}
	if (j == n) w[n]=n;
}

#define NR_SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define NR_M 7
#define NR_STACK 50

void 
sort2 (unsigned long n, COMP_PRECISION arr[], COMP_PRECISION brr[])
{
	unsigned long i,ir=n,j,k,l=1;
	int *istack,jstack=0;
	COMP_PRECISION a,b,temp;

	istack=ivector(1,NR_STACK);
	for (;;) {
		if (ir-l < NR_M) {
			for (j=l+1;j<=ir;j++) {
				a=arr[j];
				b=brr[j];
				for (i=j-1;i>=1;i--) {
					if (arr[i] <= a) break;
					arr[i+1]=arr[i];
					brr[i+1]=brr[i];
				}
				arr[i+1]=a;
				brr[i+1]=b;
			}
			if (!jstack) {
				free_ivector(istack,1,NR_STACK);
				return;
			}
			ir=istack[jstack];
			l=istack[jstack-1];
			jstack -= 2;
		} else {
			k=(l+ir) >> 1;
			NR_SWAP(arr[k],arr[l+1])
			NR_SWAP(brr[k],brr[l+1])
			if (arr[l+1] > arr[ir]) {
				NR_SWAP(arr[l+1],arr[ir])
				NR_SWAP(brr[l+1],brr[ir])
			}
			if (arr[l] > arr[ir]) {
				NR_SWAP(arr[l],arr[ir])
				NR_SWAP(brr[l],brr[ir])
			}
			if (arr[l+1] > arr[l]) {
				NR_SWAP(arr[l+1],arr[l])
				NR_SWAP(brr[l+1],brr[l])
			}
			i=l+1;
			j=ir;
			a=arr[l];
			b=brr[l];
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				NR_SWAP(arr[i],arr[j])
				NR_SWAP(brr[i],brr[j])
			}
			arr[l]=arr[j];
			arr[j]=a;
			brr[l]=brr[j];
			brr[j]=b;
			jstack += 2;
			if (jstack > NR_STACK) nrerror("NR_STACK too small in sort2.");
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
}
#undef NR_M
#undef NR_STACK
#undef NR_SWAP
int *
ivector (long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

void 
free_ivector (int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}
