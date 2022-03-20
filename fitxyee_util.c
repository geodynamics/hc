#include "fitxyee.h"


void fit(struct dp *data,int ndata,PRECISION *a,PRECISION *b,
	 PRECISION *siga,PRECISION *sigb,PRECISION *chi2,
	 PRECISION *q)
{
  
  int i;
  PRECISION wt,t,sxoss,sx=0.0,sy=0.0,st2=0.0,ss,sigdat;
  
  *b=0.0;

  ss=0.0;
  for (i=1;i<=ndata;i++) {
    wt=1.0/SQUARE(data[i].sigy);
    ss += wt;
    sx += data[i].x*wt;
    sy += data[i].y*wt;
  }

  
  sxoss=sx/ss;
  for (i=1;i<=ndata;i++) {
    t=(data[i].x-sxoss)/data[i].sigy;
    st2 += t*t;
    *b += t*data[i].y/data[i].sigy;
  }

  *b /= st2;
  *a=(sy-sx*(*b))/ss;
  *siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
  *sigb=sqrt(1.0/st2);
  *chi2=0.0;
  for (i=1;i<=ndata;i++)
    *chi2 += SQUARE(data[i].y-(*a)-(*b)*data[i].x);
  *q=1.0;
  sigdat=sqrt((*chi2)/(ndata-2));
  *siga *= sigdat;
  *sigb *= sigdat;
}


#define POTN 1.571000
#define BIG 1.0e30
#define PI 3.141592653589793238462643383
#define ACC NR_ACC




void fitexy(struct dp *data,int ndat,PRECISION *a,PRECISION *b,
	    PRECISION *siga,PRECISION *sigb,PRECISION *chi2,
	    PRECISION *q)
{
  int j;
  PRECISION swap,amx,amn,var[3],varx,vary,ang[7],ch[7],scale,
    bmn,bmx,d1,d2,r2,avea[3],
    dum1,dum2,dum3,dum4,dum5;
  struct fits fit[1];
  fit->xx=vector(1,ndat);
  fit->yy=vector(1,ndat);
  fit->sx=vector(1,ndat);
  fit->sy=vector(1,ndat);
  fit->ww=vector(1,ndat);
  /* compute variance */
  avevar(data,ndat,avea,var);
  varx = var[1];vary = var[2];
  /* reassign and rescale */
  scale=sqrt(varx/vary);
  fit->nn=ndat;
  for (j=1;j<=ndat;j++) {
    //fprintf(stderr,"%i %g %g %g %g\n",j,data[j].x,data[j].y, data[j].sigx, data[j].sigy);
    fit->xx[j]=data[j].x;
    fit->yy[j]=data[j].y*scale;
    fit->sx[j]=data[j].sigx;
    fit->sy[j]=data[j].sigy*scale;
    fit->ww[j]=sqrt(SQUARE(fit->sx[j])+SQUARE(fit->sy[j]));
  }
  fitline(fit->xx,fit->yy,fit->nn,fit->ww,1,&dum1,b,&dum2,&dum3,&dum4,&dum5);
  fit->offs=ang[1]=0.0;
  ang[2]=atan(*b);
  ang[4]=0.0;
  ang[5]=ang[2];
  ang[6]=POTN;
  for (j=4;j<=6;j++) ch[j]=chixy(ang[j],fit);
  mnbrak(&ang[1],&ang[2],&ang[3],&ch[1],&ch[2],&ch[3],(PRECISION (*)(void))chixy,fit);
  *chi2=brent(ang[1],ang[2],ang[3],(PRECISION (*)(void))chixy,fit,ACC,b);
  *chi2=chixy(*b,fit);
  *a=fit->aa;
  *q=gammq(0.5*(fit->nn-2),*chi2*0.5);
  for (r2=0.0,j=1;j<=fit->nn;j++) r2 += fit->ww[j];
  r2=1.0/r2;
  bmx=BIG;
  bmn=BIG;
  fit->offs=(*chi2)+1.0;
  for (j=1;j<=6;j++) {
    if (ch[j] > fit->offs) {
      d1=fabs(ang[j]-(*b));
      while (d1 >= PI) d1 -= PI;
      d2=PI-d1;
      if (ang[j] < *b) {
	swap=d1;
	d1=d2;
	d2=swap;
      }
      if (d1 < bmx) bmx=d1;
      if (d2 < bmn) bmn=d2;
    }
  }
  if (bmx < BIG) {
    bmx=zbrent((PRECISION (*)(void))chixy,fit,*b,*b+bmx,ACC)-(*b);
    amx=fit->aa-(*a);
    bmn=zbrent((PRECISION (*)(void))chixy,fit,*b,*b-bmn,ACC)-(*b);
    amn=fit->aa-(*a);
    *sigb=sqrt(0.5*(bmx*bmx+bmn*bmn))/(scale*SQUARE(cos(*b)));
    *siga=sqrt(0.5*(amx*amx+amn*amn)+r2)/scale;
  } else (*sigb)=(*siga)=BIG;
  *a /= scale;
  *b=tan(*b)/scale;
  free_vector(fit->ww,1,ndat);
  free_vector(fit->sy,1,ndat);
  free_vector(fit->sx,1,ndat);
  free_vector(fit->yy,1,ndat);
  free_vector(fit->xx,1,ndat);
}
#undef POTN
#undef BIG
#undef PI
#undef ACC


#define ITMAX 10000
#define CGOLD 0.3819660
#define ZEPS 1.0e-14
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

PRECISION brent(ax,bx,cx,f,fit,tol,xmin)
PRECISION (*f)(),*xmin,ax,bx,cx,tol;
struct fits *fit;
{
	int iter;
	PRECISION a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	PRECISION e=0.0;

	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=(*f)(x,fit);
	for (iter=1;iter<=ITMAX;iter++) {
	  xm=0.5*(a+b);
	  tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
	  if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
	    *xmin=x;
	    return fx;
	  }
	  if (fabs(e) > tol1) {
	    r=(x-w)*(fx-fv);
	    q=(x-v)*(fx-fw);
	    p=(x-v)*q-(x-w)*r;
	    q=2.0*(q-r);
	    if (q > 0.0) p = -p;
	    q=fabs(q);
	    etemp=e;
	    e=d;
	    if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	      d=CGOLD*(e=(x >= xm ? a-x : b-x));
	    else {
	      d=p/q;
	      u=x+d;
	      if (u-a < tol2 || b-u < tol2)
		d=SIGN(tol1,xm-x);
	    }
	  } else {
	    d=CGOLD*(e=(x >= xm ? a-x : b-x));
	  }
	  u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
	  fu=(*f)(u,fit);
	  if (fu <= fx) {
	    if (u >= x) a=x; else b=x;
	    SHFT(v,w,x,u)
	      SHFT(fv,fw,fx,fu)
	      } else {
		if (u < x) a=u; else b=u;
		if (fu <= fw || w == x) {
		  v=w;
		  w=u;
		  fv=fw;
		  fw=fu;
		} else if (fu <= fv || v == x || v == w) {
		  v=u;
		  fv=fu;
		}
	      }
	}
	nrerror("Too many iterations in brent");
	*xmin=x;
	return fx;
}
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT


#define BIG 1.0e30


PRECISION chixy(bang,fit)
PRECISION bang;
struct fits *fit;
{
  int j;
  PRECISION ans,avex=0.0,avey=0.0,sumw=0.0,b;
  
  b=tan(bang);
  for (j=1;j<=fit->nn;j++) {
    fit->ww[j] = SQUARE(b*fit->sx[j])+SQUARE(fit->sy[j]);
    sumw += (fit->ww[j] = (fit->ww[j] == 0.0 ? BIG : 1.0/fit->ww[j]));
    avex += fit->ww[j]*fit->xx[j];
		avey += fit->ww[j]*fit->yy[j];
  }
  if (sumw == 0.0) sumw = BIG;
  avex /= sumw;
  avey /= sumw;
  fit->aa=avey-b*avex;
  for (ans = -(fit->offs),j=1;j<=fit->nn;j++)
    ans += fit->ww[j]*SQUARE(fit->yy[j]-fit->aa-b*fit->xx[j]);
  return ans;
}
#undef BIG

PRECISION gammq(a,x)
PRECISION a,x;
{
	void gcf(),gser();
	void nrerror();
	PRECISION gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine gammq");
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return 1.0-gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return gammcf;
	}
}
#define ITMAX 10000
#define EPS NR_EPS
#define FPMIN 1.0e-30

void gcf(gammcf,a,x,gln)
PRECISION *gammcf,*gln,a,x;
{
	PRECISION gammln();
	void nrerror();
	int i;
	PRECISION an,b,c,d,del,h;

	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > ITMAX) {
	  fprintf(stderr,"%g %g %g\n",a,x,*gln);
	  nrerror("a too large, ITMAX too small in gcf");
	}
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}
#undef ITMAX
#undef EPS
#undef FPMIN


#define ITMAX 1000
#define EPS 5.0e-15

void gser(gamser,a,x,gln)
PRECISION *gamser,*gln,a,x;
{
	PRECISION gammln();
	void nrerror();
	int n;
	PRECISION sum,del,ap;

	*gln=gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) nrerror("x less than 0 in routine gser");
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		nrerror("a too large, ITMAX too small in routine gser");
		return;
	}
}
#undef ITMAX
#undef EPS
/* CAUTION: This is the traditional K&R C (only) version of the Numerical
   Recipes utility file nrutil.c.  Do not confuse this file with the
   same-named file nrutil.c that is supplied in the same subdirectory or
   archive as the header file nrutil.h.  *That* file contains both ANSI and
   traditional K&R versions, along with #ifdef macros to select the
   correct version.  *This* file contains only traditional K&R.           */

#include <stdio.h>
#define NR_END 1
#define FREE_ARG char*

void nrerror(error_text)
char error_text[];
/* Numerical Recipes standard error handler */
{


	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

PRECISION *vector(nl,nh)
long nh,nl;
/* allocate a PRECISION vector with subscript range v[nl..nh] */
{
	PRECISION *v;

	v=(PRECISION *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(PRECISION)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

int *ivector(nl,nh)
long nh,nl;
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

unsigned char *cvector(nl,nh)
long nh,nl;
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
	unsigned char *v;

	v=(unsigned char *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
	if (!v) nrerror("allocation failure in cvector()");
	return v-nl+NR_END;
}

unsigned long *lvector(nl,nh)
long nh,nl;
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
	unsigned long *v;

	v=(unsigned long *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(long)));
	if (!v) nrerror("allocation failure in lvector()");
	return v-nl+NR_END;
}

PRECISION *dvector(nl,nh)
long nh,nl;
/* allocate a PRECISION vector with subscript range v[nl..nh] */
{
	PRECISION *v;

	v=(PRECISION *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(PRECISION)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}

PRECISION **matrix(nrl,nrh,ncl,nch)
long nch,ncl,nrh,nrl;
/* allocate a PRECISION matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	PRECISION **m;

	/* allocate pointers to rows */
	m=(PRECISION **) malloc((unsigned int)((nrow+NR_END)*sizeof(PRECISION*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(PRECISION *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(PRECISION)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

PRECISION **dmatrix(nrl,nrh,ncl,nch)
long nch,ncl,nrh,nrl;
/* allocate a PRECISION matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	PRECISION **m;

	/* allocate pointers to rows */
	m=(PRECISION **) malloc((unsigned int)((nrow+NR_END)*sizeof(PRECISION*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(PRECISION *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(PRECISION)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

int **imatrix(nrl,nrh,ncl,nch)
long nch,ncl,nrh,nrl;
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((unsigned int)((nrow+NR_END)*sizeof(int*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

PRECISION **submatrix(a,oldrl,oldrh,oldcl,oldch,newrl,newcl)
PRECISION **a;
long newcl,newrl,oldch,oldcl,oldrh,oldrl;
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
	long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
	PRECISION **m;

	/* allocate array of pointers to rows */
	m=(PRECISION **) malloc((unsigned int) ((nrow+NR_END)*sizeof(PRECISION*)));
	if (!m) nrerror("allocation failure in submatrix()");
	m += NR_END;
	m -= newrl;

	/* set pointers to rows */
	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

PRECISION **convert_matrix(a,nrl,nrh,ncl,nch)
PRECISION *a;
long nch,ncl,nrh,nrl;
/* allocate a PRECISION matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	PRECISION **m;

	/* allocate pointers to rows */
	m=(PRECISION **) malloc((unsigned int) ((nrow+NR_END)*sizeof(PRECISION*)));
	if (!m) nrerror("allocation failure in convert_matrix()");
	m += NR_END;
	m -= nrl;

	/* set pointers to rows */
	m[nrl]=a-ncl;
	for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
}

PRECISION ***f3tensor(nrl,nrh,ncl,nch,ndl,ndh)
long nch,ncl,ndh,ndl,nrh,nrl;
/* allocate a PRECISION 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	PRECISION ***t;

	/* allocate pointers to pointers to rows */
	t=(PRECISION ***) malloc((unsigned int)((nrow+NR_END)*sizeof(PRECISION**)));
	if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(PRECISION **) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(PRECISION*)));
	if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(PRECISION *) malloc((unsigned int)((nrow*ncol*ndep+NR_END)*sizeof(PRECISION)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

void free_vector(v,nl,nh)
PRECISION *v;
long nh,nl;
/* free a PRECISION vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(v,nl,nh)
int *v;
long nh,nl;
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_cvector(v,nl,nh)
long nh,nl;
unsigned char *v;
/* free an unsigned char vector allocated with cvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_lvector(v,nl,nh)
long nh,nl;
unsigned long *v;
/* free an unsigned long vector allocated with lvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_dvector(v,nl,nh)
PRECISION *v;
long nh,nl;
/* free a PRECISION vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_matrix(m,nrl,nrh,ncl,nch)
PRECISION **m;
long nch,ncl,nrh,nrl;
/* free a PRECISION matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_dmatrix(m,nrl,nrh,ncl,nch)
PRECISION **m;
long nch,ncl,nrh,nrl;
/* free a PRECISION matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_imatrix(m,nrl,nrh,ncl,nch)
int **m;
long nch,ncl,nrh,nrl;
/* free an int matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_submatrix(b,nrl,nrh,ncl,nch)
PRECISION **b;
long nch,ncl,nrh,nrl;
/* free a submatrix allocated by submatrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_convert_matrix(b,nrl,nrh,ncl,nch)
PRECISION **b;
long nch,ncl,nrh,nrl;
/* free a matrix allocated by convert_matrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_f3tensor(t,nrl,nrh,ncl,nch,ndl,ndh)
PRECISION ***t;
long nch,ncl,ndh,ndl,nrh,nrl;
/* free a PRECISION  f3tensor allocated by f3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}
void avevar(struct dp *data,unsigned long n,
	    PRECISION *ave,PRECISION *var)
{
  unsigned long j;
  PRECISION s[3],ep[3];

  for (ave[1]=0.0,j=1;j<=n;j++) ave[1] += data[j].x;
  for (ave[2]=0.0,j=1;j<=n;j++) ave[2] += data[j].y;
  ave[1] /= n;
  ave[2] /= n;
  var[1]=var[2]=ep[1]=ep[2]=0.0;
  for (j=1;j<=n;j++) {
    s[1]=data[j].x-ave[1];
    s[2]=data[j].y-ave[2];
    ep[1] += s[1];
    ep[2] += s[2];
    var[1] += s[1]*s[1];
    var[2] += s[2]*s[2];
  }
  var[1]=(var[1]-ep[1]*ep[1]/n)/(n-1);
  var[2]=(var[2]-ep[2]*ep[2]/n)/(n-1);
}


void fitline(PRECISION *x,PRECISION *y,int ndata,
	     PRECISION *sig,int mwt,
	     PRECISION *a,PRECISION *b,
	     PRECISION *siga,PRECISION *sigb,
	     PRECISION *chi2,PRECISION *q)
{
  int i;
  PRECISION wt,t,sxoss,sx=0.0,sy=0.0,st2=0.0,ss,sigdat;
  
  *b=0.0;
  if (mwt) {
    ss=0.0;
    for (i=1;i<=ndata;i++) {
      wt=1.0/SQUARE(sig[i]);
      ss += wt;
      sx += x[i]*wt;
      sy += y[i]*wt;
    }
  } else {
    for (i=1;i<=ndata;i++) {
      sx += x[i];
      sy += y[i];
    }
    ss=ndata;
  }
  sxoss=sx/ss;
  if (mwt) {
    for (i=1;i<=ndata;i++) {
      t=(x[i]-sxoss)/sig[i];
      st2 += t*t;
      *b += t*y[i]/sig[i];
    }
  } else {
    for (i=1;i<=ndata;i++) {
      t=x[i]-sxoss;
      st2 += t*t;
      *b += t*y[i];
    }
  }
  *b /= st2;
  *a=(sy-sx*(*b))/ss;
  *siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
  *sigb=sqrt(1.0/st2);
  *chi2=0.0;
  if (mwt == 0) {
    for (i=1;i<=ndata;i++)
      *chi2 += SQUARE(y[i]-(*a)-(*b)*x[i]);
    *q=1.0;
    sigdat=sqrt((*chi2)/(ndata-2));
    *siga *= sigdat;
    *sigb *= sigdat;
  } else {
    for (i=1;i<=ndata;i++)
      *chi2 += SQUARE((y[i]-(*a)-(*b)*x[i])/sig[i]);
    *q=gammq(0.5*(ndata-2),0.5*(*chi2));
  }
}

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

void mnbrak(ax,bx,cx,fa,fb,fc,func,fit)
PRECISION (*func)(),*ax,*bx,*cx,*fa,*fb,*fc;
struct fits *fit;
{
  PRECISION ulim,u,r,q,fu,dum;
  
  *fa=(*func)(*ax,fit);
  *fb=(*func)(*bx,fit);
  if (*fb > *fa) {
    SHFT(dum,*ax,*bx,dum)
      SHFT(dum,*fb,*fa,dum)
      }
  *cx=(*bx)+GOLD*(*bx-*ax);
  *fc=(*func)(*cx,fit);
  while (*fb > *fc) {
    r=(*bx-*ax)*(*fb-*fc);
    q=(*bx-*cx)*(*fb-*fa);
    u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
      (2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
    ulim=(*bx)+GLIMIT*(*cx-*bx);
    if ((*bx-u)*(u-*cx) > 0.0) {
      fu=(*func)(u,fit);
      if (fu < *fc) {
	*ax=(*bx);
	*bx=u;
	*fa=(*fb);
	*fb=fu;
	return;
      } else if (fu > *fb) {
	*cx=u;
	*fc=fu;
	return;
      }
      u=(*cx)+GOLD*(*cx-*bx);
      fu=(*func)(u,fit);
    } else if ((*cx-u)*(u-ulim) > 0.0) {
      fu=(*func)(u,fit);
      if (fu < *fc) {
	SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
	  SHFT(*fb,*fc,fu,(*func)(u,fit))
	  }
    } else if ((u-ulim)*(ulim-*cx) >= 0.0) {
      u=ulim;
      fu=(*func)(u,fit);
    } else {
      u=(*cx)+GOLD*(*cx-*bx);
      fu=(*func)(u,fit);
    }
    SHFT(*ax,*bx,*cx,u)
      SHFT(*fa,*fb,*fc,fu)
      }
}
#undef GOLD
#undef GLIMIT
#undef TINY
#undef SHFT
#define ITMAX 1000
#define EPS 5.0e-15

PRECISION zbrent(func,fit,x1,x2,tol)
PRECISION (*func)(),tol,x1,x2;
struct fits *fit;
{
  int iter;
  PRECISION a=x1,b=x2,c=x2,d,e,min1,min2;
  PRECISION fa=(*func)(a,fit),fb=(*func)(b,fit),
    fc,p,q,r,s,tol1,xm;
  
  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
    nrerror("Root must be bracketed in zbrent");
  fc=fb;
  for (iter=1;iter<=ITMAX;iter++) {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c=a;
      fc=fa;
      e=d=b-a;
    }
    if (fabs(fc) < fabs(fb)) {
      a=b;
      b=c;
      c=a;
      fa=fb;
      fb=fc;
      fc=fa;
    }
    tol1=2.0*EPS*fabs(b)+0.5*tol;
    xm=0.5*(c-b);
    if (fabs(xm) <= tol1 || fb == 0.0) return b;
    if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
      s=fb/fa;
      if (a == c) {
	p=2.0*xm*s;
	q=1.0-s;
      } else {
	q=fa/fc;
	r=fb/fc;
	p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
	q=(q-1.0)*(r-1.0)*(s-1.0);
      }
      if (p > 0.0) q = -q;
      p=fabs(p);
      min1=3.0*xm*q-fabs(tol1*q);
      min2=fabs(e*q);
      if (2.0*p < (min1 < min2 ? min1 : min2)) {
	e=d;
	d=p/q;
      } else {
	d=xm;
	e=d;
      }
    } else {
      d=xm;
      e=d;
    }
    a=b;
    fa=fb;
    if (fabs(d) > tol1)
      b += d;
    else
      b += SIGN(tol1,xm);
    fb=(*func)(b,fit);
  }
  nrerror("Maximum number of iterations exceeded in zbrent");
  return 0.0;
}
#undef ITMAX
#undef EPS
PRECISION gammln(xx)
PRECISION xx;
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

int comparef(struct dp *a,struct dp *b)
{
  if(a->x < b->x)
    return -1;
  if(a->x == b->x)
    return 0;
  else
    return 1;
}
