#include "fitxyee.h"
/*

  fit a line through x y data with uncertainties in both x and y

  COPYRIGHT NUMERICAL RECIPES IN C, p.668, do not distribute without
  permission


  minor modifications:

  - using data structure instead of x,y,sigx,sigy
  - removed global variables and put those into fit structure

*/

/* numerical recipes routines from here on */
void nr_fit(struct nr_dp *data,int ndata,HC_PREC *a,HC_PREC *b,
	 HC_PREC *siga,HC_PREC *sigb,HC_PREC *chi2,
	 HC_PREC *q)
{
  
  int i;
  HC_PREC wt,t,sxoss,sx=0.0,sy=0.0,st2=0.0,ss,sigdat;
  
  *b=0.0;

  ss=0.0;
  for (i=1;i<=ndata;i++) {
    wt=1.0/NR_SQUARE(data[i].sigy);
    ss += wt;
    sx += data[i].x*wt;
    sy += data[i].y*wt;
  }

  
  sxoss=sx/ss;
  for (i=1;i <= ndata;i++) {
    t=(data[i].x-sxoss)/data[i].sigy;
    st2 += t*t;
    *b += t*data[i].y/data[i].sigy;
  }

  *b /= st2;
  *a=(sy-sx*(*b))/ss;
  *siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
  *sigb=sqrt(1.0/st2);
  *chi2=0.0;
  for (i=1;i <= ndata;i++)
    *chi2 += NR_SQUARE(data[i].y-(*a)-(*b)*data[i].x);
  *q=1.0;
  sigdat=sqrt((*chi2)/(ndata-2));
  *siga *= sigdat;
  *sigb *= sigdat;
}


#define POTN 1.571000
#define BIG 1.0e30



void nr_fitexy(struct nr_dp *data,int ndat,HC_PREC *a,HC_PREC *b,
	       HC_PREC *siga,HC_PREC *sigb,HC_PREC *chi2,
	       HC_PREC *q)
{
  int j;
  HC_PREC swap,amx,amn,var[3],varx,vary,ang[7],ch[7],scale,
    bmn,bmx,d1,d2,r2,avea[3],
    dum1,dum2,dum3,dum4,dum5;
  struct nr_fits fit[1];
  fit->xx=nr_vector(1,ndat);
  fit->yy=nr_vector(1,ndat);
  fit->sx=nr_vector(1,ndat);
  fit->sy=nr_vector(1,ndat);
  fit->ww=nr_vector(1,ndat);
  /* compute variance */
  nr_avevar(data,ndat,avea,var);
  varx = var[1];vary = var[2];
  /* reassign and rescale */
  scale = sqrt(varx/vary);
  fit->nn=ndat;
  for (j=1;j <= ndat;j++) {
    //fprintf(stderr,"%i %g %g %g %g\n",j,data[j].x,data[j].y, data[j].sigx, data[j].sigy);
    fit->xx[j]=data[j].x;
    fit->yy[j]=data[j].y*scale;
    fit->sx[j]=data[j].sigx;
    fit->sy[j]=data[j].sigy*scale;
    fit->ww[j]=sqrt(NR_SQUARE(fit->sx[j])+NR_SQUARE(fit->sy[j]));
  }
  nr_fitline(fit->xx,fit->yy,fit->nn,fit->ww,1,&dum1,b,&dum2,&dum3,&dum4,&dum5);
  fit->offs=ang[1]=0.0;
  ang[2]=atan(*b);
  ang[4]=0.0;
  ang[5]=ang[2];
  ang[6]=POTN;
  for (j=4;j<=6;j++)
    ch[j]=nr_chixy(ang[j],fit);
  nr_mnbrak(&ang[1],&ang[2],&ang[3],&ch[1],&ch[2],&ch[3],(HC_PREC (*)(void))nr_chixy,fit);
  *chi2=nr_brent(ang[1],ang[2],ang[3],(HC_PREC (*)(void))nr_chixy,fit,NR_ACC,b);
  *chi2=nr_chixy(*b,fit);
  *a=fit->aa;
  *q=nr_gammq(0.5*(fit->nn-2),*chi2*0.5);
  for (r2=0.0,j=1;j<=fit->nn;j++)
    r2 += fit->ww[j];
  r2=1.0/r2;
  bmx=BIG;
  bmn=BIG;
  fit->offs=(*chi2)+1.0;
  for (j=1;j<=6;j++) {
    if (ch[j] > fit->offs) {
      d1=fabs(ang[j]-(*b));
      while (d1 >= HC_PI) d1 -= HC_PI;
      d2=HC_PI-d1;
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
    bmx=nr_zbrent((HC_PREC (*)(void))nr_chixy,fit,*b,*b+bmx,NR_ACC)-(*b);
    amx=fit->aa-(*a);
    bmn=nr_zbrent((HC_PREC (*)(void))nr_chixy,fit,*b,*b-bmn,NR_ACC)-(*b);
    amn=fit->aa-(*a);
    *sigb=sqrt(0.5*(bmx*bmx+bmn*bmn))/(scale*NR_SQUARE(cos(*b)));
    *siga=sqrt(0.5*(amx*amx+amn*amn)+r2)/scale;
  } else (*sigb)=(*siga)=BIG;
  *a /= scale;
  *b=tan(*b)/scale;
  nr_free_vector(fit->ww,1,ndat);
  nr_free_vector(fit->sy,1,ndat);
  nr_free_vector(fit->sx,1,ndat);
  nr_free_vector(fit->yy,1,ndat);
  nr_free_vector(fit->xx,1,ndat);
}
#undef POTN
#undef BIG



#define ITMAX 10000
#define CGOLD 0.3819660
#define ZEPS 1.0e-14
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

HC_PREC nr_brent(ax,bx,cx,f,fit,tol,xmin)
HC_PREC (*f)(),*xmin,ax,bx,cx,tol;
struct nr_fits *fit;
{
	int iter;
	HC_PREC a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	HC_PREC e=0.0;

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
		d=NR_SIGN(tol1,xm-x);
	    }
	  } else {
	    d=CGOLD*(e=(x >= xm ? a-x : b-x));
	  }
	  u=(fabs(d) >= tol1 ? x+d : x+NR_SIGN(tol1,d));
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
	nr_error("Too many iterations in brent");
	*xmin=x;
	return fx;
}
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT


#define BIG 1.0e30


HC_PREC nr_chixy(bang,fit)
HC_PREC bang;
struct nr_fits *fit;
{
  int j;
  HC_PREC ans,avex=0.0,avey=0.0,sumw=0.0,b;
  
  b=tan(bang);
  for (j=1;j<=fit->nn;j++) {
    fit->ww[j] = NR_SQUARE(b*fit->sx[j])+NR_SQUARE(fit->sy[j]);
    sumw += (fit->ww[j] = (fit->ww[j] == 0.0 ? BIG : 1.0/fit->ww[j]));
    avex += fit->ww[j]*fit->xx[j];
		avey += fit->ww[j]*fit->yy[j];
  }
  if (sumw == 0.0) sumw = BIG;
  avex /= sumw;
  avey /= sumw;
  fit->aa=avey-b*avex;
  for (ans = -(fit->offs),j=1;j<=fit->nn;j++)
    ans += fit->ww[j]*NR_SQUARE(fit->yy[j]-fit->aa-b*fit->xx[j]);
  return ans;
}
#undef BIG

HC_PREC nr_gammq(a,x)
HC_PREC a,x;
{
	void nr_gcf(),nr_gser();
	void nr_error();
	HC_PREC gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) nr_error("Invalid arguments in routine gammq");
	if (x < (a+1.0)) {
	  nr_gser(&gamser,a,x,&gln);
	  return 1.0-gamser;
	} else {
	  nr_gcf(&gammcf,a,x,&gln);
	  return gammcf;
	}
}
#define ITMAX 10000
#define FPMIN 1.0e-30

void nr_gcf(gammcf,a,x,gln)
HC_PREC *gammcf,*gln,a,x;
{
	HC_PREC nr_gammln();
	void nr_error();
	int i;
	HC_PREC an,b,c,d,del,h;

	*gln=nr_gammln(a);
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
		if (fabs(del-1.0) < NR_EPS) break;
	}
	if (i > ITMAX) {
	  fprintf(stderr,"%g %g %g\n",a,x,*gln);
	  nr_error("a too large, ITMAX too small in gcf");
	}
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}
#undef ITMAX
#undef FPMIN


#define ITMAX 1000

void nr_gser(gamser,a,x,gln)
HC_PREC *gamser,*gln,a,x;
{
	HC_PREC nr_gammln();
	void nr_error();
	int n;
	HC_PREC sum,del,ap;

	*gln=nr_gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) nr_error("x less than 0 in routine gser");
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*NR_EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		nr_error("a too large, ITMAX too small in routine gser");
		return;
	}
}
#undef ITMAX

/* CAUTION: This is the traditional K&R C (only) version of the Numerical
   Recipes utility file nrutil.c.  Do not confuse this file with the
   same-named file nrutil.c that is supplied in the same subdirectory or
   archive as the header file nrutil.h.  *That* file contains both ANSI and
   traditional K&R versions, along with #ifdef macros to select the
   correct version.  *This* file contains only traditional K&R.           */

#include <stdio.h>
#define NR_END 1
#define FREE_ARG char*

void nr_error(error_text)
char error_text[];
/* Numerical Recipes standard error handler */
{


	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

HC_PREC *nr_vector(nl,nh)
long nh,nl;
/* allocate a HC_PREC vector with subscript range v[nl..nh] */
{
	HC_PREC *v;

	v=(HC_PREC *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(HC_PREC)));
	if (!v) nr_error("allocation failure in vector()");
	return v-nl+NR_END;
}

int *nr_ivector(nl,nh)
long nh,nl;
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nr_error("allocation failure in ivector()");
	return v-nl+NR_END;
}

unsigned char *nr_cvector(nl,nh)
long nh,nl;
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
	unsigned char *v;

	v=(unsigned char *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
	if (!v) nr_error("allocation failure in cvector()");
	return v-nl+NR_END;
}

unsigned long *nr_lvector(nl,nh)
long nh,nl;
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
	unsigned long *v;

	v=(unsigned long *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(long)));
	if (!v) nr_error("allocation failure in lvector()");
	return v-nl+NR_END;
}

HC_PREC *nr_dvector(nl,nh)
long nh,nl;
/* allocate a HC_PREC vector with subscript range v[nl..nh] */
{
	HC_PREC *v;

	v=(HC_PREC *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(HC_PREC)));
	if (!v) nr_error("allocation failure in dvector()");
	return v-nl+NR_END;
}

HC_PREC **nr_matrix(nrl,nrh,ncl,nch)
long nch,ncl,nrh,nrl;
/* allocate a HC_PREC matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	HC_PREC **m;

	/* allocate pointers to rows */
	m=(HC_PREC **) malloc((unsigned int)((nrow+NR_END)*sizeof(HC_PREC*)));
	if (!m) nr_error("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(HC_PREC *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(HC_PREC)));
	if (!m[nrl]) nr_error("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

HC_PREC **nr_dmatrix(nrl,nrh,ncl,nch)
long nch,ncl,nrh,nrl;
/* allocate a HC_PREC matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	HC_PREC **m;

	/* allocate pointers to rows */
	m=(HC_PREC **) malloc((unsigned int)((nrow+NR_END)*sizeof(HC_PREC*)));
	if (!m) nr_error("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(HC_PREC *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(HC_PREC)));
	if (!m[nrl]) nr_error("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

int **nr_imatrix(nrl,nrh,ncl,nch)
long nch,ncl,nrh,nrl;
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((unsigned int)((nrow+NR_END)*sizeof(int*)));
	if (!m) nr_error("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) nr_error("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

HC_PREC **nr_submatrix(a,oldrl,oldrh,oldcl,oldch,newrl,newcl)
HC_PREC **a;
long newcl,newrl,oldch,oldcl,oldrh,oldrl;
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
	long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
	HC_PREC **m;

	/* allocate array of pointers to rows */
	m=(HC_PREC **) malloc((unsigned int) ((nrow+NR_END)*sizeof(HC_PREC*)));
	if (!m) nr_error("allocation failure in submatrix()");
	m += NR_END;
	m -= newrl;

	/* set pointers to rows */
	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

HC_PREC **nr_convert_matrix(a,nrl,nrh,ncl,nch)
HC_PREC *a;
long nch,ncl,nrh,nrl;
/* allocate a HC_PREC matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	HC_PREC **m;

	/* allocate pointers to rows */
	m=(HC_PREC **) malloc((unsigned int) ((nrow+NR_END)*sizeof(HC_PREC*)));
	if (!m) nr_error("allocation failure in convert_matrix()");
	m += NR_END;
	m -= nrl;

	/* set pointers to rows */
	m[nrl]=a-ncl;
	for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
}

HC_PREC ***nr_f3tensor(nrl,nrh,ncl,nch,ndl,ndh)
long nch,ncl,ndh,ndl,nrh,nrl;
/* allocate a HC_PREC 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	HC_PREC ***t;

	/* allocate pointers to pointers to rows */
	t=(HC_PREC ***) malloc((unsigned int)((nrow+NR_END)*sizeof(HC_PREC**)));
	if (!t) nr_error("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(HC_PREC **) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(HC_PREC*)));
	if (!t[nrl]) nr_error("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(HC_PREC *) malloc((unsigned int)((nrow*ncol*ndep+NR_END)*sizeof(HC_PREC)));
	if (!t[nrl][ncl]) nr_error("allocation failure 3 in f3tensor()");
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

void nr_free_vector(v,nl,nh)
HC_PREC *v;
long nh,nl;
/* free a HC_PREC vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void nr_free_ivector(v,nl,nh)
int *v;
long nh,nl;
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void nr_free_cvector(v,nl,nh)
long nh,nl;
unsigned char *v;
/* free an unsigned char vector allocated with cvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void nr_free_lvector(v,nl,nh)
long nh,nl;
unsigned long *v;
/* free an unsigned long vector allocated with lvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void nr_free_dvector(v,nl,nh)
HC_PREC *v;
long nh,nl;
/* free a HC_PREC vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void nr_free_matrix(m,nrl,nrh,ncl,nch)
HC_PREC **m;
long nch,ncl,nrh,nrl;
/* free a HC_PREC matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void nr_free_dmatrix(m,nrl,nrh,ncl,nch)
HC_PREC **m;
long nch,ncl,nrh,nrl;
/* free a HC_PREC matrix allocated by dmatrix() */
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

void nr_free_submatrix(b,nrl,nrh,ncl,nch)
HC_PREC **b;
long nch,ncl,nrh,nrl;
/* free a submatrix allocated by submatrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void nr_free_convert_matrix(b,nrl,nrh,ncl,nch)
HC_PREC **b;
long nch,ncl,nrh,nrl;
/* free a matrix allocated by convert_matrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void nr_free_f3tensor(t,nrl,nrh,ncl,nch,ndl,ndh)
HC_PREC ***t;
long nch,ncl,ndh,ndl,nrh,nrl;
/* free a HC_PREC  f3tensor allocated by f3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}
void nr_avevar(struct nr_dp *data,unsigned long n,
	    HC_PREC *ave,HC_PREC *var)
{
  unsigned long j;
  HC_PREC s[3],ep[3];

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


void nr_fitline(HC_PREC *x,HC_PREC *y,int ndata,
	     HC_PREC *sig,int mwt,
	     HC_PREC *a,HC_PREC *b,
	     HC_PREC *siga,HC_PREC *sigb,
	     HC_PREC *chi2,HC_PREC *q)
{
  int i;
  HC_PREC wt,t,sxoss,sx=0.0,sy=0.0,st2=0.0,ss,sigdat;
  
  *b=0.0;
  if (mwt) {
    ss=0.0;
    for (i=1;i<=ndata;i++) {
      wt=1.0/NR_SQUARE(sig[i]);
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
      *chi2 += NR_SQUARE(y[i]-(*a)-(*b)*x[i]);
    *q=1.0;
    sigdat=sqrt((*chi2)/(ndata-2));
    *siga *= sigdat;
    *sigb *= sigdat;
  } else {
    for (i=1;i<=ndata;i++)
      *chi2 += NR_SQUARE((y[i]-(*a)-(*b)*x[i])/sig[i]);
    *q=nr_gammq(0.5*(ndata-2),0.5*(*chi2));
  }
}

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

void nr_mnbrak(ax,bx,cx,fa,fb,fc,func,fit)
  HC_PREC (*func)(),*ax,*bx,*cx,*fa,*fb,*fc;
struct nr_fits *fit;
{
  HC_PREC ulim,u,r,q,fu,dum;
  
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
      (2.0*NR_SIGN(NR_FMAX(fabs(q-r),TINY),q-r));
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

HC_PREC nr_zbrent(func,fit,x1,x2,tol)
HC_PREC (*func)(),tol,x1,x2;
struct nr_fits *fit;
{
  int iter;
  HC_PREC a=x1,b=x2,c=x2,d,e,min1,min2;
  HC_PREC fa=(*func)(a,fit),fb=(*func)(b,fit),
    fc,p,q,r,s,tol1,xm;
  
  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
    nr_error("Root must be bracketed in zbrent");
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
    tol1=2.0*NR_EPS*fabs(b)+0.5*tol;
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
      b += NR_SIGN(tol1,xm);
    fb=(*func)(b,fit);
  }
  nr_error("Maximum number of iterations exceeded in zbrent");
  return 0.0;
}
#undef ITMAX

HC_PREC nr_gammln(xx)
HC_PREC xx;
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

int nr_comparef(struct nr_dp *a,struct nr_dp *b)
{
  if(a->x < b->x)
    return -1;
  if(a->x == b->x)
    return 0;
  else
    return 1;
}
