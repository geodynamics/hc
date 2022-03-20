/* 


read in two grid files and compute the correlation of points within a
certain radius around each location

this assumes spherical grids by default, but see switch

grdgrd2correlation: usage:
	grdgrd2correlation file1.grd file2.grd [radius, 500 km] [remove_lin, 0] [only_print_stat, 1] [is_global, 1] [lin_reg_mode, 2]

grdgrd2correlation: output:
	lon lat correlation slope number_of_points diff_from_global


for regional grids, this assumes 0...360 convention


*/
#include "hc_ggrd.h"
#include "hc.h"
#include "fitxyee.h"		/* part of numerical recipes */

void calc_mean(struct nr_dp *, int ,double *);
void calc_std(struct nr_dp *, int ,double *,double *);

double correlation(struct nr_dp *, double * , int );
void calc_azi_r_lonlat(double ,double ,double ,double ,
		       double *,double *);
double mod(double , double );	/*  */

int main(int argc, char **argv)
{
  struct ggrd_gt g1[1],g2[1];
  char char_dummy,gmt_flags[100];
  int n,nphi,i,nnan,ret1,ret2;

  int is_global = 1;
  static int use_dist_weighting = 1;
  double x,y,x1,y1,dx,dy,rad_km,xloc,yloc,dri,corr,mean_ratio,
    xmin,xmax,ymin,ymax,tmpx,tmpy,dxmin,dxmax,
    weight_sum,weight;
  double mean[2],std[2],r,a,b,siga,sigb,chi2,q,dr,
    dphi,phi,aglobal,bglobal,cloc[2];
  struct nr_dp *data;
  int lin_reg_mode = 2;		/* 1: errors in x 2: errors in x and y */
  int remove = 0;		/* match the layer correlation? */
  int pstat = 0;		/* exit after global stats */
  FILE *out;

  //double x_dump = -10, y_dump = 40; /* for sample dumps, if wanted */
  double x_dump = -999, y_dump = -999;
  char out_name[500];
  
  sprintf(out_name,"tmp.%g.%g.sample.dump",x_dump,y_dump);
  
  data=(struct nr_dp *)malloc(sizeof(struct nr_dp));

  rad_km = 500;
  if(argc < 3){
    fprintf(stderr,"%s: usage:\n\t%s file1.grd file2.grd [radius, %g km] [remove_lin, %i] [only_print_stat, %i] [is_global, %i] [lin_reg_mode, %i]\n\n",
	    argv[0],argv[0],rad_km,remove,pstat,is_global,lin_reg_mode);
    fprintf(stderr,"%s: output:\n\tlon lat correlation slope number_of_points diff_from_global\n\n",argv[0]);
    exit(-1);
  }
  if(argc > 3 )
    sscanf(argv[3],"%lf",&rad_km);
  if(argc > 4 )
    sscanf(argv[4],"%i",&remove);
  if(argc > 5 )
    sscanf(argv[5],"%i",&pstat);
  if(argc > 6 )
    sscanf(argv[6],"%i",&is_global);
  if(argc > 7 )
    sscanf(argv[7],"%i",&lin_reg_mode);



  if(is_global)
    sprintf(gmt_flags,"-fg");
  else
    sprintf(gmt_flags,"-f0x,1y");
  
  fprintf(stderr,"%s: assuming grids %s and %s are %s, computing for radius %g km, remove_lin: %i, linreg mode: %i\n",
	  argv[0],argv[1],argv[2],(is_global)?("global"):("regional"),rad_km,remove,
	  lin_reg_mode);
  
  if(ggrd_grdtrack_init_general(FALSE,argv[1],&char_dummy,
				gmt_flags,g1,TRUE,FALSE,FALSE)){
    fprintf(stderr,"%s: error reading %s\n",argv[0],argv[1]);
    exit(-1);
  }
  if(ggrd_grdtrack_init_general(FALSE,argv[2],&char_dummy,
				gmt_flags,g2,TRUE,FALSE,FALSE)){
    fprintf(stderr,"%s: error reading %s\n",argv[0],argv[2]);
    exit(-1);
  }
  fprintf(stderr,"%s: read data ok\n",argv[0]);
  /* 

     compute total correlation and best-fit slopes from roughly even
     area sampling - output in on different dampling
     
  */
  if(is_global){
    ymin = -89.75;ymax = 89.75;
    //ymin = -89.5;ymax = 89.5;
    xmin = 0;  xmax = 360;
    dy = .25;
    //dy = .5;
  }else{
    ymin = MAX(g1->south,g2->south);
    ymax = MIN(g1->north,g2->north);
    xmin = MAX(g1->west,g2->west);
    xmax = MIN(g1->east,g2->east);
    //dy = .2;			/* output sampling is different */
    dy = .25;
  }



  dxmin = 1e20;dxmax = -1e20;
  fprintf(stderr,"%s: sampling %g-%g-%g %g-%g-%g\n",argv[0],xmin,dy/cos(((ymin+ymax)/2)/GGRD_PIF),xmax,ymin,dy,ymax);
  for(n=nnan=0,y1 = ymax+1e-7,y = ymin;y <= y1;y += dy){
    dx = dy/cos(y/GGRD_PIF);	/* adjust for sphericity */
    if(dx < dxmin)
      dxmin = dx;
    if(dx > dxmax)
      dxmax = dx;
    for(x=xmin,x1=xmax-dx+1e-7;x <= x1;x += dx){

      if((!ggrd_grdtrack_interpolate_xy(x,y,g1,&tmpx,FALSE))||
	 (!ggrd_grdtrack_interpolate_xy(x,y,g2,&tmpy,FALSE))){
	fprintf(stderr,"%s: interpolation error for lon %g lat %g\n",argv[0],x,y);
	exit(-1);
      }

      if(finite(tmpx) && finite(tmpy)){ /* only use non NaN */
	data[n].x = tmpx;
	data[n].y = tmpy;
	/* uncertainties? */
	data[n].sigx = data[n].sigy = 1.0;
	data=(struct nr_dp *)realloc(data,sizeof(struct nr_dp)*(n+2));
	n++;
      }else{
	nnan++;
      }
    }
  }
  fprintf(stderr,"%s: scanned region -R%g/%g/%g/%g, -I(%g-%g)/%g, nsample: %i nnan: %i\n",
	  argv[0],xmin,xmax,ymin, ymax,dxmin,dxmax,dy,n,nnan);
 
  /* 

  compute global values

  */
  calc_mean(data,n,mean);
  calc_std(data,n,std,mean);
  corr = correlation(data,mean,n);	/* compute correlation */
  //fprintf(stderr,"%s: mean: %g %g std %g %g corr %g\n",argv[0],mean[0],mean[1],std[0],std[1],corr);
  switch(lin_reg_mode){			/* fit a linear relation ymod = a + b * x */
  case 1:
    /* best fit slope, only error in y */
    nr_fit((data-1),n,&aglobal, &bglobal,&siga,&sigb,&chi2,&q);
    break;
  case 2:
    /* best fit slope */
    nr_fitexy((data-1),n,&aglobal,&bglobal,&siga,&sigb,&chi2,&q);
    break;
  default:
    fprintf(stderr,"%s: linear regression mode %i undefined\n",argv[0],lin_reg_mode);
    exit(-1);
    break;
  }


  fprintf(stderr,"%s: first  grid: mean: %11g std: %11g\n",argv[0],mean[0],std[0]);
  fprintf(stderr,"%s: second grid: mean: %11g std: %11g\n",argv[0],mean[1],std[1]);
  fprintf(stderr,"%s: correlation: %11g best-fit: offset: %11g slope: %11g (%s)\n",
	  argv[0],corr,aglobal,bglobal,(lin_reg_mode==1)?("y err"):("x and y err"));
  if(pstat){
    printf("%g %.5e %.5e\n",corr,aglobal,bglobal);
    exit(-1);
  }
  /* 
     
  regional correlation

  */
  if(remove){
    fprintf(stderr,"%s: removing trend by correcting second grid by %g + %g * x\n",argv[0],aglobal,bglobal);
  }

  dri = 10;			/* spacing of circle sampling */
  dr = rad_km / dri;
  
  if(is_global){		/* output dampling */
    ymin = -89;ymax = 89;
    xmin = 0;  xmax = 359;
    dy = dx = 1;
  }else{
    ymin = MAX(g1->south,g2->south);
    ymax = MIN(g1->north,g2->north);
    xmin = MAX(g1->west,g2->west);
    xmax = MIN(g1->east,g2->east);
    dy = dx = .2;
  }

  for(y1 = ymax+1e-7,y = ymin;y <= y1;y += dy){
    for(x=xmin,x1=xmax+1e-7;x <= x1;x += dx){

      /* start local loop */
      n = 0;
      mean_ratio = 0.0;
      weight_sum = 0.0;
      data=(struct nr_dp *)realloc(data,sizeof(struct nr_dp));
      /* set up a circular sampling region around x,y */
      for(r=dr;r <= rad_km+1e-7;r += dr){
	dphi = dr / r * GGRD_TWOPI/dri;
	nphi = (int)(GGRD_TWOPI/dphi);
	dphi = GGRD_TWOPI/(double)nphi;
	for(phi = 0;phi < GGRD_TWOPI;phi+=dphi){
	  /* 
	     compute new location 
	  */
	  calc_azi_r_lonlat(x,y,r,phi,&xloc,&yloc);

	  /* 
	     assign to data array 
	  */
	  ret1=ggrd_grdtrack_interpolate_xy(xloc,yloc,g1,&tmpx,FALSE);
	  ret2=ggrd_grdtrack_interpolate_xy(xloc,yloc,g2,&tmpy,FALSE);
	  //fprintf(stderr,"%i %i\n",ret1,ret2);
	  if(finite(tmpx) && finite(tmpy)){
	    /*  */
	    data[n].x = tmpx;
	    data[n].y = tmpy;
	    if(use_dist_weighting){
	      weight = 1/r;
	      data[n].sigx = data[n].sigy = 1/weight;
	      mean_ratio += log(tmpy/tmpx) * weight;
	    }else{
	      weight = 1;
	      mean_ratio += log(tmpy/tmpx) * weight;
	      data[n].sigx = data[n].sigy = 1/weight;
	    }
	    if(remove){
	      /* remove global trend from y? */
	      data[n].y -= aglobal +  data[n].x * bglobal;
	    }
	    /*  */
	    data=(struct nr_dp *)
	      realloc(data,sizeof(struct nr_dp)*(n+2));
	    //fprintf(stderr,"%11g %11g %11g %11g\n",xloc,yloc,data[n].x,data[n].y);
	    weight_sum += weight;
	    n++;
	  }
	} /* phi loop */
      }  /* dr loop */
      //fprintf(stderr,"%g %g %i - %g %g\n",x,y,n,dr,rad_km);
      if(x == x_dump && y == y_dump){
	fprintf(stderr,"%s: writing samples to %s\n",argv[0],out_name);
	out = fopen(out_name,"w");
	for(i=0;i<n;i++)
	  fprintf(out,"%g %g %g %g\n",data[i].x,data[i].y,data[i].sigx,data[i].sigy);
	fclose(out);
      }
      if(0){
	for(i=0;i<n;i++)
	  fprintf(stderr,"%g %g %g %g\n",data[i].x,data[i].y,data[i].sigx,data[i].sigy);
	exit(-1);
      }

      mean_ratio = exp(mean_ratio/weight_sum);
      if(n > 5){
	/* compute means for both */
	calc_mean(data,n,mean);	
	calc_std(data,n,std,mean);
	/* 
	   central values 
	*/
	calc_azi_r_lonlat(x,y,0,phi,&xloc,&yloc);
	ret1=ggrd_grdtrack_interpolate_xy(xloc,yloc,g1,(cloc),  FALSE);
	ret2=ggrd_grdtrack_interpolate_xy(xloc,yloc,g2,(cloc+1),FALSE);
	/* 
	   output is 
	   
	   lon lat correlation slope number_of_points diff_from_global mean_ratio
	   
	*/
	//fprintf(stderr,"%g %g %g %g\n",aglobal, cloc[0],bglobal,cloc[1]);
	if(finite(cloc[0]) && finite(cloc[1])){
	  if(std[0] > 1e-5 && std[1] > 1e-5){
	    /* compute local correlation */
	    corr = correlation(data,mean,n);
	    /* best fit with errors in x and y  */
	    nr_fitexy((data-1),n,&a,&b,&siga,&sigb,&chi2,&q);
	    fprintf(stdout,"%11g %11g %11g %11g %5i %11g %11g\n",
		    x,y,corr,b,n,(aglobal + cloc[0] * bglobal)-cloc[1],
		    mean_ratio);
	  }else{
	    fprintf(stdout,"%11g %11g NaN NaN %5i %11g NaN\n",
		    x,y,n,(aglobal + cloc[0] * bglobal)-cloc[1]);


	  }
	}
      }
    }
  }
    


}
/* linear correlation coefficient */
double correlation(struct nr_dp *data, double *mean, int n)
{
  int i,nuse;
  double s1,s2,s3,dx,dy;
  s1 = s2 = s3 = 0.0;
  for(i=0;i<n;i++){
    if(finite(data[i].x) && finite(data[i].y)){
      dx = data[i].x - mean[0];
      dy = data[i].y - mean[1];
      s1 += (dx*dy);
      s2 += (dx*dx);
      s3 += (dy*dy);
      nuse++;
    }
  }
  return s1/(sqrt(s2)*sqrt(s3));
}


void calc_mean(struct nr_dp *data, int n, double *mean)
{
  int i,nuse;
  mean[0] = mean[1] = 0.0;
  for(i=nuse=0;i < n;i++){
    if(finite(data[i].x) && finite(data[i].y)){
      //fprintf(stderr,"%g %g\n",data[i].x,data[i].y);
      mean[0] += data[i].x;
      mean[1] += data[i].y;
      nuse++;
    }
  }
  if(nuse){
    mean[0] /= (double)nuse;
    mean[1] /= (double)nuse;
  }
}

void calc_std(struct nr_dp *data, int n, double *std, double *mean)
{
  int i,nuse;
  double tmp;
  std[0] = std[1] = 0.0;
  for(i=nuse=0;i < n;i++){
    if(finite(data[i].x) && finite(data[i].y)){
      tmp = (data[i].x-mean[0]);
      std[0] += tmp*tmp;
      tmp = (data[i].y-mean[1]);
      std[1] += tmp*tmp;
      nuse++;
    }
  }
  if(nuse){
    std[0] /= (double)nuse;
    std[1] /= (double)nuse;
    std[0] = sqrt(std[0]);
    std[1] = sqrt(std[1]);
  }
}
/* 

given lon (x) and lat (y) in degrees, r in km and phi in radians,
compute the location (lon/lat, deg) of a point with azimuth phi and
distance r

*/
void calc_azi_r_lonlat(double x,double y,double r,double phi,
		       double *xloc,double *yloc)
{
  double dlon,lon1,lat1,lat,lon,d;
  lon1 = x / GGRD_PIF;		/* lon in rad */
  lat1 = y / GGRD_PIF;		/* lat in rad */
  d = (r/6371) ;	/* distance in rad */

  lat = asin(sin(lat1) * cos(d) + cos(lat1)*sin(d)*cos(phi));
  dlon=atan2(sin(phi)*sin(d)*cos(lat1),cos(d)-sin(lat1)*sin(lat));
  lon = mod( lon1 - dlon +GGRD_PI,GGRD_TWOPI) - GGRD_PI;

  *xloc = lon * GGRD_PIF;
  if(*xloc < 0)
    *xloc += 360;
  *yloc = lat * GGRD_PIF;

}
double mod(double y, double x)
{
  return  y - x*floor(y/x);
}
