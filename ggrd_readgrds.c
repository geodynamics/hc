#include "hc.h"
#include "gmt.h"
/*

  read in velocities from a set of GMT grd files named 
  
  vr.i.grd, vt.i.grd, and vp.i.grd
  
  or 1/vr.i.grd 1/... through n/vr.i.grd ..
  
  where i runs from 1 to N, and N is the number of lines in the file dfilename
  which has the depth of each layer in positive numbers in units of km

  and the 1/ ... n/ directory mode is chosen if velocity fields at different
  times are specified

  the following variables refer to the model structure, eg. vr means MDP->vr

  on return, vr, vt, and vp will hold the velocities in r, theta, and phi direction
  on n[R] layers with n[HC_PHI] and n[HC_THETA] points in longitudinal and 
  latitudinal direction, resp

  n[5]: n[R], n[HC_PHI], n[HC_THETA], n[TPPROD] which is n[HC_PHI]*n[HC_THETA]
  and n[NRNTNP] which is n[R]*n[HC_THETA]*n[R]

  r will hold the radial coordinate of each layer in ascending order

  all theta/phi coordinates of the grids will be from 0+dtheta/2 .. Pi-dtheta/2 
   and 0 .. 2Pi-dphi
  
  dtheta = Pi/n[HC_THETA]
  dphi =  2Pi/n[HC_PHI]
  
  velscale is the scaling velocity (divide vel by this number), output
  is in cm/yr NOT deg


  $Id: ggrd_readgrds.c,v 1.5 2006/01/22 01:11:34 becker Exp becker $

  input_mode determines if we read GMT grd files or our own double prec binary 
  format

*/

/* init a v structure  */
void ggrd_init_vstruc(struct ggrd_vel *v)
{
  v->read_gmt = TRUE;		/* read GMT by default */
  v->amode = GGRD_NORMAL;
  v->init = FALSE;
  v->history = FALSE;
  v->vr = v->vt = v->vp = NULL;
  v->thist.init = FALSE;
  v->velscale =  1.0; 
  v->rcmb = GGRD_RCMB_ND;
}

/* 
   read velocities 
*/
void ggrd_read_vel_grids(struct ggrd_vel *v, /* velocity structure,
						should be initialized first
						
					     */
			 GGRD_CPREC scale, /* divide all velocities by this 
					      factor */
			 hc_boolean verbose, /* verbosity level */
			 hc_boolean zero_boundary_vr, /* zero out top and 
							 bottom layers 
							 radial velocity 
						      */
			 char *prefix /* start filenames with this
					 prefix */
			 )
{
  FILE *in,*out;
  int i,j,k,l,level,os,os1,dummy[4]={0,0,0,0},ivt,*index;
  hc_boolean 
    init = FALSE,
    wraparound = FALSE,
    pixelreg = FALSE,
    weighted = TRUE;
  char sname[GGRD_STRLEN],suffix[50],loc_prefix[50],
    vsfile_loc[GGRD_STRLEN],tfilename[GGRD_STRLEN];
  float *fgrd;
  double *dgrd;
  GGRD_CPREC minphi,mintheta,omaxphi,maxtheta,std[4],rms[4],
    mean[4],ddummy,*weights,theta,tmp=0.0;
  /* gmt  */
  struct GRD_HEADER header[1];

  in = out = NULL;
  fgrd = NULL;dgrd = NULL;
  minphi = FLT_MAX;
  omaxphi=FLT_MIN;
  mintheta=FLT_MAX;
  maxtheta = FLT_MIN;
  weights=NULL;

  v->velscale = scale;
  if(fabs(v->velscale) < HC_EPS_PREC){
    fprintf(stderr,"ggrd_read_vel_grids: error: velocity scale is zero\n");
    exit(-1);
  }
  if(!v->init){
    //

    // read time intervals for velocities from file
    sprintf(tfilename,"%s%s",prefix,GGRD_THFILE);
    ggrd_read_time_intervals(&v->thist,tfilename,
			     v->history,verbose);
    //
    // read depth layers on which velocities are specified from files
    // this also creates a sorting array
    //
    sprintf(vsfile_loc,"%s%s",prefix,GGRD_DFILE);
    ggrd_read_depth_levels(v,&index,vsfile_loc,verbose);
    /*
      
      read the velocities in binary format, either GMT grd or double bin
      
    */
    if(v->n[HC_R] < 1)
      GGRD_PE("ggrd_read_vel_grids: error: should have more than one layer, check depth file");
    if(v->read_gmt){
      /* prepare filenames  */
      if(verbose)
	fprintf(stderr,"ggrd_read_vel_grids: reading grd files\n");
      strcpy(suffix,"grd");
    }else{
      if(verbose)
	fprintf(stderr,"ggrd_read_vel_grids: reading bin files\n");
      strcpy(suffix,"bin");
    }
    if(v->amode == GGRD_ONLY_VEL_STATS){
      sprintf(vsfile_loc,"%s.%s",prefix,GGRD_VSFILE);
      fprintf(stderr,"ggrd_read_vel_grids: writing z rms_vr rms_vt rms_vp rms_vh to %s\n",
	      vsfile_loc);
      out = hc_open(vsfile_loc,"w","ggrd_read_vel_grids");
    }
    for(ivt=0;ivt < v->thist.nvtimes;ivt++){
      if(v->history)
	fprintf(stderr,"ggrd_read_vel_grids: reading velocities for time [%12g, %12g] from %3i/\n",
		v->thist.vtimes[ivt*3],v->thist.vtimes[ivt*3+2],ivt+1);
      for(i=0;i < v->n[HC_R];i++){
	//
	// determine number of grd file based on resorted arrays
	//
	level = index[i]+1;// level numbers should go from 1 .. N 
	for(j=0;j<3;j++){
	  if(v->history)
	    sprintf(loc_prefix,"%i/",ivt+1);
	  else
	    sprintf(loc_prefix,"./");
	  // filenames
	  if(j==0)
	    sprintf(sname,"%s%svr.%i.%s",
		    prefix,loc_prefix,level,suffix);
	  else if(j==1)
	    sprintf(sname,"%s%svt.%i.%s",
		    prefix,loc_prefix,level,suffix);
	  else
	    sprintf(sname,"%s%svp.%i.%s",
		    prefix,loc_prefix,level,suffix);
	  if(v->read_gmt){
	    if(GMT_cdf_read_grd_info (sname,header) == -1){
	      fprintf(stderr,"ggrd_read_vel_grids: error opening GMT grd file %s\n",sname);
	      exit(-1);
	    }
	  }else{
	    in = hc_open(sname,"r","ggrd_read_vel_grids");
	    // read header type of information
	    header->node_offset=FALSE;
	    fread(&header->x_min, sizeof(double), 1, in);
	    fread(&header->x_max, sizeof(double), 1, in);
	    fread(&header->y_min, sizeof(double), 1, in);
	    fread(&header->y_max, sizeof(double), 1, in);
	    fread(&header->x_inc, sizeof(double), 1, in);
	    fread(&header->y_inc, sizeof(double), 1, in);
	    fread(&header->nx, sizeof(int), 1, in);
	    fread(&header->ny, sizeof(int), 1, in);
	  }
	  if(!init){
	    /* 
	       obtain grid dimensions and check if they are the way we
	       like it, ie.  lon lat such that
	       
	       0 <= phi <= 2pi-dphi and 0+dtheta/2<=theta<=Pi-dtheta/2 */
	    pixelreg=(header->node_offset ? TRUE : FALSE);
	    minphi=  LON2PHI(header->x_min+(pixelreg?header->x_inc/2.0:0.0));	
	    omaxphi= LON2PHI(header->x_max-(pixelreg?header->x_inc/2.0:0.0));
	    maxtheta=LAT2THETA(header->y_min+(pixelreg?header->y_inc/2.0:0.0));
	    mintheta=LAT2THETA(header->y_max-(pixelreg?header->y_inc/2.0:0.0));
	    v->dphi=  DEG2RAD( header->x_inc);
	    v->dtheta=DEG2RAD( header->y_inc);
	    if(HC_DIFFERENT(minphi,0.0) || 
	       HC_DIFFERENT(mintheta,v->dtheta*0.5) || 
	       HC_DIFFERENT(maxtheta,GGRD_PI - v->dtheta*0.5) || 
	       (HC_DIFFERENT(omaxphi,TWOPI) && 
		HC_DIFFERENT(omaxphi,TWOPI - v->dphi))){
	      fprintf(stderr,"ggrd_read_vel_grids: expecting 0/360(or %g)/%g/%g range, problem with %s\n",
		      360-RAD2DEG(v->dphi),-90+RAD2DEG(v->dtheta*0.5),
		      90-RAD2DEG(v->dtheta*0.5),sname);
	      fprintf(stderr,"ggrd_read_vel_grids: expected range in radians: t: %g/%g p: %g/%g\n",
		      mintheta,maxtheta,minphi,omaxphi);
	      fprintf(stderr,"ggrd_read_vel_grids: expected range in degrees: %g/%g/%g/%g\n",
		      PHI2LON(minphi),PHI2LON(omaxphi),
		      THETA2LAT(maxtheta),THETA2LAT(mintheta));
	      
	      fprintf(stderr,"ggrd_read_vel_grids: xy extreme: %g %g %g %g\n",
		      header->x_min,header->x_max,header->y_min,header->y_max);
	      exit(-1);
	    }
	    //
	    // check if we should throw away double entries at 0 and 360
	    if(!HC_DIFFERENT(omaxphi,TWOPI)){
	      v->n[HC_PHI] = header->nx - 1;
	      wraparound = TRUE;
	    }else{
	      v->n[HC_PHI] = header->nx;
	      wraparound = FALSE;
	    }
	    v->n[HC_THETA] = header->ny;
	    if(HC_DIFFERENT(v->dtheta,GGRD_PI /
			 ((GGRD_CPREC)(v->n[HC_THETA])))||
	       HC_DIFFERENT(v->dphi,TWOPI/
			 ((GGRD_CPREC)(v->n[HC_PHI])))){
	      fprintf(stderr,"ggrd_read_vel_grids: spacing error: ndx/dx phi: %g/%g theta: %g/%g\n",
		      TWOPI/v->n[HC_PHI],v->dphi,
		      GGRD_PI/v->n[HC_THETA],v->dtheta);
	      exit(-1);
	    }
	    //
	    // set auxiliary grid dimensions
	    //
	    v->n[HC_TPPROD] = v->n[HC_THETA]  * v->n[HC_PHI];// ny * nx
	    v->n[HC_NRNTNP] = v->n[HC_TPPROD] * v->n[HC_R];  // ny * nx * nr
	    os = v->n[HC_NRNTNP] * v->thist.nvtimes;//              ny * nx * nr *nt
	    //
	    // allocate space
	    hc_vecalloc(&v->vr,os,"ggrd_readgrds: vr");
	    hc_vecalloc(&v->vt,os,"ggrd_readgrds: vt");
	    hc_vecalloc(&v->vp,os,"ggrd_readgrds: vp");
	    if(v->read_gmt){
	      // this has to be of the original GRD file size
	      // NOT the new grid dimensions
	      fgrd = (float  *)malloc(sizeof(float)  * header->nx * header->ny);
	    }else{
	      dgrd = (double *)malloc(sizeof(double) * header->nx * header->ny);
	    }
	    if((v->read_gmt && !fgrd) ||(!v->read_gmt && !dgrd))
	      HC_MEMERROR("ggrd_read_vel_grids: velocity fields:");
	    if(weighted){
	      //
	      // need to construct 2-D array with area weights
	      //
	      hc_vecalloc(&weights,v->n[HC_TPPROD],"readgrds");
	      for(theta=mintheta,
		    k=0;k < v->n[HC_THETA];k++,theta += v->dtheta){
		tmp = sin(theta);
		for(l=0;l < v->n[HC_PHI];l++)
		  weights[k*v->n[HC_PHI]+l] = tmp;
	      }
	    }
	    if(verbose)
	      fprintf(stderr,"ggrd_read_vel_grids: x: %g/%g/%g nx: %i y: %g/%g/%g ny: %i wrap: %i v_c: %g\n",
		      PHI2LON(minphi),PHI2LON(omaxphi),
		      RAD2DEG(v->dphi),v->n[HC_PHI],
		      THETA2LAT(maxtheta),THETA2LAT(mintheta),
		      RAD2DEG(v->dtheta),v->n[HC_THETA],wraparound,
		      v->velscale);
	    init = TRUE;
	  }else{
	    if(HC_DIFFERENT(minphi,LON2PHI(header->x_min+(pixelreg?header->x_inc/2.0:0.0)))||
	       HC_DIFFERENT(omaxphi,LON2PHI(header->x_max-(pixelreg?header->x_inc/2.0:0.0)))||
	       HC_DIFFERENT(maxtheta,LAT2THETA(header->y_min+(pixelreg?header->y_inc/2.0:0.0)))||
	       HC_DIFFERENT(mintheta,LAT2THETA(header->y_max-(pixelreg?header->y_inc/2.0:0.0)))||
	       HC_DIFFERENT(v->dphi,DEG2RAD(header->x_inc))||
	       HC_DIFFERENT(v->dtheta,DEG2RAD( header->y_inc))){
	      fprintf(stderr,"ggrd_read_vel_grids: grd files have different size, grd: %s\n",
		      sname);
	      exit(-1);
	    }
	  }
	  if(v->read_gmt)
#ifdef USE_GMT4
	    // read the netcdf GRD file
	    GMT_cdf_read_grd (sname,header,fgrd, 0.0, 0.0, 0.0, 0.0, 
			      dummy, 0, NC_FLOAT);
#else
	    GMT_cdf_read_grd (sname,header,fgrd, 0.0, 0.0, 0.0, 0.0, 
			      dummy, 0);
#endif
	  else{
	    fread(dgrd,sizeof(double),header->nx*header->ny,in);
	    fclose(in);
	  }
	  //
	  // leave velocities in cm/yr
	  //
	  // AND: leave those pointer calculations here, since we
	  // do not initially have the size of the arrays
	  //
	  os1  = v->n[HC_NRNTNP] * ivt;
	  os1 += v->n[HC_TPPROD] * i;
	  // these should theoretically be == zero
	  if(j == HC_R){
	    //
	    // vr
	    //
	    if(zero_boundary_vr &&
	       (1.0 - v->rlevels[i] < HC_EPS_PREC)){
	      if(verbose)
		fprintf(stderr,"ggrd_read_vel_grids: WARNING: assuming level %3i is at surface and setting vr to zero\n",
			level);
	      ggrd_resort_and_check((v->vr+os1),fgrd,dgrd,v->n[HC_PHI],
				    v->n[HC_THETA],wraparound,1.0/v->velscale,
				    v->read_gmt,TRUE,0.0);
	    }else if((zero_boundary_vr)&&(v->rlevels[i] < v->rcmb)){
	      if(verbose)
		fprintf(stderr,"ggrd_read_vel_grids: WARNING: assuming level %3i is at CMB     and setting vr to zero\n",
			level);
	      ggrd_resort_and_check((v->vr+os1),fgrd,dgrd,v->n[HC_PHI],
				    v->n[HC_THETA],wraparound,1.0/v->velscale,
				    v->read_gmt,TRUE,0.0);
	    }else
	      ggrd_resort_and_check((v->vr+os1),fgrd,dgrd,v->n[HC_PHI],
				    v->n[HC_THETA],wraparound,1.0/v->velscale,
			       v->read_gmt,FALSE,ddummy);
	    hc_calc_mean_and_stddev((v->vr+os1),&ddummy,v->n[HC_TPPROD],
				    (mean+j),(std+j),(rms+j),FALSE,weighted,weights);
	  }else if(j == HC_THETA){
	    //
	    // vtheta
	    //
	    ggrd_resort_and_check((v->vt+os1),fgrd,dgrd,v->n[HC_PHI],
				  v->n[HC_THETA],wraparound,1.0/v->velscale,
				  v->read_gmt,FALSE,ddummy);
	    hc_calc_mean_and_stddev((v->vt+os1),&ddummy,v->n[HC_TPPROD],
				    (mean+j),(std+j),(rms+j),FALSE,weighted,weights);
	  }else{
	    //
	    // vphi
	    //
	    if(j != HC_PHI)
	      GGRD_PE("ggrd_read_vel_grds: index error");
	    ggrd_resort_and_check((v->vp+os1),fgrd,dgrd,v->n[HC_PHI],v->n[HC_THETA],
				  wraparound,1.0/v->velscale,
				  v->read_gmt,FALSE,ddummy);
	    hc_calc_mean_and_stddev((v->vp+os1),&ddummy,v->n[HC_TPPROD],
				    (mean+j),(std+j),(rms+j),FALSE,weighted,weights);
	    //
	    // and horizontal stats, put those in the 4th element 
	    // of mean
	    //
	    hc_calc_mean_and_stddev((v->vp+os1),(v->vt+os1),v->n[HC_TPPROD],
				    (mean+3),(std+3),(rms+3),TRUE,weighted,weights);
	  }
	}
	if(verbose)
	  fprintf(stderr,"ggrd_read_depth_levels: %13s: l: %3i i: %3i r: %9.7f z: %9.2f %s mean/RMS: vr: %9.2e/%9.2e vt: %9.2e/%9.2e vp: %9.2e/%9.2e\n",
		  sname,level,i,v->rlevels[i],
		  HC_Z_DEPTH(v->rlevels[i]),
		  (weighted?"weighted":"unweighted"),mean[HC_R]*v->velscale,
		  rms[HC_R]*v->velscale,mean[HC_THETA]*v->velscale,
		  rms[HC_THETA]*v->velscale,mean[HC_PHI]*v->velscale,
		  rms[HC_PHI]*v->velscale);
	if(v->amode == GGRD_ONLY_VEL_STATS)// velocity statistics output
	  fprintf(out,"%14.5e %14.5e %14.5e %14.5e %14.5e %5i %13.5f\n",
		  HC_Z_DEPTH(v->rlevels[i]),rms[HC_R]*v->velscale,
		  rms[HC_THETA]*v->velscale,rms[HC_PHI]*v->velscale,
		  rms[3]*v->velscale,ivt+1,
		  ((v->history)?(v->thist.vtimes[ivt*3+1]):(0.0)));
      }
    }
    /* free sorting array */
    free(index);
    if(v->read_gmt)
      free(fgrd);
    else
      free(dgrd);
    if(weighted)
      free(weights);
    if(v->amode == GGRD_ONLY_VEL_STATS){
      fclose(out);
      fprintf(stderr,"ggrd_read_vel_grids: exiting after printing vel stats\n");
      exit(0);
    }
    v->init = TRUE;
  }else{
    GGRD_PE("ggrd_read_vel_grds: error, already initialized");
  }
}

/* 

deal with some of the GMT vs. other array issues and handle NaNs

*/
void ggrd_resort_and_check(GGRD_CPREC *a,float *fb,double *db,
			   int m, int n,hc_boolean wrap,
			   GGRD_CPREC factor,hc_boolean read_gmt,
			   hc_boolean set_to_constant,
			   GGRD_CPREC constant)
{
  int i,j,nm,os1,os2,boff;
  static hc_boolean warned = FALSE;
  nm = m*n;
  if(read_gmt){
    // check for NaNs
    for(i=0;i < nm;i++)
      if(!finite(fb[i])){
	fb[i] = 0.0;
	if(!warned){
	  fprintf(stderr,"WARNING: at least one NaN entry in the data has been replaced with zero\n");
	  warned=TRUE;
	}
      }
  }else{
    // check for NaNs
    for(i=0;i<nm;i++)
      if(!finite(db[i])){
	db[i]=0.0;
	if(!warned){
	  fprintf(stderr,"WARNING: at least one NaN entry in the data has been replaced with zero\n");
	  warned=TRUE;
	}
      }
  }
  if(read_gmt){
    // see if we should average the 0 and 360 entries in b 
    // which might thus also be of dimension (m+1)*n, really
    if(wrap){
      boff = m+1;
      for(i=os1=os2=0;i<n;i++,os1+=m,os2+=boff){
	a[os1] = ((GGRD_CPREC)((fb[os2] + fb[os2+m])/2.0))*factor;
	for(j=1;j<m;j++)
	  a[os1+j] = ((GGRD_CPREC)fb[os2+j])*factor;
      }
    }else{
      for(i=os1=0;i<n;i++,os1+=m)
	for(j=0;j<m;j++)
	  a[os1+j] = ((GGRD_CPREC)fb[os1+j])*factor;
    }
  }else{// our own format, use doubles
    if(wrap){
      boff = m+1;
      for(i=os1=os2=0;i<n;i++,os1+=m,os2+=boff){
	a[os1] = ((GGRD_CPREC)((db[os2] + db[os2+m])/2.0))*factor;
	for(j=1;j<m;j++)
	  a[os1+j] = ((GGRD_CPREC)db[os2+j])*factor;
      }
    }else{
      for(i=os1=0;i<n;i++,os1+=m)
	for(j=0;j<m;j++)
	  a[os1+j] = ((GGRD_CPREC)db[os1+j])*factor;
    }
  }
  if(set_to_constant){
#ifdef DEBUG
    fprintf(stderr,"ggrd_resort_and_check: WARNING: setting this field to constant: %g\n",
	    constant);
#endif
    for(i=os1=0;i<n;i++,os1+=m)
      for(j=0;j<m;j++)
	a[os1+j] = constant;
  }
}
/* 
     
   read in depth levels from file and assign to r vector
   and create sorting index
   
*/
void ggrd_read_depth_levels(struct ggrd_vel *v,
			    int **index,char *filename,
			    hc_boolean verbose)
{
  FILE *in;
  int i;
  GGRD_CPREC *rnew;
  hc_boolean warned = FALSE;
  in=hc_open(filename,"r","ggrd_read_depth_levels");
  /* set counters */
  v->n[HC_R]=0;
  v->rlevels=(GGRD_CPREC *)malloc(sizeof(GGRD_CPREC));
  if(!v->rlevels)
    HC_MEMERROR("ggrd_read_depth_levels");
  while(fscanf(in,HC_FLT_FORMAT,(v->rlevels + v->n[HC_R]))==1){
    if(v->n[HC_R] > 1)		/* test, if sorted */
      if(fabs(v->rlevels[v->n[HC_R]] - v->rlevels[v->n[HC_R]-1]) < 1e-7)
	GGRD_PE("ggrd_read_depth_levels: error: two radii are at same level");
    if(v->rlevels[v->n[HC_R]] < 0){
      /* flip sign */
      v->rlevels[v->n[HC_R]] = -v->rlevels[v->n[HC_R]];
      if((!warned) && (verbose)){
	fprintf(stderr,"ggrd_read_depth_levels: WARNING: flipping sign of depth levels in %s\n",
		GGRD_DFILE);
	warned=TRUE;
      }
    }
    /* radius of levels */
    v->rlevels[v->n[HC_R]] = HC_ND_RADIUS(v->rlevels[v->n[HC_R]]);
    if((v->rlevels[v->n[HC_R]] > 1)||
       (v->rlevels[v->n[HC_R]] < GGRD_RCMB_ND)){
      // check for above surface or below CMB
      fprintf(stderr,"ggrd_read_depth_levels: radius %g out of range\n",v->rlevels[v->n[HC_R]]);
      exit(-1);
    }
    v->n[HC_R]++;
    v->rlevels=(GGRD_CPREC *)realloc(v->rlevels,sizeof(GGRD_CPREC)*
				       (v->n[HC_R]+1));
    if(!v->rlevels)
      HC_MEMERROR("ggrd_read_depth_levels");
  }
  fclose(in);
  // sort and create index
  *index=(int *)malloc(sizeof(int)*v->n[HC_R]);
  if(! *index)
    HC_MEMERROR("ggrd_read_depth_levels");
  hc_indexx(v->n[HC_R],(v->rlevels-1),(*index-1));
  // reassign
  rnew=(GGRD_CPREC *)malloc(sizeof(GGRD_CPREC)*v->n[HC_R]);
  for(i=0;i < v->n[HC_R];i++)
    rnew[i] = v->rlevels[(*index)[i]];
  for(i=0;i < v->n[HC_R];i++)
    v->rlevels[i] = rnew[i];
  free(rnew);
  if(verbose)
    fprintf(stderr,"ggrd_read_depth_levels: read %i levels from %s, r_min: %g r_max: %g \n",
	    v->n[HC_R],GGRD_DFILE,v->rlevels[0],v->rlevels[v->n[HC_R]-1]);
}
