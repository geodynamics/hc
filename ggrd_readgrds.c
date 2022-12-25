#include "hc.h"
#include "hc_ggrd.h"
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
int finite(double );		/* why? */
    
/* init a v structure  */
void ggrd_init_vstruc(struct ggrd_master *ggrd)
{
  /* directly velocity related */
  ggrd->v.read_gmt = TRUE;		/* read GMT by default */
  ggrd->v.init = FALSE;
  ggrd->v.history = FALSE;

  ggrd->v.vr = ggrd->v.vt = ggrd->v.vp = NULL;
  ggrd->v.velscale =  1.0; 
  ggrd->v.rcmb = GGRD_RCMB_ND;
  /* rlevels */
  ggrd->v.rlevels = NULL;
  ggrd->v.rl_warned = FALSE;
  /* interpolation structure */
  ggrd->v.vd.init = FALSE;
  ggrd->v.vd.reduce_r_stencil = FALSE;
  ggrd->v.vd.z_warned = FALSE; 
  ggrd->v.vd.w_warned = FALSE;

}

/* 
   read velocities, those are more restricted with regard to the geographic region 
   (-R0/359/-89.5/89.5 scheme for -I1/1)


   if v->use_age is set:

   initialize seafloor ages specified at the beginning of each plate
   tectonic stage. this needs nr_stages+1 files, and later ages will
   be linearly interpolated


*/
int ggrd_read_vel_grids(struct ggrd_master *ggrd, /* ggrd master structure
						     should be initialized first
						  */
			GGRD_CPREC scale, /* divide all velocities by this 
					     factor */
			hc_boolean verbose, /* verbosity level */
			hc_boolean zero_boundary_vr, /* zero out top and 
							bottom layers 
							radial velocity 
						     */
			char *prefix, /* start filenames with this
					prefix */
			ggrd_boolean use_nearneighbor)
{
  FILE *in,*out;
  int i,j,k,l,level,os,os1,ivt,*index,rcheck;

  //int dummy[4]={0,0,0,0};	/* GMT  < 4.5.1 */
  GMT_LONG dummy[4]={0,0,0,0};	/* GMT >= 4.5.1 */

  hc_boolean 
    init = FALSE,
    wraparound = FALSE,
    pixelreg = FALSE,
    weighted = TRUE;
  char sname[GGRD_CHAR_LENGTH],suffix[50],loc_prefix[50],*char_dummy=NULL,
    vsfile_loc[GGRD_CHAR_LENGTH],tfilename[GGRD_CHAR_LENGTH];
  float *fgrd;
  double *dgrd;
  GGRD_CPREC minphi,mintheta,omaxphi,maxtheta,std[4],rms[4],
    mean[4],ddummy,*weights,theta,tmp=0.0;
  /* gmt  */
  struct GRD_HEADER header[1];
#ifndef USE_GMT3
  GMT_io_init ();/* Init the table i/o structure */
  GMT_grdio_init();
  GMT_program = "g";
  GMT_make_fnan (GMT_f_NaN);
  GMT_make_dnan (GMT_d_NaN);
#endif

  in = out = NULL;
  fgrd = NULL;dgrd = NULL;
  minphi   = HC_FLT_MAX;
  omaxphi  = HC_FLT_MIN;
  mintheta = HC_FLT_MAX;
  maxtheta = HC_FLT_MIN;
  weights=NULL;

  ggrd->v.velscale = scale;
  if(fabs(ggrd->v.velscale) < HC_EPS_PREC){
    fprintf(stderr,"ggrd_read_vel_grids: error: velocity scale is zero\n");
    return(-1);
  }
  if(!ggrd->v.init){
    //

    // read time intervals for velocities from file
    sprintf(tfilename,"%s%s",prefix,GGRD_THFILE);
    /* 
       if ggrd->v.history is set, will look for different time intervals 
    */
    ggrd_init_thist_from_file(&ggrd->time_hist,tfilename,ggrd->v.history,verbose);
    if(ggrd->age_control){
      /* 

      initialize seafloor ages specified at the beginning of each 
      plate tectonic stage. this needs nr_stages+1 files, and 
      later ages will be linearly interpolated

      */
      
      if(!ggrd->v.history){
	fprintf(stderr,"ggrd_read_vel_grids: error: for ages, need history input\n");
	return(-1);
      }
      if(verbose)
	fprintf(stderr,"ggrd_read_vel_grids: expecting %i (nt) + 1 age grids\n",
		ggrd->time_hist.nvtimes);
      ggrd->nage = ggrd->time_hist.nvtimes + 1;
      
      /* important to use calloc so that some flags are set to zero */
      ggrd->ages = (struct  ggrd_gt *)calloc(ggrd->nage,
					  sizeof(struct ggrd_gt));
      ggrd->age_time = (GGRD_CPREC *)malloc(ggrd->nage*sizeof(GGRD_CPREC));
      if(!ggrd->ages || ! ggrd->age_time){
	fprintf(stderr,"ggrd_read_vel_grids: memory error\n");
	return -5;
      }
      /* 
	 read in the age grids 
      */
      for(ivt=0;ivt < ggrd->nage;ivt++){
	ggrd->ages[ivt].bandlim = ggrd->age_bandlim;
	sprintf(tfilename,"%s%i/age.grd",prefix,ivt+1);
	if(ggrd_grdtrack_init_general(FALSE,tfilename,char_dummy, /* load file */
				      "-Lx",(ggrd->ages+ivt),verbose,
				      FALSE,use_nearneighbor)){
	  fprintf(stderr,"ggrd_read_vel_grids: file error\n");
	  return -10;
	}
	if(ivt < ggrd->nage-1)	/* assign beginning of stage as time 
				   for seafloor age */
	  ggrd->age_time[ivt] = ggrd->time_hist.vtimes[ivt*3];
	else			/* end of last stage */
	  ggrd->age_time[ivt] = ggrd->time_hist.vtimes[(ivt-1)*3+2];
	if(verbose)
	  fprintf(stderr,"ggrd_read_vel_grids: read %s for seafloor age at time %g\n",
		  tfilename,ggrd->age_time[ivt]);
      }
      /* end age init */
    }
    //
    // read depth layers on which velocities are specified from files
    // this also creates a sorting array
    //
    sprintf(vsfile_loc,"%s%s",prefix,GGRD_DFILE);
    ggrd_read_depth_levels(ggrd,&index,vsfile_loc,verbose);
    /*
      
      read the velocities in binary format, either GMT grd or double bin
      
    */
    if(ggrd->v.n[HC_R] < 1)
      GGRD_PE("ggrd_read_vel_grids: error: should have more than one layer, check depth file");
    if(ggrd->v.read_gmt){
      /* prepare filenames  */
      if(verbose)
	fprintf(stderr,"ggrd_read_vel_grids: reading grd files\n");
      strcpy(suffix,"grd");
    }else{
      if(verbose)
	fprintf(stderr,"ggrd_read_vel_grids: reading bin files\n");
      strcpy(suffix,"bin");
    }
    if(ggrd->amode == GGRD_ONLY_VEL_STATS){
      sprintf(vsfile_loc,"%s.%s",prefix,GGRD_VSFILE);
      fprintf(stderr,"ggrd_read_vel_grids: writing z rms_vr rms_vt rms_vp rms_vh to %s\n",
	      vsfile_loc);
      out = ggrd_open(vsfile_loc,"w","ggrd_read_vel_grids");
    }
    for(ivt=0;ivt < ggrd->time_hist.nvtimes;ivt++){
      if((ggrd->v.history)&&(verbose))
	fprintf(stderr,"ggrd_read_vel_grids: reading velocities for time [%12g, %12g] from %3i/\n",
		ggrd->time_hist.vtimes[ivt*3],
		ggrd->time_hist.vtimes[ivt*3+2],ivt+1);
      for(i=0;i < ggrd->v.n[HC_R];i++){
	//
	// determine number of grd file based on resorted arrays
	//
	level = index[i]+1;// level numbers should go from 1 .. N 
	for(j=0;j<3;j++){
	  if(ggrd->v.history)
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
	  if(ggrd->v.read_gmt){
#ifdef USE_GMT3	    	/* old */
	    if(GMT_cdf_read_grd_info (sname,header) == -1){
	      fprintf(stderr,"ggrd_read_vel_grids: error opening GMT grd file %s\n",sname);
	      return(-2);
	    }
#else  /* new */
	    if(GMT_read_grd_info (sname,header) == -1){
	      fprintf(stderr,"ggrd_read_vel_grids: error opening GMT grd file %s\n",sname);
	      return(-2);
	    }
#endif
	  }else{
	    in = ggrd_open(sname,"r","ggrd_read_vel_grids");
	    //
	    // read header type of information
	    //
	    header->node_offset=FALSE;
	    rcheck =fread(&header->x_min, sizeof(double), 1, in);
	    rcheck+=fread(&header->x_max, sizeof(double), 1, in);
	    rcheck+=fread(&header->y_min, sizeof(double), 1, in);
	    rcheck+=fread(&header->y_max, sizeof(double), 1, in);
	    rcheck+=fread(&header->x_inc, sizeof(double), 1, in);
	    rcheck+=fread(&header->y_inc, sizeof(double), 1, in);
	    rcheck+=fread(&header->nx, sizeof(int), 1, in);
	    rcheck+=fread(&header->ny, sizeof(int), 1, in);
	    if(rcheck != 8){
	      fprintf(stderr,"ggrd_read_vel_grids: error reading header values\n");
	      return(-4);
	    }
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
	    ggrd->v.dphi=  DEG2RAD( header->x_inc);
	    ggrd->v.dtheta=DEG2RAD( header->y_inc);
	    if(HC_DIFFERENT(minphi,0.0) || 
	       HC_DIFFERENT(mintheta,ggrd->v.dtheta*0.5) || 
	       HC_DIFFERENT(maxtheta,GGRD_PI - ggrd->v.dtheta*0.5) || 
	       (HC_DIFFERENT(omaxphi,GGRD_TWOPI) && 
		HC_DIFFERENT(omaxphi,GGRD_TWOPI - ggrd->v.dphi))){
	      fprintf(stderr,"ggrd_read_vel_grids: expecting 0/360(or %g)/%g/%g range, problem with %s\n",
		      360-RAD2DEG(ggrd->v.dphi),-90+RAD2DEG(ggrd->v.dtheta*0.5),
		      90-RAD2DEG(ggrd->v.dtheta*0.5),sname);
	      fprintf(stderr,"ggrd_read_vel_grids: expected range in radians: t: %g/%g p: %g/%g\n",
		      mintheta,maxtheta,minphi,omaxphi);
	      fprintf(stderr,"ggrd_read_vel_grids: expected range in degrees: %g/%g/%g/%g\n",
		      PHI2LON(minphi),PHI2LON(omaxphi),
		      THETA2LAT(maxtheta),THETA2LAT(mintheta));
	      
	      fprintf(stderr,"ggrd_read_vel_grids: xy extreme: %g %g %g %g\n",
		      header->x_min,header->x_max,header->y_min,header->y_max);
	      return(-2);
	    }
	    //
	    // check if we should throw away double entries at 0 and 360
	    if(!HC_DIFFERENT(omaxphi,GGRD_TWOPI)){
	      ggrd->v.n[HC_PHI] = header->nx - 1;
	      wraparound = TRUE;
	    }else{
	      ggrd->v.n[HC_PHI] = header->nx;
	      wraparound = FALSE;
	    }
	    ggrd->v.n[HC_THETA] = header->ny;
	    if(HC_DIFFERENT(ggrd->v.dtheta,GGRD_PI /
			 ((GGRD_CPREC)(ggrd->v.n[HC_THETA])))||
	       HC_DIFFERENT(ggrd->v.dphi,GGRD_TWOPI/
			 ((GGRD_CPREC)(ggrd->v.n[HC_PHI])))){
	      fprintf(stderr,"ggrd_read_vel_grids: spacing error: ndx/dx phi: %g/%g theta: %g/%g\n",
		      GGRD_TWOPI/ggrd->v.n[HC_PHI],ggrd->v.dphi,
		      GGRD_PI/ggrd->v.n[HC_THETA],ggrd->v.dtheta);
	      return(-3);
	    }
	    //
	    // set auxiliary grid dimensions
	    //
	    ggrd->v.n[HC_TPPROD] = ggrd->v.n[HC_THETA]  * ggrd->v.n[HC_PHI];// ny * nx
	    ggrd->v.n[HC_NRNTNP] = ggrd->v.n[HC_TPPROD] * ggrd->v.n[HC_R];  // ny * nx * nr
	    os = ggrd->v.n[HC_NRNTNP] * ggrd->time_hist.nvtimes;//              ny * nx * nr *nt
	    //
	    // allocate space
	    ggrd_vecalloc(&ggrd->v.vr,os,"ggrd_readgrds: vr");
	    ggrd_vecalloc(&ggrd->v.vt,os,"ggrd_readgrds: vt");
	    ggrd_vecalloc(&ggrd->v.vp,os,"ggrd_readgrds: vp");
	    if(ggrd->v.read_gmt){
	      // this has to be of the original GRD file size
	      // NOT the new grid dimensions
	      fgrd = (float  *)malloc(sizeof(float)  * header->nx * header->ny);
	    }else{
	      dgrd = (double *)malloc(sizeof(double) * header->nx * header->ny);
	    }
	    if((ggrd->v.read_gmt && !fgrd) ||(!ggrd->v.read_gmt && !dgrd))
	      HC_MEMERROR("ggrd_read_vel_grids: velocity fields:");
	    if(weighted){
	      //
	      // need to construct 2-D array with area weights
	      //
	      ggrd_vecalloc(&weights,ggrd->v.n[HC_TPPROD],"readgrds");
	      for(theta=mintheta,
		    k=0;k < ggrd->v.n[HC_THETA];k++,theta += ggrd->v.dtheta){
		tmp = sin(theta);
		for(l=0;l < ggrd->v.n[HC_PHI];l++)
		  weights[k*ggrd->v.n[HC_PHI]+l] = tmp;
	      }
	    }
	    if(verbose)
	      fprintf(stderr,"ggrd_read_vel_grids: x: %g/%g/%g nx: %i y: %g/%g/%g ny: %i wrap: %i v_c: %g\n",
		      PHI2LON(minphi),PHI2LON(omaxphi),
		      RAD2DEG(ggrd->v.dphi),ggrd->v.n[HC_PHI],
		      THETA2LAT(maxtheta),THETA2LAT(mintheta),
		      RAD2DEG(ggrd->v.dtheta),ggrd->v.n[HC_THETA],wraparound,
		      ggrd->v.velscale);
	    init = TRUE;
	  }else{
	    if(HC_DIFFERENT(minphi,LON2PHI(header->x_min+(pixelreg?header->x_inc/2.0:0.0)))||
	       HC_DIFFERENT(omaxphi,LON2PHI(header->x_max-(pixelreg?header->x_inc/2.0:0.0)))||
	       HC_DIFFERENT(maxtheta,LAT2THETA(header->y_min+(pixelreg?header->y_inc/2.0:0.0)))||
	       HC_DIFFERENT(mintheta,LAT2THETA(header->y_max-(pixelreg?header->y_inc/2.0:0.0)))||
	       HC_DIFFERENT(ggrd->v.dphi,DEG2RAD(header->x_inc))||
	       HC_DIFFERENT(ggrd->v.dtheta,DEG2RAD( header->y_inc))){
	      fprintf(stderr,"ggrd_read_vel_grids: grd files have different size, grd: %s\n",
		      sname);
	      exit(-1);
	    }
	  }
	  if(ggrd->v.read_gmt){
#ifndef USE_GMT3
	    GMT_read_grd (sname,header,fgrd, 0.0, 0.0, 0.0, 0.0, 
			  dummy,0);
#else
	    GMT_cdf_read_grd (sname,header,fgrd, 0.0, 0.0, 0.0, 0.0, 
			      dummy, 0);
#endif
	  }else{
	    rcheck=fread(dgrd,sizeof(double),header->nx*header->ny,in);
	    if(rcheck!=header->nx*header->ny){
	      fprintf(stderr,"ggrd_read_vel_grids: error reading data values\n");
	      return(-5);
	    }
	    fclose(in);
	  }
	  //
	  // leave velocities in cm/yr
	  //
	  // AND: leave those pointer calculations here, since we
	  // do not initially have the size of the arrays
	  //
	  os1  = ggrd->v.n[HC_NRNTNP] * ivt;
	  os1 += ggrd->v.n[HC_TPPROD] * i;
	  // these should theoretically be == zero
	  if(j == HC_R){
	    //
	    // vr
	    //
	    if(zero_boundary_vr &&
	       (1.0 - ggrd->v.rlevels[i] < HC_EPS_PREC)){
	      if(verbose)
		fprintf(stderr,"ggrd_read_vel_grids: WARNING: assuming level %3i is at surface and setting vr to zero\n",
			level);
	      ggrd_resort_and_check((ggrd->v.vr+os1),fgrd,dgrd,ggrd->v.n[HC_PHI],
				    ggrd->v.n[HC_THETA],wraparound,1.0/ggrd->v.velscale,
				    ggrd->v.read_gmt,TRUE,0.0,&ggrd->v.vd.w_warned);
	    }else if((zero_boundary_vr)&&(ggrd->v.rlevels[i] < ggrd->v.rcmb)){
	      if(verbose)
		fprintf(stderr,"ggrd_read_vel_grids: WARNING: assuming level %3i is at CMB     and setting vr to zero\n",
			level);
	      ggrd_resort_and_check((ggrd->v.vr+os1),fgrd,dgrd,ggrd->v.n[HC_PHI],
				    ggrd->v.n[HC_THETA],wraparound,1.0/ggrd->v.velscale,
				    ggrd->v.read_gmt,TRUE,0.0,&ggrd->v.vd.w_warned);
	    }else
	      ggrd_resort_and_check((ggrd->v.vr+os1),fgrd,dgrd,ggrd->v.n[HC_PHI],
				    ggrd->v.n[HC_THETA],wraparound,1.0/ggrd->v.velscale,
			       ggrd->v.read_gmt,FALSE,ddummy,&ggrd->v.vd.w_warned);
	    ggrd_calc_mean_and_stddev((ggrd->v.vr+os1),&ddummy,ggrd->v.n[HC_TPPROD],
				    (mean+j),(std+j),(rms+j),FALSE,weighted,weights);
	  }else if(j == HC_THETA){
	    //
	    // vtheta
	    //
	    ggrd_resort_and_check((ggrd->v.vt+os1),fgrd,dgrd,ggrd->v.n[HC_PHI],
				  ggrd->v.n[HC_THETA],wraparound,1.0/ggrd->v.velscale,
				  ggrd->v.read_gmt,FALSE,ddummy,&ggrd->v.vd.w_warned);
	    ggrd_calc_mean_and_stddev((ggrd->v.vt+os1),&ddummy,ggrd->v.n[HC_TPPROD],
				    (mean+j),(std+j),(rms+j),FALSE,weighted,weights);
	  }else{
	    //
	    // vphi
	    //
	    if(j != HC_PHI)
	      GGRD_PE("ggrd_read_vel_grds: index error");
	    ggrd_resort_and_check((ggrd->v.vp+os1),fgrd,dgrd,ggrd->v.n[HC_PHI],ggrd->v.n[HC_THETA],
				  wraparound,1.0/ggrd->v.velscale,
				  ggrd->v.read_gmt,FALSE,ddummy,&ggrd->v.vd.w_warned);
	    ggrd_calc_mean_and_stddev((ggrd->v.vp+os1),&ddummy,ggrd->v.n[HC_TPPROD],
				    (mean+j),(std+j),(rms+j),FALSE,weighted,weights);
	    //
	    // and horizontal stats, put those in the 4th element 
	    // of mean
	    //
	    ggrd_calc_mean_and_stddev((ggrd->v.vp+os1),(ggrd->v.vt+os1),ggrd->v.n[HC_TPPROD],
				    (mean+3),(std+3),(rms+3),TRUE,weighted,weights);
	  }
	}
	if(verbose)
	  fprintf(stderr,"ggrd_read_depth_levels: %13s: l: %3i i: %3i r: %9.7f z: %9.2f %s mean/RMS: vr: %9.2e/%9.2e vt: %9.2e/%9.2e vp: %9.2e/%9.2e\n",
		  sname,level,i,ggrd->v.rlevels[i],
		  HC_Z_DEPTH(ggrd->v.rlevels[i]),
		  (weighted?"weighted":"unweighted"),mean[HC_R]*ggrd->v.velscale,
		  rms[HC_R]*ggrd->v.velscale,mean[HC_THETA]*ggrd->v.velscale,
		  rms[HC_THETA]*ggrd->v.velscale,mean[HC_PHI]*ggrd->v.velscale,
		  rms[HC_PHI]*ggrd->v.velscale);
	if(ggrd->amode == GGRD_ONLY_VEL_STATS)// velocity statistics output
	  fprintf(out,"%14.5e %14.5e %14.5e %14.5e %14.5e %5i %13.5f\n",
		  HC_Z_DEPTH(ggrd->v.rlevels[i]),rms[HC_R]*ggrd->v.velscale,
		  rms[HC_THETA]*ggrd->v.velscale,rms[HC_PHI]*ggrd->v.velscale,
		  rms[3]*ggrd->v.velscale,ivt+1,
		  ((ggrd->v.history)?(ggrd->time_hist.vtimes[ivt*3+1]):(0.0)));
      }
    }
    /* free sorting array */
    free(index);
    if(ggrd->v.read_gmt)
      free(fgrd);
    else
      free(dgrd);
    if(weighted)
      free(weights);
    if(ggrd->amode == GGRD_ONLY_VEL_STATS){
      fclose(out);
      fprintf(stderr,"ggrd_read_vel_grids: exiting after printing vel stats\n");
      return(0);
    }
    ggrd->v.init = TRUE;
  }else{
    GGRD_PE("ggrd_read_vel_grds: error, already initialized");
  }
  return 0;
}

/* 

deal with some of the GMT vs. other array issues and handle NaNs

*/
void ggrd_resort_and_check(GGRD_CPREC *a,float *fb,double *db,
			   int m, int n,hc_boolean wrap,
			   GGRD_CPREC factor,hc_boolean read_gmt,
			   hc_boolean set_to_constant,
			   GGRD_CPREC constant,
			   ggrd_boolean *warned)
{
  int i,j,nm,os1,os2,boff;
  nm = m*n;
  if(read_gmt){
    // check for NaNs
    for(i=0;i < nm;i++)
      if(!finite((double)fb[i])){
	fb[i] = 0.0;
	if(!(*warned)){
	  fprintf(stderr,"WARNING: at least one NaN entry in the data has been replaced with zero\n");
	  *warned=TRUE;
	}
      }
  }else{
    // check for NaNs
    for(i=0;i<nm;i++)
      if(!finite((double)db[i])){
	db[i]=0.0;
	if(!(*warned)){
	  fprintf(stderr,"WARNING: at least one NaN entry in the data has been replaced with zero\n");
	  *warned=TRUE;
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
#ifdef HC_DEBUG
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
void ggrd_read_depth_levels(struct ggrd_master *ggrd,
			    int **index,char *filename,
			    hc_boolean verbose)
{
  FILE *in;
  int i;
  GGRD_CPREC *rnew;

  in = ggrd_open(filename,"r","ggrd_read_depth_levels");
  /* set counters */
  ggrd->v.n[HC_R]=0;
  ggrd->v.rlevels=(GGRD_CPREC *)realloc(ggrd->v.rlevels,sizeof(GGRD_CPREC));
  if(!ggrd->v.rlevels)
    HC_MEMERROR("ggrd_read_depth_levels");
  while(fscanf(in,"%lf",(ggrd->v.rlevels + ggrd->v.n[HC_R]))==1){
    if(ggrd->v.n[HC_R] > 1)		/* test, if sorted */
      if(fabs(ggrd->v.rlevels[ggrd->v.n[HC_R]] - ggrd->v.rlevels[ggrd->v.n[HC_R]-1]) < 1e-7)
	GGRD_PE("ggrd_read_depth_levels: error: two radii are at same level");
    if(ggrd->v.rlevels[ggrd->v.n[HC_R]] < 0){
      /* flip sign */
      ggrd->v.rlevels[ggrd->v.n[HC_R]] = -ggrd->v.rlevels[ggrd->v.n[HC_R]];
      if((!ggrd->v.rl_warned) && (verbose)){
	fprintf(stderr,"ggrd_read_depth_levels: WARNING: flipping sign of depth levels in %s\n",
		GGRD_DFILE);
	ggrd->v.rl_warned = TRUE;
      }
    }
    /* radius of levels */
    ggrd->v.rlevels[ggrd->v.n[HC_R]] = HC_ND_RADIUS(ggrd->v.rlevels[ggrd->v.n[HC_R]]);
    if((ggrd->v.rlevels[ggrd->v.n[HC_R]] > 1)||
       (ggrd->v.rlevels[ggrd->v.n[HC_R]] < GGRD_RCMB_ND)){
      // check for above surface or below CMB
      fprintf(stderr,"ggrd_read_depth_levels: radius %g out of range\n",ggrd->v.rlevels[ggrd->v.n[HC_R]]);
      exit(-1);
    }
    ggrd->v.n[HC_R]++;
    ggrd->v.rlevels=(GGRD_CPREC *)realloc(ggrd->v.rlevels,sizeof(GGRD_CPREC)*
				       (ggrd->v.n[HC_R]+1));
    if(!ggrd->v.rlevels)
      HC_MEMERROR("ggrd_read_depth_levels");
  }
  fclose(in);
  // sort and create index
  *index=(int *)malloc(sizeof(int)*ggrd->v.n[HC_R]);
  if(! *index)
    HC_MEMERROR("ggrd_read_depth_levels");
  ggrd_indexx(ggrd->v.n[HC_R],(ggrd->v.rlevels-1),(*index-1));
  // reassign
  rnew=(GGRD_CPREC *)malloc(sizeof(GGRD_CPREC)*ggrd->v.n[HC_R]);
  for(i=0;i < ggrd->v.n[HC_R];i++)
    rnew[i] = ggrd->v.rlevels[(*index)[i]];
  for(i=0;i < ggrd->v.n[HC_R];i++)
    ggrd->v.rlevels[i] = rnew[i];
  free(rnew);
  if(verbose)
    fprintf(stderr,"ggrd_read_depth_levels: read %i levels from %s, r_min: %g r_max: %g \n",
	    ggrd->v.n[HC_R],GGRD_DFILE,ggrd->v.rlevels[0],ggrd->v.rlevels[ggrd->v.n[HC_R]-1]);
}
