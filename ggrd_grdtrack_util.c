/* 
   
   subroutines for 3-D interpolation of scalar data in GMT grd
   filesbased on GMT3.4.3 grdtrack


   $Id: ggrd_grdtrack_util.c,v 1.6 2006/02/16 02:18:03 twb Exp twb $


   original comments for grdtrack from GMT at bottom of file
 
*/

#include "hc_ggrd.h"

#ifndef ONEEIGHTYOVERPI
#define ONEEIGHTYOVERPI  57.295779513082320876798154814105
#endif
#include <math.h>
#include <string.h>
#include <math.h>

#ifndef irint
#define irint(x) ((int)rint(x))
#endif


void ggrd_init_master(struct ggrd_master *ggrd)
{
  ggrd->mat_control = ggrd->mat_control_init = FALSE;
  ggrd->ray_control = ggrd->ray_control_init = FALSE;
  ggrd->vtop_control = ggrd->vtop_control_init = FALSE;
  ggrd->age_control = ggrd->age_control_init = FALSE;
  ggrd->nage = 0;
  ggrd->age_bandlim = 200.;
  ggrd->sf_init = FALSE;
  ggrd->time_hist.init = FALSE;
  ggrd->temp.init = ggrd->use_temp = FALSE;
  ggrd->comp.init = ggrd->use_comp = FALSE;
  ggrd->time_hist.vstage_transition = 0.1; /* in Ma, transition */
  ggrd->time_hist.interpol_time_lin = FALSE;
  ggrd->time_hist.called = FALSE;
  /* 3-D velocity settings  */
  ggrd_init_vstruc(ggrd);
}


/* 

wrapper

*/



/* for debugging, mostly */

void ggrd_grdinfo(char *filename)
{
  struct ggrd_gt g[1];
  char cdummy='c';
  ggrd_grdtrack_init_general(FALSE,filename,&cdummy,
			     "",g,2,FALSE,FALSE);
  fprintf(stderr,"ggrd_grdinfo: %s W: %g E: %g S: %g N: %g\n",
	  filename,g->grd->x_min,g->grd->x_max,g->grd->y_min,g->grd->y_max);
  
}

/* 
   init structure and files for the general case, this is a wrapper
   for what's below
   
   is_three: TRUE/FALSE for 3-D/2-D
   grdfile: filename for 2-D and prefix for 3-D
   depth_file: file with depth layers for 3-D
   edge_info_string: as in GMT 
   g: ggrd_gt structure, pass allocated
   use_nearneighbor: FALSE, will use bilinear interpolation TRUE: use nearneighbor
   
   returns error code, 0 for regular execution

*/
int ggrd_grdtrack_init_general(ggrd_boolean is_three,
			       char *grdfile, char *depth_file,
			       char *gmt_edgeinfo_string,
			       struct ggrd_gt *g,
			       ggrd_boolean verbose,
			       ggrd_boolean change_z_sign,
			       ggrd_boolean use_nearneighbor)
{
  /* this is a switch and can be left in */
  //int pad[4];			/* GMT < 4.5.1 */
  GMT_LONG pad[4];			/* GMT >= 4.5.1 */
  int i,j;
  float zavg,tmp;
#ifndef USE_GMT3
  GMT_LONG interpolant;
  if(use_nearneighbor)
    interpolant = BCR_NEARNEIGHBOR; /* no interpolation */
  else
    interpolant = BCR_BILINEAR; /* bilinear is default for
				   interpolation */
#else  /* GMT 3 */
  ggrd_boolean interpolant = TRUE;
  if(use_nearneighbor){
    fprintf(stderr,"ggrd_grdtrack_init_general: nearneighbor not implemented for GMT3\n");
    return 3;
  }
#endif
  /* clear all entries */
  g->east = g->west = g->south = g->north = 0.0;
  g->is_three = is_three;
  g->init = FALSE;
  if(ggrd_grdtrack_init(&g->west,&g->east,
			&g->south,&g->north,&g->f,&g->mm,
			grdfile,&g->grd,&g->edgeinfo,
			gmt_edgeinfo_string,&g->geographic_in,
			pad,is_three,depth_file,&g->z,&g->nz,
			interpolant,verbose,change_z_sign,
			g->loc_bcr))
    return 2;
  /* 

  check bandlimited maximum with positivity constraint

  */
  g->fmaxlim = (float *)malloc(sizeof(float)*g->nz);
  for(i=0;i < g->nz;i++){	/* loop through layers */
    //g->fmaxlim[i] = g->grd[i].z_min;
    g->fmaxlim[i] = 0.0;
    for(j=0;j < g->mm;j++){	/* loop trough entries */
      tmp = fabs(g->f[i*g->mm+j]);
	//if((g->f[i*g->mm+j] < g->bandlim) &&(g->f[i*g->mm+j] > g->fmaxlim[i]))
      if((tmp < g->bandlim) && tmp > g->fmaxlim[i])
	g->fmaxlim[i] = tmp;
    }
    /* min: g->grd[i].z_min
       max: g->grd[i].z_max,
       bandlim_max: g->fmaxlim[i] */
    //fprintf(stderr,"%g %g %g %g\n", g->grd[i].z_min,g->grd[i].z_max,g->bandlim,g->fmaxlim[i]);
  }
  if(is_three){
    /* 
       check how the depth levels are specified for debugging 
    */
    zavg = 0.0;
    for(i=0;i < g->nz;i++)
      zavg += g->z[i];
    if(zavg > 0)
      g->zlevels_are_negative = FALSE;
    else
      g->zlevels_are_negative = TRUE;
  }else{
    g->zlevels_are_negative = FALSE;
  }
  //  if(change_z_sign)		/* reverse logic */
  //g->zlevels_are_negative = (g->zlevels_are_negative)?(FALSE):(TRUE);
  if(verbose){
    fprintf(stderr,"ggrd_grdtrack_init_general: initialized from %s, %s, bcflag: %s.\n",
	    grdfile,(is_three)?("3-D"):("1-D"),	
	    gmt_edgeinfo_string);
    if(is_three){
      fprintf(stderr,"ggrd_grdtrack_init_general: depth file %s, %i levels, %s.\n",
	      depth_file,g->nz,
	      (g->zlevels_are_negative)?("negative z levels"):("positive z levels"));
    }
  }
  g->init = TRUE;
  return 0;
}


/* 

take log10, 10^x and/or scale complete grid. log10 and 10^x applies first

 */
int ggrd_grdtrack_rescale(struct ggrd_gt *g,
			  ggrd_boolean take_log10, /* take log10() */
			  ggrd_boolean take_power10, /* take 10^() */
			  ggrd_boolean rescale, /* rescale? */
			  double scale	/* factor for rescaling */)
{
  int i,j,k;
  if(!g->init){
    fprintf(stderr,"ggrd_grdtrack_rescale: error: ggrd not initialized\n");
    return 1;
  }
  for(i=0;i < g->nz;i++){	/* loop through depths */
    k = i*g->mm;
    for(j=0;j < g->mm;j++,k++){
      if(take_log10){
	g->f[k] = log10(g->f[k]);
      }
      if(take_power10){
	g->f[k] = pow(10.0,g->f[k]);
      }
      if(rescale)
	g->f[k] *= scale;
    }
  }
  return 0;
}


/* 

for 3-D spherical

   interpolation wrapper, uses r, theta, phi input. return value and TRUE if success,
   undefined and FALSE else
*/
ggrd_boolean ggrd_grdtrack_interpolate_rtp(double r,double t,double p,
					   struct ggrd_gt *g,
					   double *value,
					   ggrd_boolean verbose,
					   ggrd_boolean shift_to_pos_lon,
					   double radius_planet_in_km)
{
  double x[3];
  ggrd_boolean result;
  if(!g->init){			/* this will not necessarily work */
    fprintf(stderr,"ggrd_grdtrack_interpolate_rtp: error, g structure not initialized\n");
    return FALSE;
  }
  if(!g->is_three){
    fprintf(stderr,"ggrd_grdtrack_interpolate_rtp: error, g structure is not 3-D\n");
    return FALSE;
  }
  /* 
     convert coordinates to lon / lat / z
  */
  x[0] = p * ONEEIGHTYOVERPI; /* lon */
  if(shift_to_pos_lon){
    /* make sure we are in 0 ... 360 system ? */
    if(x[0]<0)
      x[0]+=360.0;
    if(x[0]>=360)
      x[0]-=360.0;
  }
  x[1] = 90.0 - t * ONEEIGHTYOVERPI; /* lat */

  x[2] = (1.0-r) * radius_planet_in_km;	/* depth in [km] */

  if(g->zlevels_are_negative)	/* adjust for depth */
    x[2] = -x[2];

  result = ggrd_grdtrack_interpolate(x,TRUE,g->grd,g->f,
				     g->edgeinfo,g->mm,g->z,
				     g->nz,value,verbose,
				     g->loc_bcr);
  return result;
}

/* 

this is almost redundant, use lon lat in degrees and z in [km] depth

*/
ggrd_boolean ggrd_grdtrack_interpolate_lonlatz(double lon,double lat,double z,
					       struct ggrd_gt *g,
					       double *value,
					       ggrd_boolean verbose)
{
  double x[3];
  ggrd_boolean result;
  if(!g->init){			/* this will not necessarily work */
    fprintf(stderr,"ggrd_grdtrack_interpolate_rtp: error, g structure not initialized\n");
    return FALSE;
  }
  if(!g->is_three){
    fprintf(stderr,"ggrd_grdtrack_interpolate_rtp: error, g structure is not 3-D\n");
    return FALSE;
  }
  /* 
     convert coordinates to lon / lat / z
  */
  x[0] = lon;
  x[1] = lat;
  x[2] = z;
  if(g->zlevels_are_negative)	/* adjust for depth */
    x[2] = -x[2];
  
  result = ggrd_grdtrack_interpolate(x,TRUE,g->grd,g->f,
				     g->edgeinfo,g->mm,g->z,
				     g->nz,value,verbose,
				     g->loc_bcr);
  return result;
}


/* 
   for 3-D  cartesian

   interpolation wrapper, uses x, y, z input. return value and TRUE if success,
   undefined and FALSE else
   this mean lon lat z
*/
ggrd_boolean ggrd_grdtrack_interpolate_xyz(double x,double y,
					   double z,
					   struct ggrd_gt *g,
					   double *value,
					   ggrd_boolean verbose)
{
  double xloc[3];
  ggrd_boolean result;
  if(!g->init){			/* this will not necessarily work */
    fprintf(stderr,"ggrd_grdtrack_interpolate_xyz: error, g structure not initialized\n");
    return FALSE;
  }
  if(!g->is_three){
    fprintf(stderr,"ggrd_grdtrack_interpolate_xyz: error, g structure is not 3-D\n");
    return FALSE;
  }
  /* 
     convert coordinates
  */
  xloc[0] = x; /* lon, x */
  xloc[1] = y; /* lat, y */
  xloc[2] = z;	/* depth, z*/
  if(g->zlevels_are_negative)	/* adjust for depth */
    xloc[2] = -xloc[2];
  result = ggrd_grdtrack_interpolate(xloc,TRUE,g->grd,g->f,
				     g->edgeinfo,g->mm,g->z,
				     g->nz,value,verbose,
				     g->loc_bcr);
  return result;
}

/* 

for 2-D spherical 

interpolation wrapper, uses theta, phi input. 
return value and TRUE if success,
undefined and FALSE else

*/
ggrd_boolean ggrd_grdtrack_interpolate_tp(double t,double p,
					  struct ggrd_gt *g,
					  double *value,
					  ggrd_boolean verbose,
					  ggrd_boolean shift_to_pos_lon)
{
  double x[3];
  ggrd_boolean result;
  if(!g->init){			/* this will not necessarily work */
    fprintf(stderr,"ggrd_grdtrack_interpolate_tp: error, g structure not initialized\n");
    return FALSE;
  }
  if(g->is_three){
    fprintf(stderr,"ggrd_grdtrack_interpolate_tp: error, g structure is not 2-D\n");
    return FALSE;
  }
  /* 
     convert coordinates
  */
  x[0] = p * ONEEIGHTYOVERPI; /* lon */
  if(shift_to_pos_lon){
    if(x[0] < 0)
      x[0] += 360.0;
    if(x[0] >= 360.0)
      x[0]-=360.0;
  }
  x[1] = 90.0 - t * ONEEIGHTYOVERPI; /* lat */
  x[2] = 1.0;
  result = ggrd_grdtrack_interpolate(x,FALSE,g->grd,g->f,
				     g->edgeinfo,g->mm,g->z,g->nz,
				     value,verbose,g->loc_bcr);
  return result;
}

/* 
   for 2-D cartesian

   interpolation wrapper, uses x,y input
   return value and TRUE if success,
   undefined and FALSE else
*/
ggrd_boolean ggrd_grdtrack_interpolate_xy(double xin,double yin,
					   struct ggrd_gt *g,
					   double *value,
					   ggrd_boolean verbose)
{
  double x[3];
  ggrd_boolean result;
  if(!g->init){			/* this will not necessarily work */
    fprintf(stderr,"ggrd_grdtrack_interpolate_xy: error, g structure not initialized\n");
    return FALSE;
  }
  if(g->is_three){
    fprintf(stderr,"ggrd_grdtrack_interpolate_xy: error, g structure is not 2-D\n");
    return FALSE;
  }
  x[0] = xin;
  x[1] = yin;
  x[2] = 0.0;
  result = ggrd_grdtrack_interpolate(x,FALSE,g->grd,g->f,g->edgeinfo,
				     g->mm,g->z,g->nz,value,verbose,
				     g->loc_bcr);
  //fprintf(stderr,"%g %g %g\n",x[0],x[1],*value);
  return result;
}

/* 

free structure

*/
void ggrd_grdtrack_free_gstruc(struct ggrd_gt *g)
{
  free(g->grd);
  free(g->edgeinfo);
  free(g->f);
  if(g->is_three)
    free(g->z);
}
/* 

given a location vector in spherical theta, phi system (xp[3]) and a
Cartesian rotation vector wx, wy, wz (omega[3]), find the spherical
velocities vr, vtheta,vphi

*/
void ggrd_find_spherical_vel_from_rigid_cart_rot(double *vr,
						 double *vtheta,
						 double *vphi,
						 double *xp,
						 double *omega)
{
  double vp[3],polar_base[3][3],ct,cp,st,sp,xc[3],tmp,vc[3];
  int i;
  /* cos and sin theta and phi */
  ct=cos(xp[1]);cp=cos(xp[2]);
  st=sin(xp[1]);sp=sin(xp[2]);
  /* convert location to Cartesian */
  tmp =  st * xp[0];
  xc[0]= tmp * cos(xp[2]);	/* x */
  xc[1]= tmp * sin(xp[2]);	/* y */
  xc[2]= ct * xp[0];		/* z */
  /* v = omega \cross r */
  vc[0] = omega[1]*xc[2] - omega[2]*xc[1];
  vc[1] = omega[2]*xc[0] - omega[0]*xc[2];
  vc[2] = omega[0]*xc[1] - omega[1]*xc[0];
  /* get basis */
  polar_base[0][0]= st * cp;polar_base[0][1]= st * sp;polar_base[0][2]= ct;
  polar_base[1][0]= ct * cp;polar_base[1][1]= ct * sp;polar_base[1][2]= -st;
  polar_base[2][0]= -sp;polar_base[2][1]= cp;polar_base[2][2]= 0.0;
  /* convert */
  for(i=0;i<3;i++){
    vp[i]  = polar_base[i][0] * vc[0];
    vp[i] += polar_base[i][1] * vc[1];
    vp[i] += polar_base[i][2] * vc[2];
  }
  vr[0] = vp[0];
  vtheta[0] = vp[1];
  vphi[0]=vp[2];
  


}
						 
/* 

initialize

*/
#ifndef USE_GMT3
int ggrd_grdtrack_init(double *west, double *east,double *south, double *north, 
			/* geographic bounds,
			   set all to zero to 
			   get the whole range from the
			   input grid files
			*/
		       float **f,	/* data, pass as empty */
		       int *mm,  /* size of data */
		       char *grdfile,	/* name, or prefix, of grd file with scalars */
		       struct GRD_HEADER **grd,	/* pass as empty */
		       struct GMT_EDGEINFO **edgeinfo, /* pass as empty */
		       char *edgeinfo_string, /* -fg/ -L type flags from GMT, can be empty */
		       ggrd_boolean *geographic_in, /* this is normally TRUE */
		       
		       //int *pad,	/* [4] array with padding (output) GMT<4.5.1*/
		       GMT_LONG *pad,
		       ggrd_boolean three_d, char *dfile, 	/* depth file name */
		       float **z,	/* layers, pass as NULL */
		       int *nz,		/* number of layers */
		       GMT_LONG interpolant, /* linear/cubic? */
		       ggrd_boolean verbose,
		       ggrd_boolean change_depth_sign, /* change the
							  sign of the
							  depth
							  levels to go from depth (>0) to z (<0) */
		       struct GMT_BCR *loc_bcr)
#else
int ggrd_grdtrack_init(double *west, double *east,
		       double *south, double *north, 
		       float **f,int *mm,char *grdfile,
		       struct GRD_HEADER **grd,
		       struct GMT_EDGEINFO **edgeinfo,
		       char *edgeinfo_string, 
		       ggrd_boolean *geographic_in,
		       int *pad,ggrd_boolean three_d, 
		       char *dfile, float **z,int *nz,		
		       ggrd_boolean interpolant,
		       ggrd_boolean verbose,
		       ggrd_boolean change_depth_sign,
		       struct BCR *loc_bcr)
#endif
{
  FILE *din;
  float dz1,dz2;
  struct GRD_HEADER ogrd;
  int i,one_or_zero,nx,ny,mx,my;
  char filename[BUFSIZ*2],*cdummy;
  static int gmt_init = FALSE;
  /* 
     deal with edgeinfo 
  */
  *edgeinfo = (struct GMT_EDGEINFO *)
    GMT_memory (VNULL, (size_t)1, sizeof(struct GMT_EDGEINFO), "ggrd_grdtrack_init");
  /* init with nonsense to avoid compiler warning */
  ogrd.x_min = ogrd.y_min =ogrd.x_max = ogrd.y_max = -100;
  ogrd.x_inc = ogrd.y_inc = -1;
  ogrd.node_offset = 0;ogrd.nx = ogrd.ny = -1;
#ifndef USE_GMT3

  if(!gmt_init){
    /* this should be OK as is. init only once globally */
    GMT_program = "ggrd";
    GMT_make_fnan (GMT_f_NaN);
    GMT_make_dnan (GMT_d_NaN);
    GMT_io_init ();/* Init the table i/o structure */
    GMT_grdio_init();
    if(strcmp(edgeinfo_string,"-fg")==0){
      GMT_io.in_col_type[GMT_X] = GMT_io.out_col_type[GMT_X] = GMT_IS_LON;
      GMT_io.in_col_type[GMT_Y] = GMT_io.out_col_type[GMT_Y] = GMT_IS_LAT;
    }
    if(strcmp(edgeinfo_string,"-fx")==0){
      GMT_io.in_col_type[GMT_X] = GMT_io.out_col_type[GMT_X] = GMT_IS_LON;
    }
    if(strcmp(edgeinfo_string,"-fy")==0){
      GMT_io.in_col_type[GMT_Y] = GMT_io.out_col_type[GMT_Y] = GMT_IS_LAT;
    }
    gmt_init = TRUE;
  }

#endif
  /* 
     init first edgeinfo (period/global?)
  */
  GMT_boundcond_init (*edgeinfo);
  /* check if geographic */
  if (strlen(edgeinfo_string)>2){ /* the boundary flag was set */
    /* parse */
    GMT_boundcond_parse (*edgeinfo, (edgeinfo_string+2));
    if ((*edgeinfo)->gn)
      *geographic_in = 1;
    else if((*edgeinfo)->nxp == -1)
      *geographic_in = 2;
    else
      *geographic_in = 0;
  }else{
    *geographic_in = 0;
  }
  if(verbose >= 2)
    if(*geographic_in)
      fprintf(stderr,"ggrd_grdtrack_init: detected geographic region from geostring: %s\n",
	      edgeinfo_string);
  
  *z = (float *) GMT_memory 
    (VNULL, (size_t)1, sizeof(float), "ggrd_grdtrack_init");
 
  if(three_d){
    /*
      
    three D part first
    
    */
    /* 
       init the layers
    */
    din = fopen(dfile,"r");
    if(!din){
      fprintf(stderr,"ggrd_grdtrack_init: could not open depth file %s\n",
	      dfile);
      return 1;
    }
    /* read in the layers */
    *nz = 0;
    dz1 = -1;
    while(fscanf(din,"%f",(*z+ (*nz))) == 1){ 
      if(change_depth_sign)
	*(*z+ (*nz)) = -(*(*z+ (*nz)));
      /* read in each depth layer */
      *z = (float *) GMT_memory ((void *)(*z), (size_t)((*nz)+2), sizeof(float), "ggrd_grdtrack_init");
      if(*nz > 0){		/* check for increasing layers */
	if(dz1 < 0){	
	  /* init first interval */
	  dz1 = *(*z+(*nz)) - *(*z+(*nz)-1);
	  dz2 = dz1;
	}else{
	  /* later intervals */
	  dz2 = *(*z+(*nz)) - *(*z+(*nz)-1);
	}
	if(dz2 <= 0.0){		/* check for monotonic increase */
	  fprintf(stderr,"%s: error: levels in %s have to increase monotonically: n: %i dz; %g\n",
		  "ggrd_grdtrack_init",dfile,*nz,dz2);
	  return 2;
	}
      }
      *nz += 1;
    }
    fclose(din);
    /* end layer init"ggrd_grdtrack_initialization */
    if(*nz < 2){
      fprintf(stderr,"%s: error: need at least two layers in %s\n",
	      "ggrd_grdtrack_init", dfile);
      return 3;
    }
    if(verbose)
      fprintf(stderr,"%s: read %i levels from %s between zmin: %g and zmax: %g\n",
	      "ggrd_grdtrack_init",*nz,dfile,*(*z+0),*(*z+(*nz)-1));
  }else{
    *nz = 1;
    *(*z) = 0.0;
    if(verbose >= 2)
      fprintf(stderr,"ggrd_grdtrack_init: single level at z: %g\n",*(*z));
  }
  /* 
     get nz grd and edgeinfo structures 
  */
  *grd = (struct GRD_HEADER *)
    GMT_memory (NULL, (size_t)(*nz), sizeof(struct GRD_HEADER), "ggrd_grdtrack_init");
  *edgeinfo = (struct GMT_EDGEINFO *)
    GMT_memory (*edgeinfo, (size_t)(*nz), sizeof(struct GMT_EDGEINFO), "ggrd_grdtrack_init");
  if(verbose >= 2)
    fprintf(stderr,"ggrd_grdtrack_init: mem alloc ok\n");
#ifndef USE_GMT3  
  /* init the header */
  GMT_grd_init (*grd,0,&cdummy,FALSE);
#endif
  if(*nz == 1){
    if(verbose >= 2)
      
#ifdef USE_GMT3		/* old */
      fprintf(stderr,"ggrd_grdtrack_init: opening single file %s, GMT3 mode\n",grdfile);
    if (GMT_cdf_read_grd_info (grdfile,(*grd))) {
      fprintf (stderr, "%s: error opening file %s\n", 
	       "ggrd_grdtrack_init", grdfile);
      return 4;
    }
   
#else  /* >=4.1.2 */
    if(verbose >= 2)
      fprintf(stderr,"ggrd_grdtrack_init: opening single file %s, GMT4 mode\n",grdfile);
    if(GMT_read_grd_info (grdfile,*grd)){
      fprintf (stderr, "%s: error opening file %s for header\n", 
	       "ggrd_grdtrack_init", grdfile);
      return 4;
    }
#endif 
  }else{
    /* loop through headers for testing purposess */
    for(i=0;i<(*nz);i++){
      sprintf(filename,"%s.%i.grd",grdfile,i+1);
#ifdef USE_GMT3
      if (GMT_cdf_read_grd_info (filename, (*grd+i))) {
	fprintf (stderr, "%s: error opening file %s (-D option was used)\n", 
		 "ggrd_grdtrack_init", filename);
	return 6;
      }
#else  /* gmt 4 */
      if (GMT_read_grd_info (filename,(*grd+i))) {
	fprintf (stderr, "%s: error opening file %s (-D option was used)\n", 
		 "ggrd_grdtrack_init", filename);
	return 6;
      }
#endif
      if(i == 0){
	/* save the first grid parameters */
	ogrd.x_min = (*grd)[0].x_min;
	ogrd.y_min = (*grd)[0].y_min;
	ogrd.x_max = (*grd)[0].x_max;
	ogrd.y_max = (*grd)[0].y_max;
	ogrd.x_inc = (*grd)[0].x_inc;
	ogrd.y_inc = (*grd)[0].y_inc;
	ogrd.node_offset = (*grd)[0].node_offset;
	ogrd.nx = (*grd)[0].nx;
	ogrd.ny = (*grd)[0].ny;
	/* 
	   
	make sure we are in 0 ... 360 system

	*/
	if((ogrd.x_min < 0)||(ogrd.x_max<0)){
	  fprintf(stderr,"%s: WARNING: geographic grids should be in 0..360 lon system (found %g - %g)\n",
		  "ggrd_grdtrack_init",ogrd.x_min,ogrd.x_max);
	}
      }else{
	/* test */
	if((fabs(ogrd.x_min -  (*grd)[i].x_min)>5e-7)||
	   (fabs(ogrd.y_min -  (*grd)[i].y_min)>5e-7)||
	   (fabs(ogrd.x_max -  (*grd)[i].x_max)>5e-7)||
	   (fabs(ogrd.y_max -  (*grd)[i].y_max)>5e-7)||
	   (fabs(ogrd.x_inc -  (*grd)[i].x_inc)>5e-7)||
	   (fabs(ogrd.y_inc -  (*grd)[i].y_inc)>5e-7)||
	   (fabs(ogrd.nx    -  (*grd)[i].nx)>5e-7)||
	   (fabs(ogrd.ny    -  (*grd)[i].ny)>5e-7)||
	   (fabs(ogrd.node_offset - (*grd)[i].node_offset)>5e-7)){
	  fprintf(stderr,"%s: error: grid %i out of %i has different dimensions or setting from first\n",
		 "ggrd_grdtrack_init",i+1,(*nz));
	  return 8;
	}
      }
    }
  }
  if(verbose > 2)
    fprintf(stderr,"ggrd_grdtrack_init: read %i headers OK, grids appear to be same size\n",*nz);
  if (fabs(*west - (*east)) < 5e-7) {	/* No subset asked for , west same as east*/
    *west =  (*grd)[0].x_min;
    *east =  (*grd)[0].x_max;
    *south = (*grd)[0].y_min;
    *north = (*grd)[0].y_max;
  }
  one_or_zero = ((*grd)[0].node_offset) ? 0 : 1;
  nx = irint ( (*east - *west) / (*grd)[0].x_inc) + one_or_zero;
  ny = irint ( (*north - *south) / (*grd)[0].y_inc) + one_or_zero;
  /* real size of data */
  //nn = nx * ny;

  /* padded */
  mx = nx + 4;
  my = ny + 4;
  /* 
     get space for all layers
  */
  *mm = mx * my;

  *f = (float *) calloc((*mm) * (*nz) ,sizeof (float));
  if(!(*f)){
    fprintf(stderr,"ggrd_grdtrack_init: f memory error, mm: %i (%i by %i) by nz: %i \n",*mm,mx,my, *nz);
    return 9;
  }
  if(verbose >= 2){
    fprintf(stderr,"ggrd_grdtrack_init: mem alloc 2 ok, %g %g %g %g %i %i\n",
	    *west,*east,*south,*north,nx,ny);
  }
  /* 
     pad on sides 
  */
  pad[0] = pad[1] = pad[2] = pad[3] = 2;
  for(i=0;i < (*nz);i++){
    /* 
       loop through layers
    */
    if(i != 0)			/* copy first edgeinfo over */
      memcpy((*edgeinfo+i),(*edgeinfo),sizeof(struct GMT_EDGEINFO));
    if((*nz) == 1){
      sprintf(filename,"%s",grdfile);
    }else{			/* construct full filename */
      sprintf(filename,"%s.%i.grd",grdfile,i+1);
    }
    if (verbose) 
      fprintf(stderr,"ggrd_grdtrack_init: reading grd file %s (%g - %g (%i) %g - %g (%i); geo: %i flag: %s\n",
	      filename,*west,*east,nx,*south,*north,ny,
	      *geographic_in,edgeinfo_string);

    /* 
       read the grd files
    */
#ifndef USE_GMT3
    /* GMT 4 */
    if (GMT_read_grd (filename,(*grd+i), (*f+i* (*mm)), 
		      *west, *east, *south, *north, 
		      pad, FALSE)) {
      fprintf (stderr, "%s: error reading file %s\n", "ggrd_grdtrack_init", grdfile);
      return 10;
    }
    //fprintf(stderr,"%g %g %i %i %i %i\n",(*grd)->z_scale_factor,(*grd)->z_add_offset,nx,ny,mx,my);
#else
    /* old GMT */
    if (GMT_cdf_read_grd (filename, (*grd+i), (*f+i* (*mm)), 
			  *west, *east, *south, *north, 
			  pad, FALSE)) {
      fprintf (stderr, "%s: error reading file %s\n", "ggrd_grdtrack_init", grdfile);
      return 10;
    }
#endif
    /* 
       prepare the boundaries 
    */
    GMT_boundcond_param_prep ((*grd+i), (*edgeinfo+i));
    if(i == 0){
      
      /* 
	 Initialize bcr structure, this can be the same for 
	 all grids as long as they have the same dimensions

      */
#ifndef USE_GMT3
      GMT_bcr_init ((*grd+i), pad, interpolant,1.0,loc_bcr);
#else
      my_GMT_bcr_init ((*grd+i), pad, interpolant,loc_bcr);
#endif
     }
    /* Set boundary conditions  */
    GMT_boundcond_set ((*grd+i), (*edgeinfo+i), pad, 
		       (*f+i*(*mm)));
  } /* end layer loop */
  if(verbose){
    ggrd_print_layer_avg(*f,*z,mx,my,*nz,stderr,pad);
  }
  return 0;

}

void ggrd_print_layer_avg(float *x,float *z,int nx, int ny, 
			  int m,FILE *out,
			  GMT_LONG *pad) /* >= 4.5.1 */
			  //int *pad)
{
  int i,j,k,yl,xl,l,nxny,nxnyr;
  float *tmp;
  nxny = nx*ny;		/* size with padding */
  if(pad[0]+pad[1]+pad[2]+pad[3] == 0){
    for(i=0;i < m;i++){
      fprintf(stderr,"ggrd_grdtrack_init: layer %3i at depth %11g, mean: %11g rms: %11g\n",
	      i+1,z[i],ggrd_gt_mean((x+i*nxny),nxny),
	      ggrd_gt_rms((x+i*nxny),nxny));
    }
  }else{

    nxnyr = (nx-pad[0]-pad[1]) * (ny-pad[2]-pad[3]); /* actual data */
    tmp = (float *)malloc(sizeof(float)*nxnyr);
    if(!tmp)GGRD_MEMERROR("ggrd_print_layer_avg");
    xl = nx - pad[1];
    yl = ny - pad[2];
    for(i=0;i < m;i++){		/* loop through depths */
      for(l=0,j=pad[3];j < yl;j++)
	for(k=pad[0];k < xl;k++,l++)
	  tmp[l] = x[i*nxny + j * nx + k];
      fprintf(stderr,"ggrd_grdtrack_init: layer %3i at depth %11g, mean: %11g rms: %11g\n",
	      i+1,z[i],ggrd_gt_mean(tmp,nxnyr),
	      ggrd_gt_rms(tmp,nxnyr));
    }
    free(tmp);
  }
}


/* 

interpolate value 

 */
#ifndef USE_GMT3
ggrd_boolean ggrd_grdtrack_interpolate(double *in, /* lon/lat/z [2/3] in degrees/km */
				       ggrd_boolean three_d, /* use 3-D inetrpolation or 2-D? */
				       struct GRD_HEADER *grd, /* grd information */
				       float *f,	/* data array */
				       struct GMT_EDGEINFO *edgeinfo, /* edge information */
				       int mm, /* nx * ny */
				       float *z, /* depth layers */
				       int nz,	/* number of depth layers */
					double *value, /* output value */
				       ggrd_boolean verbose,
				       struct GMT_BCR *loc_bcr)
#else
ggrd_boolean ggrd_grdtrack_interpolate(double *in, /* lon/lat/z [2/3] in degrees/km */
				       ggrd_boolean three_d, /* use 3-D inetrpolation or 2-D? */
				       struct GRD_HEADER *grd, /* grd information */
				       float *f,	/* data array */
				       struct GMT_EDGEINFO *edgeinfo, /* edge information */
				       int mm, /* nx * ny */
				       float *z, /* depth layers */
				       int nz,	/* number of depth layers */
				       double *value, /* output value */
				       ggrd_boolean verbose,
				       struct BCR *loc_bcr
				       )
#endif
{
  int i1,i2;
  double fac1,fac2,val1,val2;
  static ggrd_boolean zwarned;	/* this should move, leave for now */
  /* If point is outside grd area, 
     shift it using periodicity or skip if not periodic. */

  *value = NAN;			/* use NAN as default so that error
				   returns don't have some number in
				   case user forgets to check */
  /* check if in bounds */
  while ( (in[1] < grd[0].y_min) && (edgeinfo[0].nyp > 0) ) 
    in[1] += (grd[0].y_inc * edgeinfo[0].nyp);
  if (in[1] < grd[0].y_min){

    return FALSE;
  }  
  while ( (in[1] > grd[0].y_max) && (edgeinfo[0].nyp > 0) )
    in[1] -= (grd[0].y_inc * edgeinfo[0].nyp);
  if (in[1] > grd[0].y_max) {
    return FALSE;
  }
  while ( (in[0] < grd[0].x_min) && (edgeinfo[0].nxp > 0) ) 
    in[0] += (grd[0].x_inc * edgeinfo[0].nxp);
  if (in[0] < grd[0].x_min) {
    return FALSE;
  }
  while ( (in[0] > grd[0].x_max) && (edgeinfo[0].nxp > 0) ) 
    in[0] -= (grd[0].x_inc * edgeinfo[0].nxp);
  if (in[0] > grd[0].x_max) {
    return FALSE;
  }
  /* 
     interpolate 
  */
  if(three_d){
    ggrd_gt_interpolate_z(in[2],z,nz,&i1,&i2,&fac1,&fac2,verbose,&zwarned);
    /* 
       we need these calls to reset the bcr.i and bcr.j counters 
       otherwise the interpolation routine would assume we have the same 
       grid
       
       TO DO:
       
       now, we still need the same grid dimensions, else a separate bcr
       variable has to be introduced
       
       
       
    */
#ifndef USE_GMT3
    val1 = GMT_get_bcr_z((grd+i1), in[0], in[1], (f+i1*mm), (edgeinfo+i1),loc_bcr);
    val2 = GMT_get_bcr_z((grd+i2), in[0], in[1], (f+i2*mm), (edgeinfo+i2),loc_bcr);
#else
    ggrd_global_bcr_assign(loc_bcr);
    val1 = GMT_get_bcr_z((grd+i1), in[0], in[1], (f+i1*mm), (edgeinfo+i1));
    ggrd_global_bcr_assign(loc_bcr);
    val2 = GMT_get_bcr_z((grd+i2), in[0], in[1], (f+i2*mm), (edgeinfo+i2));
#endif
    /*      fprintf(stderr,"z(%3i/%3i): %11g z: %11g z(%3i/%3i): %11g f1: %11g f2: %11g v1: %11g v2: %11g rms: %11g %11g\n",   */
    /* 	      i1+1,nz,z[i1],in[2],i2+1,nz,z[i2],fac1,fac2,  */
    /*        	      val1,val2,rms((f+i1*mm),mm),rms((f+i2*mm),mm));   */
    *value  = fac1 * val1;
    *value += fac2 * val2;
  }else{
    /* single layer */
#ifndef USE_GMT3
    *value = GMT_get_bcr_z(grd, in[0], in[1], f, edgeinfo,loc_bcr);
#else
    ggrd_global_bcr_assign(loc_bcr);
    *value = GMT_get_bcr_z(grd, in[0], in[1], f, edgeinfo);
#endif
  }
  if(verbose)
    fprintf(stderr,"ggrd_interpolate: lon: %g lat: %g val: %g\n",in[0],in[1],*value);
  return TRUE;
}
/*
  
  read in times for time history of velocities, if needed

  if read_thistory is TRUE:

  the input file for the tectonic stages is in format

  t_start^1 t_stop^1
  ....
  t_start^nvtimes t_stop^nvtimes
  
  expecting ascending time intervals, which should have smaller time first and have no gaps

  on return, the vtimes vector is in the format 

     t_left^1 t_mid^1 t_right^1 
     t_left^2 ...
     ...
     t_left^nvtimes t_mid^nvtimes t_right^nvtimes 
     

  else, will only init for one step
  
  returns error code 
*/
int ggrd_init_thist_from_file(struct ggrd_t *thist,
			      char *input_file,
			      ggrd_boolean read_thistory,
			      ggrd_boolean verbose)
{
  FILE *in;
  double ta,tb;
  ggrd_boolean opened_file = FALSE;
  if(thist->init){
    fprintf(stderr,"ggrd_read_time_intervals: error: already initialized\n");
    return 1;
  }
  ggrd_vecalloc(&thist->vtimes,3,"rti: 1");
  
  if(read_thistory){
    in = fopen(input_file,"r");
    if(!in){
      if(verbose)
	fprintf(stderr,"ggrd_read_time_intervals: WARNING: could not open file %s\n",
		input_file);
      opened_file = FALSE;
    }else{
      opened_file = TRUE;
    }
  }
  if(opened_file){
    thist->nvtimes = thist->nvtimes3 = 0;
    while(fscanf(in,"%lf %lf",&ta,&tb) == 2){
      thist->vtimes[thist->nvtimes3] = ta;
      thist->vtimes[thist->nvtimes3+2] = tb;
      if(thist->nvtimes > 0){
	if((*(thist->vtimes+thist->nvtimes3+2) < *(thist->vtimes+thist->nvtimes3))||
	   (*(thist->vtimes+thist->nvtimes3) < *(thist->vtimes+(thist->nvtimes-1)*3))||
	   (*(thist->vtimes+thist->nvtimes3+2) < *(thist->vtimes+(thist->nvtimes-1)*3+2))||
	   (fabs(*(thist->vtimes+(thist->nvtimes - 1)*3+2) - *(thist->vtimes+thist->nvtimes3))>5e-7)){
	  GGRD_PE("ggrd_read_time_intervals: error, expecting ascending time intervals");
	  GGRD_PE("ggrd_read_time_intervals: which should have smaller time first and have no gaps");
	}
      }
      // compute mid point
      *(thist->vtimes + thist->nvtimes3+1) =  
	(*(thist->vtimes+thist->nvtimes3) + *(thist->vtimes + thist->nvtimes3+2))/2.0;
      thist->nvtimes += 1;
      thist->nvtimes3 += 3;
      ggrd_vecrealloc(&thist->vtimes,thist->nvtimes3+3,"rti: 2");
    }
    thist->tmin = *(thist->vtimes+0);
    thist->tmax = *(thist->vtimes+ (thist->nvtimes-1) * 3 +2);
    if(!(thist->nvtimes)){
      fprintf(stderr,"ggrd_read_time_intervals: error, no times read from %s\n",
	      input_file);
      return 3;
    }else{
      if(verbose){
	fprintf(stderr,"ggrd_read_time_intervals: read %i time intervals from %s\n",
		thist->nvtimes,input_file);
	fprintf(stderr,"ggrd_read_time_intervals: t_min: %g t_max: %g\n",
		thist->tmin,thist->tmax);
      }
    }
    fclose(in);

  }else{
    /* 
       only one time step, or no file
    */
    thist->nvtimes = 1;
    thist->nvtimes3 = thist->nvtimes * 3;
    *(thist->vtimes+0) = *(thist->vtimes+1) = 
      *(thist->vtimes+2) = thist->tmin = thist->tmax= 0.0;
    if(verbose)
      fprintf(stderr,"ggrd_read_time_intervals: only one timestep / constant fields\n");
  }
  thist->called = FALSE;
  thist->init = TRUE;
  return 0;
}

/* 

use linear interpolation for z direction

*/
void ggrd_gt_interpolate_z(double z,float *za,int nz,
			   int *i1, int *i2, 
			   double *fac1, double *fac2,
			   ggrd_boolean verbose,
			   ggrd_boolean *zwarned)
{
  int nzm1;
  if((!(*zwarned))&&verbose){
    if((z<za[0])||(z>za[nz-1])){
      fprintf(stderr,"interpolate_z: WARNING: at least one z value extrapolated\n");
      fprintf(stderr,"interpolate_z: zmin: %g z: %g zmax: %g\n",
	      za[0],z,za[nz-1]);
      
      *zwarned = TRUE;
    }
  }
  nzm1 = nz-1;
  *i2 = 0;
  while((za[*i2]<z)&&(*i2 < nzm1))
    *i2 += 1;
  if(*i2 == 0)
    *i2 += 1;
  *i1 = *i2 - 1;

  *fac2 = ((z - (double)za[*i1])/((double)za[*i2]-(double)za[*i1]));
  *fac1 = 1.0 - *fac2;
}

//
//     produce inetrpolation weights 
//     given a vtimes(nvtimes * 3) vector with 
//             t_left t_mid t_right 
//     entries
//
// this routine will smooth transitions between tectonic stages
//
//  INPUT:
//
//  time: time to be interpolated to
//  dxlimit: transition width time

//
//
//        
//  OUTPUT:
//  i1, i2: indices of the two time intervals
//  f1, f2: weights of these intervals 
//
//
//     WARNING: if the routine gets called with the same time twice, will
//     not recalculate the weights
//
// will return 1 if OK, none 1 if error
//
/* 
   if thist.interpol_time_lin is set, will interpolate linearly from

   t1left
   t2mid
   t3mid
   t4mid
   ...



 */

int ggrd_interpol_time(GGRD_CPREC time,struct ggrd_t *thist,
		       int *i1,int *i2,GGRD_CPREC *f1,GGRD_CPREC *f2)
{

  GGRD_CPREC tloc,xll,dxlimit;
  int i22,i;
  /* respect stages with some smoothing */
  dxlimit = thist->vstage_transition;

  if(!thist->init){
    fprintf(stderr,"ggrd_interpol_time: error: thist is not init\n");
    return 0;
  }
  tloc = time;
  //
  //     special case: only one interval
  //
  if(thist->nvtimes == 1){
    *i1 = 0;
    *i2 = 0;
    *f1 = 1.0;
    *f2 = 0.0;
  }else{ // more than one stage
    if(!thist->called){
      /* 
	 init loop

      */
      if(thist->interpol_time_lin){
	/* linear approach */
	ggrd_vecalloc(&thist->tl,thist->nvtimes,"tinterp 1");
	thist->tl[0] = thist->vtimes[0]; /* left */
	for(i=1;i<thist->nvtimes;i++)
	  thist->tl[i] =  thist->vtimes[i*3+1]; /* center */
      }else{

	/* check if dxlimit is larger than actual stage intervals */
	for(i=0;i < thist->nvtimes;i++){
	  xll =  thist->vtimes[i*3+2] - thist->vtimes[i*3];
	  if(dxlimit > xll){
	    fprintf(stderr,"inter_vel: adjusting transition width to %g from time intervals\n",
		    xll);
	    dxlimit = xll;
	  }
	}
	thist->ntlim = thist->nvtimes-1;
	/* init some more factors */
	thist->xllimit=  -dxlimit/2.0;
	thist->xrlimit=   dxlimit/2.0;
      }
      thist->time_old = thist->vtimes[0] - 100; /* to make sure we get a
					    first fix */
      thist->called = TRUE;
    } /* end init */
    
    if(fabs(tloc - thist->time_old)>1e-7){       
      /* 
	 need to get new factors 
      */
      if(thist->interpol_time_lin){
	thist->iright_loc=1; 
	while((tloc >  thist->tl[thist->iright_loc])&&(thist->iright_loc < thist->nvtimes-1))
	  thist->iright_loc++;
	thist->ileft_loc=thist->iright_loc-1; 
	thist->f2_loc = (tloc-thist->tl[thist->ileft_loc])/(thist->tl[thist->iright_loc] - thist->tl[thist->ileft_loc]);
	thist->f1_loc = 1-thist->f2_loc;
      }else{
	//
	//     obtain the time intervals to the left and right of time
	//     
	if(tloc < thist->vtimes[0]){
	  if(fabs(tloc - thist->vtimes[0]) < 0.1){
	    tloc = thist->vtimes[0];
	  }else{
	    fprintf(stderr,"inter_vel: error, time: %g  too small\n",tloc);
	    fprintf(stderr,"inter_vel: while vtimes are given from %g to %g\n",
		    thist->vtimes[0],thist->vtimes[thist->nvtimes3-1]);
	    return 0;
	  }
	}
	if(tloc > thist->vtimes[thist->nvtimes3-1]){
	  if(fabs(tloc - thist->vtimes[thist->nvtimes3-1]) < 0.1){
	    tloc = thist->vtimes[thist->nvtimes3-1];
	  }else{
	    fprintf(stderr,"inter_vel: error, time: %g too large\n",tloc);
	    fprintf(stderr,"inter_vel: while vtimes are given from %g to %g\n",
		    thist->vtimes[0],thist->vtimes[thist->nvtimes3-1]);
	    return 0;
	  }
	}
	thist->iright_loc=0;  // right interval index
	i22=1;         // right interval midpoint
	//     find the right interval such that its midpoint is larger or equal than time
	while((tloc > thist->vtimes[i22]) && 
	      (thist->iright_loc < thist->ntlim)){
	  thist->iright_loc++;
	  i22 += 3;
	}
	if(thist->iright_loc == 0){
	  thist->iright_loc = 1;
	  i22 = 4;
	}
	thist->ileft_loc = thist->iright_loc - 1; // left interval index
	//
	//     distance from right boundary of left interval 
	//     (=left boundary of right interval) 
	// normalized version
	//xll = 2.0 * (tloc - thist->vtimes[i22-1])/(thist->vtimes[i22] - thist->vtimes[ileft_loc*3-2]);
	xll = tloc - thist->vtimes[i22-1]; /* real time version */
	//
	//     this will have xll go from -0.5 to 0.5 around the transition between stages
	//     which is at xl1=0
	//
	//     vf1_loc and vf2_loc are the weights for velocities within the left and right
	//     intervals respectively
	//
	if(xll < thist->xllimit){ // xllimit should be 1-dx, dx~0.1
	  thist->f1_loc = 1.0;
	  thist->f2_loc = 0.0;
	}else{
	  if(xll > thist->xrlimit){ // xrlimit should be 1+dx, dx~0.1
	    thist->f1_loc = 0.0;
	    thist->f2_loc = 1.0;
	  }else {            // in between 
	    xll =     (xll-thist->xllimit)/dxlimit; // normalize by transition width
	    thist->f2_loc = ((1.0 - cos(xll * GGRD_PI))/2.0); // this goes from 0 to 1
	    //     weight for left velocities
	    thist->f1_loc = 1.0 - thist->f2_loc;
	  }
	}
	//fprintf(stderr,"%g %g %g %i %i %g %g\n",time,f1_loc,f2_loc,ileft_loc,iright_loc,thist->vtimes[3*ileft_loc+2],thist->vtimes[3*iright_loc]);

      }	/* end the stage branch */
      thist->time_old = tloc;
    }
    //     assign the local variables to the output variables
    //     
    *i1 = thist->ileft_loc;
    *i2 = thist->iright_loc;
    *f1 = thist->f1_loc;
    *f2 = thist->f2_loc;
  } // end ntimes>1 part
  return 1;
}

/* 

interpolate seafloor age at location vt,vp and time age

this is a scalar interpolation that needs special care

*/
int interpolate_seafloor_ages(GGRD_CPREC xt, GGRD_CPREC xp,
			      GGRD_CPREC age,struct ggrd_master *ggrd, 
			      GGRD_CPREC *seafloor_age)
{
  int left, right,i;
  GGRD_CPREC f1,f2;
  double a1,a2;
  static ggrd_boolean shift_to_pos_lon = FALSE;

  if(!ggrd->sf_init){
    /* 

    init the constants
    
    */
    ggrd->sf_ntlim = ggrd->nage - 1;
    if(!ggrd->nage){
      fprintf(stderr,"interpolate_seafloor_ages: not initialized?!\n");
      return -2;
    }
    for(i=1;i < ggrd->nage;i++)
      if(ggrd->age_time[i] < ggrd->age_time[i-1]){
	fprintf(stderr,"interpolate_seafloor_ages: error: times need to be sorted monotnically increasing\n");
	return -3;
      }
    ggrd->sf_old_age =  ggrd->age_time[0] - 1000; /* so that we get interpolation
					       factors  */
    ggrd->sf_init = TRUE;
  } 

  /* check range */
  if((age < ggrd->age_time[0]) || (age > ggrd->age_time[ggrd->sf_ntlim])){
    fprintf(stderr,"interpolate_seafloor_ages: age: %g out of bounds [%g;%g]\n",
	    age,ggrd->age_time[0],ggrd->age_time[ggrd->sf_ntlim]);
    return -3;
  }

  if(fabs(age- ggrd->sf_old_age) > 1e-8){
    /* 
       time interpolation 
    */
    right = 0;

    while((right < ggrd->sf_ntlim) && (ggrd->age_time[right] < age))
      right++;
    if(right == 0)
      right++;

    left = right - 1;
    f2 = (age - ggrd->age_time[left])/
      (ggrd->age_time[right]-ggrd->age_time[left]);
    f1 = 1.0-f2;
    //fprintf(stderr,"sai: %g %g %g\t%g %g\n",ggrd->age_time[left],age,ggrd->age_time[right],f1,f2);
    ggrd->sf_old_age = age;
    ggrd->sf_old_left = left;ggrd->sf_old_right = right;
    ggrd->sf_old_f1 = f1;ggrd->sf_old_f2 = f2;
  }else{			/* reuse time interpolation */
    left = ggrd->sf_old_left;right = ggrd->sf_old_right;
    f1 = ggrd->sf_old_f1;f2 = ggrd->sf_old_f2;
  }
  if(!ggrd_grdtrack_interpolate_tp((double)xt,(double)xp,
				   (ggrd->ages+left),&a1,FALSE,shift_to_pos_lon)){
    fprintf(stderr,"interpolate_seafloor_ages: interpolation error left\n");
    return -6;
  }
  if(a1 > ggrd->ages[left].fmaxlim[0]) /* limit to bandlim */
    a1 = ggrd->ages[left].bandlim;
  if(!ggrd_grdtrack_interpolate_tp((double)xt,(double)xp,
				   (ggrd->ages+right),&a2,FALSE,shift_to_pos_lon)){
    fprintf(stderr,"interpolate_seafloor_ages: interpolation error right\n");
    return -6;
  }
  if(a2 > ggrd->ages[right].fmaxlim[0]) 
    a2 = ggrd->ages[right].bandlim;
  *seafloor_age = (GGRD_CPREC)(f1 * a1 + f2 * a2);
  if(*seafloor_age < 0)
    *seafloor_age = 0.0;	/* and >= 0 */
  return 0;
}


/* 
   open a file safely and give an error message if there was
   a problem
*/
FILE *ggrd_open(char *name, char *mode, char *program)
{
  FILE *in;
  if((in=fopen(name,mode)) == NULL){
    fprintf(stderr,"%s: error: can not open file %s for mode %s access\n",
	    program,name,mode);
    exit(-1);
  }
  return in;
}

/* general floating point vector allocation */
void ggrd_vecalloc(double **x,int n,char *message)
{
  *x = (double *)malloc(sizeof(double)*(size_t)n);
  if(! (*x)){
    fprintf(stderr,"mem error: %s\n",message);
    exit(-1);
  }
}

/* general version */
void ggrd_vecrealloc(double **x,int n,char *message)
{
  *x = (double *)realloc(*x,sizeof(double)*(size_t)n);
  if(!(*x)){
    fprintf(stderr,"mem error: %s\n",message);
    exit(-1);
  }
}

/*

  calculate mean and standard deviation of x_i
  if hypoth is set, calc mean and stddev of sqrt(x_i^2 + y_i^2)
  if weighted is set, uses weights[n] for weighting the mean and so on

  this way of calculating the stddev is inaccurate but fast

*/
void ggrd_calc_mean_and_stddev(GGRD_CPREC *x, GGRD_CPREC *y,int n,GGRD_CPREC *mean,
			       GGRD_CPREC *stddev,GGRD_CPREC *rms, 
			       ggrd_boolean hypoth, ggrd_boolean weighted,GGRD_CPREC *weight)
{
  GGRD_CPREC sum1=0.0,sum2=0.0,tmp,ws;
  int i,nuse;
  if(n <= 1){
    fprintf(stderr,"ggrd_calc_mean_and_stddev: error: n: %i\n",n);
    exit(-1);
  }
  ws=0.0;
  if(hypoth){// sqrt(x^2+y+2)
    for(i=nuse=0;i < n;i++){
      if(finite(x[i])&&finite(y[i])){
	if(weighted){
	  tmp = hypot(x[i],y[i]) * weight[i];
	  ws += weight[i];
	}else{
	  tmp = hypot(x[i],y[i]);
	  ws += 1.0;
	}
	sum1 += tmp;sum2 += tmp * tmp;
	nuse++;
      }
    }
  }else{
    for(i=nuse=0;i < n;i++){
      if(finite(x[i])){
	if(weighted){
	  tmp  = x[i] * weight[i];
	  ws += weight[i];
	}else{
	  tmp = x[i];
	  ws += 1.0;
	}
	sum1 += tmp;sum2 += tmp*tmp;
	nuse++;
      }
    }
  }
  // standard deviation
  tmp = (ws * sum2 - sum1 * sum1) / (ws*(ws-1.0));
  if(tmp > 0)
    *stddev = sqrt(tmp);
  else
    *stddev = 0.0;
  *rms  = sqrt(sum2 / ws);// RMS 
  *mean = sum1 / ws;      // mean 
}
/*

  sort array by order, returns index

  GIVE INPUT INDEX SHIFTED
  (ie.    x-1)
  output will be 0...n-1

  based on numerical recipes

  
*/


#define GGRD_INDEXX_SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define GGRD_INDEXX_M 7
#define GGRD_INDEXX_NSTACK 5000

void ggrd_indexx(int n,GGRD_CPREC *arr, int *indx)
{
  int i,indxt,ir=n,itemp,j,k,l=1;
  int jstack=0,*istack;
  GGRD_CPREC a;
  //
  istack=(int *)malloc(sizeof(int)*(GGRD_INDEXX_NSTACK+1));
  if(!istack)
    GGRD_MEMERROR("indexx");
  //
  for (j=1;j<=n;j++) 
    indx[j]=j;
  for (;;) {
    if (ir-l < GGRD_INDEXX_M) {
      for (j=l+1;j<=ir;j++) {
	indxt=indx[j];
	a=arr[indxt];
	for (i=j-1;i>=1;i--) {
	  if (arr[indx[i]] <= a) break;
	  indx[i+1]=indx[i];
	}
	indx[i+1]=indxt;
      }
      if (jstack == 0) break;
      ir=istack[jstack--];
      l=istack[jstack--];
    } else {
      k=(l+ir) >> 1;
      GGRD_INDEXX_SWAP(indx[k],indx[l+1]);
      if (arr[indx[l+1]] > arr[indx[ir]]) {
	GGRD_INDEXX_SWAP(indx[l+1],indx[ir]);
      }
      if (arr[indx[l]] > arr[indx[ir]]) {
	GGRD_INDEXX_SWAP(indx[l],indx[ir]);
      }
      if (arr[indx[l+1]] > arr[indx[l]]) {
	GGRD_INDEXX_SWAP(indx[l+1],indx[l]);
      }
      i=l+1;
      j=ir;
      indxt=indx[l];
      a=arr[indxt];
      for (;;) {
	do i++; while (arr[indx[i]] < a);
	do j--; while (arr[indx[j]] > a);
	if (j < i) break;
	GGRD_INDEXX_SWAP(indx[i],indx[j]);
      }
      indx[l]=indx[j];
      indx[j]=indxt;
      jstack += 2;
      if (jstack > GGRD_INDEXX_NSTACK) {
	fprintf(stderr,
		"indexx: GGRD_INDEXX_NSTACK too small");
	exit(-1);
      }
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
  free(istack);
  // go back to normal numbering
  for(i=1;i<=n;i++)
    indx[i]--;
}
#undef GGRD_INDEXX_M
#undef GGRD_INDEXX_NSTACK
#undef GGRD_INDEXX_SWAP



/* 


WARNING: those are in accurate and do not take the padding into account!


 */
/* compute simple RMS */
float ggrd_gt_rms(float *x,int n)
{
  int i,nuse;
  float rms = 0.0;
  for(i=nuse=0;i<n;i++){
    if(finite(x[i])){
      rms += x[i]*x[i];
      nuse++;
    }
  }
  return sqrt(rms/(float)nuse);
}
/* compute simple mean */
float ggrd_gt_mean(float *x,int n)
{
  int i,nuse;
  float mean = 0.0;
  for(i=nuse=0;i<n;i++){
    if(finite(x[i])){
      mean += x[i];
      nuse++;
    }
  }
  return mean/(float)nuse;
}
#ifdef USE_GMT3
/* 

this is aweful, but works?


*/
void ggrd_global_bcr_assign(struct BCR *loc_bcr)
{
  /* copy to global */
  memcpy((void *)(&bcr),(void *)loc_bcr,sizeof(struct BCR));
}


void my_GMT_bcr_init (struct GRD_HEADER *grd, int *pad, 
		      int interpolant,struct BCR *loc_bcr)
{

  GMT_bcr_init(grd,pad,interpolant);
  /* assign to local bcr */
  memcpy((void *)loc_bcr,(void *)(&bcr),sizeof(struct BCR));
}



#endif
/*


*	Id: grdtrack.c,v 1.5.4.6 2003/04/17 22:39:25 pwessel Exp 
*
*	Copyright (c) 1991-2003 by P. Wessel and W. H. F. Smith
*	See COPYING file for copying and redistribution conditions.
*
*	This program is free software; you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation; version 2 of the License.
*
*	This program is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*
*	Contact info: gmt.soest.hawaii.edu
*--------------------------------------------------------------------

 * grdtrack reads a xyfile, opens the 2d binary gridded grdfile, 
 * and samples the dataset at the xy positions with a bilinear or bicubic
 * interpolant.  This new data is added to the input as an extra column
 * and printed to standard output.  In order to evaluate derivatives along
 * the edges of the grdfile region, we assume natural bicubic spline
 * boundary conditions (d2z/dn2 = 0, n being the normal to the edge;
 * d2z/dxdy = 0 in the corners).  Rectangles of size x_inc by y_inc are 
 * mapped to [0,1] x [0,1] by affine transformation, and the interpolation
 * done on the normalized rectangle.
 *
 * Author:	Walter H F Smith
 * Date:	23-SEP-1993
 * 
 * Based on the original grdtrack, which had this authorship/date/history:
 *
 * Author:	Paul Wessel
 * Date:	29-JUN-1988
 * Revised:	5-JAN-1990	PW: Updated to v.2.0
 *		4-AUG-1993	PW: Added -Q
 *		14-AUG-1998	PW: GMT 3.1
 *  Modified:	10 Jul 2000 3.3.5  by PW to allow plain -L to indicate geographic coordinates
 * Version:	3.4.3
 */
