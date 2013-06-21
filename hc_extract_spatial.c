#include "hc.h"
/* 


extract part of a solution of a HC run and convert to spatial


*/

int main(int argc, char **argv)
{
  int ilayer,nvsol,ndsol=0,mode,shps,loop,i1,i2,nlat,nlon,
    ivec,lc,ndata,ndata_all,ndata_d,npoints,i,j,
    poff,shps_read=0,shps_read_d=0;
  FILE *in;
  struct sh_lms *vsol=NULL,*dsol=NULL;
  struct hcs *model;
  HC_PREC zlabel;
  hc_boolean binary_in = TRUE, verbose = FALSE,read_dsol=FALSE;
  HC_PREC *data,*plm=NULL,*xpos,*xvec,lon,lat,theta,phi,xtmp[3],pvec[3],
    *xscalar;
  HC_PREC polar_base[9];
  hc_struc_init(&model);
  /* 
     deal with parameters
  */
  ilayer = 0;
  mode = 1;
  switch(argc){
  case 3:
    sscanf(argv[2],"%i",&ilayer);
    break;
  case 4:
    sscanf(argv[2],"%i",&ilayer);
    sscanf(argv[3],"%i",&mode);
    break;
  case 5:
    sscanf(argv[2],"%i",&ilayer);
    sscanf(argv[3],"%i",&mode);
    read_dsol = TRUE;
    break;
  default:
    fprintf(stderr,"%s: usage\n%s sol.file layer [mode,%i] [scalar.sol]\n\n",
	    argv[0],argv[0],mode);
    fprintf(stderr,"extracts spatial solution (velocity or stress, v) from output file sol.file\n");
    fprintf(stderr,"         if scalar.sol argument is given, will also read in a scalar for VTK output\n");
    fprintf(stderr,"layer: 1...nset\n");
    fprintf(stderr,"\tif ilayer= 1..nset, will print one layer\n");
    fprintf(stderr,"\t          -1, will select nset (the top layer)\n");
    fprintf(stderr,"\t          -2, will print all layers\n");
    fprintf(stderr,"mode: 1...4\n");
    fprintf(stderr,"\tif mode = 1, will print lon lat z v_r \n");
    fprintf(stderr,"\t          2, will print lon lat z v_theta v_phi \n");
    fprintf(stderr,"\t          3, will print lon lat z v_r v_theta v_phi\n");
    fprintf(stderr,"\t          4, will print the depth levels of all layers\n");
    fprintf(stderr,"\t          5, compute all depth levels (set ilayer=-2) and write VTK file, ASCII\n");
    fprintf(stderr,"\t          6, compute all depth levels (set ilayer=-2) and write VTK file, BINARY\n");
    exit(-1);
    break;
  }
  if((mode == 4)||(mode==5)||(mode==6))
    ilayer = -2;
  /* 
     read in velocity/traction solution
  */
  in = ggrd_open(argv[1],"r","hc_extract_spatial");
  shps_read = hc_read_sh_solution(model,&vsol,in,binary_in,verbose);
  fclose(in);
  nvsol = model->nradp2 * shps_read;
  /* 
     deal with selection
  */
  loop = 0;
  if(ilayer == -1)
    ilayer = model->nradp2;
  else if(ilayer == -2){
    ilayer = model->nradp2;
    loop =1;
  }
  if((ilayer < 1)||(ilayer > model->nradp2)){
    fprintf(stderr,"%s: ilayer (%i) out of range, use 1 ... %i\n",
	    argv[0],ilayer,model->nradp2);
    exit(-1);
  }
  /* set up layer bounds */
  if(loop){
    i1=0;i2=model->nradp2-1;
  }else{
    i1=ilayer-1;i2 = i1;
  }
  /* detect number of expansions */
  if(mode == 1){
    shps = 1;			/* r */
  }else if(mode == 2){
    shps = 2;			/* theta,phi */
  }else if((mode == 3)||(mode == 5)||(mode==6)){
    shps = 3;			/* r,theta,phi */
  }else{
    shps = 1;
  }
  if(shps > shps_read){
    fprintf(stderr,"%s: solution file only had %i expansions, mode %i requests %i\n",
	    argv[0],shps_read,mode,shps);
    exit(-1);
  }
  /* 
     density solution or other scalar
  */
  if(read_dsol){
    if((mode != 5)&&(mode != 6))
      HC_ERROR("hc_extract_spatial","error, only mode 5 and  can handle scalar input");
    in = ggrd_open(argv[4],"r","hc_extract_spatial");
    shps_read_d = hc_read_sh_solution(model,&dsol,in,binary_in,
				    verbose);
    fclose(in);
    ndsol = model->nradp2 * shps_read_d;
  }
  /* 
     
     room for spatial expansion 

  */
  npoints = (vsol+i1*shps_read)->npoints;
  if((vsol+i1*shps_read)->type != SH_RICK)
    HC_ERROR("sh_extract_spatial","SH_RICK type required");
  /* geographic set up */
  nlat = (vsol+i1*shps_read)->rick.nlat;
  nlon = (vsol+i1*shps_read)->rick.nlon;
  
  ndata =     npoints * shps ;
  ndata_d =   npoints * shps_read_d;
  ndata_all = npoints * (shps + shps_read_d);

  if((mode == 5)||(mode==6)){			/* save all layers */
    hc_vecalloc(&data,model->nradp2 * ndata_all,"hc_extract_spatial");
  }else
    hc_vecalloc(&data, ndata_all,"hc_extract_spatial");
  for(lc=0,ilayer=i1;ilayer <= i2;ilayer++,lc++){
    /* 
       output 
    */

    zlabel = HC_Z_DEPTH(model->r[ilayer]);
    switch(mode){
    case 1:
      /*  */
      if(verbose)
	fprintf(stderr,"%s: printing v_r at layer %i (depth: %g)\n",argv[0],ilayer,
		(double)zlabel);

      ivec=FALSE;sh_compute_spatial((vsol+ilayer*shps_read),ivec,TRUE,&plm,data,verbose);
      sh_print_spatial_data_to_stream((vsol+ilayer*shps_read),shps,data,TRUE,zlabel,stdout);
      break;
    case 2:
      /*  */
      if(verbose)
	fprintf(stderr,"%s: printing v_theta v_phi SHE at layer %i (depth: %g)\n",argv[0],ilayer,(double)zlabel);
      ivec=TRUE;sh_compute_spatial((vsol+ilayer*shps_read+1),ivec,TRUE,&plm,data,verbose);
      sh_print_spatial_data_to_stream((vsol+ilayer*shps_read+1),shps,data,TRUE,zlabel,stdout);
      break;
    case 3:
      if(verbose)
	fprintf(stderr,"%s: printing v_r v_theta v_phi SHE at layer %i (depth: %g)\n",argv[0],ilayer,(double)zlabel);
      ivec=FALSE;sh_compute_spatial((vsol+ilayer*shps_read),  ivec,TRUE,&plm,data,verbose); /* radial */
      ivec=TRUE; sh_compute_spatial((vsol+ilayer*shps_read+1),ivec,TRUE,&plm,(data+npoints),verbose); /* theta,phi */
      sh_print_spatial_data_to_stream((vsol+ilayer*shps_read),shps,data,TRUE,zlabel,stdout);
      break;
    case 4:
      fprintf(stdout,"%5i %11g\n",ilayer,(double)HC_Z_DEPTH(model->r[ilayer]));
      break;
    case 5:			/* compute all and store */
    case 6:
      ivec=FALSE;sh_compute_spatial((vsol+ilayer*shps_read),  ivec,TRUE,&plm,(data+lc*ndata_all),verbose); /* radial */
      ivec=TRUE; sh_compute_spatial((vsol+ilayer*shps_read+1),ivec,TRUE,&plm,(data+lc*ndata_all+npoints),verbose); /* theta,phi */
      if(read_dsol){
	if(!shps_read_d)
	  HC_ERROR("sh_extract_spatial","logic error");
	ivec=FALSE;sh_compute_spatial((dsol+ilayer*shps_read_d),ivec,TRUE,&plm,(data+lc*ndata_all+npoints*shps),verbose); /* radial */
      }
      break;
    default:
      fprintf(stderr,"%s: error, mode %i undefined\n",argv[0],mode);
      exit(-1);
      break;
    }
  }
  /* clear and exit */
  sh_free_expansion(vsol,nvsol);
  if(read_dsol)
    sh_free_expansion(dsol,ndsol);
  free(plm);
  /*  */
  if((mode == 5)||(mode==6)){
    /* 
       print the already stored properties
    */
    if(shps != 3)HC_ERROR("hc_extract_spatial","shps has to be 3 for mode 5 and 6");
    /* convert */
    hc_vecalloc(&xpos,model->nradp2 * ndata,"hc_extract_spatial"); 
    hc_vecalloc(&xvec,model->nradp2 * ndata,"hc_extract_spatial");
    if(read_dsol)
      hc_vecalloc(&xscalar,model->nradp2 * ndata_d,"hc_extract_spatial");
    for(i=0;i < npoints;i++){	/* loop through all points */
      /* lon lat coordinates */
      sh_get_coordinates((vsol+i1*3),i,&lon,&lat);
      theta = LAT2THETA(lat);phi = LON2PHI(lon);
      xtmp[0] = xtmp[1] = sin(theta);
      xtmp[0] *= cos(phi);	/* x */
      xtmp[1] *= sin(phi);	/* y */
      xtmp[2] = cos(theta);	/* z */
      /* for conversion */
      calc_polar_base_at_theta_phi(theta,phi,polar_base);
      for(ilayer=0;ilayer < model->nradp2;ilayer++){
	/* this is the slow data storage loop but it avoids
	   recomputing the polar basis vector */
	poff = ilayer * ndata + i*shps;	/* point offset */
	for(j=0;j < 3;j++){
	  xpos[poff+j]   = xtmp[j] * model->r[ilayer]; /* cartesian coordinates */
	}
	/* data are stored a bit weirdly, this makes for lots of
	   jumping around in memory ... */
	pvec[0] = data[ilayer*ndata_all +           i];
	pvec[1] = data[ilayer*ndata_all + npoints  +i];
	pvec[2] = data[ilayer*ndata_all + npoints*2+i];
	lonlatpv2cv_with_base(pvec,polar_base,(xvec+poff));
	/* assign scalar fata if any */
	for(j=0;j < shps_read_d;j++)
	  xscalar[j * model->nradp2 * ndata_d + ilayer * npoints  + i] = 
	    data[ilayer * ndata_all + npoints*(shps+j) + i];
      }
    }
    free(data);
    /* print in VTK format */
    hc_print_vtk(stdout,xpos,xvec,npoints,model->nradp2,(mode==6),
		 shps_read_d,xscalar,nlon,nlat);
    free(xvec);free(xpos);
    if(shps_read_d)
      free(xscalar);
    
  }else{
    free(data);
  }

  return 0;
}
