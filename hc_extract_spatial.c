#include "hc.h"
/* 


extract part of a solution of a HC run and convert to spatial



*/

int main(int argc, char **argv)
{
  int ilayer,nsol,mode,shps,loop,i1,i2,ivec,lc,ndata,npoints,i,j,
    poff,ilayer3;
  FILE *in;
  struct sh_lms *sol=NULL;
  struct hcs *model;
  HC_PREC zlabel;
  hc_boolean binary_in = TRUE, verbose = FALSE;
  float *data,*plm=NULL,*xpos,*xvec,lon,lat,theta,phi,xtmp[3],
    pvec[3];
  double polar_base[9];
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
  default:
    fprintf(stderr,"%s: usage\n%s sol.file layer [mode,%i] \n\n",
	    argv[0],argv[0],mode);
    fprintf(stderr,"extracts spatial solution x (velocity or stress) from HC run\n");
    fprintf(stderr,"layer: 1...nset\n");
    fprintf(stderr,"\tif ilayer= 1..nset, will print one layer\n");
    fprintf(stderr,"\t          -1, will select nset (the top layer)\n");
    fprintf(stderr,"\t          -2, will print all layers\n");
    fprintf(stderr,"mode: 1...4\n");
    fprintf(stderr,"\tif mode = 1, will print lon lat z v_r \n");
    fprintf(stderr,"\t          2, will print lon lat z v_theta v_phi \n");
    fprintf(stderr,"\t          3, will print lon lat z v_r v_theta v_phi\n");
    fprintf(stderr,"\t          4, will print the depth levels of all layers\n");
    fprintf(stderr,"\t          5, compute all depth levels and write VTK file, ASCII\n");
    fprintf(stderr,"\t          6, compute all depth levels and write VTK file, BINARY\n");
    exit(-1);
    break;
  }
  if((mode == 4)||(mode==5)||(mode==6))
    ilayer = -2;
  /* 
     read in solution
  */
  in = ggrd_open(argv[1],"r","hc_extract_sh_layer");
  hc_read_sh_solution(model,&sol,in,binary_in,verbose);
  fclose(in);
  nsol = model->nradp2 * 3;
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
  /* room for spatial expansion */
  npoints = (sol+i1*3)->npoints;
  ndata = npoints * shps;
  if((mode == 5)||(mode==6))			/* save all layers */
    hc_svecalloc(&data,model->nradp2 * ndata,"hc_extract_spatial");
  else
    hc_svecalloc(&data,ndata,"hc_extract_spatial");
  for(lc=0,ilayer=i1,ilayer3=ilayer*3;ilayer <= i2;ilayer++,lc++,ilayer3+=3){
    /* 
       output 
    */

    zlabel = HC_Z_DEPTH(model->r[ilayer]);
    switch(mode){
    case 1:
      /*  */
      if(verbose)
	fprintf(stderr,"%s: printing v_r at layer %i (depth: %g)\n",argv[0],ilayer,zlabel);

      ivec=FALSE;sh_compute_spatial((sol+ilayer3),ivec,TRUE,&plm,data,verbose);
      sh_print_spatial_data_to_file((sol+ilayer3),shps,data,TRUE,zlabel,stdout);
      break;
    case 2:
      /*  */
      if(verbose)
	fprintf(stderr,"%s: printing v_theta v_phi SHE at layer %i (depth: %g)\n",argv[0],ilayer,zlabel);
      ivec=TRUE;sh_compute_spatial((sol+ilayer3+1),ivec,TRUE,&plm,data,verbose);
      sh_print_spatial_data_to_file((sol+ilayer3+1),shps,data,TRUE,zlabel,stdout);
      break;
    case 3:
      if(verbose)
	fprintf(stderr,"%s: printing v_r v_theta v_phi SHE at layer %i (depth: %g)\n",argv[0],ilayer,zlabel);
      ivec=FALSE;sh_compute_spatial((sol+ilayer3),  ivec,TRUE,&plm,data,verbose); /* radial */
      ivec=TRUE; sh_compute_spatial((sol+ilayer3+1),ivec,TRUE,&plm,(data+npoints),verbose); /* theta,phi */
      sh_print_spatial_data_to_file((sol+ilayer3),shps,data,TRUE,zlabel,stdout);
      break;
    case 4:
      fprintf(stdout,"%5i %11g\n",ilayer,HC_Z_DEPTH(model->r[ilayer]));
      break;
    case 5:			/* compute all and store */
    case 6:
      ivec=FALSE;sh_compute_spatial((sol+ilayer3),  ivec,TRUE,&plm,(data+lc*ndata),verbose); /* radial */
      ivec=TRUE; sh_compute_spatial((sol+ilayer3+1),ivec,TRUE,&plm,(data+lc*ndata+npoints),verbose); /* theta,phi */
      break;
    default:
      fprintf(stderr,"%s: error, mode %i undefined\n",argv[0],mode);
      exit(-1);
      break;
    }
  }
  /* clear and exit */
  sh_free_expansion(sol,nsol);
  free(plm);
  /*  */
  if((mode == 5)||(mode==6)){
    if(shps != 3)HC_ERROR("hc_extract_spatial","shps has to be 3 for mode 5 and 6");
    /* convert */
    hc_svecalloc(&xpos,model->nradp2 * ndata,"hc_extract_spatial"); /* ndata = npoints * 3 */
    hc_svecalloc(&xvec,model->nradp2 * ndata,"hc_extract_spatial");
    for(i=0;i < npoints;i++){	/* loop through all points */
      /* lon lat coordinates */
      sh_get_coordinates((sol+i1*3),i,&lon,&lat);
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
	for(j=0;j<3;j++){
	  xpos[poff+j]   = xtmp[j] * model->r[ilayer]; /* cartesian coordinates */
	}
	/* data are stored a bit weirdly, this makes for lots of
	   jumping around in memory ... */
	pvec[0] = data[ilayer*ndata+i];
	pvec[1] = data[ilayer*ndata+npoints+i];
	pvec[2] = data[ilayer*ndata+npoints*2+i];
	lonlatpv2cv_with_base(pvec,polar_base,(xvec+poff));
      }
    }
    free(data);
    /* print in VTK format */
    hc_print_vtk(stdout,xpos,xvec,npoints,model->nradp2,(mode==6));
    free(xvec);free(xpos);
  }else{
    free(data);
  }

  return 0;
}
