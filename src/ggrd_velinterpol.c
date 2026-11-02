//
//
//    routines that deal with velocity interpolation in time and space
//
//    $Id: ggrd_velinterpol.c,v 1.4 2006/01/22 01:11:34 becker Exp $
//
//    find the velocities (and their first derivatives, if icalc_der is set)
//     with respect to r,theta,phi
//
//     dvr(0,1,2,3):     v_r, d_r v_r, d_theta v_r, d_phi v_r
//     dvtheta(0,1,2,3): v_theta, d_r v_theta, d_theta v_theta, d_phi v_theta
//     dvphi(0,1,2,3):   v_phi, d_r v_phi, d_theta v_phi, d_phi v_phi
//    
//     if icalc_der is TRUE, will compute derivatives, else only velocities
//     
//     order is the order of interpolation, e.g. 3
//
//
//     input are velocity fields vr,vtheta,vphi given laterally in dtheta/dphi
//     spaced ntheta*nphi layers at nr radial levels at rlevels[nr]
//
//
//     returns normalized velocities, scaled by v->scale
//
//     dtrange determines the time range used for transitioning
//
//
#include "hc.h"
#ifndef NO_GMT
#include "hc_ggrd.h"

int ggrd_find_vel_and_der(GGRD_CPREC *xloc,
			  GGRD_CPREC time,GGRD_CPREC dtrange,
			  struct ggrd_master *ggrd,int order,
			  hc_boolean icalc_der,
			  hc_boolean verbose,
			  GGRD_CPREC *dvr,GGRD_CPREC *dvtheta,
			  GGRD_CPREC *dvphi)
{
  GGRD_CPREC rnet,vrloc,vphiloc,vthetaloc;
  int i,k,j,m,iorder,idindex,ilim,ishift,igrid[3][GGRD_MAX_ORDERP1],index,lorder;
  GGRD_CPREC grid[GGRD_MAX_ORDERP1*3];
  struct wgt{			/* weights from polynomial interpolation */
    GGRD_CPREC w[GGRD_MAX_ORDERP1][GGRD_MAX_IORDERP1];
  };
  struct wgt weights[3];
   
  
  if(!ggrd->v.vd.init){
    //     
    //     do some checks if called for the first time
    //     
    if(order > GGRD_MAX_ORDER){
      fprintf(stderr,"ggrd_find_vel_and_der: error: order %i too large, max is %i\n",
	      order,GGRD_MAX_ORDER);
      return(-1);
    }
    ggrd->v.vd.old_order = order;
    ggrd->v.vd.orderp1 = order+1;
    if(ggrd->v.n[HC_R] < order+1){
      if(verbose)
	fprintf(stderr,"ggrd_find_vel_and_der: WARNING: reducing r stencil to nl-1: %i\n",
		ggrd->v.n[HC_R] -1 );
      ggrd->v.vd.reduce_r_stencil = TRUE;
    }
    if(ggrd->v.n[HC_PHI] < order+1){
      fprintf(stderr,"ggrd_find_vel_and_der: need at least four lon levels\n");
      fprintf(stderr,"ggrd_find_vel_and_der: using polynomial interpolation\n");
      return(-1);
    }
    if(ggrd->v.n[HC_THETA] < order+1){
      fprintf(stderr,"ggrd_find_vel_and_der: need at least four lat levels\n");
      fprintf(stderr,"ggrd_find_vel_and_der: using polynomial interpolation\n");
      return(-1);
    }
    /* test levels */
    for(i=1;i < ggrd->v.n[HC_R];i++){
      if(ggrd->v.rlevels[i] <= ggrd->v.rlevels[i-1]){
	fprintf(stderr,"ggrd_find_vel_and_der: error:\n");
	fprintf(stderr,"ggrd_find_vel_and_der: rlevels have to be ascending\n");
	fprintf(stderr,"ggrd_find_vel_and_der: i: %i r(i): %g r(i-1): %g\n",i,
		ggrd->v.rlevels[i],ggrd->v.rlevels[i-1]);
	return(-1);
      }
    }

    for(i=0;i < 3;i++){
      if(i==HC_R)
	lorder = (ggrd->v.vd.reduce_r_stencil)?(ggrd->v.n[HC_R]-1):(order);
      else
	lorder = order;
      if(lorder < 0){
	fprintf(stderr,"ggrd_find_vel_and_der: error: (reduced) order smaller tan zero: %i\n",
		lorder);
	return(-1);
      }
      /* loop through dimensions */
      /* all directions will be interpolated up to the 
	 default order + 1 */
      ggrd->v.vd.istencil[i] = lorder + 1;
      //     these are offsets for the interpolation routine
      ggrd->v.vd.isshift[i] = (int)(ggrd->v.vd.istencil[i]/2.0);
    }
    //     
    //     this is for derivatives, initialize once as zeroes
    //
    for(i=0;i < 1+GGRD_MAX_IORDER*3;i++)
      ggrd->v.vd.ider[i] = 0.0;
    ggrd->v.vd.init = TRUE;
  } /* end of init loop */
  if(order != ggrd->v.vd.old_order){
    fprintf(stderr,"ggrd_find_vel_and_der: error: order (%i) shouldn't change, old: %i\n",
	    order,ggrd->v.vd.old_order);
    return(-1);
  }
  //     
  //     for the summing up routine of weights
  //
  if(icalc_der){        // compute derivatives
    /* first derivatives */
    iorder = 1;
    idindex = 1 + iorder*3;
  }else{
    /* no derivatives */
    idindex = 1;
    iorder = 0;
  }
  if(iorder > GGRD_MAX_IORDER){
    fprintf(stderr,"ggrd_find_vel_and_der: error: dorder: %i max is %i\n",
	    iorder,GGRD_MAX_IORDER);
    return(-1);
  }
  /* 

  begin interpolation

  */

  /* 
     check range 
  */
  if(xloc[HC_PHI]<0)
    xloc[HC_PHI] += GGRD_TWOPI;
  if(xloc[HC_PHI]>GGRD_TWOPI)
    xloc[HC_PHI] -= GGRD_TWOPI;
  if((xloc[HC_R]<0) || (xloc[HC_R]>1) || (xloc[HC_THETA]<0) ||
     (xloc[HC_THETA] > GGRD_PI) || (xloc[HC_PHI] < 0) || 
     (xloc[HC_PHI] > GGRD_TWOPI)){
    fprintf(stderr,"ggrd_find_vel_and_der: coordinate x{%g, %g, %g} (lon: %g, lat: %g, z: %g) out of range\n",
	    xloc[HC_R],xloc[HC_THETA],xloc[HC_PHI],
	    PHI2LON(xloc[HC_PHI]),THETA2LAT(xloc[HC_THETA]),
	    HC_Z_DEPTH(xloc[HC_R]));
    return(-1);
  }


  //     
  //     RADIAL COMPONENT
  //     
  //     default order-1 interpolation for radial direction
  ggrd->v.vd.ixtracer[HC_R] = -1;
  ilim = ggrd->v.n[HC_R] - 1;
  i=0;
  while ((ggrd->v.vd.ixtracer[HC_R] == -1) && (i < ilim)){
    if(xloc[HC_R] <= ggrd->v.rlevels[i]){
      ggrd->v.vd.ixtracer[HC_R] = i;
      break;
    }
    i++;
  }
  if(ggrd->v.vd.ixtracer[HC_R] == -1){   // no depth levels found, tracer is above surface
    ggrd->v.vd.ixtracer[HC_R] = ilim;          // assign last layer, x_r should be corrected by the RK routines
  }
  //     
  //     pick indices of grid points for the radial stencil
  //     
  for(i=0;i < ggrd->v.vd.istencil[HC_R];i++)
    igrid[HC_R][i] = ggrd->v.vd.ixtracer[HC_R] - ggrd->v.vd.isshift[HC_R] + i;
  //     
  //     make sure all grid points exist 
  //     
  ishift = igrid[HC_R][0];
  if(ishift < 0)
    for(i=0;i < ggrd->v.vd.istencil[HC_R];i++)
      igrid[HC_R][i] -= ishift;
  //     same for upper limit
  ishift = igrid[HC_R][ggrd->v.vd.istencil[HC_R]-1] - ilim;
  if(ishift > 0)
    for(i=0;i < ggrd->v.vd.istencil[HC_R];i++)
      igrid[HC_R][i] -= ishift;
  //     find values of r for each grid point
  for(j=HC_R*(ggrd->v.vd.orderp1),i=0;i < ggrd->v.vd.istencil[HC_R];i++)
    grid[j+i] = ggrd->v.rlevels[igrid[HC_R][i]];


  //     
  //     THETA COMPONENT
  //
  ilim = ggrd->v.n[HC_THETA]-1;
  ggrd->v.vd.ixtracer[HC_THETA] = (int)(xloc[HC_THETA]/ggrd->v.dtheta);
  //     pick grid points
  for(i=0;i < ggrd->v.vd.istencil[HC_THETA];i++)
    igrid[HC_THETA][i] = ggrd->v.vd.ixtracer[HC_THETA] - ggrd->v.vd.isshift[HC_THETA]+i;
  //     
  //     adust grid points to avoid wrap-around
  //     
  ishift = igrid[HC_THETA][0];
  if(ishift < 0){
    for(i=0;i < ggrd->v.vd.istencil[HC_THETA];i++)
      igrid[HC_THETA][i] -= ishift;
  }
  //     same for upper limit
  ishift = igrid[HC_THETA][ggrd->v.vd.istencil[HC_THETA]-1] - ilim;
  if(ishift > 0)
    for(i=0;i < ggrd->v.vd.istencil[HC_THETA];i++)
      igrid[HC_THETA][i] -= ishift;
  //
  // find values of theta: since given on dtheta/2 .... pi-dtheta/2
  // theta_i = (i+0.5)*dtheta, i=0,1,...
  //
  for(j=HC_THETA*(ggrd->v.vd.orderp1),i=0;i < ggrd->v.vd.istencil[HC_THETA];i++)
    grid[j+i] = (igrid[HC_THETA][i] + 0.5) * ggrd->v.dtheta;
  //     
  //     now for phi
  //     
  ilim = ggrd->v.n[HC_PHI]-1;
  ggrd->v.vd.ixtracer[HC_PHI] = (int)(xloc[HC_PHI]/ggrd->v.dphi+.5);
  //pick grid points
  for(i=0;i < ggrd->v.vd.istencil[HC_PHI];i++){
    igrid[HC_PHI][i] = ggrd->v.vd.ixtracer[HC_PHI] - ggrd->v.vd.isshift[HC_PHI] + i;
    //     
    //     wrap around 
    //     
    if(igrid[HC_PHI][i] > ilim)
      igrid[HC_PHI][i] -= ggrd->v.n[HC_PHI];
    if(igrid[HC_PHI][i] < 0)
      igrid[HC_PHI][i] += ggrd->v.n[HC_PHI];
  }
  //
  // find values of phi. phi_i = i * dphi
  //
  for(j=HC_PHI*(ggrd->v.vd.orderp1),i=0;i < ggrd->v.vd.istencil[HC_PHI];i++)
    grid[j+i] = igrid[HC_PHI][i] * ggrd->v.dphi;
  

#ifdef HC_DEBUG
  //     
  //     check if all indices are ok
  //     
  for(i=0;i < ggrd->v.vd.istencil[HC_R];i++){
    if((igrid[HC_R][i]< 0)||(igrid[HC_R][i] >= ggrd->v.n[HC_R])){
      fprintf(stderr,"ggrd_find_vel_and_der: row %i r index %i error\n",i,igrid[HC_R][i]);
      return(-1);
    }
  }
  for(i=0;i < ggrd->v.vd.istencil[HC_THETA];i++){
    if((igrid[HC_THETA][i] < 0) || (igrid[HC_THETA][i] >= ggrd->v.n[HC_THETA])){
      fprintf(stderr,"ggrd_find_vel_and_der: row %i theta index %i error \n",i,igrid[HC_THETA][i]);
      return(-1);
    }
  }
  for(i=0;i < ggrd->v.vd.istencil[HC_PHI];i++){
    if((igrid[HC_PHI][i] < 0)||(igrid[HC_PHI][i] >= ggrd->v.n[HC_PHI])){
      fprintf(stderr,"ggrd_find_vel_and_der: row %i phi index %i error\n",
	      i,igrid[HC_PHI][i]);
      return(-1);
    }
  }
  if(idindex > 4){
    fprintf(stderr,"ggrd_find_vel_and_der: second derivatives not implemented\n");
    return(-1);
  }

  if(verbose >= 2){		/* debugging output */
    fprintf(stderr,"ggrd_velinterpol: x={%g, %g, %g} [%i (%i), %i (%i), %i(%i)]\n",
	    xloc[HC_R],xloc[HC_THETA],xloc[HC_PHI],ggrd->v.vd.ixtracer[HC_R],ggrd->v.n[HC_R],
	    ggrd->v.vd.ixtracer[HC_THETA],ggrd->v.n[HC_THETA],
	    ggrd->v.vd.ixtracer[HC_PHI],ggrd->v.n[HC_PHI]);
    for(i=0;i < 3;i++){
      rnet = GGRD_TWOPI;
      fprintf(stderr,"ggrd_velinterpol: dim: %i:",i);
      for(j=0;j < ggrd->v.vd.istencil[i];j++){
	fprintf(stderr,"%.5f (%3i) ",grid[i*(ggrd->v.vd.orderp1)+j],igrid[i][j]);
	/* find min distance to stencil point */
	vrloc = fabs(grid[i*(ggrd->v.vd.orderp1)+j]-xloc[i]);
	if(vrloc < rnet){rnet=vrloc;k=j;}
      }
      fprintf(stderr,"\tms: %.4f(%i)\n",
	      (GGRD_CPREC)k/(GGRD_CPREC)(ggrd->v.vd.istencil[i]-1),ggrd->v.vd.istencil[i]);
    }
  }
#endif 
  //     
  //     POLYNOMIAL
  //     
  //     compute all the weights for each stencil
  //
  for(i=0;i < 3;i++){		/* loop through spatial dimension */
    ggrd_weights(xloc[i],(grid+i*(ggrd->v.vd.orderp1)),ggrd->v.vd.istencil[i],iorder,weights[i].w);
  }
  //     
  //     first calculate velocities only (idindex=0) or vel and  derivatives 
  //     of velocity (e.g. v_(r,r)) if needed (idindex=3)
  // 

  /* 
     first velocities
  */
  dvr[0] = dvtheta[0] = dvphi[0] = 0.0;
  for(i=0;i < ggrd->v.vd.istencil[HC_R];i++){   // radial 
    for(j=0; j < ggrd->v.vd.istencil[HC_THETA];j++){ // theta 
      for(k=0; k < ggrd->v.vd.istencil[HC_PHI];k++){ // phi
	rnet  = weights[HC_R].w[i][0];
	rnet *= weights[HC_THETA].w[j][0];
	rnet *= weights[HC_PHI].w[k][0];
	index = igrid[HC_R][i] * ggrd->v.n[HC_TPPROD] + 
	  igrid[HC_THETA][j] * ggrd->v.n[HC_PHI] + 
	  igrid[HC_PHI][k];
	ggrd_get_velocities(&vrloc,&vthetaloc,&vphiloc,index,
			    ggrd,time,dtrange);
	dvr[0]     += rnet * vrloc;
	dvtheta[0] += rnet * vthetaloc;
	dvphi[0]   += rnet * vphiloc;
      }
    }
  }


  for(m=1;m < idindex;m++){           //m=0 -> no derivative
    dvr[m]=0.0;            //m=_R_ -> derivative wrt r
    dvtheta[m]=0.0;        //m=THETA -> derivative wrt theta
    dvphi[m]=0.0;         //m=PHI -> derivative wrt phi
    //     
    //     this is the derivative yes/no array
    //     set once, and switch off again below//
    //     
    ggrd->v.vd.ider[m] = 1;
    for(i=0;i < ggrd->v.vd.istencil[HC_R];i++){   
      for(j=0; j < ggrd->v.vd.istencil[HC_THETA];j++){ 
	for(k=0; k < ggrd->v.vd.istencil[HC_PHI];k++){ 
	  rnet  = weights[HC_R].w[i][ggrd->v.vd.ider[HC_R+1]];
	  rnet *= weights[HC_THETA].w[j][ggrd->v.vd.ider[HC_THETA+1]];
	  rnet *= weights[HC_PHI].w[k][ggrd->v.vd.ider[HC_PHI+1]];
	  index = igrid[HC_R][i] * ggrd->v.n[HC_TPPROD] + 
	    igrid[HC_THETA][j] * ggrd->v.n[HC_PHI] + 
	    igrid[HC_PHI][k];
	  ggrd_get_velocities(&vrloc,&vthetaloc,&vphiloc,index,
			      ggrd,time,dtrange);
	  dvr[m]     += rnet * vrloc;
	  dvtheta[m] += rnet * vthetaloc;
	  dvphi[m]   += rnet * vphiloc;
	}
      }
    }
    ggrd->v.vd.ider[m]=0;              // reset derivative switxh  to zero
  }
  /* succesful return */
  return 0;
}
//
//     get_velocities
//
//     
//     obtain the time-interpolated velocities while index
//     specifies the 3-D position in the arrays that are
//     vr(nrntnp*nvtimes) long. the vtimes array is nvtimes*3 and has 
//     t_left t_mid t_right for each interval in a row
//
//
//     dtrange: time range used to transition between plate tectonic stages
//
void ggrd_get_velocities(GGRD_CPREC *vrloc,GGRD_CPREC *vthetaloc,
			 GGRD_CPREC *vphiloc,
			 int index, struct ggrd_master *ggrd,
			 GGRD_CPREC time,GGRD_CPREC dtrange)
{
  int index1,i1,i2;
  GGRD_CPREC vf1,vf2;
  if((index < 0) || (index >= ggrd->v.n[HC_NRNTNP]) ){
    HC_ERROR("ggrd_get_velocities","index out of bounds");
    exit(-1);
  }
  if(ggrd->time_hist.nvtimes == 1){
    // only one time-step, steady-state calculation
    *vrloc=      ggrd->v.vr[index];
    *vthetaloc = ggrd->v.vt[index];
    *vphiloc=    ggrd->v.vp[index];
  } else {
    /* transition between stages */
    ggrd->time_hist.vstage_transition = dtrange;
    /* 
       false by default, could change that  
    */
    //ggrd->time_hist.interpol_time_lin = FALSE;
    
    // interpolate in time
    if(!ggrd_interpol_time(time,&ggrd->time_hist,&i1,&i2,&vf1,&vf2))
      exit(-1);
    if(fabs(vf1) > 1e-7){
      index1 = i1 * ggrd->v.n[HC_NRNTNP] + index;
      *vrloc=      ggrd->v.vr[index1] * vf1 ;
      *vthetaloc = ggrd->v.vt[index1] * vf1; 
      *vphiloc=    ggrd->v.vp[index1] * vf1 ;
    }else{
      *vrloc = *vthetaloc = *vphiloc = 0.0;
    }
    if(fabs(vf2) > 1e-7){
      index1 = i2 * ggrd->v.n[HC_NRNTNP] + index;
      *vrloc     += ggrd->v.vr[index1] * vf2;
      *vthetaloc += ggrd->v.vt[index1] * vf2;
      *vphiloc   += ggrd->v.vp[index1] * vf2;
    }
  }
}



void ggrd_weights(GGRD_CPREC xi,GGRD_CPREC *x,
		  int n,int m,
		  GGRD_CPREC c[GGRD_MAX_ORDERP1][GGRD_MAX_IORDERP1])
{
  //
  // copied from Fornberg(1996),p.168
  // calculates weights for 1-d interpolations
  //
  // INPUT PARAMETERS:
  // xi: point at which approximations are to be accurate
  // x : xcoords for grid points, array dimensioned to x(0:n)
  // n : # of grid points 
  // m : highest order of derivative to be approximated
  //
  // OUTPUT PARAMETER:
  // c:  weights, array dimensioned  c(0:n, 0:m)
  //     the element c(j,k) contains the weight to be applied
  //     at x(j) when the kth derivative is approximated by
  //     a stencil extending over x(0), x(1),...,x(n)
  //*********************************************

  GGRD_CPREC c1,c2,c3,c4,c5;
  int i,j,k,mn,os;
  if((n > GGRD_MAX_ORDERP1)||(m > GGRD_MAX_IORDER)){
    /* check limits */
    fprintf(stderr,"ggrd_weights: n(order+1): %i (max: %i) m(der order): %i (max: %i) out of bounds\n",
	    n,GGRD_MAX_ORDERP1,m,GGRD_MAX_IORDER);
    exit(-1);
  }
  c1 = 1.0;
  c4 = x[0] - xi;
  for(k=0;k <= m;k++)
    for(j=0;j < n;j++)
      c[j][k] = 0.0;
  
  c[0][0] = 1.0;
  for(i=1;i < n;i++){
    mn = HC_MIN(i,m);
    c2 = 1.0;
    c5 = c4;
    c4 = x[i] - xi;
    os = i - 1;
    for(j=0;j <= os;j++){
      c3 = x[i] - x[j];
      c2 *= c3;
      for(k=mn;k >= 1;k--)
	c[i][k] = c1*(((GGRD_CPREC)k)*c[i-1][k-1] - c5 * c[i-1][k])/c2;
      c[i][0] = -c1 * c5 * c[i-1][0]/c2;
      for(k=mn;k >= 1;k--)
	c[j][k] = (c4*c[j][k] - ((GGRD_CPREC)k) * c[j][k-1])/c3;
      c[j][0] = c4 * c[j][0]/c3;
    }
    c1 = c2;
  }
}

#endif
