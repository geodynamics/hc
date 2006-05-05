/*
  =====================================================================

                        CitcomS Version 1.0
                  ---------------------------------

                              Authors:
           Louis Moresi, Shijie Zhong, Lijie Han, Eh Tan,
           Clint Conrad, Michael Gurnis, and Eun-seo Choi
          (c) California Institute of Technology 1994-2003

    Agreeing to a non-commercial license is required to use this program.
            Please check http://geoframework.org for details.
               Free for non-commercial academic use ONLY.
      This program is distributed WITHOUT ANY WARRANTY whatsoever.

      =====================================================================

      Copyright January 2003, by the California Institute of Technology.
      ALL RIGHTS RESERVED. United States Government Sponsorship Acknowledged.
      
      =====================================================================


*/

#include "element_definitions.h"


#include "global_defs.h"
#include "proto.h"
#ifdef USE_GGRD
#include "hc.h"			/* Hager & O'Connell routines */
#endif

/*
minor changes in this file by TWB

$Id: Boundary_conditions.c,v 1.12 2004/12/03 02:38:25 becker Exp becker $

see the rlog of RCS
*/

/* ========================================== */

void velocity_boundary_conditions(struct All_variables *E)
{
  int node,d,j,noz,lv;
  if((E->control.vbcs_file)&&(E->mesh.topvbc != 1))
    myerror("logic error: reading vel from file but top v boundary is free-slip",E);
  if((E->mesh.topvbc==1)&&(!E->control.vbcs_file)&&
     (hypot(E->control.VBXtopval,E->control.VBYtopval)>1e-7)){
    fprintf(stderr,"WARNING: non-zero fixed surface velocities\n");
  }

  
  if(E->control.vsideflow){
    /* sideflow velocity boundary conditions */
    if(E->mesh.botvbc != 1)
      myerror("sideflow VBC but bottom free slip",E);
  }
  if(E->control.vbcs_file && E->control.vsideflow)
    myerror("sideflow BC and VBC from file not synchronized yet",E);
  for(lv=E->mesh.gridmax;lv>=E->mesh.gridmin;lv--)
    for (j=1;j<=E->sphere.caps_per_proc;j++)     {
      noz = E->lmesh.NOZ[lv];
      /* 
	 top velocity boundary conditions 
      */
      if(E->mesh.topvbc != 1) {	
	/* specified stress (free-slip for stress=0) */
	horizontal_bc(E,E->sphere.cap[j].VB,noz,1,0.0,VBX,0,lv,j);	 
	horizontal_bc(E,E->sphere.cap[j].VB,noz,3,0.0,VBZ,1,lv,j);
	horizontal_bc(E,E->sphere.cap[j].VB,noz,2,0.0,VBY,0,lv,j);	 
	horizontal_bc(E,E->sphere.cap[j].VB,noz,1,E->control.VBXtopval,SBX,1,lv,j);
	horizontal_bc(E,E->sphere.cap[j].VB,noz,3,0.0,                 SBZ,0,lv,j);
	horizontal_bc(E,E->sphere.cap[j].VB,noz,2,E->control.VBYtopval,SBY,1,lv,j);
      }else{	
	/* top specified velocities */
        horizontal_bc(E,E->sphere.cap[j].VB,noz,1,E->control.VBXtopval,VBX,1,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,noz,3,0.0,                 VBZ,1,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,noz,2,E->control.VBYtopval,VBY,1,lv,j); 
        horizontal_bc(E,E->sphere.cap[j].VB,noz,1,0.0,SBX,0,lv,j);	 
        horizontal_bc(E,E->sphere.cap[j].VB,noz,3,0.0,SBZ,0,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,noz,2,0.0,SBY,0,lv,j);	 
        if(E->control.vbcs_file)   {
	  /* read in the velocity boundary */
          read_velocity_boundary_from_file(E);
	}
	/* else{
	   /* commented out 2/26/00 CPC - remove program specification 
	     of velocity BC
	     renew_top_velocity_boundary(E);
	     }
	*/
      }
      /* 
	 bottom velocity boundary condition 
      */

      if(E->mesh.botvbc != 1) {	/* specified stress */
        horizontal_bc(E,E->sphere.cap[j].VB,1,1,0.0,VBX,0,lv,j); 
        horizontal_bc(E,E->sphere.cap[j].VB,1,3,0.0,VBZ,1,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,1,2,0.0,VBY,0,lv,j);	 
        horizontal_bc(E,E->sphere.cap[j].VB,1,1,E->control.VBXbotval,SBX,1,lv,j); 
        horizontal_bc(E,E->sphere.cap[j].VB,1,3,0.0,                 SBZ,0,lv,j);  
        horizontal_bc(E,E->sphere.cap[j].VB,1,2,E->control.VBYbotval,SBY,1,lv,j); 
      }else{	/* bottom specified velocities */
        horizontal_bc(E,E->sphere.cap[j].VB,1,1,E->control.VBXbotval,VBX,1,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,1,3,0.0,                 VBZ,1,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,1,2,E->control.VBYbotval,VBY,1,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,1,1,0.0,SBX,0,lv,j);	 
        horizontal_bc(E,E->sphere.cap[j].VB,1,3,0.0,SBZ,0,lv,j); 
        horizontal_bc(E,E->sphere.cap[j].VB,1,2,0.0,SBY,0,lv,j);	 
      }
    }    /* end for j and lv */
  if(E->control.vsideflow){
    if((E->mesh.botvbc != 1)||(E->mesh.topvbc != 1)){
      myerror("sideflow VBCs set but either top or bottom of box free-slip",E);
    }
    /* set all sides no-slip and read in velocities */
    velocity_side_bc(E);
  }else{
    /* set other boundary conditions on side of box free slip */
    velocity_refl_vert_bc(E);
  }
  if(E->control.verbose)
    for (j=1;j<=E->sphere.caps_per_proc;j++)
      for (node=1;node<=E->lmesh.nno;node++)
	fprintf(E->fp_out,"m=%d VB== %d %g %g %g flag %u %u %u\n",
		j,node,
		E->sphere.cap[j].VB[1][node],
		E->sphere.cap[j].VB[2][node],
		E->sphere.cap[j].VB[3][node],
		E->node[j][node]&VBX,
		E->node[j][node]&VBY,
		E->node[j][node]&VBZ);
  
  /* If any imposed internal velocity structure it goes here */

 
  return; 
}

/* ========================================== */

void temperature_boundary_conditions(struct All_variables *E)
{
  int j,lev,noz;
  
  lev = E->mesh.levmax;
  if(E->parallel.me == 0)
    fprintf(stderr,"temperature_boundary_conditions: mesh.levmax: %i mesh.gridmax: %i\n",
	    E->mesh.levmax, E->mesh.gridmax);
  temperature_refl_vert_bc(E);

  for (j=1;j<=E->sphere.caps_per_proc;j++)    {
    noz = E->lmesh.noz;
    /* top */
    if(E->mesh.toptbc == 1)    { /* temperature */
      horizontal_bc(E,E->sphere.cap[j].TB,noz,3,E->control.TBCtopval,TBZ,1,lev,j);
      horizontal_bc(E,E->sphere.cap[j].TB,noz,3,E->control.TBCtopval,FBZ,0,lev,j);
    }
    else   {			/* flux */
      horizontal_bc(E,E->sphere.cap[j].TB,noz,3,E->control.TBCtopval,TBZ,0,lev,j);
      horizontal_bc(E,E->sphere.cap[j].TB,noz,3,E->control.TBCtopval,FBZ,1,lev,j);
    }
    /* bottom */
    if(E->mesh.bottbc == 1)    { /* temperature */
      horizontal_bc(E,E->sphere.cap[j].TB,1,3,E->control.TBCbotval,TBZ,1,lev,j);
      horizontal_bc(E,E->sphere.cap[j].TB,1,3,E->control.TBCbotval,FBZ,0,lev,j);
      }
    else        {		/* heat flux */
      horizontal_bc(E,E->sphere.cap[j].TB,1,3,E->control.TBCbotval,TBZ,0,lev,j);
      horizontal_bc(E,E->sphere.cap[j].TB,1,3,E->control.TBCbotval,FBZ,1,lev,j); 
      }
    if((E->control.temperature_bound_adj) || 
       (E->control.lith_age_time))  {
      /* 
	 read in ages to make age < 0 flag available for boundary 
	 condition routine temperature_lith_adj TWB
      */
      get_lith_ages(E);
      /* 
	 set the regions in which to use lithosphere files to determine temperature 
	 note that this is called if the lithosphere age in inputted every time step
	 OR it is only maintained in the boundary regions 
      */
      temperature_lith_adj(E,lev);
    }
    temperatures_conform_bcs(E);
  }     /* end for j sphere cap loop */

  return; 
}

/* ========================================== 
   
this deals with the side boundary conditions, needed for the 
regional case

*/

void velocity_refl_vert_bc(E)
     struct All_variables *E;
{
  int m,i,j,ii,jj;
  int node1,node2;
  int level,nox,noy,noz;
  const int dims=E->mesh.nsd;

  /*  for two YOZ planes   */


  if ((E->parallel.me_locl[1]==0) || 
      (E->parallel.me_locl[1]==E->parallel.nprocxl-1)){
    /* for the 0YZ planes */
    for (m=1;m<=E->sphere.caps_per_proc;m++)  
      for(j=1;j<=E->lmesh.noy;j++)
	for(i=1;i<=E->lmesh.noz;i++)  {
	  node1 = i + (j-1)*E->lmesh.noz*E->lmesh.nox;
	  node2 = node1 + (E->lmesh.nox-1)*E->lmesh.noz;

	  ii = i + E->lmesh.nzs - 1;
	  if (E->parallel.me_locl[1]==0 )  { /* left */
	    /* no flow in x direction */
	    E->sphere.cap[m].VB[1][node1] = 0.0;
	    if((ii != 1) && (ii != E->mesh.noz))
	      /* no vertical flow */
              E->sphere.cap[m].VB[3][node1] = 0.0; 
	  }
	  if (E->parallel.me_locl[1]==E->parallel.nprocxl-1)  {
	    /* right */
	    /* no flow in x direction */
	    E->sphere.cap[m].VB[1][node2] = 0.0;
	    if((ii != 1) && (ii != E->mesh.noz))
	      /* no vertical flow */
              E->sphere.cap[m].VB[3][node2] = 0.0; 
	  }
        }      /* end loop for i and j */
  }

  /*  for two XOZ  planes  */
  if (E->parallel.me_locl[2]==0) /* top */
    for (m=1;m<=E->sphere.caps_per_proc;m++)  
      for(j=1;j<=E->lmesh.nox;j++)
        for(i=1;i<=E->lmesh.noz;i++)       {
          node1 = i + (j-1)*E->lmesh.noz;
          ii = i + E->lmesh.nzs - 1;
	  /* no flow in y */
          E->sphere.cap[m].VB[2][node1] = 0.0;
          if((ii != 1) && (ii != E->mesh.noz))
	    /* no flow in z */
            E->sphere.cap[m].VB[3][node1] = 0.0;
	}    /* end of loop i & j */
  if (E->parallel.me_locl[2]==E->parallel.nprocyl-1) /* bottom */
    for (m=1;m<=E->sphere.caps_per_proc;m++)  
      for(j=1;j<=E->lmesh.nox;j++)
        for(i=1;i<=E->lmesh.noz;i++)       {
          node2 = (E->lmesh.noy-1)*E->lmesh.noz*E->lmesh.nox + i + (j-1)*E->lmesh.noz;
          ii = i + E->lmesh.nzs - 1;
	  /* no flow in y */
          E->sphere.cap[m].VB[2][node2] = 0.0;
          if((ii != 1) && (ii != E->mesh.noz))
	    /* no flow in z */
            E->sphere.cap[m].VB[3][node2] = 0.0;
	}    /* end of loop i & j */

  
  /* 

     all vbc's apply at all levels  

  */
  for(level=E->mesh.levmax;level>=E->mesh.levmin;level--) {
    
    if ( (E->control.CONJ_GRAD && (level==E->mesh.levmax)) ||
	 (E->control.NMULTIGRID))  {
      noz = E->lmesh.NOZ[level] ;
      noy = E->lmesh.NOY[level] ;
      nox = E->lmesh.NOX[level] ;
      
      for (m=1;m<=E->sphere.caps_per_proc;m++)  { 
	if ((E->parallel.me_locl[1]==0) || 
	    (E->parallel.me_locl[1]==E->parallel.nprocxl-1)) {
	  /* YZ planes (no flow in x) */
	  for(j=1;j<=noy;j++)
	    for(i=1;i<=noz;i++) {
	      node1 = i + (j-1)*noz*nox;
	      node2 = node1 + (nox-1)*noz;
	      ii = i + E->lmesh.NZS[level] - 1;
	      if (E->parallel.me_locl[1]==0 )  { /* left */
		/* no flow in x */
		E->NODE[level][m][node1] = E->NODE[level][m][node1] | VBX;
		E->NODE[level][m][node1] = E->NODE[level][m][node1] & (~SBX);
		if((ii!=1) && (ii!=E->mesh.NOZ[level])) {
		  /* free slip in y and z */
		  E->NODE[level][m][node1] = E->NODE[level][m][node1] & (~VBY);
		  E->NODE[level][m][node1] = E->NODE[level][m][node1] | SBY;
		  E->NODE[level][m][node1] = E->NODE[level][m][node1] & (~ VBZ);
		  E->NODE[level][m][node1] = E->NODE[level][m][node1] | SBZ;
		}
	      }
	      if (E->parallel.me_locl[1]==E->parallel.nprocxl-1)  { /* right */
		/* no flow in x */
		E->NODE[level][m][node2] = E->NODE[level][m][node2] | VBX;
		E->NODE[level][m][node2] = E->NODE[level][m][node2] & (~SBX);
		if((ii!=1) && (ii!=E->mesh.NOZ[level])) {
		  /* free slip in y and z */
		  E->NODE[level][m][node2] = E->NODE[level][m][node2] & (~VBY);
		  E->NODE[level][m][node2] = E->NODE[level][m][node2] | SBY;
		  E->NODE[level][m][node2] = E->NODE[level][m][node2] & (~ VBZ);
		  E->NODE[level][m][node2] = E->NODE[level][m][node2] | SBZ;
		}
	      }
	    }   /* end for loop i & j */
	  
	}

	/* X-Z planes (no flow in y)  */
	if (E->parallel.me_locl[2]==0) /* back */
	  for(j=1;j<=nox;j++)
	    for(i=1;i<=noz;i++) {
	      node1 = i + (j-1)*noz;
	      ii = i + E->lmesh.NZS[level] - 1;
	      jj = j + E->lmesh.NXS[level] - 1;
	      /* no flow in y */
	      E->NODE[level][m][node1] = E->NODE[level][m][node1] | VBY;
	      E->NODE[level][m][node1] = E->NODE[level][m][node1] & (~SBY);
	      if((ii!= 1) && (ii != E->mesh.NOZ[level]))  {
		/* free slip in z */
                E->NODE[level][m][node1] = E->NODE[level][m][node1] & (~VBZ);
                E->NODE[level][m][node1] = E->NODE[level][m][node1] | SBZ;
	      }
	      if((jj!=1) && (jj!=E->mesh.NOX[level]) && (ii!=1) && (ii!=E->mesh.NOZ[level])){
		/* free slip in x */
                E->NODE[level][m][node1] = E->NODE[level][m][node1] & (~VBX);
                E->NODE[level][m][node1] = E->NODE[level][m][node1] | SBX;
	      }
	    }    /* end for loop i & j  */
	
	if (E->parallel.me_locl[2]==E->parallel.nprocyl-1) /* front */
	  for(j=1;j<=nox;j++)
	    for(i=1;i<=noz;i++)       {
	      node2 = (noy-1)*noz*nox + i + (j-1)*noz;
	      ii = i + E->lmesh.NZS[level] - 1;
	      jj = j + E->lmesh.NXS[level] - 1;
	      /* no slip in y */
	      E->NODE[level][m][node2] = E->NODE[level][m][node2] | VBY;
	      E->NODE[level][m][node2] = E->NODE[level][m][node2] & (~SBY);
	      if((ii!= 1) && (ii != E->mesh.NOZ[level]))  {
		/* free slip in z */
                E->NODE[level][m][node2] = E->NODE[level][m][node2] & (~VBZ);
                E->NODE[level][m][node2] = E->NODE[level][m][node2] | SBZ;
	      }
	      if((jj!=1) && (jj!=E->mesh.NOX[level]) && (ii!=1) && (ii!=E->mesh.NOZ[level])){
		/* free slip in x */
                E->NODE[level][m][node2] = E->NODE[level][m][node2] & (~VBX);
                E->NODE[level][m][node2] = E->NODE[level][m][node2] | SBX;
	      }
            }
	
      }       /* end for m  */
    }
  }       /*  end for loop level  */
  
  return;
}

/* ========================================== 
   
this assigns velocity boundary conditions on all sides of the box


*/

void velocity_side_bc(E)
     struct All_variables *E;
{
  int m,i,j;
  int node;
  int level,nox,noy,noz,nxnz;
  int elx,elz,ely;
  const int dims=E->mesh.nsd;
  float age;
  
  
  for(level=E->mesh.levmax;level>=E->mesh.levmin;level--) {
    
    noz = E->lmesh.NOZ[level] ;
    noy = E->lmesh.NOY[level] ;
    nox = E->lmesh.NOX[level] ;

    elx=E->lmesh.ELX[level];
    elz=E->lmesh.ELZ[level];
    ely=E->lmesh.ELY[level];


    
    nxnz = nox * noz;
    /* 
       
    index(i_x, j_y, k_z) = k_z + (i_x-1)*n_z + (j_y-1)*n_x*n_z 
    
    i=1,...,n
    */
    for (m=1;m <= E->sphere.caps_per_proc;m++)  { 
      if (E->parallel.me_locl[1]==0){
	/* 
	   YZ (phi, r) North, x=1
	*/
	for(j=1;j<=noy;j++)
	  for(i=1;i<=noz;i++) {
	    node = i + (j-1)*nxnz;
	    vsideflow_assign(node,level,m,0,age,TRUE,TRUE,E,FALSE);
	  }
      }
      if(E->parallel.me_locl[1]==E->parallel.nprocxl-1){
	/* 
	   YZ (phi, r) South, x=nx
	*/
	for(j=1;j<=noy;j++)
	  for(i=1;i<=noz;i++) {
	    node = i + (j-1)*nxnz + (nox-1)*noz;
	    vsideflow_assign(node,level,m,1,age,TRUE,TRUE,E,FALSE);
	  }
      }
      if(E->parallel.me_locl[2]==0){
	/* 
	   XZ (theta, r) West, y=1
	*/
	for(j=1;j<=nox;j++)
	  for(i=1;i<=noz;i++) {
	    node = i + (j-1)*noz;
	    vsideflow_assign(node,level,m,2,age,TRUE,TRUE,E,FALSE);
	  }
      }
      if(E->parallel.me_locl[2]==E->parallel.nprocyl-1){
	/* 
	   
	XZ (theta, r) East,y=ny
	
	*/
	for(j=1;j<=nox;j++)
	  for(i=1;i<=noz;i++) {
	    node = i + (j-1)*noz + (noy-1)*nxnz;
	    vsideflow_assign(node,level,m,3,age,TRUE,TRUE,E,FALSE);
	  }
      }
      if (E->parallel.me_locl[3]==0){
	/* 
	   XY (theta, phi) bottom, z=1
	*/
	for(j=1;j<=noy;j++)
	  for(i=1;i<=nox;i++) {
	    node = (i-1)*noz + (j-1)*nxnz +1;
	    vsideflow_assign(node,level,m,4,age,TRUE,TRUE,E,FALSE);
	  }
      }
      if(E->parallel.me_locl[3]==E->parallel.nproczl-1){
	/* 
	   XY top, z=nz
	*/
	for(j=1;j<=noy;j++)
	  for(i=1;i<=nox;i++) {
	    node = (i-1) * noz + (j-1)*nxnz + noz;
	    vsideflow_assign(node,level,m,5,age,TRUE,TRUE,E,FALSE);
	  }
	}
    }
    /* print stats */
    vsideflow_assign(node,level,m,5,age,FALSE,FALSE,E,TRUE);
  }
  return;
}
/* 
   assign side flow boundary conditions 
   code: code of plane (see above)
*/

void vsideflow_assign( node, level, m, code, age, assign_bc, assign_val, E,
		       print_stats)
     int node,level,m,code;
     float age;
     boolean assign_bc,assign_val,print_stats;
     struct All_variables *E;
{
#ifdef USE_GGRD
  HC_CPREC xloc[3],time,vr[1],vt[1],vp[1];
  static int order = 3;		/* order of velocity interpolation */
  static boolean calc_derivatives = FALSE, /* don't need derivatives of velocities */
     init = FALSE;
  static HC_CPREC vscale;
  /* for stats */
  static float vmean[6][4],vmin[6][4],vmax[6][4];
  static int vn[6];
  static boolean do_stats = TRUE;
  float vadd,vtot,area,r1r2,dp,dt;
#endif
  /*  */
  int i,j;
#ifndef USE_GGRD
  myerror("vsideflow_assign relies on ggrd",E);
#else
  if(!init){
    /* all processors load velocities into memory */
    /* 
       initialize velocity structure
    */
    E->control.ggrd_v = (struct ggrd_vel *)calloc(1,sizeof(struct ggrd_vel));
    ggrd_init_vstruc(E->control.ggrd_v);
    /* 
       read in scaled velocities 
    */
    vscale = (HC_CPREC) E->data.timedir/E->data.scalev;
    if(E->parallel.me==0)
      fprintf(stderr,"initializing sideflow BCs from files in %s, scale %g\n",
	      E->control.vsideflow_dir,vscale);
    ggrd_read_vel_grids(E->control.ggrd_v,1.0/vscale,
			FALSE,TRUE,E->control.vsideflow_dir);
    /* init means and extrena */
    for(i=0;i < 6;i++){
      for(j=1;j <= 3;j++){
	vmean[i][j] =  0.0;
	vmin[i][j] = 1e20;
	vmax[i][j] = -1e20;
      }
      vn[i] = 0;
    }
    init = TRUE;
  }
  if(assign_bc){
    /* 
       prescribed flow in all components 
    */
    E->NODE[level][m][node] = E->NODE[level][m][node] | VBX;
    E->NODE[level][m][node] = E->NODE[level][m][node] & (~SBX);
    E->NODE[level][m][node] = E->NODE[level][m][node] | VBY;
    E->NODE[level][m][node] = E->NODE[level][m][node] & (~SBY);
    E->NODE[level][m][node] = E->NODE[level][m][node] | VBZ;
    E->NODE[level][m][node] = E->NODE[level][m][node] & (~SBZ);
  }
  if(assign_val){
    /* 
       find velocities 
    */
    /* coordinated of node r, theta, phi, code */
    xloc[HC_R] = (HC_CPREC)E->sx[m][3][node]; /* r */
    xloc[HC_THETA] = (HC_CPREC)E->sx[m][1][node]; /* theta */
    xloc[HC_PHI] = (HC_CPREC)E->sx[m][2][node]; /* phi */
    time = (HC_CPREC)age;
    /* interpolate */
    ggrd_find_vel_and_der(xloc,time,E->control.ggrd_v,order,calc_derivatives,
			  FALSE,vr,vt,vp);
    /*
      fprintf(stderr,"%g %g %g\t%g %g %g\n",PHI2LON(xloc[HC_PHI]),
      90-THETA2LAT(xloc[HC_THETA]),-(1-xloc[HC_R])*6371,
      vp[0]/vscale,-vt[0]/vscale,vr[0]/vscale);
    */
    /* assign */
    E->sphere.cap[m].VB[1][node] = (float)vt[0];
    E->sphere.cap[m].VB[2][node] = (float)vp[0];
    E->sphere.cap[m].VB[3][node] = (float)vr[0];
    if(do_stats){
      /* 
	 do some stats 
      */
      if((code >= 6)||(code<0)){
	fprintf(stderr,"vsideflow_assign: p: %2i: code error: %i\n",
		E->parallel.me,code);
	exit(-1);
      }
      /* mean velocities on sides */
      vn[code]++;
      for(i=1;i <= 3;i++){
 	vmean[code][i] +=  E->sphere.cap[m].VB[i][node]; 
	if(vmin[code][i] > E->sphere.cap[m].VB[i][node]) 
	  vmin[code][i] =  E->sphere.cap[m].VB[i][node]; 
	if(vmax[code][i] < E->sphere.cap[m].VB[i][node]) 
	  vmax[code][i] =  E->sphere.cap[m].VB[i][node]; 
      }
    }
  }
  if(print_stats){
    vtot = 0.0;
    for(code=0;code < 6;code++){
      if(vn[code]){
	fprintf(stderr,"vsideflow_assign: p: %2i: s%i: ",E->parallel.me,code);
	for(j=1;j <= 3;j++)
	  /* output of min/mean/max */
	  fprintf(stderr,"v%i %9.2e/%9.2e/%9.2e ",
		  j,vmin[code][j],vmean[code][j]/(float)vn[code],
		  vmax[code][j]);
	/* sums of velocity.normal on plane */
	r1r2 = (E->sphere.ro * E->sphere.ro - E->sphere.ri * E->sphere.ri)/
	  ((float)E->parallel.nproczl);
	dp = (E->control.fi_max - E->control.fi_min)/((float)E->parallel.nprocyl);
	dt = (E->control.theta_max-E->control.theta_min)/((float)E->parallel.nprocxl);
	switch(code){
	case 0: 		/* North */
	  vadd = -vmean[code][1];
	  area = sin(E->control.theta_min)*r1r2*dp/2.0;
	  break;
	case 1: 		/* South */
	  vadd =  vmean[code][1];
	  area = sin(E->control.theta_max)*r1r2*dp/2.0;
	  break;
	case 2: 		/* West */
	  vadd = -vmean[code][2];
	  area = dt/2.0*r1r2;
	  break;
	case 3: 		/* East */
	  vadd =  vmean[code][2];
	  area = dt/2.0*r1r2;
	  break;
	case 4: 		/* bottom */
	  vadd =  vmean[code][3];
	  area = E->sphere.ri * E->sphere.ri * dp / ((float)E->parallel.nprocxl) * 
	    (cos(E->control.theta_min)-cos(E->control.theta_max));
	  break;
	case 5: 		/* top */
	  vadd = -vmean[code][3];
	  area = E->sphere.ro * E->sphere.ro * dp / ((float)E->parallel.nprocxl) * 
	    (cos(E->control.theta_min)-cos(E->control.theta_max));
	  break;
	}
	fprintf(stderr,", Sum(n.v): %11.4e A: %11.4e SA: %11.4e\n",vadd, area,area*vadd);
	vtot += vadd * area;
      }else{
	fprintf(stderr,"vsideflow_assign: p: %2i: s%i: nothing assigned\n",
		E->parallel.me,code);
      }
    }
    fprintf(stderr,"vsideflow_assign: total vnsum: %12.5e\n",vtot);
  }
#endif
}
/* 
   this does: ?

*/
void temperature_refl_vert_bc(E)
     struct All_variables *E;
{
  int i,j,m;
  int node1,node2;
  const int dims=E->mesh.nsd;

 /* Temps and bc-values  at top level only */

  if (E->parallel.me_locl[1]==0 || E->parallel.me_locl[1]==E->parallel.nprocxl-1)
    for(m=1;m<=E->sphere.caps_per_proc;m++)
      for(j=1;j<=E->lmesh.noy;j++)
	for(i=1;i<=E->lmesh.noz;i++) {
	  node1 = i + (j-1)*E->lmesh.noz*E->lmesh.nox;
	  node2 = node1 + (E->lmesh.nox-1)*E->lmesh.noz;
	  if (E->parallel.me_locl[1]==0 )                   {
	    E->node[m][node1] = E->node[m][node1] & (~TBX);
	    E->node[m][node1] = E->node[m][node1] | FBX;
	    E->sphere.cap[m].TB[1][node1] = 0.0;
	  }
	  if (E->parallel.me_locl[1]==E->parallel.nprocxl-1)   {
	    E->node[m][node2] = E->node[m][node2] & (~TBX);
	    E->node[m][node2] = E->node[m][node2] | FBX;
	    E->sphere.cap[m].TB[1][node2] = 0.0;
	  }
        }       /* end for loop i & j */
  
    if (E->parallel.me_locl[2]==0)
      for(m=1;m<=E->sphere.caps_per_proc;m++)
	for(j=1;j<=E->lmesh.nox;j++)
	  for(i=1;i<=E->lmesh.noz;i++) {
	    node1 = i + (j-1)*E->lmesh.noz;
	    E->node[m][node1] = E->node[m][node1] & (~TBY);
	    E->node[m][node1] = E->node[m][node1] | FBY;
	    E->sphere.cap[m].TB[2][node1] = 0.0;
	  }

    if (E->parallel.me_locl[2]==E->parallel.nprocyl-1)
      for(m=1;m<=E->sphere.caps_per_proc;m++)
	for(j=1;j<=E->lmesh.nox;j++)
	  for(i=1;i<=E->lmesh.noz;i++) {
	    node2 = i +(j-1)*E->lmesh.noz + (E->lmesh.noy-1)*E->lmesh.noz*E->lmesh.nox;
	    E->node[m][node2] = E->node[m][node2] & (~TBY);
	    E->node[m][node2] = E->node[m][node2] | FBY;
	    E->sphere.cap[m].TB[3][node2] = 0.0;
          }    /* end loop for i and j */

    return;
}


/*  =========================================================  */
    
/* 

this actually assigns the horizontal boundary condition

*/
void horizontal_bc(E,BC,ROW,dirn,value,mask,onoff,level,m)
     struct All_variables *E;
     float *BC[];
     int ROW;
     int dirn;
     float value;
     unsigned int mask;
     char onoff;
     int level,m;

{
  int i,j,node,rowl;
  const int dims=E->mesh.nsd;

  /* safety feature */
  if(dirn > E->mesh.nsd) 
    return;

  if (ROW==1) 
    rowl = 1;
  else 
    rowl = E->lmesh.NOZ[level];
   
  if ( ((ROW == 1) && (E->parallel.me_locl[3] == 0)) ||
       ((ROW == E->lmesh.NOZ[level]) && (E->parallel.me_locl[3] == E->parallel.nproczl-1 ))) {
    
    /* turn bc marker to zero */
    if (onoff == 0){
      /* switches BC off */
      for(j=1;j<=E->lmesh.NOY[level];j++)
    	for(i=1;i<=E->lmesh.NOX[level];i++)     {
    	  node = rowl+(i-1)*E->lmesh.NOZ[level]+(j-1)*E->lmesh.NOX[level]*E->lmesh.NOZ[level];
    	  E->NODE[level][m][node] = E->NODE[level][m][node] & (~ mask);
	}        /* end for loop i & j */
    }else{
      /* switch bc on */
      for(j=1;j<=E->lmesh.NOY[level];j++)
        for(i=1;i<=E->lmesh.NOX[level];i++)       {
    	  node = rowl+(i-1)*E->lmesh.NOZ[level]+(j-1)*E->lmesh.NOX[level]*E->lmesh.NOZ[level];
    	  E->NODE[level][m][node] = E->NODE[level][m][node] | (mask);
	  
    	  if(level==E->mesh.levmax)   /* NB */
    	    BC[dirn][node] = value;   
    	  }     /* end for loop i & j */
      }

    }             /* end for if ROW */
    
  return;
}


void velocity_apply_periodic_bcs(E)
    struct All_variables *E;
{
  int n1,n2,level;
  int i,j,ii,jj;
  const int dims=E->mesh.nsd;

  fprintf(E->fp,"Periodic boundary conditions\n");

  return;
  }

void temperature_apply_periodic_bcs(E)
    struct All_variables *E;
{
 const int dims=E->mesh.nsd;

 fprintf(E->fp,"pERIodic temperature boundary conditions\n");
   
  return;
  }



void strip_bcs_from_residual(E,Res,level)
    struct All_variables *E;
    double **Res;
    int level;
{
    int m,i;
     
  for (m=1;m<=E->sphere.caps_per_proc;m++)
    if (E->num_zero_resid[level][m])
      for(i=1;i<=E->num_zero_resid[level][m];i++)
         Res[m][E->zero_resid[level][m][i]] = 0.0;

    return;
    }


void get_bcs_id_for_residual(E,level,m)
    struct All_variables *E;
    int level,m;
  {

    int i,j;
     
    const int dims=E->mesh.nsd,dofs=E->mesh.dof;
    const int nno=E->lmesh.NNO[level];
    const int neq=E->lmesh.NEQ[level];
    const int addi_dof=additional_dof[dims];

   j = 0;
   for(i=1;i<=nno;i++) { 
      if ( (E->NODE[level][m][i] & VBX) != 0 )  {
	j++;
        E->zero_resid[level][m][j] = E->ID[level][m][i].doff[1];
	}
      if ( (E->NODE[level][m][i] & VBY) != 0 )  {
	j++;
        E->zero_resid[level][m][j] = E->ID[level][m][i].doff[2];
	}
      if ( (E->NODE[level][m][i] & VBZ) != 0 )  {
	j++;
        E->zero_resid[level][m][j] = E->ID[level][m][i].doff[3];
	}
      }

    E->num_zero_resid[level][m] = j;
 
    return;
}
    
/* 

set the fixed temperature boundary conditions (at each timestep)
with possible adjustments at the sides

*/
void temperatures_conform_bcs(E)
     struct All_variables *E;
{
  int m,j,nno,node,nox,noz,noy,gnox,gnoy,gnoz,nodeg,ii,i,k;
  unsigned int type;
  float ttt2,ttt3,fff2,fff3;
  float r1,t1,f1,t0,tmp;
  float e_4;
  FILE *fp1, *fp2;
  float rlith_lim1,rlith_lim2;
  int output;
  static boolean been_here = FALSE;
  
  /* for side boundary adjustment */
  rlith_lim1 = E->sphere.ro - E->control.depth_bound_adj;
  /* for lith age adjustment */
  rlith_lim2 = E->sphere.ro - E->control.lith_age_depth;

  e_4=1.e-4;
  nno=E->lmesh.nno;
  output = 0;
  /* 
     boundaries for side frames around box to fix 
     temperatures 
  */
  ttt2=E->control.theta_min + E->control.width_bound_adj; /* north */
  ttt3=E->control.theta_max - E->control.width_bound_adj; /* south */
  fff2=E->control.fi_min + E->control.width_bound_adj; /* west */
  fff3=E->control.fi_max - E->control.width_bound_adj; /* east */
  
  gnox=E->mesh.nox;
  gnoy=E->mesh.noy;
  gnoz=E->mesh.noz;
  nox=E->lmesh.nox;
  noy=E->lmesh.noy;
  noz=E->lmesh.noz;
  /* 
     read in lithospheric ages 
  */
  if((!been_here) || (E->control.lith_age_time)){
    get_lith_ages(E);
    been_here = TRUE;
  }
  /* 

     NOW SET THE TEMPERATURES IN THE BOUNDARY REGIONS 
     
     starts at > 0 now TWB

  */
  if((E->monitor.solution_cycles > 0) && 
     (E->control.temperature_bound_adj))   {
    if(!E->control.age_init)
      myerror("temperatures_conform_bcs: ages not initialized",E);
    for(m=1;m<=E->sphere.caps_per_proc;m++)
      for(i=1;i<=noy;i++)  
	for(j=1;j<=nox;j++) 
	  for(k=1;k<=noz;k++)  {
	    nodeg=E->lmesh.nxs-1+j+(E->lmesh.nys+i-2)*gnox;
	    node=k+(j-1)*noz+(i-1)*nox*noz;
	    /* coordinates */
	    t1=E->sx[m][1][node];
	    f1=E->sx[m][2][node];
	    r1=E->sx[m][3][node];

	    /* if NOT right on the bottom or top boundary */
	    if((fabs(r1-E->sphere.ro)>=e_4) && 
	       (fabs(r1-E->sphere.ri)>=e_4))  { 

	      if( ((t1 <= ttt2) && (r1 >= rlith_lim1)) ||
		  ((t1 >= ttt3) && (r1 >= rlith_lim1)) ) {

		/* if < (width) from x bounds 
		   AND (depth) from top */
		t0 = half_space_temp(E,r1,E->age_t[nodeg],
				     E->sphere.cap[m].TB[1][node],FALSE);
		E->sphere.cap[m].TB[1][node]=t0;
		E->sphere.cap[m].TB[2][node]=t0;
		E->sphere.cap[m].TB[3][node]=t0;
	      }

	      if( ((f1 <= fff2) || (f1 >= fff3))){
		/* if < (width) from y bounds */
		t0 = half_space_temp(E, r1, E->age_t[nodeg],
				     E->sphere.cap[m].TB[1][node],FALSE);
		E->sphere.cap[m].TB[1][node]=t0;
		E->sphere.cap[m].TB[2][node]=t0;
		E->sphere.cap[m].TB[3][node]=t0;
	      }
	    }
	  }     /* end k   */
      
  }  /*  end of solution cycles  && temperature_bound_adj */

  /* 
     NOW SET THE TEMPERATURES IN THE LITHOSPHERE IF CHANGING EVERY TIME STEP 
     changed this to start from > 0 TWB
  */
  if((E->monitor.solution_cycles > 0) && (E->control.lith_age_time))   {
    if(E->parallel.me==0)
      fprintf(stderr,"temperature_conform_bcs: re-adjusting lith_age T, step: %i, age: %g\n",
	      E->monitor.solution_cycles,find_age_in_MY(E));
    for(m=1;m <= E->sphere.caps_per_proc;m++)
      for(i=1;i <= noy;i++)  
	for(j=1;j <= nox;j++) 
	  for(k=1;k <= noz;k++)  {
	    nodeg=E->lmesh.nxs-1+j+(E->lmesh.nys+i-2)*gnox;
	    node=k+(j-1)*noz+(i-1)*nox*noz;
	    /* radial coordinate */
	    r1 = E->sx[m][3][node];
	    if((fabs(r1 - E->sphere.ro)>=e_4) && 
	       (fabs(r1 - E->sphere.ri)>=e_4))  {
	      /* if NOT right on the radial boundaries */
	      if((E->age_t[nodeg] >= 0) && (r1 >= rlith_lim2)){
		t0 = half_space_temp(E,r1,E->age_t[nodeg],
				     E->sphere.cap[m].TB[1][node],FALSE);

		tmp = E->control.mantle_temp * 
			erf(0.5 / sqrt(E->age_t[nodeg]) * 
			    (E->sphere.ro-r1)* E->data.radius_km /E->data.layer_km);
		if(fabs(tmp - t0) > 1e-6){
		  fprintf(stderr,"error: t_age %g and t_new %g different for age %g\n",
			  tmp, t0, E->age_t[nodeg]);
		  exit(-1);
		}

		E->sphere.cap[m].TB[1][node] = t0;
		E->sphere.cap[m].TB[2][node] = t0;
		E->sphere.cap[m].TB[3][node] = t0;
	      }
	    }
	  }     /* end k   */
  }   /*  end of solution cycles  && lith_age_time */
    
  /* end control.lith_age=true */
    

  /* 

  
  check for the top/bottom boundary condition temperatures


  */
  for(j=1;j<=E->sphere.caps_per_proc;j++)
    for(node=1;node<=E->lmesh.nno;node++)  {

      type = (E->node[j][node] & (TBX | TBZ | TBY));

      switch (type) {
      case 0:  /* no match, next node */
	break;
      case TBX:
	E->T[j][node] = E->sphere.cap[j].TB[1][node];
	break;
      case TBZ:
	E->T[j][node] = E->sphere.cap[j].TB[3][node];
	break;
      case TBY:
	E->T[j][node] = E->sphere.cap[j].TB[2][node];
	break;
      case (TBX | TBZ):     /* clashes ! */
	E->T[j][node] = 0.5 * (E->sphere.cap[j].TB[1][node] + E->sphere.cap[j].TB[3][node]);
	break;
      case (TBX | TBY):     /* clashes ! */
	E->T[j][node] = 0.5 * (E->sphere.cap[j].TB[1][node] + E->sphere.cap[j].TB[2][node]);
	break;
      case (TBZ | TBY):     /* clashes ! */
	E->T[j][node] = 0.5 * (E->sphere.cap[j].TB[3][node] + E->sphere.cap[j].TB[2][node]);
	break;
      case (TBZ | TBY | TBX):     /* clashes ! */
	E->T[j][node] = 0.3333333 * (E->sphere.cap[j].TB[1][node] + E->sphere.cap[j].TB[2][node] + E->sphere.cap[j].TB[3][node]);
	break;
      } /* end switch */

        
    } /* next node */

  return;

}


void velocities_conform_bcs(E,U)
    struct All_variables *E;
    double **U;
{ 
    int node,d,m;

    const unsigned int typex = VBX;
    const unsigned int typez = VBZ;
    const unsigned int typey = VBY;
    const int addi_dof = additional_dof[E->mesh.nsd];

    const int dofs = E->mesh.dof;
    const int nno = E->lmesh.nno;

    for(m=1;m<=E->sphere.caps_per_proc;m++)   {
      for(node=1;node<=nno;node++) {

        if (E->node[m][node] & typex)  
	      U[m][E->id[m][node].doff[1]] = E->sphere.cap[m].VB[1][node]; 
 	if (E->node[m][node] & typey)  
	      U[m][E->id[m][node].doff[2]] = E->sphere.cap[m].VB[2][node]; 
	if (E->node[m][node] & typez)  
	      U[m][E->id[m][node].doff[3]] = E->sphere.cap[m].VB[3][node]; 
        } 
      } 

    return;
}

void temperature_lith_adj(E,lv)
    struct All_variables *E;
    int lv;
 {
   int i,j,k,node,nno,nodeg,m,itmp;
   int gnox,nox,noy,noz;
   float ttt2,ttt3,fff2,fff3;
   FILE *fp1;
   char output[CHARBUF_SIZE];
   float rlith_lim1,rlith_lim2;
  
   rlith_lim1 = E->sphere.ro - E->control.depth_bound_adj;
   rlith_lim2 = E->sphere.ro - E->control.lith_age_depth;

   nno=E->lmesh.nno;
     
   gnox=E->mesh.nox;

   nox=E->lmesh.nox;
   noy=E->lmesh.noy;
   noz=E->lmesh.noz;


   ttt2=E->control.theta_min + E->control.width_bound_adj;
   ttt3=E->control.theta_max - E->control.width_bound_adj;
   fff2=E->control.fi_min + E->control.width_bound_adj;
   fff3=E->control.fi_max - E->control.width_bound_adj;
   fprintf(stderr,"temperature_lith_adj: (0) lv: %i gm: %i\n",lv,E->mesh.gridmax);

   /* NOTE: To start, the relevent bits of "node" are zero. Thus, they only
      get set to TBX/TBY/TBZ if the node is in one of the bounding regions.
      Also note that right now, no matter which bounding region you are in,
      all three get set to true. CPC 6/20/00 */
    
   if (E->control.temperature_bound_adj == 1) {
     /* 
	this was if(lv=E->mesh.gridmax) before. it doesn't make sense to me that 
	lv should be set here, since it's not passed on anywhere TWB
     */
     if(lv == E->mesh.gridmax) { 
       if(E->parallel.me == 0)
	 fprintf(stderr,"temperature_lith_adj: (1) setting BC type for temperature_bound_adj, step %i\n",
		 E->monitor.solution_cycles);
       for(m=1;m <= E->sphere.caps_per_proc;m++){
	 for(i=1;i<=noy;i++) {
	   for(j=1;j<=nox;j++){ 
	     for(k=1;k<=noz;k++)  {
	       nodeg=E->lmesh.nxs-1+j+(E->lmesh.nys+i-2)*gnox;
	       node=k+(j-1)*noz+(i-1)*nox*noz;
	       if( ((E->sx[m][1][node]<=ttt2) && 
		    (E->sx[m][3][node]>=rlith_lim1)) || 
		   ((E->sx[m][1][node]>=ttt3) && 
		    (E->sx[m][3][node]>=rlith_lim1)) ) 
		 /* if < (width) from x bounds AND (depth) from top */
		 {
		   E->node[m][node]=E->node[m][node] | TBX;
		   E->node[m][node]=E->node[m][node] & (~FBX);
		   E->node[m][node]=E->node[m][node] | TBY;
		   E->node[m][node]=E->node[m][node] & (~FBY);
		   E->node[m][node]=E->node[m][node] | TBZ;
		   E->node[m][node]=E->node[m][node] & (~FBZ);
		 }
	       
	       if( ((E->sx[m][2][node]<=fff2) && 
		    (E->sx[m][3][node]>=rlith_lim1)) ) 
		 /* if fi is < (width) from side AND z is < (depth) from top */
		 {
		   E->node[m][node]=E->node[m][node] | TBX;
		   E->node[m][node]=E->node[m][node] & (~FBX);
		   E->node[m][node]=E->node[m][node] | TBY;
		   E->node[m][node]=E->node[m][node] & (~FBY);
		   E->node[m][node]=E->node[m][node] | TBZ;
		   E->node[m][node]=E->node[m][node] & (~FBZ);
		 }
	       
	       if( ((E->sx[m][2][node]>=fff3) && 
		    (E->sx[m][3][node]>=rlith_lim1)) ) 
		 /* if fi is < (width) from side AND z is < (depth) from top */
		 {
		   E->node[m][node]=E->node[m][node] | TBX;
		   E->node[m][node]=E->node[m][node] & (~FBX);
		   E->node[m][node]=E->node[m][node] | TBY;
		   E->node[m][node]=E->node[m][node] & (~FBY);
		   E->node[m][node]=E->node[m][node] | TBZ;
		   E->node[m][node]=E->node[m][node] & (~FBZ);
		 }
	     }
	   }
	 }
       }
     } /* end lv = gridmax */
   } /* end E->control.temperature_bound_adj==1 */

   if(E->control.lith_age_time == 1) {
     /* 
	this was if(lv=E->mesh.gridmax) before. it doesn't make sense to me that 
	lv should be set here, since it's not passed on anywhere TWB
     */
     if(lv == E->mesh.gridmax){
       if(E->parallel.me == 0)
	 fprintf(stderr,"temperature_lith_adj: (2) setting BC type for lith_age_time, step %i\n",
		 E->monitor.solution_cycles);
       if(!E->control.age_init)
	 myerror("temperatures_lith_adj: ages not initialized",E);
       itmp = nox*noz;
       for(m=1;m <= E->sphere.caps_per_proc;m++){
	 for(i=1;i <= noy;i++) {
	   for(j=1;j <= nox;j++){ 
	     for(k=1;k <= noz;k++)  {
	       nodeg=E->lmesh.nxs-1+j+(E->lmesh.nys+i-2)*gnox;
	       node=k+(j-1)*noz+(i-1)*itmp;
	       if((E->sx[m][3][node] >= rlith_lim2) &&
		  /* if closer than (lith_age_depth) from top */
		  (E->age_t[nodeg] >= 0.0)){
		 /* and age >= 0  */
		 E->node[m][node]=E->node[m][node] | TBX;
		 E->node[m][node]=E->node[m][node] & (~FBX);
		 E->node[m][node]=E->node[m][node] | TBY;
		 E->node[m][node]=E->node[m][node] & (~FBY);
		 E->node[m][node]=E->node[m][node] | TBZ;
		 E->node[m][node]=E->node[m][node] & (~FBZ);
	       }
	     }
	   }
	 }
       }
     }
   } /* end E->control.lith_age_time==1 */

   return;
}

/* 
   
this routine was previously used to specify surface velocity 
boundary conditions

typically not called anymore

*/
void renew_top_velocity_boundary(E)
     struct All_variables *E;
{
  int i,k,lev;
  int nox,noz,noy,nodel;
  float fxx10,fxx20,fyy1,fyy2,fxx0,fxx,fyy;
  float vxx1,vxx2,vxx,vvo,vvc;
  float fslope,vslope;
  static float fxx1,fxx2;
  
  FILE *fp;
  char output_file[CHARBUF_SIZE];


  myerror("renew_top_velocity_boundary: TEST FIRST. exiting",E);

  nox=E->lmesh.nox;
  noz=E->lmesh.noz;
  noy=E->lmesh.noy;
  lev=E->mesh.levmax;
  
  fxx10=1.0;
  fyy1=0.76;
  fxx20=1.0;   /* (fxx1,fyy1), (fxx2,fyy2) the initial coordinates of the trench position */
  fyy2=0.81;
  
  vxx1=-2.*2.018e0;
  
  vvo=6.*2.018e0;   
  vvc=-2.*2.018e0;     /* vvo--oceanic plate velocity; vvc--continental plate velocity      */
  
  if(E->advection.timesteps>1)  {
    fxx1=fxx1+E->advection.timestep*vxx1;
    fxx2=fxx2+E->advection.timestep*vxx1;
  }else  {
    fxx1=fxx10;
    fxx2=fxx20;
  }
  
  fprintf(stderr,"%f %f\n",fxx1,fxx2);
   
   if (E->parallel.me_locl[3] == E->parallel.nproczl-1 ) {
     for(k=1;k<=noy;k++)
       for(i=1;i<=nox;i++)   {
	 nodel = (k-1)*nox*noz + (i-1)*noz+noz;
	 fyy=E->SX[lev][1][1][nodel];
	 if (fyy < fyy1 || fyy >fyy2 )   {
	   E->sphere.cap[1].VB[1][nodel]=0.0;
	   E->sphere.cap[1].VB[2][nodel]=-vvc;
	   E->sphere.cap[1].VB[3][nodel]=0.0;
	 }    /* the region outside of the domain bounded by the trench length  */
	 else if (fyy>=fyy1 && fyy <=fyy2)  {
	   if ((E->SX[lev][1][2][nodel]>=0.0) && (E->SX[lev][1][2][nodel]<= fxx1)) {
	     E->sphere.cap[1].VB[1][nodel]=0.0;
	     E->sphere.cap[1].VB[2][nodel]=vvo;
	     E->sphere.cap[1].VB[3][nodel]=0.0;
	   }
	   else if ( E->SX[lev][1][2][nodel]>fxx1 && E->SX[lev][1][2][nodel]<fxx2) {
	     E->sphere.cap[1].VB[1][nodel]=0.0;
	     E->sphere.cap[1].VB[2][nodel]=vxx1;
	     E->sphere.cap[1].VB[3][nodel]=0.0;
	   }
	   else if ( E->SX[lev][1][2][nodel]>=fxx2) {
	     E->sphere.cap[1].VB[1][nodel]=0.0;
	     E->sphere.cap[1].VB[2][nodel]=vvc;
	     E->sphere.cap[1].VB[3][nodel]=0.0;
	   }
	 }   /* end of else if (fyy>=fyy1 && fyy <=fyy2)  */
	 
       }  /* end if for(i=1;i<nox;i++)  */
   }    /* end of E->parallel.me_locl[3]   */
   
   return;
 }



