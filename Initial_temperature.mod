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


minor changes by TWB

$Id: Initial_temperature.c,v 1.19 2004/12/03 02:38:25 becker Exp becker $



=====================================================================
*/

#include <math.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "proto.h"
/*
minor changes in this file by TWB

$Id: Initial_temperature.c,v 1.19 2004/12/03 02:38:25 becker Exp becker $

see the rlog of RCS
*/
#include <stdlib.h> /* for "system" command */
#include <strings.h>
#ifdef USE_GGRD
#include "hc.h"			/* Hager & O'Connell routines */
#endif


void convection_initial_temperature(E)
     struct All_variables *E;
{

    int i,j,k,p,node,ii,jj,m,mm,ll;
    int noz2;
    char input_file[1000],input_s[1000];
    float lscale,ficenter,rcenter,dist2,lscalex,rscale,age;

    double con,t1,f1,r1,tgrad,tbot,tadd;
    float rad,beta,notusedhere; 
    float a,b,c,d,e,f,g; 
    FILE *fp,*fp1,*fp2;

    float v1,v2,v3,tt;
    float e_4;
    /* for other boundary conditions, i moved this here
       to avoid C++ usage of variable declaration, TWB
    */
    int number_of_perturbations;
    int perturb_ll[32], perturb_mm[32], load_depth[32];
    float perturb_mag[32];
    float tlen, flen, rlen;
    
    boolean rezip=FALSE;


    const int dims=E->mesh.nsd;
    rad = 180.0/M_PI;
    e_4=1.e-4;


    noz2=(E->mesh.noz-1)/2+1;  

    /* 
       do we need to read temperature/density from file?

       TWB

    */
    input_boolean("ggrd_tinit",&(E->control.ggrd_tinit),"off",E->parallel.me); /* switch */
    input_double("ggrd_tinit_scale",&(E->control.ggrd_tinit_scale),"1.0",E->parallel.me); /* scale */
    input_string("ggrd_tinit_gfile",E->control.ggrd_tinit_gfile,"",E->parallel.me); /* grids */
    input_string("ggrd_tinit_dfile",E->control.ggrd_tinit_dfile,"",E->parallel.me); /* depth.dat latyers */

    if ((E->control.restart || E->control.post_p))    {
      if(E->control.ggrd_tinit){
	fprintf(stderr,"p %i: ggrd_tinit and (restart || post_p) set at same time?!\n",
		E->parallel.me);
	exit(-1);
      }
      /* 

      USED IF RESTARTING FROM A PREVIOUS RUN. CPC 1/28/00 
      

      */
      /* reset the solution cycles TWB */
      E->monitor.solution_cycles = E->monitor.solution_cycles_init;
      /* 
	 old velocity file 
      */
      sprintf(input_file,"%s.velo.%d.%d",
	      E->control.old_P_file,E->parallel.me,
	      E->monitor.solution_cycles_init);
      /* 
	 open file or zipped version of it
      */
      rezip = open_file_zipped(input_file,&fp,E);
      if(fscanf(fp,"%i %i %f",&ll,&mm,&notusedhere) != 3)
	myerror("read error restart velocities:1 ",E);

      
      for(m=1;m<=E->sphere.caps_per_proc;m++)  {
	if(fscanf(fp,"%i %i",&ll,&mm) != 2)
	  myerror("read error restart velocities:2",E);
	if(mm != E->lmesh.nno)
	  myerror("mismatch of lmesh.nno in input file",E);
	for(i=1;i<=E->lmesh.nno;i++)  {
	  /* 
	     read in velocities and temperature from restart file
	  */
	  if(fscanf(fp,"%f %f %f %f",&v1,&v2,&v3,&tt)!=4)
	    myerror("read error restart velocities: 3",E);
	  E->T[m][i] = tt;	/* temperature */
	  if(E->control.restart_with_vp){
	    /* 

	    use velocities only when pressure are also read in 
	    TWB

	    */
	    E->sphere.cap[m].V[1][i] = v1; /* v_theta */
	    E->sphere.cap[m].V[2][i] = v2; /* v_phi */
	    E->sphere.cap[m].V[3][i] = v3; /* v_r */
	  }
	  /* maybe make 0 <= T <= 1 here again? TWB */
	}
      }	/* end m loop */
      fclose (fp);
      if(rezip)
	compress_file(input_file);
      my_report(E,"read temperatures for restart");
      if(E->control.restart_with_vp){
	/* 
	   
	read pressures at nodes

	*/
	my_report(E,"reading pressures and velocities for restart");
	record(E,"reading pressures and velocities for restart");
	sprintf(input_file,"%s.p.%d.%d",
	      E->control.old_P_file,E->parallel.me,
	      E->monitor.solution_cycles_init);
	/* 
	   open file or zipped version of it
	*/
	rezip = open_file_zipped(input_file,&fp,E);
	if(fscanf(fp,"%i %i %f",&ll,&mm,&notusedhere)!=3)
	  myerror("read error restart pressures:1b ",E);
	for(m=1;m<=E->sphere.caps_per_proc;m++)  {
	  if(fscanf(fp,"%i %i",&ll,&mm) != 2)
	    myerror("read error restart pressures:2b",E);
	  if(mm != E->lmesh.nno)
	    myerror("mismatch of lmesh.nno in input file",E);
	  for(i=1;i<=E->lmesh.nno;i++)  {
	    /* 
	       read in velocities and temperature from restart file
	    */
	    if(fscanf(fp,"%f",&v1)!=1)
	      myerror("read error restart pressures: 3",E);
	    E->NP[m][i] = v1;	/* nodal pressure */
	  }	


	}	/* end m loop */
	
	fclose (fp);
	if(rezip)
	  compress_file(input_file);
	/* 
	   assign the pressure to the element centres TWB
	*/
	p_to_centres(E,E->NP,E->P,E->mesh.levmax);
      }
      /* 
	 use this section to include imposed ages at the surface CPC 4/27/00 
	 
	 in general:

	 lith_age=1: only half-space cooling according to age
	             

      */
      if(E->control.lith_age >= 1)   {
	if(E->control.lith_age_time == 0){
	  report(E,"WARNING: restart: no new age dependent structure");
	  record(E,"WARNING: restart: no new age dependent structure");
	}else{
	  /* 
	     
	  opening lithosphere age info every timestep 
	  
	  */
	  get_lith_ages(E);
	  /* assign age based temperature based on a restart */
	  assign_age_based_temperature(E,TRUE,e_4);
	}
      }
      /* 
	 END RESTART LOOP
      */
    } else { 
      
      /* 

      NOT USING A RESTART FILE OR POSTPROCESSING, GENERAL INITIAL 
      TEMPERATURE BOUNDARY CONDITION

      */
      if(E->control.lith_age >= 1 )   {
	if(E->control.ggrd_tinit){
	  fprintf(stderr,"p %i: ggrd_tinit and lith_age at same time?!\n",
		  E->parallel.me);
	  exit(-1);
	}
	/* 

	lithospheric age control
	used if the lithosphere age is given in an input file. CPC 1/28/00 

	*/
	if(E->control.lith_age_time == 1)   { 
	  /* 
	     we are opening lithosphere age info at every timestep -
	     the naming is different
	  */
	  age=find_age_in_MY(E);
	  sprintf(input_file,"%s%0.0f",E->control.lith_age_file,
		  age);
	  if(E->parallel.me==0)  {
	    fprintf(E->fp,"%s %s\n","Initial Lithosphere age info:",input_file);
	  }
	} else {     
	  /* 
	     just open lithosphere age info here and only once
	  */
	  sprintf(input_file,"%s",E->control.lith_age_file);
	}
	fp1=fopen(input_file,"r");
	if (fp1 == NULL) {
          fprintf(E->fp,"(Convection.c #3b) Cannot open %s\n",input_file);
          exit(8);
	}
	/* 
	   read age from file 
	*/
	read_gsnode_age_from_file(E,fp1,input_file);
        fclose(fp1);
	if(E->parallel.me==0)  
	  fprintf(stderr,"convection_initial_temperature: read %s, max age: %g\n",
		  input_file,E->init_max_age*E->data.scalet_Ma);
	/* 
	   assign the temperature based on the surface age field, set the 
	   restart flag to FALSE since we are truely initializing the fields
	*/
	assign_age_based_temperature(E,FALSE,e_4);
	
      	/* end lith_age specified */
      }else{
	if(E->control.ggrd_tinit){
#ifdef USE_GGRD
	  if(E->parallel.me==0)  
	    fprintf(stderr,"convection_initial_temperature: using GMT grd files for temperatures\n");
	  /* 
	     
	  read in tempeatures/density from GMT grd files
	  

	  */
	  /* initialize the GMT grid files */
	  ggrd_grdtrack_init_general(TRUE,E->control.ggrd_tinit_gfile, 
				     E->control.ggrd_tinit_dfile,"",
				     &E->control.ggrd_tinit_d,TRUE);
	  /* 
	     
	  interpolate densities to temperature

	  */
	  for(m=1;m <= E->sphere.caps_per_proc;m++)
	    for(i=1;i<=E->lmesh.noy;i++)  
	      for(j=1;j<=E->lmesh.nox;j++) 
		for(k=1;k<=E->lmesh.noz;k++)  {
		  ii = k + E->lmesh.nzs - 1;
		  node=k+(j-1)*E->lmesh.noz+
		    (i-1)*E->lmesh.nox*E->lmesh.noz;
		  ggrd_grdtrack_interpolate_rtp((double)E->sx[m][3][node],
						(double)E->sx[m][1][node],
						(double)E->sx[m][2][node],
						E->control.ggrd_tinit_d,&tadd);
		  /* assign temperature */
		  E->T[m][node] = 0.5 + tadd * E->control.ggrd_tinit_scale;
		}
	  /* free the structure, not needed anymore */
	  ggrd_grdtrack_free_gstruc(E->control.ggrd_tinit_d);

	  /* 
	     end temperature/density from GMT grd init
	  */
#else
	  fprintf(stderr,"p %i: error, need to use GGRD\n",E->parallel.me);
	  exit(-1);

#endif
	}else{

	  /* 
	     
	  no lith age specified, no restart, no T from grd files
	  
	  assign a linear gradient 

	  */
	  if((E->mesh.toptbc==0)&&(E->mesh.bottbc==0)) 
	    myerror("initial_temperature: can not deal with two heat flux boundary conditions",E);
	  if(E->mesh.toptbc==0)
	    myerror("initial_temperature: top heat flux BC not implemented yet",E);
	  if(E->mesh.bottbc == 1){
	    /* bottom has specified temperature */
	    tbot =  E->control.TBCbotval;
	  }else{
	    /* 
	       bottom has specified heat flux
	       start with unity bottom temperature
	    */
	    tbot = 1.0;
	  }
	  tgrad = (E->control.TBCtopval + tbot)/ (E->sphere.ro - E->sphere.ri);
	  for(m=1;m<=E->sphere.caps_per_proc;m++)
	    for(i=1;i<=E->lmesh.noy;i++)  
	      for(j=1;j<=E->lmesh.nox;j++) 
		for(k=1;k<=E->lmesh.noz;k++)  {
		  ii = k + E->lmesh.nzs - 1;
		  /*  local node number */
		  node=k+(j-1)*E->lmesh.noz+
		    (i-1)*E->lmesh.nox*E->lmesh.noz;
		  r1 =  E->sx[m][3][node] - E->sphere.ri;
		  E->T[m][node] = tbot - tgrad * r1 ;
		}
	}
      }
      /* 

	 read in the number of additional perturbations 
	 
      */
      input_int("num_perturbations",&number_of_perturbations,"0,0,32",m);
      if (number_of_perturbations > 0) {
	if(E->control.ggrd_tinit || (E->control.lith_age >= 1)){
	  fprintf(stderr,"p %i: WARNING: ggrd_tinit or lith_age but also perturbation?\n",
		  E->parallel.me);
	}

	/* 
	   ADD TEMPERATURE PERTURBATIONS

	   this was changed as an addition here
	   
	*/
        m = E->parallel.me;
	tlen = M_PI / (E->control.theta_max - E->control.theta_min);
	flen = M_PI / (E->control.fi_max - E->control.fi_min);
	
	/* This part put a temperature anomaly at depth where the global 
	   node number is equal to load_depth. The horizontal pattern of
	   the anomaly is given by spherical harmonic ll & mm. */
	if (! input_float_vector("perturbmag",number_of_perturbations,perturb_mag,m) ) {
	  fprintf(stderr,"Missing input parameter: 'perturbmag'\n");
	  parallel_process_termination();
	}
	if (! input_int_vector("perturbm",number_of_perturbations,perturb_mm,m) ) {
	  fprintf(stderr,"Missing input parameter: 'perturbm'\n");
	  parallel_process_termination();
	}
	if (! input_int_vector("perturbl",number_of_perturbations,perturb_ll,m) ) {
	  fprintf(stderr,"Missing input parameter: 'perturbml'\n");
	  parallel_process_termination();
	}
	/* 	  if (! input_int_vector("perturblayer",number_of_perturbations,load_depth,m) ) { */
	/* 	    fprintf(stderr,"Missing input parameter: 'perturblayer'\n"); */
	/* 	    parallel_process_termination(); */
	/* 	  } */
	for(m=1;m<=E->sphere.caps_per_proc;m++)
	  for(i=1;i<=E->lmesh.noy;i++)  
	    for(j=1;j<=E->lmesh.nox;j++) 
	      for(k=1;k<=E->lmesh.noz;k++)  {
		ii = k + E->lmesh.nzs - 1;
		/*  local node number */
		node=k+(j-1)*E->lmesh.noz+
		  (i-1)*E->lmesh.nox*E->lmesh.noz;
		t1 = (E->sx[m][1][node] - E->control.theta_min) * tlen;
		f1 = (E->sx[m][2][node] - E->control.fi_min) * flen;
		r1 =  E->sx[m][3][node] - E->sphere.ri;
		/* the background gradient is now taken care of above TWB */
		for (p=0; p < number_of_perturbations; p++) {
		  mm = perturb_mm[p];
		  ll = perturb_ll[p];
		  con = perturb_mag[p];
		  /* calculate perturbation as addition */
		  tadd = con*cos(ll*f1)*cos(mm*t1)*
		    sin(M_PI*r1/(E->sphere.ro - E->sphere.ri));
		  E->T[m][node] += tadd;
		  /* limit 0 ... T ... 1 */
		  E->T[m][node] = max(min(E->T[m][node], 1.0), 0.0);
		}
	      }

      }	/* end for added perturbations */
      
    }   /* end for regular init (no restart)  */
    /* 
       assign T as boundary condition
    */
    assign_T_as_TBC(E,e_4);
    /* 
       make sure all the top/bottom boundary conditions are matched
    */
    temperatures_conform_bcs(E);
    
    if (E->control.verbose)  {
      /* 
	 
       output 
      
      */
      report(E,"writing initial temperatures to logfile");
      fprintf(E->fp_out,"output_temperature\n");
      for(m=1;m<=E->sphere.caps_per_proc;m++)        {
	fprintf(E->fp_out,"for cap %d\n",E->sphere.capid[m]);
	for (j=1;j<=E->lmesh.nno;j++)
	  fprintf(E->fp_out,"X = %.6e Z = %.6e Y = %.6e T[%06d] = %.6e \n",
		  E->sx[m][1][j],
		  E->sx[m][2][j],
		  E->sx[m][3][j],j,
		  E->T[m][j]);
      }
      fflush(E->fp_out);
    }
 

    return; 
}


/* 

assign a temperature based on age

if this is called from a restart file, 
this will reassign the temperatures within the lithosphere


*/
/*  #define DEBUG */
void assign_age_based_temperature(E, from_restart, e_4)
     int from_restart;
     struct All_variables *E;
     double e_4;
{

  int i,j,k,m,nodes,node;
  float r1;

 
  if(E->parallel.me==0)
    if(from_restart == 1)
      fprintf(stderr,"assign_age_based_temperature: initializing from restart for lith_age=%i\n",
	      E->control.lith_age);
    else
      fprintf(stderr,"assign_age_based_temperature: initializing for lith_age=%i\n",
	      E->control.lith_age);
  if(!E->control.age_init)
    myerror("assign_age_based_temperature: error: ages not initialized",E);

  for(m=1;m<=E->sphere.caps_per_proc;m++)
    for(i=1;i<=E->lmesh.noy;i++)  
      for(j=1;j<=E->lmesh.nox;j++) 
	for(k=1;k<=E->lmesh.noz;k++)  {
	  /* surface node number */
	  nodes = E->lmesh.nxs-1+j+(E->lmesh.nys+i-2)*
	    E->mesh.nox;
	  /* local node number */
	  node = k + (j-1) * E->lmesh.noz + 
	    (i-1) * E->lmesh.nox * E->lmesh.noz;
	  /* assign r/r_E coordinate */
	  r1 = E->sx[m][3][node];
	  if(from_restart)
	    E->T[m][node] = half_space_temp(E,r1,E->age_t[nodes],E->T[m][node],FALSE);
	  else
	    E->T[m][node] = half_space_temp(E,r1,E->age_t[nodes],E->T[m][node],TRUE);
#ifdef DEBUG
	  fprintf(stderr,"lon: %11g lat: %11g z: %11g age: %11g T: %11g\n",
		  RAD2DEG(E->sx[m][2][node]),90.-RAD2DEG(E->sx[m][1][node]),
		  (E->sphere.ro-E->sx[m][3][node])*
		  E->data.radius_km,
		  E->age_t[nodes]*E->data.scalet_Ma,
		  E->T[m][node]);
#endif
		  
	}
}

/* 

assign half space cooling temperature, 
if r1 > E->sphere.ro-E->control.lith_age_depth, and age >= 0 

assumes 0 ... r1 ... 1

compute the non-dim temperature in a half-space at 
radius (normalized by radius_earth) r1 and 
non-dim age (age/t_c) 

for other cases, there are other options. e.g:

is nd_age is negative, will assume continental region and 
assign different values



*/
float half_space_temp(E, r1, nd_age, orig_temp, initial_assign)
     struct All_variables *E;
     float r1,nd_age,orig_temp;
     boolean initial_assign;
{
  float temp,zp,z,x;
  if(nd_age >= 0){		
    /* 
       real age, underneath oceanic plate 
    */
    if(r1 >= (E->sphere.ro - E->control.lith_age_depth) ) { 
      /* 
	 within lithospheric limits, assign half space cooling
	 at all times!
      */
      z = (E->sphere.ro-r1)* E->data.radius_km; /* depth in km */
      zp = z/E->data.layer_km;	/* normalized depth */
      if(nd_age == 0.0){
	/* 
	   zero age, at ridge
	*/
	if(fabs(z) < 1e-4)		/* surface */
	  temp = 0.0;
	else
	  temp = E->control.mantle_temp;
      }else{
	/* 
	   away from ridge, normal half-space cooling
	*/
	x = zp * 0.5 / sqrt(nd_age);
	temp = E->control.mantle_temp * erf(x);
      }
    }else{
      /* 
	 below lithosphere, there will a uniform temperature initially
      */
      if(initial_assign){
	/* initial temperature assignment  */
	temp = E->control.mantle_temp;
      }else{
	/* called from restart or at later timestep */
	temp = orig_temp;
      }
    }
  }else{
    /* 
       underneath continental plate 
    */
    if(initial_assign){
      /*
	if age is smaller than zero,
	assign a temperature gradient
	depending on the depth, if first init
      */
      /* E->T[m][node] = E->control.age_tgrad * (E->sphere.ro - r1)/E->sphere.ri; */
      temp = E->control.mantle_temp;
      /*	      fprintf(stderr,"r: %g grad: %g T: %g\n",
		      (E->sphere.ro - r1), E->control.age_tgrad,E->T[m][node] ); */
    }else{
      temp = orig_temp;
    }
  }
  return temp;
}
/* 
   assign the global array temperature to the caps as thermal boundary conditions
   only in interior

*/
void assign_T_as_TBC(E,e_4)
     struct All_variables *E;
     float e_4;
{
  int i,j,k,m,node;
  float r1;
  /* assign to cap  */
  for(m=1;m<=E->sphere.caps_per_proc;m++)
    for(i=1;i<=E->lmesh.noy;i++)  
      for(j=1;j<=E->lmesh.nox;j++) 
	for(k=1;k<=E->lmesh.noz;k++)  {
	  /* local node number */
	  node=k+(j-1)*E->lmesh.noz+
	    (i-1)*E->lmesh.nox*E->lmesh.noz;
	  r1=E->sx[m][3][node];

	  if(fabs(r1-E->sphere.ro)>=e_4 &&   fabs(r1-E->sphere.ri)>=e_4)  {
	    E->sphere.cap[m].TB[1][node] = E->T[m][node];
	    E->sphere.cap[m].TB[2][node] = E->T[m][node];
	    E->sphere.cap[m].TB[3][node] = E->T[m][node];
	  }
	}
}

/* 
   
routine dealing with reading ages from file,
either only initially, or for each timestep
calls other routines for reading, mainly establishes the 
existance of the age_t array

*/
void get_lith_ages(E)
     struct All_variables *E;
{
  static int been_here = 0;
  static int local_solution_cycles = -1;
  static float old_time = -6e10; /* time at which ages were last read */
  int gnox,gnoy,output;
  FILE *fp1;
  
  char output_file[CHARBUF_SIZE];
  output=0;
  gnox=E->mesh.nox;
  gnoy=E->mesh.noy;
  
  if(E->control.lith_age >= 1)   { /* if specifying lithosphere age */
    if(E->control.lith_age_time == 1)   {  
      /* 
	 
      to open files every timestep 

      */
      if(been_here == 0) {
	/* first call */
	E->age_t=(float *) safe_malloc((gnox*gnoy+1)*sizeof(float));
	local_solution_cycles = E->monitor.solution_cycles-1;
	been_here++;

      }
      if (local_solution_cycles < E->monitor.solution_cycles) {
	output = 1;
	local_solution_cycles++; /* update so that output only happens once */
      }
      /* 
	 mode 2 for reading ages, this will interpolate 
      */
      if((old_time == -6e10)||(find_age_in_MY(E) != old_time)){
	read_input_files_for_timesteps(E,2,output);
	old_time = find_age_in_MY(E);
 	E->control.age_init = TRUE;
      }
      /* end E->control.lith_age_time == true */
    } else {  
      /* 
	 otherwise, just open for the first timestep 
      */
      /* NOTE: This is only used if we are adjusting the boundaries */
      if(been_here == 0) {
	E->age_t=(float*) safe_malloc((gnox*gnoy+1)*sizeof(float));
	been_here++;
      } /* end been_here */
      if((old_time == -6e10) || (find_age_in_MY(E) != old_time)){
	
	sprintf(output_file,"%s",E->control.lith_age_file);
	fp1=fopen(output_file,"r");
	if (fp1 == NULL) {
	  sprintf(output_file,"get_lith_ages: Can't open file %s",E->control.lith_age_file);
	  myerror(output_file,E);
	}
	/* read in the ages from file */
	read_gsnode_age_from_file(E, fp1, output_file);
	fclose(fp1);
	old_time = find_age_in_MY(E); 
	if(E->parallel.me==0)
	  fprintf(stderr,"get_lith_ages: read %s, max age: %g Ma at time: %g Ma\n",
		  output_file,E->init_max_age*E->data.scalet_Ma,
		  old_time);
	E->control.age_init = TRUE;
       }


    } /* end E->control.lith_age_time == false */ 
  } /* end lith_age >= 1 */
}  
/* 


read in ages in millions of years from fp1 
using all global surface nodes, and assign to  
E->age_t[node], which will be non-dimensionalized age, ie.
we scale the age after reading it in

*/
void read_gsnode_age_from_file(E, fp1, input_file)
     struct All_variables *E;
     FILE *fp1;
     char *input_file;
{
  int node,i,j;
  char mstring[1000];
  for(E->init_max_age=-1e20,
	i=1;i<=E->mesh.noy;i++) {
    for(j=1;j<=E->mesh.nox;j++) {
      /* global surface nodes */
      node=j+(i-1)*E->mesh.nox;
      if(fscanf(fp1,"%f",&(E->age_t[node])) != 1){
	sprintf(mstring,"error reading from %s, line %i",
		input_file,node);
	myerror(mstring,E);
      }
      /* normalize by diffusive time scale  */
      E->age_t[node] /= E->data.scalet_Ma;
      /* find maximum age */
      if(E->age_t[node] > E->init_max_age)
	E->init_max_age = E->age_t[node] ;
    }
  }
 }
