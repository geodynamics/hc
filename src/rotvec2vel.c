#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
/* 

utility program, here for historical reasons


 */
#define DEF_FIXED -1

//
// reads in rotation vector file and xy code file with plate codes
//
// output is lon lat vp vt
//
// these are the coded plates as in the lon lat code tripel that is 
// read from stdin
//

// some factors
#define PIOVERONEEIGHTY 0.0174532925199433
#define PIHALF 1.5707963267949
#define PI 3.14159265358979
//#define VELFACTOR  (6371009.0/1.0e6*1.0e2) 
//#define PIOVERONEEIGHTY_TIMES_VELFACTOR (PIOVERONEEIGHTY*VELFACTOR) 
#define PIOVERONEEIGHTY_TIMES_VELFACTOR 11.1195083724191
FILE *rv_myopen(const char *,const char *);


int 
main (argc, argv)
int argc;
char *argv[];
{
  int i,j,k,nplt,nrp,fixed_plate,*stats,minplate,maxplate,normalize=0,unity_vec=0,
    *assigned,hit,*name,code,code3;
  double *rotvel,*points,x,y,z,stheta,cphi,ctheta,sphi,vx,vy,vz,vt,vp,
    redvec[3],length;
  char *allplates,pname[4],tmpstring[800];
  FILE *vel;
  int model_code = 1,coded_plates,code_length;
  if(argc == 2)
    fixed_plate=DEF_FIXED;
  if(argc > 2)
    sscanf(argv[2],"%i",&fixed_plate);
  if(argc == 4){
    if(strcmp(argv[3],"bird02")==0)
      model_code = 2;
    if(strcmp(argv[3],"morvel")==0)
      model_code = 3;
   if(strcmp(argv[3],"geodvel")==0)
      model_code = 4;
  }
  if(argc==1){
    fprintf(stderr,"%s rotvector_file [fixed_plate] [model, nuvel]\n",argv[0]);
    fprintf(stderr,"\t reads points in (x y plate_code) format from stdin\n");
    fprintf(stderr,"\t and calculates velocities (vp,vt) from rotvector_file\n");
    fprintf(stderr,"\t rotvectorfile is in\n");
    fprintf(stderr,"\t   platename_1 wx_1 wy_1 wz_1\n");
    fprintf(stderr,"\t   platename_2 wx_2 wy_2 wz_2...\n");
    fprintf(stderr,"\t format, wi_j in input is in [degrees/Myr]\n");
    fprintf(stderr,"\t plate_code is a number from 1...n in the rotvectorfile\n");
    fprintf(stderr,"\t fixed_plate is the reference plate number (1,2,...,n) in the poles list,\n");
    fprintf(stderr,"\t  if set to -1, no change in motions.\n");
    fprintf(stderr,"\t  if set to -2, normalize all omega vectors to 1 deg/Myr length\n");
    fprintf(stderr,"\t  if set to -3, use 1 deg/Myr w_x vector, rest zero\n");
    fprintf(stderr,"\t  if set to -4, use 1 deg/Myr w_y vector, rest zero\n");
    fprintf(stderr,"\t  if set to -5, use 1 deg/Myr w_z vector, rest zero\n");
    fprintf(stderr,"\t  By default set to %i.\n",DEF_FIXED);
    fprintf(stderr,"\t output is in\n\tlon lat v_p v_t\n\tformat in cm/yr\n");
    fprintf(stderr,"\t model can be nuvel, bird02, morvel, or geodvel\n");
    exit(-1);
  }
  if(model_code == 1){
    /* nuvel  */
    fprintf(stderr,"%s: init for NUVEL\n",argv[0]);
    code_length = 3;
    coded_plates = 14;
    allplates = (char *)malloc(sizeof(char)*(coded_plates*(code_length+2)));
    sprintf(allplates,"%s","ANT AUS AFR PAC EUR NAM NAZ COC CAR ARA PHI SAM IND JDF ");
  }else if(model_code ==2){
    /* bird02 */
    fprintf(stderr,"%s: init for Bird 2002\n",argv[0]);
    code_length = 2;
    coded_plates = 52;
    allplates = (char *)malloc(sizeof(char)*(coded_plates*(code_length+2)));
    sprintf(allplates,"%s","AF AM AN AP AR AS AT AU BH BR BS BU CA CL CO CR EA EU FT GP IN JF JZ KE MA MN MO MS NA NB ND NH NI NZ OK ON PA PM PS RI SA SB SC SL SO SS SU SW TI TO WL YA ");
  }else if(model_code==3){
    /* morvel */
    fprintf(stderr,"%s: init for MORVEL-NNR56 (Argus et al., 2011)\n",argv[0]);
    code_length = 2;
    coded_plates = 56;
    allplates = (char *)malloc(sizeof(char)*(coded_plates*(code_length+2)));
    sprintf(allplates,"%s","am an AP ar AS AT au BH BR BS BU ca CL co cp CR EA eu FT GP in jf JZ KE lw MA MN MO mq MS na nb NB ND NH NI nz OK ON pa PM ps ri sa SB sc SL sm sr SS su sw TI TO WL yz ");

  }else if(model_code == 4){
    /* geodvel  */
    fprintf(stderr,"%s: init for GEODVEL\n",argv[0]);
    code_length = 3;
    coded_plates = 10;
    allplates = (char *)malloc(sizeof(char)*(coded_plates*(code_length+2)));
    sprintf(allplates,"%s","ANT ARA AUS IND NAM NAZ NUB PAC SAM SOM ");
  }else{
    fprintf(stderr,"%s: error, model %i not defined\n",argv[0],model_code);
    exit(-1);
  }
  assigned = (int *)calloc(coded_plates,sizeof(int));
  name = (int *)calloc(coded_plates,sizeof(int));
  rotvel=(double *)malloc(sizeof(double)*3*coded_plates);
  stats=(int *)calloc(coded_plates,sizeof(int));

  /*  */

  vel=rv_myopen(argv[1],"r");

  fprintf(stderr,"%s: reading from %s, fixed plate: %i\n",argv[0],argv[1],fixed_plate);
  nplt=0;
  // skip initial comment lines
  while(fscanf(vel,"%[^\n]%*1c",tmpstring)!=EOF){
    if(tmpstring[0]!='#' && tmpstring[0]!=' '){
      if(sscanf(tmpstring,"%s %lf %lf %lf",pname,&vx,&vy,&vz)==4){
	for(hit=i=0;i < coded_plates;i++)
	  if(strncmp(pname,(allplates+(code_length+1)*i),code_length)==0){
	    if(assigned[i]){
	      fprintf(stderr,"%s: plate %s was assigned already, is it twice in list?\n",
		      argv[0],pname);
	      exit(-1);
	    }
	    assigned[i]=hit=1;
	    k=i*3;
	    *(rotvel+k++)=vx;
	    *(rotvel+k++)=vy;
	    *(rotvel+k)=vz;

	    name[nplt]=i;
	  }
	if(!hit){
	  fprintf(stderr,"%s: could not find internal code for input plate code %s\n",
		  argv[0],pname);
	  fprintf(stderr,"%s: will ignore this code\n",argv[0]);
	  fprintf(stderr,"%s: internal codes are:\n",argv[0]);
	  for(j=0;j < coded_plates;j++){
	    strncpy(pname,(allplates+j*(code_length+1)),code_length);
	    fprintf(stderr,"%s: plate nr %i, code: %s\n",argv[0],j+1,pname);
	  }
	}
	nplt++;
      }
    }
  }
  fclose(vel);
  fprintf(stderr,"%s: read in %i poles\n",argv[0],nplt);

  // assign zero rotation vectors to all plates we did not find in the list
  for(i=0,k=0;i < coded_plates;i++,k+=3)
    if(!assigned[i])
      for(j=0;j<3;j++)
	rotvel[k+j]=0.0;
  
  if(fixed_plate > 0){
    fprintf(stderr,"%s: input plate number %i(%i) will be fixed\n",
	    argv[0],fixed_plate,name[fixed_plate-1]);
    if(fixed_plate > nplt){
      fprintf(stderr,"%s: fixed plate number (%i) greater than plate numbers (1...%i)\n\n",
	      argv[0],fixed_plate,nplt);
    }
    for(i=0;i < 3;i++)
      redvec[i]= *(rotvel+name[fixed_plate-1]*3+i);
  }else{
    for(i=0;i < 3;i++)
      redvec[i]=0.0;
    if(fixed_plate == -2){
      normalize=1;
      fprintf(stderr,"%s: normalizing all vectors\n",argv[0]);
    }else if(fixed_plate < -2 && fixed_plate > -6){
      unity_vec = -fixed_plate - 2;
      fprintf(stderr,"%s: using unity for %ith component of omega, zero else\n",
	      argv[0],unity_vec);
    }
  }
  //
  // convert rotation poles
  // 
  for(k=i=0;i < coded_plates;i++,k+=3){
    if(assigned[i]){
      strncpy(pname,(allplates+(code_length+1)*i),code_length);
      // original vector
      fprintf(stderr,"%s: plate %3i %3s wx: %12g wy: %12g wz: %12g",
	      argv[0],i+1,pname,*(rotvel+k),*(rotvel+k+1),*(rotvel+k+2));
      for(j=0;j < 3;j++)	/* correction */
	*(rotvel+k+j) -= redvec[j];

      if(i == name[fixed_plate-1])
	fprintf(stderr," [deg/Myr] (fixed)\n");
      else
	fprintf(stderr," [deg/Myr]\n");
     
      // normalize
      for(length=0.0,j=0;j<3;j++)
	length += *(rotvel+k+j)* *(rotvel+k+j);
      length=sqrt(length);
      if(normalize)
	for(j=0;j < 3;j++)
	  *(rotvel+k+j) /= length;
      // or assign unit length
      if(unity_vec){
	for(j=0;j < 3;j++){
	  if(unity_vec-1==j){
	    *(rotvel+k+j)=1.0;
	  }else{
	    *(rotvel+k+j)=0.0;
	  }
	}
      }
      for(length=0.0,j=0;j < 3;j++)
	length+= *(rotvel+k+j)* *(rotvel+k+j);
      length=sqrt(length);
      // converted 
      fprintf(stderr,"%s:               wx: %12g wy: %12g wz: %12g [deg/Myr]\n",
	      argv[0],*(rotvel+k),*(rotvel+k+1),*(rotvel+k+2));
      fprintf(stderr,"%s:           length: %12g [deg/Myr]\n",argv[0],length);

      // rescale the omega vector such that it will yield cm/yr when
      // multiplied with r vector of unit length
      *(rotvel+k)  *= PIOVERONEEIGHTY_TIMES_VELFACTOR; 
      *(rotvel+k+1)*= PIOVERONEEIGHTY_TIMES_VELFACTOR; 
      *(rotvel+k+2)*= PIOVERONEEIGHTY_TIMES_VELFACTOR;
    }
  }
  minplate=1e6;
  maxplate=-1;
  //
  // read in lon lat code tripels
  //
  nrp=1;
  points=(double *)malloc(sizeof(double)*3*nrp);
  k=0;
  while(fscanf(stdin,"%lf %lf %lf",(points+k),(points+k+1),(points+k+2))==3){
    if(*(points+k)<0.0)
      *(points+k) += 360.0;
    *(points+k)  *=PIOVERONEEIGHTY;
    *(points+k+1)*=PIOVERONEEIGHTY;
    /* fix poles */
    *(points+k+1)= PIHALF- *(points+k+1);
    if(fabs(points[k+1]) < 1e-5)
      points[k+1] = 1e-5;
    if(fabs(points[k+1] -PIHALF) < 1e-5)
      points[k+1] = PIHALF - 1e-5;

    // change code here fore more efficient assignment
    *(points+k+2)-=1;
    if(*(points+k+2)<minplate)
      minplate= *(points+k+2);
    if(*(points+k+2)>maxplate)
      maxplate= *(points+k+2);
    nrp++;
    k += 3;
    if((points=((double *)realloc(points,sizeof(double)*3*nrp)))==NULL){
      fprintf(stderr,"%s: memerror, too many points for x array, %i\n",argv[0],nrp);
      exit(-1);
    }
  }
  nrp--;

  fprintf(stderr,"%s: read %i points for velocities, calculating...\n",argv[0],nrp);
  fprintf(stderr,"%s: input plate codes run from %3i to %3i\n",argv[0],minplate+1,maxplate+1);
  if(minplate<0){
    fprintf(stderr,"%s: expect code to start from 1\n",argv[0]);exit(-1);
  }
  if(maxplate>nplt-1 || maxplate-minplate>nplt){
    fprintf(stderr,"%s: could only read %i plate codes from rotation vector file\n",
	    argv[0],nplt);
    fprintf(stderr,"%s: will set velocities for all other codes to zero\n",argv[0]);
  }

  fprintf(stderr,"%s: now output...\n",argv[0]);

#define WX ( *(rotvel+code3)   )
#define WY ( *(rotvel+code3+1) )
#define WZ ( *(rotvel+code3+2) )

  for(i=k=0;i<nrp;i++,k+=3){
    code=(int)*(points+k+2);
    if(assigned[code]){
      code3 = code * 3;
      // we have included the radius factor in the omega vector,
      // hence all locations are calculated with radius unity
      x=(cphi=cos(*(points+k)))*(stheta=sin(*(points+k+1)));
      y=(sphi=sin(*(points+k)))*stheta;
      z=(ctheta=cos(*(points+k+1)));
      vx = WY * z - WZ * y;
      vy = WZ * x - WX * z;
      vz = WX * y - WY * x;
      *(stats+code)+=1;
      /* v_theta goes southward */
      vt= ctheta*cphi*vx + ctheta*sphi*vy - stheta*vz;
      /* v_phi goes eastward */
      vp= -sphi*vx + cphi*vy;
      fprintf(stdout,"%12g %12g %20g %20g\n",
	      (1.0/PIOVERONEEIGHTY)* *(points+k),
	      90.0-(1.0/PIOVERONEEIGHTY)* *(points+k+1),
	      vp,vt);
    }else{
      ;
      //fprintf(stdout,"%12g %12g %20g %20g\n",
      //	      (1.0/PIOVERONEEIGHTY)* *(points+k),
      //      90.0-(1.0/PIOVERONEEIGHTY)* *(points+k+1),0.,0.);
    }
  }
  for(i=0;i < coded_plates;i++)
    if(assigned[i]){
      strncpy(pname,(allplates+(code_length+1)*i),3);
      fprintf(stderr,"%s: %3s (%3i) (%9.2e %9.2e %9.2e) assigned %7i/%7i, %15.3f%% \n",
	      argv[0],pname,i+1,
	      *(rotvel+i*3)/PIOVERONEEIGHTY_TIMES_VELFACTOR,
	      *(rotvel+i*3+1)/PIOVERONEEIGHTY_TIMES_VELFACTOR,
	      *(rotvel+i*3+2)/PIOVERONEEIGHTY_TIMES_VELFACTOR,
	      *(stats+i),nrp,(double)(*(stats+i))/(double)nrp*100.0);
    }
  fprintf(stderr,"%s: done\n",argv[0]);
  return 0;
}
FILE *
rv_myopen (name, modus)
const char *name;
const char *modus;
{
  FILE *tmp;
  if( (tmp = (FILE *)fopen(name,modus)) == NULL)
    {	fprintf(stderr,"Error opening \"%s\" for \"%s\".\n",name,modus);
	fprintf(stderr,"Exiting.\n");
	exit(-1);};
  return ((FILE *)tmp);
}	
