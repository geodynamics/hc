#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "hcplates.h"

/*

*/

int main(int argc, char **argv)
{
	struct hc_plates_params plates[1];
	struct hc_plates_arrays A[1];
	char file1[HC_CHAR_LENGTH], file2[HC_CHAR_LENGTH], file3[HC_CHAR_LENGTH], file4[HC_CHAR_LENGTH];
	char file5[HC_CHAR_LENGTH];

	printf("Init parameters\n");
	hc_init_params(plates);
	
	printf("Reading command line args...\n");
	hcplates_command_line(argc,argv,plates);
	
	// Read parameter file somewhere here 
	
	fprintf(stderr,"Load file is %s\n",plates->loadfile);
	fprintf(stderr,"Unitrot file is %s\n",plates->unitrotfile);
	fprintf(stderr,"Plate map file is %s\n",plates->platemapfile);
	fprintf(stderr,"Output file is %s\n",plates->outputfile);
	fprintf(stderr,"Pole output file is %s\n",plates->polesfile);

	printf("Init arrays\n");
	hc_init_arrays(plates,A);

	printf("Initializing arrays\n");
	zero_arrays(plates,A);

	/* Hardwire filename for now */

	printf("Get loads\n");
	//char filename[]=plates->loadfile;
	strncpy(file1,plates->loadfile,HC_CHAR_LENGTH);
	get_loads(plates,A,file1);

	printf("Read unit rotations coefficient file\n");
	//char file2[]="unitrot.coeffs";
	strncpy(file2,plates->unitrotfile,HC_CHAR_LENGTH);
	read_unitrots_get_coeff(plates,A,file2);

	/*Integrate torques */
	printf("Integrate torques \n");
	//char file3[]="map.plate";
	strncpy(file3,plates->platemapfile,HC_CHAR_LENGTH);
	integrate_torques(plates,A,file3);

	/* Solve for rotations */
	printf("Solving for rotations\n");
	solve_rot(plates,A);
	printf("Rotations solved!\n");
	
	/* Not doing correlations */
	
	/* Output out.forces and out.predrot */
	printf("Output routines:\n");
	//write_output(plates,A);
	printf(" Printing poles to %s \n",plates->polesfile);
	strncpy(file4,plates->polesfile,HC_CHAR_LENGTH);
	hc_poles(plates, A, file4);
	
	printf(" Printing velocity field to %s \n",plates->velgridfile);
	strncpy(file5,plates->velgridfile,HC_CHAR_LENGTH);
	hcplates_velgrid(plates,A,file5);
	
	printf("Done!\n");
	return 0;
}

	
