#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define HC_MAX_LAYERS 100 /* Hope this is enough radial layers */
#define MAX_NO_PLATE 30
#define MAX_LDIM 100
#define MAX_NPDIM 90
#define MAX_KDIM 515	/* Note: having this too high gives seg faults But we need it... weird.*/

#define HC_CHAR_LENGTH 300		/* length of char arrays */

struct hc_findplates{
	char boundary_file[HC_CHAR_LENGTH];
	char data_filename[HC_CHAR_LENGTH];
	char output_file[HC_CHAR_LENGTH];
};

/* findplate.c */
//void hc_findplates(int, char **);
void hc_findplates_command_line(int, char **,
			    struct hc_findplates *);
void hc_findplates_advance_argument(int *,int, char **);
void hc_findplates_init(struct hc_findplates *);

