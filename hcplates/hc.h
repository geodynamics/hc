/* 
  CO'N: abridged working version for plates
 
   header file for Hager & O'Connell experimental code

   Thorsten (twb@usc.edu)

   $Id: hc.h,v 1.12 2006/05/01 17:46:18 becker Exp becker $

*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

/* 

general variable type defines 

*/
#ifndef hc_boolean
#define hc_boolean unsigned short
#endif
#ifndef HC_PREC			/* 
					   precision for most C functions
					*/
#define HC_PREC double
#define HC_FLT_FORMAT "%lf"
#define HC_TWO_FLT_FORMAT "%lf %lf"
#define HC_EPS_PREC 5e-15
#endif

#ifndef HC_CPREC
#define HC_CPREC HC_PREC
#endif

#define HC_CHAR_LENGTH 300		/* length of char arrays */

#define HC_BIN_PREC float	/* precision for binary I/O */


#define HC_PI 3.1415926535897932384626433832795

#ifndef SQRT_TWO 
#define SQRT_TWO 1.41421356237309504880168872420970 
#endif
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif
#ifndef PIOVERONEEIGHTY
#define PIOVERONEEIGHTY  0.017453292519943295769236907684886 
#endif
#ifndef ONEEIGHTYOVERPI
#define ONEEIGHTYOVERPI  57.295779513082320876798154814105
#endif
#ifndef DEG2RAD
#define DEG2RAD(x) ((x) * PIOVERONEEIGHTY)
#endif
#ifndef RAD2DEG
#define RAD2DEG(x) ((x) * ONEEIGHTYOVERPI)
#endif

#ifndef FLT_MAX 
#define FLT_MAX 1e20
#endif
#ifndef FLT_MIN
#define FLT_MIN -1e20
#endif 

#define THETA2LAT(x) ( (90.0 - (x)*ONEEIGHTYOVERPI) )
#define LAT2THETA(x) ( (90.0 - (x))*PIOVERONEEIGHTY )
#define PHI2LON(x) ( RAD2DEG(x) )
#define LON2PHI(x) ( DEG2RAD(x) )
/* 

other constants

*/
// now taken from earth model
#define HC_RE_KM 6371.0		/* radius(Earth) in [km] */

//#define HC_RCMB_ND 0.546225    /* non-dim radius, ~10km above
//				  CMB  */

#define HC_TIMESCALE_YR 1e6	/* timescale [yr] */
/* 
   convert depth (>0) in km to non-dimensionalizd radius
*/
#define HC_ND_RADIUS(x) (1.0-((x)/HC_RE_KM))
/* 
   the other way around
*/
#define HC_Z_DEPTH(x) ((HC_RE_KM * (1.0-(x))))

