/* 

Earth model constants 


*/
// now taken from earth model
#define HC_RE_KM 6371.0087714		/* radius(Earth) in [km] (equivalent
					   volume radius of std ellipsoid) */


#define HC_TIMESCALE_YR 1e6	/* timescale [yr] */


#define HC_VISNOR 1e21		/* reference viscosity [Pas] */
#define HC_GACC 10.0e2		/* gravitational acceleration [cm/s2] */
#define HC_CAPITAL_G 6.6742e-11	/* gravitational constant [Nm2/kg2]*/
//6.67408 would be better, leave for backward compatibility?

#define HC_SECYR  3.1556926e7	/* seconds/year  */

/* average mantle density */
#define HC_AVG_DEN_MANTLE 4.4488 /* in g/cm^3 */
/* average core density */
#define HC_AVG_DEN_CORE 11.60101

/* 

 */
#define HC_PVEL_TSTEPS 140	/* number of timesteps for pvel init from file */

/* 

other modeling constants or initial values

*/
/* 

default constant scaling for input density file, e.g. use 0.01 for
seismic tomography models which are given in terms of % variation from
PREM

*/
#define HC_DENSITY_SCALING 0.01

/* default dln v_s/dln rho density scaling */

#define HC_D_LOG_V_D_LOG_D 0.2


/* length scale factor for velocity I/O, in units of m 
   use 0.01, if input/output is to be in cm/yr
*/
#define HC_VEL_IO_SCALE 0.01


#define HC_LMAX_DEFAULT 31

/* 
   for viscosity scans

*/

#define HC_VSCAN_VMAX 3				/* range for viscosities in log space */
#define HC_VSCAN_NLAYER_MAX 4		/* max number of layers, HAS
					   TO BE FOUR FOR NOW */

#define HC_VSCAN_DV0 0.25			/* default spacing for scan */
