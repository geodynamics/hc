#include "hc.h"



#ifndef __GMT_INCLUDED__
#define _GMT_MATH_H		/* gmt_math.h led to clashes with some other codes */
#include "gmt.h"
#include "gmt_bcr.h"
#include "gmt_math.h"
#define __GMT_INCLUDED__
#endif
#define irint(x) ((int)rint(x))

/* 
   dealing with velocity grids 

*/
#include "ggrd_grdtrack_util.h"
#include "hc_ggrd_auto_proto.h"
