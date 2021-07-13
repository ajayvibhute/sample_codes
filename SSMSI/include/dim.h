
/*******************************************************************
            Mission: ASTROSAT
            Project: SSM

          File-Name: dim.h ($Revision: 1.1 $)
#%%            Path: /data1/ravibt/ASTROSAT/SSM/etc
         Created on: 2008 Sep 26 04:56:14 PM (IST: UTC+0530)
   Last Modified on: 2009 Feb 23 11:33:58 AM (IST: UTC+0530)
             Author: Ravishankar B.T. (ravibt@isac.gov.in)
        Supervisors: Dipankar Bhattacharya (dipankar@iucaa.ernet.in)
                     Seetha S (seetha@isac.gov.in)
*******************************************************************/

/*####################################################################

$Id: dim.h,v 1.1 2008/10/27 10:07:00 ravibt Exp ravibt $
$Log: dim.h,v $
Revision 1.1  2008/10/27 10:07:00  ravibt
Initial revision

####################################################################*/

#ifndef ___DIM__H___
#define ___DIM__H___

#include <math.h>



#ifndef M_PI	/* COMPILING WITH -ansi OPTION */
#define M_PI			3.14159265358979323846
#endif
struct infoDetails
{
	int myrank;
};
/* Path for observed shadow file */
#define OBSERVED_SHADOW_PATH "../shadows/underprocess/"
#define OBSERVATION_ID "20110720:0123456789"
#define LOG_FILE "../log/ExecutionLog"

/* D2R = degree --> radian        R2D = radian --> degree */
#define D2R				(M_PI/180.0)
#define R2D				(180.0/M_PI)

#define TWOPI			(2.0*M_PI)
#define PIBYTWO			(M_PI/2.0)


/* unit = mm;   WD = WIDTH;  DET = DETECTOR */
#define MASK_WD_X		60.0
#define DET_WD_X		MASK_WD_X
#define DET_WD_Y		96.0

/* unit = mm */
#define SL_MASK_WD_Y	500.0
#define BO_MASK_WD_Y	634.0

/* unit = mm */
#define SL_HEIGHT		251.5
#define BO_HEIGHT		306.5

/* unit = degree */
#define SL_FOV_LIM_X	(atan(1.0*MASK_WD_X/SL_HEIGHT) * 180.0/M_PI)
#define BO_FOV_LIM_X	(atan(1.0*MASK_WD_X/BO_HEIGHT) * 180.0/M_PI)

/* unit = degree */
#define SL_FOV_LIM_Y	(atan(1.0*(SL_MASK_WD_Y+DET_WD_Y)/(2*SL_HEIGHT)) * 180.0/M_PI)
#define BO_FOV_LIM_Y	(atan(1.0*(BO_MASK_WD_Y+DET_WD_Y)/(2*BO_HEIGHT)) * 180.0/M_PI)

/* unit = degree */
#define SLANT_ANG		12.0

/* unit = degree */
#define CANT_ANG		45.0

#define NUM_WIRES		8
#define NUM_MASK_ELEM	63
#define NUM_DET_BINS	NUM_MASK_ELEM
#define NUM_PATTERNS	6


/* unit = mm */
#define DET_DEPTH		12.0
#define CELL_WD_Y		(1.0*DET_WD_Y/NUM_WIRES)


/* mm */
#define MASK_ORIGIN_Y	(-SL_MASK_WD_Y/2.0)
#define MASK_ORIGIN_X	(-MASK_WD_X /2.0)
#define DET_ORIGIN_X	(-DET_WD_X/2.0)
#define DET_END_X		(DET_WD_X + DET_ORIGIN_X)
#define DET_END_Y		(DET_WD_Y + DET_ORIGIN_Y)
#define DET_ORIGIN_Y	(-DET_WD_Y/2.0)

/*
  REAL DET_END_X considering that PMEW = 0.95. currently not used
  anywhere
*/
#define R_DET_END_X		(DET_END_X - (DET_WD_X-NUM_MASK_ELEM*PMEW))

/* unit = mm,  PMEW = PerMaskElementWidth */
#define PMEW			0.95
#define SL_MASK_THICK	3.0
#define BO_MASK_THICK	3.0

/* flags */
#define BOOM			0
#define SLANT1			1
#define SLANT2			2

#define NUM_CAM			3

/* unit = mm */
#define INTER_PAT_GAP	2.0
#define BO_PAT_WD_Y		((BO_MASK_WD_Y-INTER_PAT_GAP*(NUM_PATTERNS-1))/NUM_PATTERNS)
#define SL_PAT_WD_Y		((SL_MASK_WD_Y-INTER_PAT_GAP*(NUM_PATTERNS-1))/NUM_PATTERNS)

#define MAX_BIN_SIZE_FAC	1.99

/* unit = mm */
#define WINDOW_HEIGHT	6.5

/* unit = mm */
#define MASK_TILT_X_BO	0.0
#define MASK_TILT_X_SL1	0.0
#define MASK_TILT_X_SL2	0.0
#define MASK_TILT_Y_BO	0.0
#define MASK_TILT_Y_SL1	0.0
#define MASK_TILT_Y_SL2	0.0

/* unit = mm */
#define PAT_LLIM_X_BO	DET_ORIGIN_X
#define PAT_LLIM_X_SL1	DET_ORIGIN_X
#define PAT_LLIM_X_SL2	DET_ORIGIN_X

/* unit = mm */
#define PAT_TILT_LLIM_X_BO	0.0
#define PAT_TILT_LLIM_X_SL1	0.0
#define PAT_TILT_LLIM_X_SL2	0.0

#define MASK_Y_ORIGIN_BO	(-BO_MASK_WD_Y/2.0)
#define MASK_Y_ORIGIN_SL	(-SL_MASK_WD_Y/2.0)

#define PAT_TILT_LLIM_Y_BO	0.0
#define PAT_TILT_LLIM_Y_SL1	0.0
#define PAT_TILT_LLIM_Y_SL2	0.0

#define NUM_WIN_SUP_RODS	(NUM_WIRES+1)

/* mm */
#define WIN_SUP_ROD_HEIGHT		WINDOW_HEIGHT
#define WIN_SUP_ROD_DIAMETER	2.0

/* mm */
#define CAL_WIRE_DIA		1.0
#define CAL_WIRE_POS_X		0.0
#define CAL_WIRE_HEIGHT		(WIN_SUP_ROD_HEIGHT+WIN_SUP_ROD_DIAMETER)

#define MASK_PATH_FORMAT	"../config/MASK_PATTERN/mask_pattern%d.data"

/* unit = cnts/sqcm/sec */
#define CRAB_EQUI		3.75

#ifndef SUCCESS
#define SUCCESS			0
#endif

#ifndef ERROR
#define ERROR			1
#endif

#endif
