/*******************************************************************
			Mission: ASTROSAT
			Project: SSM

		  File-Name: Mask_Struct.h ($Revision: 1.2 $)
#%%			   Path: /astrosat2/ravibt/SSM/THREE_CAMERA/SIM/src
   Last Modified on: 2005 Jul 14 02:07:06 PM (IST)
			 Author: Ravi Shankar B.T. (ravibt@rri.res.in)
		 Supervisor: Dipankar Bhattacharya (dipankar@rri.res.in)
*******************************************************************/

/*

$Id: mask.h,v 1.2 2005/07/14 08:37:09 ravibt Exp $
$Log: mask.h,v $
Revision 1.2  2005/07/14 08:37:09  ravibt
Added the multiple-include-avoiding check.

Revision 1.1  2005/07/13 09:35:00  ravibt
Initial revision

*/

#ifndef __MASK_H__
#define __MASK_H__

/*
  THE FOLLOWING STRUCTURE DEFINITION (MaskPatternTypeY) IS FOR
  DEFINING THE EDGES OF THE SIX-MASK PATTERNS AND THE DIFFERENT
  MASK PLATE RIBS (SO IT IS DESCRIPTION OF THE MASK PLATE _ACROSS_
  THE MASK CODING DIRECTION).
*/
typedef struct{
	 double	*EdgeY;
	    int	*AreaY;
	    int *pat;
	    int	NumEdgesY;
} MaskPatternTypeY;

/*
  THE FOLLOWING STRUCTURE DEFINITION (MaskPatternTypeX) IS FOR
  DEFINING THE OPEN/CLOSE MASK-AREAS OF A MASK PATTERN (SO IT IS A
  DESCRIPTION OF THE MASK PLATE _ALONG_ THE MASK CODING DIRECTION).
*/
typedef struct{
	 double	*EdgeX;
	    int	*AreaX;
	    int	NumEdgesX;
} MaskPatternTypeX;

#endif
