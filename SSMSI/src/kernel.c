//This is the kernal which will execute the no of threads in  a block to calculate the svdfit
//Header 
#include"nrutil.h"
#include "svdfit.c"
#define min(a,b)	((a)<(b)?(a):(b))
#define max(a,b)	((a)>(b)?(a):(b))

__device__ void transbin_device(int icam, float xbin_det_lo, float *xbin_width, float ybin_det_lo, float *ybin_width, float thetax, float thetay, float *trans, float *ageom,int *maskpat){


	float ctx, stx, cty, sty, tx, ty,
		   maskthick_proj_x, el_extn_lo, el_extn_hi,
		   maskthick_proj_y, pat_extn_hi, pat_extn_lo,
		   xbin_det_cen, ybin_det_cen, hmasktop,
		   xshift_to_mask, yshift_to_mask, xbin_lo,
		   ybin_lo, xbin_hi, ybin_hi, xbin_mid, ybin_mid,
		   xref, yref, yrod_det, rod_extent,
		   yrod_lo[NUM_WIN_SUP_RODS], yrod_hi[NUM_WIN_SUP_RODS],
		   fpat1, ypat_lo, ylo_eff, ypat_hi, yhi_eff, totrodblock,
		   rodblock, fpat2, zcal_wire, xcal_wire,
		   calshadow_halfwidth, xcalshadow_cen,
		   xcal_lo, xcal_hi, xcal_start, xcal_end,
		   xcal_in_bin, xlo_eff, xhi_eff, xtotblock, xtotcommon,
		   xel_lo, xel_hi, blockstart, blockend, blocklen,
		   xcal_common, xfrac_block=0.0;

	int ipat1, ipat2, ielem_lo, ielem_hi, irod, iel_start, iel_end,
		iel=0;
	float mask_thickness[NUM_CAM] = {0},
			  elwidth[NUM_CAM] = {0},
			  patwidth[NUM_CAM] = {0},
			  patgap[NUM_CAM] = {0},
			  xbin_width_max[NUM_CAM] = {0},
			  ybin_width_max[NUM_CAM] = {0},
			  zmasktop_0[NUM_CAM] = {0},
			  zmasktopx_1[NUM_CAM] = {0},
			  zmasktopy_1[NUM_CAM] = {0},
			  xmin_mask_0[NUM_CAM] = {0},
			  xmin_mask_1[NUM_CAM] = {0},
			  ymin_mask_0[NUM_CAM] = {0},
			  ymin_mask_1[NUM_CAM] = {0},
			  ywinrod_0[NUM_CAM][NUM_WIN_SUP_RODS] = { {0} },
			  ywinrod_1[NUM_CAM][NUM_WIN_SUP_RODS] = { {0} },
			  winrod_dia[NUM_CAM][NUM_WIN_SUP_RODS] = { {0} },
			  zcalwire_0[NUM_CAM] = {0},
			  zcalwire_1[NUM_CAM] = {0},
			  xcalwire_0[NUM_CAM] = {0},
			  xcalwire_1[NUM_CAM] = {0},
			  calwire_dia[NUM_CAM] = {0};

	int npat[NUM_CAM] = {0}, nelem[NUM_CAM] = {0},nwinrod[NUM_CAM] = {0};

	
	/* Thickness of each mask plate */
	mask_thickness[BOOM]   = BO_MASK_THICK;
	mask_thickness[SLANT1] = SL_MASK_THICK;
	mask_thickness[SLANT2] = SL_MASK_THICK;

	/* Width of individual mask elements (x-direction) */
	elwidth[BOOM] = PMEW;
	elwidth[SLANT1] = PMEW;
	elwidth[SLANT2] = PMEW;

	/* Gap between successive patterns (closed) */
	patgap[BOOM] = INTER_PAT_GAP;
	patgap[SLANT1] = INTER_PAT_GAP;
	patgap[SLANT2] = INTER_PAT_GAP;

	/*
	  Width of each mask pattern (y-direction).  [The width in
	  x-direction is number of elements times element width]
	*/
	patwidth[BOOM] = BO_PAT_WD_Y + patgap[BOOM];
	patwidth[SLANT1] = SL_PAT_WD_Y + patgap[SLANT1];
	patwidth[SLANT2] = SL_PAT_WD_Y + patgap[SLANT2];

	zmasktop_0[BOOM] = BO_HEIGHT-WINDOW_HEIGHT;
	zmasktop_0[SLANT1] = SL_HEIGHT-WINDOW_HEIGHT;
	zmasktop_0[SLANT2] = SL_HEIGHT-WINDOW_HEIGHT;

	zmasktopx_1[BOOM] = MASK_TILT_X_BO;
	zmasktopx_1[SLANT1] = MASK_TILT_X_SL1;
	zmasktopx_1[SLANT2] = MASK_TILT_X_SL2;

	zmasktopy_1[BOOM] = MASK_TILT_Y_BO;
	zmasktopy_1[SLANT1] = MASK_TILT_Y_SL1;
	zmasktopy_1[SLANT2] = MASK_TILT_Y_SL2;

	/*
	  Left (lowest x) edge of the mask pattern.  The pattern used
	  for analysis has a virtual closed element appended to each end,
	  so the starting point is pushed back by one mask element width.

	  Tilt in the mounting of the mask pattern is accounted for by
	  using a linear coefficient: xmin_mask=xmin_mask_0+y*xmin_mask_1
	*/
	xmin_mask_0[BOOM] = PAT_LLIM_X_BO-elwidth[BOOM];
	xmin_mask_0[SLANT1] = PAT_LLIM_X_SL1-elwidth[SLANT1];
	xmin_mask_0[SLANT2] = PAT_LLIM_X_SL2-elwidth[SLANT2];

	xmin_mask_1[BOOM] = PAT_TILT_LLIM_X_BO;
	xmin_mask_1[SLANT1] = PAT_TILT_LLIM_X_SL1;
	xmin_mask_1[SLANT2] = PAT_TILT_LLIM_X_SL2;

	/*
	  Minimum value of y at which the mask pattern starts.  Linear
	  tilt is accommodated as: ymin_mask=ymin_mask_0+x*ymin_mask_1
	*/
	ymin_mask_0[BOOM] = MASK_Y_ORIGIN_BO;
	ymin_mask_0[SLANT1] = MASK_Y_ORIGIN_SL;
	ymin_mask_0[SLANT2] = MASK_Y_ORIGIN_SL;

	ymin_mask_1[BOOM] = PAT_TILT_LLIM_Y_BO;
	ymin_mask_1[SLANT1] = PAT_TILT_LLIM_Y_SL1;
	ymin_mask_1[SLANT2] = PAT_TILT_LLIM_Y_SL2;

	/*
	  True number of patterns.
	*/
	npat[BOOM] = npat[SLANT1] = npat[SLANT2] = NUM_PATTERNS;

	/*
	  True pattern length.  The maskpat array has one closed element
	  padded on either side
	*/
	nelem[BOOM] = nelem[SLANT1] = nelem[SLANT2] = NUM_MASK_ELEM;

	/*
	  At present, a common mask pattern is being considered for all
	  three SSM cameras.  Mask pattern is read in for icam=1 and is
	  assumed to be the same for the other two

	  one closed element is appended to either side of the mask pattern
	*/

	xbin_width_max[BOOM] = MAX_BIN_SIZE_FAC*elwidth[BOOM];
	xbin_width_max[SLANT1] = MAX_BIN_SIZE_FAC*elwidth[SLANT1];
	xbin_width_max[SLANT2] = MAX_BIN_SIZE_FAC*elwidth[SLANT2];

	ybin_width_max[BOOM] = MAX_BIN_SIZE_FAC*patwidth[BOOM];
	ybin_width_max[SLANT1] = MAX_BIN_SIZE_FAC*patwidth[SLANT1];
	ybin_width_max[SLANT2] = MAX_BIN_SIZE_FAC*patwidth[SLANT2];

	zcalwire_0[BOOM] = zcalwire_0[SLANT1] = zcalwire_0[SLANT2] =
	  CAL_WIRE_HEIGHT - WIN_SUP_ROD_HEIGHT + CAL_WIRE_DIA/2.0;

	zcalwire_1[BOOM] = zcalwire_1[SLANT1] = zcalwire_1[SLANT2] = 0.0;

	xcalwire_0[BOOM] = xcalwire_0[SLANT1] = xcalwire_0[SLANT2] = 0.0;

	xcalwire_1[BOOM] = xcalwire_1[SLANT1] = xcalwire_1[SLANT2] = 0.0;

	calwire_dia[BOOM] = calwire_dia[SLANT1] = calwire_dia[SLANT2] =
	  CAL_WIRE_DIA;

	int icam1, cam1, irod1;

	/*
	  Number of rods
	*/
	nwinrod[BOOM] = nwinrod[SLANT1] = nwinrod[SLANT2] = NUM_WIN_SUP_RODS;

	/*
	  Position (y-direction) of the (axis of) the rods, and the tilt
	  (linear coefficients of x): axis position = ywinrod_0+x*ywinrod_1

	  winrod_dia is the diameter of each rod
	*/
	for(icam1 = 1; icam1 <= NUM_CAM; ++icam1){ /* we have three cameras */

		switch(icam1){
			case 1: cam1 = BOOM;   break;
			case 2: cam1 = SLANT1; break;
			case 3: cam1 = SLANT2; break;
		}

		for(irod1 = 0; irod1 <= nwinrod[cam1]-1; ++irod1){
			ywinrod_0[cam1][irod1] = CELL_WD_Y*irod1 + DET_ORIGIN_Y;
			ywinrod_1[cam1][irod1] = 0.0;
			winrod_dia[cam1][irod1] = WIN_SUP_ROD_DIAMETER;
		}
	}

	/*#################################################################*/
	/*
	  First compute source-direction specific, bin-independent quantities
	*/	


	/*	
	  Projection factors in the direction of the source
	*/
	ctx=cos(thetax);
	stx=sin(thetax);
	cty=cos(thetay);
	sty=sin(thetay);
	tx=stx/ctx;
	ty=sty/cty;

	/*
	  thickness of the mask plate projected on the top surface of mask
	  (x-direction)
	*/
	maskthick_proj_x = mask_thickness[icam]*tx;

	if(maskthick_proj_x > elwidth[icam]){
		*trans = 0.0; /* theta_x outside valid range */
		/*printf("\nError: (%s:%d) ThetaX (%lf) outside valid range.\n", __FILE__, __LINE__, thetax);*/
		return;
	}

	el_extn_lo = min(maskthick_proj_x,0.0);
	el_extn_hi = max(maskthick_proj_x,0.0);

	/*
	  thickness of the mask plate projected on the top surface of mask
	  (y-direction)
	*/
	maskthick_proj_y = mask_thickness[icam]*ty;
	if (maskthick_proj_y > (patwidth[icam]-patgap[icam])){
		*trans = 0.0; /* theta_x outside valid range */
		/*printf("\nError: (%s:%d) ThetaY (%lf) outside valid range.\n", __FILE__, __LINE__, thetay);*/
		return;
	}

	pat_extn_lo = min(maskthick_proj_y,0.0);
	pat_extn_hi = max(maskthick_proj_y,0.0);

	/*#################################################################*/

	/*
	  All other quantities are bin-dependent. Compute them now.
	*/


	/*
	  Zero sized or invalid bins, return zero transmission
	*/
	if ((*xbin_width*(1-TOL) <= 0.0) || (*ybin_width*(1-TOL) <= 0.0)){
		*trans = 0.0; /* Invalid bin width */
		/*printf("\nError: (%s:%d) Either (or both) xbin_width (%lf) or ybin_width (%lf) is invalid.\n", __FILE__, __LINE__, *xbin_width, *ybin_width);*/
		return;
	}


	/*
	  Algorithm below works for binwidths less than two element widths
	  in x-direction and less than two pattern widths in the y-direction.
	  Restrict bin widths to these limits
	*/
	if(*xbin_width > xbin_width_max[icam])
		*xbin_width = xbin_width_max[icam]; /* adjusting *xbin_width */

	if (*ybin_width > ybin_width_max[icam])
		*ybin_width = ybin_width_max[icam]; /* adjusting *ybin_width */

	/* Bin centre on the detector plane */
	xbin_det_cen = xbin_det_lo + 0.5*(*xbin_width);
	ybin_det_cen = ybin_det_lo + 0.5*(*ybin_width);

	/* Geometric area of the bin  (sq. mm) */
	*ageom=*xbin_width*(*ybin_width);
	
	/*
	  Obtain the mask height at the location of the projected bin
	  allowing for linear tilt in both x and y
	*/
	hmasktop=(zmasktop_0[icam]+zmasktopx_1[icam]*xbin_det_cen
			 +zmasktopy_1[icam]*ybin_det_cen)
			 /(1.0-zmasktopx_1[icam]*tx-zmasktopy_1[icam]*ty);

	/*
	  Compute boundaries of the bin projected on the upper surface
	  of the mask plate, for parallel beam directed towards the source
	*/
	xshift_to_mask=hmasktop*tx;
	yshift_to_mask=hmasktop*ty;
	xbin_lo=xbin_det_lo+xshift_to_mask;
	ybin_lo=ybin_det_lo+yshift_to_mask;
	xbin_hi=xbin_lo+*xbin_width;
	ybin_hi=ybin_lo+*ybin_width;

	/*
	  Centre location of the projected bin
	*/
	xbin_mid=xbin_det_cen+xshift_to_mask;
	ybin_mid=ybin_det_cen+yshift_to_mask;

	/*
	  mask plate starting reference, taking into account linear x-y tilt
	*/
	xref = xmin_mask_0[icam] + xmin_mask_1[icam]*ybin_mid;
	yref = ymin_mask_0[icam] + ymin_mask_1[icam]*xbin_mid;

	/*
	  mask pattern numbers [1..npat(icam)] on which the projected bin falls
	  (in the present implementation, could at most be two patterns).
	*/
	ipat1 = (int) ((ybin_lo - yref)/patwidth[icam]);
	ipat2 = (int) ((ybin_hi - yref)/patwidth[icam]);

	/*
	  mask element numbers [2..nelem(icam)+1] within the x-bin
	  (element no. 1 and nelem(icam)+2 are outside the mask pattern
	  and are forced to be zero (closed). They are needed to take into
	  account the projection of mask plate thickness at the edges
	  of the mask pattern)
	*/
	ielem_lo = (int) ((xbin_lo-xref)/elwidth[icam]);
	ielem_hi = (int) ((xbin_hi-xref)/elwidth[icam]);

	/*
	  Projected bin falls outside the mask plate
	*/
	if((ipat2 < 0) || (ipat1 >= npat[icam])){
		*trans = 0.0; /* bin outside mask plate */
		return;
	}

	if ((ielem_hi <= 0) || (ielem_lo >= nelem[icam]+1)){ /* one elem more at end */
		*trans = 0.0; /* bin outside mask plate */
		return;
	}

	/*
	  projection of window support rods on the top surface of mask plate
	*/
	for(irod = 0; irod <= nwinrod[icam]-1; ++irod){

		/*
		  first find the y-location of the window support rod, allowing
		  for mounting tilts (linear gradient in x-direction)
		*/
		yrod_det = ywinrod_0[icam][irod]+ywinrod_1[icam][irod]*xbin_det_cen;

		/*
		  now project it on the mask: yrod_lo and yrod_hi are the two ends
		*/
		rod_extent=winrod_dia[icam][irod]/cty;
		yrod_lo[irod] = yrod_det + yshift_to_mask
						- 0.5*winrod_dia[icam][irod]*ty
						- 0.5*rod_extent;
		yrod_hi[irod]=yrod_lo[irod]+rod_extent;
	}

	/*
	  fractional contribution of each pattern
	*/
	if (ipat2 > ipat1){
		if(ipat1 < 0) fpat1 = 0.0;
		else{
			/*
			  check for blockage due to inter-pattern gap and
			  projection of mask plate thickness within the
			  portion of the bin that falls on the first pattern
			  covered.
			*/
			ypat_lo = yref + ipat1*patwidth[icam];
			ylo_eff = max((ypat_lo+pat_extn_hi), ybin_lo);
			ypat_hi = yref + (ipat1+1)*patwidth[icam];
			yhi_eff = ypat_hi-patgap[icam] + pat_extn_lo;
			fpat1 = (yhi_eff-ylo_eff)/(*ybin_width);
			if(fpat1 < 0.0) fpat1 = 0.0;
			else{
				/* check for further blockage by window support rods */
				totrodblock = 0.0;

				for(irod = 0; irod <= nwinrod[icam]-1; ++irod){
					rodblock = min(yhi_eff,yrod_hi[irod])
							   -max(ylo_eff,yrod_lo[irod]);

					if(rodblock > 0.0)
						totrodblock=totrodblock+rodblock;
				}
				fpat1=fpat1-totrodblock/(*ybin_width);

			}
		}

		if(ipat2 >= npat[icam]) fpat2 = 0.0;
		else{

			/*
			  check for blockage by pattern gap/mask thickness
			  in the second pattern covered
			*/
			ypat_lo = yref + ipat2*patwidth[icam];
			ylo_eff = max((ypat_lo+pat_extn_hi),ybin_lo);
			ypat_hi = yref + (ipat2+1)*patwidth[icam];
			yhi_eff = min((ypat_hi-patgap[icam]+pat_extn_lo),ybin_hi);
			fpat2 = (yhi_eff-ylo_eff)/(*ybin_width);
			if(fpat2 < 0.0) fpat2 = 0.0;
			else{

				/*
				  check for blockage by window support rods
				*/
				totrodblock = 0.0;
				for(irod = 0; irod <= nwinrod[icam]-1; ++irod){
					rodblock = min(yhi_eff,yrod_hi[irod])
							   -max(ylo_eff,yrod_lo[irod]);
					if (rodblock > 0.0)
						totrodblock=totrodblock+rodblock;
				}

				fpat2=fpat2-totrodblock/(*ybin_width);
			}
		}
	}
	else{
		/*
		  Bin wholly on one pattern. Set pattern 2 transmission to zero
		*/
		fpat2 = 0.0;

		/*
		  Find blockage by pattern gap and projected mask plate width
		*/
		ypat_lo = yref + ipat1*patwidth[icam];
		ylo_eff = max((ypat_lo+pat_extn_hi),ybin_lo);
		ypat_hi = yref + (ipat1+1)*patwidth[icam];
		yhi_eff = min((ypat_hi-patgap[icam]+pat_extn_lo),ybin_hi);
		fpat1 = (yhi_eff-ylo_eff)/(*ybin_width);
		if(fpat1 < 0.0) fpat1=0.0;
		else{

			/* check for blockage by window support rods */
			totrodblock = 0.0;
			for(irod = 0; irod <= nwinrod[icam]-1; ++irod){
				rodblock = min(yhi_eff,yrod_hi[irod])
						   -max(ylo_eff,yrod_lo[irod]);
				if (rodblock > 0.0)
					totrodblock = totrodblock+rodblock;
			}
			fpat1 = fpat1-totrodblock/(*ybin_width);
		}
	}

	/*
	  projection of calibration wire on mask plate top, allowing for
	  linear mounting tilt

	  for a given ybin_det_cen, first find the height from the window
	  plane and the x-position of the axis of the calibration wire
	*/
	zcal_wire = (zcalwire_0[icam]+zcalwire_1[icam]*ybin_det_cen)
				/(1.0-zcalwire_1[icam]*ty);
	xcal_wire = xcalwire_0[icam]+xcalwire_1[icam]*ybin_det_cen;

	/* then project it on the mask plane top surface */
	calshadow_halfwidth = 0.5*calwire_dia[icam]/ctx;
	xcalshadow_cen  = xcal_wire+(hmasktop-zcal_wire)*tx;

	xcal_lo = xcalshadow_cen - calshadow_halfwidth;
	xcal_hi = xcalshadow_cen + calshadow_halfwidth;

	/*
	  find the width of that part of the shadow of the calibration wire
	  which falls within the bin
	*/
	xcal_start = max(xcal_lo,xbin_lo);
	xcal_end   = min(xcal_hi,xbin_hi);
	xcal_in_bin = xcal_end-xcal_start;
	if (xcal_in_bin < 0.0) xcal_in_bin=0.0;

	/*
	  redefine effective element area to look for blockage within the
	  bin: due to projected thickness of mask plate adding to the
	  shadow, closed elements from adjacent mask elements may cast
	  shadow within the bin range.
	*/
	xlo_eff = xbin_lo - el_extn_hi;
	xhi_eff = xbin_hi - el_extn_lo;

	iel_start = (int) ((xlo_eff-xref)/elwidth[icam]);
	iel_end = (int ) ((xhi_eff-xref)/elwidth[icam]);

	/*
	  compute transmission fraction by estimating the total
	  blockage in the x direction
	*/
	*trans=0.0;
	if (fpat1 > 0.0){
		xtotblock=0.0;
		xtotcommon=0.0;
		for(iel = iel_start; iel <= iel_end; ++iel){
			/*
			  for each closed element within search range, find the
			  length of the element, after addition of mask thickness
			  projection, that falls within the bin.  Also find if any
			  part of the shadow of the calibration wire is in common
			  with the shadow of these closed mask elements.
			*/
			if(maskpat[ipat1*(NUM_MASK_ELEM+2)+iel] == 0){
				xel_lo = xref+iel*elwidth[icam];

				/*
				  if the adjacent element at lower x is open, then
				  extend the area blocked by this element by the
				  projection of mask thickness in the low-x direction

				  Note: el_extn_lo is negative or zero.

				  It is assumed that the projection of mask thickness
				  is not more than one element width.
				*/
				if (iel > 0 && maskpat[ipat1*(NUM_MASK_ELEM+2)+(iel-1)] == 1)
					xel_lo = xel_lo+el_extn_lo;

				xel_hi = xref+(iel+1)*elwidth[icam];

				/*
				  if the adjacent element at higher x is open, then
				  extend the area blocked by this element by the
				  projection of mask thickness in the high-x direction

				  Note: el_extn_hi is zero or positive.

				  It is assumed that the projection of mask thickness is
				  not more than one element width.
				*/
				if (iel <= nelem[icam] && maskpat[ipat1*(NUM_MASK_ELEM+2)+(iel+1)] == 1)
					xel_hi=xel_hi+el_extn_hi;

				/* Now estimate the blocked area within the bin */
				blockstart = max(xel_lo,xbin_lo);
				blockend   = min(xel_hi,xbin_hi);
				blocklen   = blockend-blockstart;

				if(blocklen < 0.0) blocklen=0.0;
				xtotblock = xtotblock+blocklen;

				if (xcal_in_bin > 0.0){

					/*
					  if any part of the shadow of the calibration wire
					  falls within the bin, then find how much of this
					  shadow may be in common with the shadow of closed
					  elements falling within the bin.  The remaining
					  part of the calibration wire shadow is added to
					  the total blockage.
					*/

					xcal_common = min(blockend,xcal_end)
								  - max(blockstart,xcal_start);
					if(xcal_common < 0.0) xcal_common=0.0;
					xtotcommon = xtotcommon+xcal_common;
				}
			}
		}

		/*
		  From total blockage, compute the fraction of x bin width
		  that is blocked.  The complement of this fraction, multiplied
		  by the fractional contribution of the pattern, is the
		  transmission through this part.
		*/
		xtotblock=xtotblock+xcal_in_bin-xtotcommon;
		xfrac_block=xtotblock/(*xbin_width);
		if (xfrac_block > 1.0) xfrac_block=1.0;
		if (xfrac_block < 0.0) xfrac_block=0.0;
		*trans = *trans + fpat1*(1.0-xfrac_block);
	}

	/*
	  repeat the above for the second mask pattern falling within
	  the bin, and add to the total transmission.
	*/
	if (fpat2 > 0.0){
		xtotblock = 0.0;
		xtotcommon = 0.0;

		for(iel = iel_start; iel <= iel_end; ++iel){
			/* blockage by closed elements+projected mask thickness */
			if (maskpat[ipat2*(NUM_MASK_ELEM+2)+iel] == 0){

				xel_lo = xref+iel*elwidth[icam];

				if(iel > 0 && maskpat[ipat2*(NUM_MASK_ELEM+2)+(iel-1)] == 1)
					xel_lo = xel_lo + el_extn_lo;

				xel_hi = xref+(iel+1)*elwidth[icam];

				if (iel <= nelem[icam] && maskpat[ipat2*(NUM_MASK_ELEM+2)+(iel+1)] == 1)
					xel_hi = xel_hi + el_extn_hi;

				blockstart = max(xel_lo,xbin_lo);
				blockend   = min(xel_hi,xbin_hi);
				blocklen   = blockend-blockstart;
				if(blocklen < 0.0) blocklen = 0.0;
				xtotblock = xtotblock + blocklen;
				if (xcal_in_bin > 0.0){
					/*
					  part of calibration wire shadow falling on
					  shadow of closed mask elements.
					*/
					xcal_common = min(blockend,xcal_end)
								  - max(blockstart,xcal_start);
					if (xcal_common < 0.0) xcal_common=0.0;
					xtotcommon = xtotcommon+xcal_common;
				}
			}
		}

		/* total blocked fraction in this part */
		xtotblock=xtotblock+xcal_in_bin-xtotcommon;
		xfrac_block=xtotblock/(*xbin_width);
		if (xfrac_block > 1.0) xfrac_block=1.0;
		if (xfrac_block < 0.0) xfrac_block=0.0;

		/* transmission through this part, added to total transmission */
		*trans = *trans + fpat2*(1.0-xfrac_block);
	}
	
	return;
}

__device__ void generateShadow(int icam,int *maskpat,float *bins,float thetaX,float thetaY)
{
	int counter=1;
	float i=0.0;
	float xbin_det_lo=0,ybin_det_lo=0,thetax,thetay;
	float xbin_width=1,ybin_width=1,ageom;
	float deg2rad=0.01745329252;
	float ibin=0;
	int noOfBins=nbins;
	xbin_width = DET_WD_X/noOfBins;
	ybin_width = CELL_WD_Y;	
	//converting deg to rad	
	thetax=thetaX*deg2rad;
	thetay=thetaY*deg2rad;
	for(i=0;i<=NUM_WIRES-1;i++)
	{
		ybin_det_lo = DET_ORIGIN_Y + i*ybin_width;
		for(ibin=0;ibin<=noOfBins;ibin++)
		{	
			xbin_det_lo = DET_ORIGIN_X+ibin*xbin_width;
			transbin_device(icam,xbin_det_lo,&xbin_width,ybin_det_lo,&ybin_width,thetax,thetay,&bins[counter],&ageom,maskpat);
			counter++;
		}
	}
				
				
}

__global__ void Kernel(float *x,float *y,float *thetaX,float *thetaY,float *sig,float *a,float *chisq,float *u,float *b,float *afunc,int no_another,float *x_another,int *maskpat,int noOfGrids,int icam)
{
	/* 
	calculate the index of the current matrix which will be passed for execution
	*/
	int index=blockIdx.x*blockDim.x+threadIdx.x;
	//bound check on index
	int noOfBins=rows;
	//for(int k=0;k<noOfGrids;k++)
	{
		for(int i=1;i<=noOfBins;i++)
		{
			sig[i]=1;
			
		}
	}
	if(index<noOfGrids)
	{
		
		/*
		calculate the bin values for given location thetaX and thetaY
		*/
		generateShadow(icam,maskpat,&y[index*noOfBins],thetaX[index],thetaY[index]);
		/*
		 fitting of the generated shadow with the observed shadow and return chisquare and fitted       parameters 
		*/
		svdfit(x,&y[index*noOfBins],sig,noOfBins,&a[index*colms],colms,&u[index*noOfBins*colms],&chisq[index],&b[index*noOfBins],&afunc[index*colms],no_another,&x_another[index*no_another*noOfBins]);
		
	}
	
}

////svdfit(&x[index*noOfBins],&y[index*noOfBins],&sig[index*noOfBins],noOfBins,&a[index*colms],colms,&u[index*noOfBins*colms],&v[index*colms*colms],&w[index*colms],&chisq[index],&cvm[index*colms*colms],&b[index*noOfBins],&afunc[index*noOfBins],no_another,&x_another[index*no_another*noOfBins]);
