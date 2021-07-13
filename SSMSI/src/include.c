


__device__ void kern_qe(float mfp, float thetax,float absthetax,float thetay, float *bin, int nbins){

	float  d, cellx, tanx, tany, costheta, lmx, lmy, lmz, lm, Q,
		   mfpXcx, fac, flossx, prod;
	int i;

	d = DET_DEPTH;
	cellx = DET_WD_X;

	tanx = sin(thetax)/cos(thetax);
	tany = sin(thetay)/cos(thetay);

	costheta = 1/sqrt(1 + tanx*tanx + tany*tany);

	lmx = d * tanx;

	if(thetax == 0) lmy = d * tany;
	else lmy = lmx * tany/tanx;

	if(thetay == 0){
		if(thetax == 0) lmz = d;
		else lmz = lmx/tanx;
	}
	else lmz = lmy/tany;

	lm = lmz/costheta; /* total distance of the photon track in the cell */

	Q = 1 - exp(-lm/mfp); /* Quantum efficiency */

	mfpXcx = mfp * cos(thetax); /* MeanFreePath times Cos(thetax) */
	fac = d/mfpXcx; 	
	flossx = (mfp*sin(thetax)/cellx) * (1-exp(-fac)*(1+fac));

	prod = Q * (1-flossx);

	for(i = 0; i <= nbins-1; ++i)
		bin[i] *= (float)prod;
}
__device__ void cuda_reverse_f(float *arr, int size){

	float newarr[NUM_DET_BINS];
	int i;

	/* dont want to explore C-skill to play with indices! Playing safe! */
	for(i = 0; i <= size-1; ++i)
		newarr[i] = arr[size-1-i];

	for(i = 0; i <= size-1; ++i)
		arr[i] = newarr[i];

	return;
}

__device__ void kern_srb(float mfp, float thetax, float absthetax, float thetay, float *bin, int nbins){

		float startpos, endpos, d, pbw,
			   pos[NUM_DET_BINS], srbbin[NUM_DET_BINS],
			   shift, sx, mfpXsx, x0, tshft, xm, xl, bw,
			   nm, frac, dm, srbfrac;
		int i, j;

	absthetax = fabs(thetax);

	/* very small angles make very insignificant changes in the shadow */
	if(absthetax - 0.0001 < 0) return;

	startpos = -MASK_WD_X/2.0;
	endpos = MASK_WD_X/2.0;

	d = DET_DEPTH;

	pbw = MASK_WD_X/nbins;

	for(i = 0; i <= nbins-1; ++i){
		pos[i] = startpos + pbw* (i + 0.5); /* mid-bin position values */
		srbbin[i] = 0; /* initialising srbbin */
	}

	if(thetax < 0) cuda_reverse_f(bin, nbins); /* reverse bin for neg ang */

	sx = sin(absthetax); /* Not bothered about sign any more! */
	mfpXsx = mfp * sx;

	shift = d*sx/cosf(absthetax); /* shift due to tan-thetax */

	/*
	  In the loop below, we start from lower-limit (startpos) and proceed
	  towards upper-limit (endpos), at every step considering an anode
	  length of "shift" size between x0 and xm (or smaller than
	  "shift" if we are nearing the endpos!), applying the SRB over that
	  length
	*/
	for(i = 0; i <= nbins-1; ++i){

		x0 = pos[i];

		/* Length over which the SRB will be applied */
		tshft = x0 + shift;
		xm = tshft < endpos ? tshft : endpos;

		/* dm -- DenoMinator in the SRB expression */
		dm = 1 - exp(-(xm-x0)/mfpXsx);

		for(j = i; pos[j] <= xm && j <= nbins-1; ++j){
			xl = pos[j];

			bw = xm-xl;

			/* nm - NuMerator term in SRB expression */
			nm = 1 - exp(-(pbw < bw ? pbw : bw)/mfpXsx);

			frac = nm/dm;

			srbfrac = exp(-(xl-x0)/mfpXsx)*frac;

			srbbin[j] += srbfrac * bin[i]; /* j => xl and i => x0 */
		}
	}

	if(thetax < 0)
		cuda_reverse_f(srbbin, nbins); /* negative angle, reverse the result */

	for(i = 0; i <= nbins-1; ++i)
		bin[i] = srbbin[i]; /* overwrite the input binned data */

	return;
}

int cuda_shadowgen_srb_qe_lossx(double mfp, double thetax, double thetay,float *bins,int nbins)
{
	cudaError_t cuerr;
	dim3 dimblk = dim3(NUM_WIRES);
	/* launching the kernel*/
	if((cuerr = cudaGetLastError()) != cudaSuccess){
		fprintf(stderr, "Error: (pre srb launch) \"%s\"\n", cudaGetErrorString(cuerr));
		return ERROR;
	}

	kern_srb2 <<<1, dimblk>>> ((float)mfp, (float) thetax,0,(float) thetay, bins,nbins);

	if((cuerr = cudaGetLastError()) != cudaSuccess){
		fprintf(stderr, "Error: (post srb launch) \"%s\"\n", cudaGetErrorString(cuerr));
		exit(0);
		return ERROR;
	}
return SUCCESS;
}
