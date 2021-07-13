// This function is used to find the covariance matix 

__device__ void svdvar(float *v, int ma, float w[], float *cvm)
{
	int k,j,i;
	float sum,wti[1024];

	
	for (i=1;i<=ma;i++) {
		wti[i]=0.0;
		if (w[i]) wti[i]=1.0/(w[i]*w[i]);
	}
	for (i=1;i<=ma;i++) {
		for (j=1;j<=i;j++) {
			for (sum=0.0,k=1;k<=ma;k++) sum += v[i*ma+k]*v[j*ma+k]*wti[k];
			cvm[j*ma+i]=cvm[i*ma+j]=sum;
		}
	}
	
}

