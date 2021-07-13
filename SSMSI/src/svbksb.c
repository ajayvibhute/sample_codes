/*
	This function is used for the backsubstituation and find the coefficients
*/
__device__ void svbksb(float *u, float *w, float *v, int m, int n, float *b, float *x)
{
	int jj,j,i;
	float s,tmp[1024];

	/*
		Back substitution of the matrices to find the un
	*/
	for (j=1;j<=n;j++) {
		s=0.0;
		if (w[j]) {
			for (i=1;i<=m;i++) s += u[i*n+j]*b[i];
			s /= w[j];
		}
		tmp[j]=s;
	}
	for (j=1;j<=n;j++) {
		s=0.0;
		for (jj=1;jj<=n;jj++) s += v[j*n+jj]*tmp[jj];
		x[j]=s;
	}
	
}

