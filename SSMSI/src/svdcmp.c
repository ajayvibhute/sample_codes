
/*
	Given a matrix a[1..m][1..n] this function computes it's singular value decomposition.
        A=U.W.V'. The matrix U replaces a on output. The diagonal matrix of singular values W is
        output as a vector w[1..n]. The matrix V is output as v[1..n][1..n]
	
*/

__device__ float pythag(float a, float b);
__device__ void svdcmp(float *a, int m, int n, float *w, float *v)
{

	int flag,i,its,j,jj,k,l,nm;
	float anorm,c,f,g,h,s,scale,x,y,z,rv1[1024];


	g=scale=anorm=0.0;
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		/*
			House holder reduction to bidiagonal form
		*/
		if (i <= m) { 
			for (k=i;k<=m;k++) scale += fabs(a[k*n+i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					a[k*n+i] /= scale;
					s += a[k*n+i]*a[k*n+i];
				}
				f=a[i*n+i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i*n+i]=f-g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=i;k<=m;k++) s += a[k*n+i]*a[k*n+j];
					f=s/h;
					for (k=i;k<=m;k++) a[k*n+j] += f*a[k*n+i];
				}
				for (k=i;k<=m;k++) a[k*n+i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(a[i*n+k]);
			if (scale) {
				for (k=l;k<=n;k++) {
					a[i*n+k] /= scale;
					s += a[i*n+k]*a[i*n+k];
				}
				f=a[i*n+l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i*n+l]=f-g;
				for (k=l;k<=n;k++) rv1[k]=a[i*n+k]/h;
				for (j=l;j<=m;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[j*n+k]*a[i*n+k];
					for (k=l;k<=n;k++) a[j*n+k] += s*rv1[k];
				}
				for (k=l;k<=n;k++) a[i*n+k] *= scale;
			}
		}
		anorm=FMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	/*
		Accumulation of right hand transformations	
	*/
	for (i=n;i>=1;i--) {
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++)/* double division  to avoid possible underflow */
					v[j*n+i]=(a[i*n+j]/a[i*n+l])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[i*n+k]*v[k*n+j];
					for (k=l;k<=n;k++) v[k*n+j] += s*v[k*n+i];
				}
			}
			for (j=l;j<=n;j++) v[i*n+j]=v[j*n+i]=0.0;
		}
		v[i*n+i]=1.0;
		g=rv1[i];
		l=i;
	}
	/*
		Accumulation of left hand transformations	
	*/
	for (i=IMIN(m,n);i>=1;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<=n;j++) a[i*n+j]=0.0;
		if (g) {
			g=1.0/g;
			for (j=l;j<=n;j++) {
				for (s=0.0,k=l;k<=m;k++) s += a[k*n+i]*a[k*n+j];
				f=(s/a[i*n+i])*g;
				for (k=i;k<=m;k++) a[k*n+j] += f*a[k*n+i];
			}
			for (j=i;j<=m;j++) a[j*n+i] *= g;
		} else for (j=i;j<=m;j++) a[j*n+i]=0.0;
		++a[i*n+i];
	}
	for (k=n;k>=1;k--) {
		/*
			Diagonalization of the bidiagonal form: Loop over singular values, and over a allowed iterations 
		*/
		for (its=1;its<=30;its++) {
			flag=1;
			/*
				Test for splitting  Note that rv1[1] is always 0
			*/
			for (l=k;l>=1;l--) {
				nm=l-1;
				if ((float)(fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((float)(fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag) { /* Cancellation of rv1[1], if rv1[1] >1*/
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((float)(fabs(f)+anorm) == anorm) break;
					g=w[i];
					h=pythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=1;j<=m;j++) {
						y=a[j*n+nm];
						z=a[j*n+i];
						a[j*n+nm]=y*c+z*s;
						a[j*n+i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {/* Convergence Singular value is made non negative */
				if (z < 0.0) {
					w[k] = -z;
					for (j=1;j<=n;j++) v[j*n+k] = -v[j*n+k];
				}
				break;
			}
			//if (its == 30) nrerror("no convergence in 30 svdcmp iterations");
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			/*
				Next QR transformation 
			*/
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj*n+j];
					z=v[jj*n+i];
					v[jj*n+j]=x*c+z*s;
					v[jj*n+i]=z*c-x*s;
				}
				z=pythag(f,h);
				/* Rotation can be arbitrary if z=0 */
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=1;jj<=m;jj++) {
					y=a[jj*n+j];
					z=a[jj*n+i];
					a[jj*n+j]=y*c+z*s;
					a[jj*n+i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	
}

