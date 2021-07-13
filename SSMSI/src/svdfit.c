//This code is borrowed from Numerical Recipes and modified to work on GPUs

//This function calculate the coeffiecient by using another two functions svdcmp and svbksb.
//This function also calculate the chisquare value for given matrix u
//Header files
#include "pythag.c"
#include "svdcmp.c"
#include "svbksb.c"
#include "svdvar.c"
#include"fpoly.c"
#define tol 1.0e-5
/*
  ndata=No of Bins
  ma=No of Colms
*/
/*
	Given a set of data points x[1..ndata], y[1..ndata] with individual standard deviations sig[1..ndata]
	use  χ2 minimization to determine the coefficient a[1..ma] of the fitting function
	Here we solve the fitting equation using singular value decomposition of the ndata by ma matrix
*/
//declaring the device functions
__device__ void svbksb(float *u, float *w, float *v, int m, int n, float *b,float x[]);
__device__ void svdcmp(float *a, int m, int n, float *w, float *v);
__device__ void svdfit(float *x, float *y, float *sig,const int ndata, float *a,const int ma,float *u,float *chisq,float *b,float *afunc,int no_another,float *x_another)
{
	float *w=NULL,*v=NULL;

	w=(float*)malloc(sizeof(float)*ma);
	v=(float*)malloc(sizeof(float)*(ma+1)*(ma+1));	
	float wmax=0,tmp=0,thresh=0,sum=0;
	int j=0,i=0;	

	//Generate the matrix u using the fpoly function which generates the functions of x.
	for (i=1;i<=ndata;i++) 
	{
		/* Accumulate coefficient of the fitting matrix	*/
		if(no_another==0)
		{	
			fpoly(y[i],afunc,ma);		
		}
		else
		{
			fpoly(y[i],afunc,ma,no_another,&x_another[(i-1)*no_another]);
		}
		tmp=1.0/sig[i];
		for (j=1;j<=ma;j++) u[i*ma+j]=afunc[j]*tmp;
		b[i]=x[i]*tmp;
	}
	/*svdcmp decomposes the matrix u into three matrices using Singular Value Decomposition.*/
	svdcmp(&u[0],ndata,ma,w,v);
	//Set the threadshold values
	wmax=0.0;
	for (j=1;j<=ma;j++)
		if (w[j] > wmax) wmax=w[j];
	thresh=tol*wmax;
	for (j=1;j<=ma;j++)
		if (w[j] < thresh) w[j]=0.0;

	//svd back substitution for calculating coefficients of matrix(generated from given set of data)
	svbksb(&u[0],w,v,ndata,ma,b,a);

	//Calculate the chisquare value of the given matrix(generated from given set of data)
	*chisq=0.0;
	for (i=1;i<=ndata;i++) {
		if(no_another==0)
		{
			fpoly(y[i],afunc,ma);		
		}
		else
		{
			fpoly(x[i],afunc,ma,no_another,&x_another[(i-1)*no_another]);
		}
		for (sum=0.0,j=1;j<=ma;j++) sum += a[j]*afunc[j];
		*chisq += (tmp=(x[i]-sum)/sig[i],tmp*tmp);
	}
	//calculate the coverience matrix for the given matrix	
	//svdvar(v,ma,w,cvm);
	
	free(w);
	free(v);
}
