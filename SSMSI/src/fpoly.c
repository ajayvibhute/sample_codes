/*
	This function genrate the functions of x for the given x value for the generation of matrix.
        The generated matrix will be used for the svd decomposition. 
     
*/

__device__ void fpoly(float x, float *p, int np)
{
	int j;

	p[1]=1.0;
	for (j=2;j<=np;j++) p[j]=p[j-1]*x;
}
/*
	 This function generate np-no_other functions of x and merge it with no_another functions of x
*/
__device__ void fpoly(float x, float *p, int np,int no_another,float *x_another)
{
	
	np-=no_another;
	p[1]=1.0;
	for (int j=2;j<=np;j++) p[j]=p[j-1]*x;
//extending the another function of x with the generated function of x
	for(int j=np+1,i=1;i<=no_another;j++,i++)
	{
		p[j]=x_another[i];
	}

	
}
