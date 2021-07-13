//This function will return minimum no between given two numbers
__device__ float IMIN(float a,float b)
{
		if(a<b)
			return a;
		else
			return b;
}
//This function will return the larger no 
__device__ float FMAX(float a,float b)
{
		if(a>b)
			return a;
		else 
			return b;
}
// this will return the absolute value of x
__device__ float SIGN(float a,float b)
{
		if(b>=0.0)
			return fabs(a);
		else
			return (-fabs(a));
}
// This will return the square of given no.
__device__ float SQR(float a)
{
	return a*a;
}

__device__ void nrerror(char error_text[])
{
	
}
