
/*
	Computes (a*a+b*b)^0.5 withous destructive underflow or over flow 
*/
__device__ float pythag(float a,float b)
{
	float absa=0.0,absb=0.0;
	absa=fabs(a);
	absb=fabs(b);
	if(absa>absb)
	{
		return absa*sqrtf(1.0+SQR(absb/absa));
	}
	else
	{
		return (absb == 0.0 ? 0.0 : absb*sqrtf(1.0+SQR(absa/absb)));
	}
}
