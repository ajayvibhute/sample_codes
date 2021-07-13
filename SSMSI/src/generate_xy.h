// This function is generates x and y coordinates and save it to the file named SVD_Input.txt
#include"transbin.c"
#include "dim.h"
#define nbins 240

/* This function will genrate a image for given thetaX and thetaY which will be used for fitting.
   This shadow will be generated only in simulation mode in actual case 
   this shadow will be given by base station
 */

void generateShadow(float thetaX1,float thetaY1, float intensity,int flag,float *y,int icam)
{
	
	int counter=1;
	float i=0;
	float xbin_det_lo=0,ybin_det_lo=0,thetax1,thetay1;
	float xbin_width=1,ybin_width=1,ageom;
	float deg2rad=0.01745329252;
	float ibin=0;
	int k=0;
	float first;
	transbininit();

	xbin_width = DET_WD_X/nbins;
	ybin_width = CELL_WD_Y;	
	
	thetax1=thetaX1*deg2rad;
	thetay1=thetaY1*deg2rad;
	float error;
	if(flag==0)
	{

		//for(k=0;k<3;k++)
		{
	
		for(i = 0; i <= NUM_WIRES-1; ++i)
		{
					
			ybin_det_lo = DET_ORIGIN_Y + i*ybin_width;
			for(ibin=0;ibin<=nbins;ibin++)
			{	
				xbin_det_lo = DET_ORIGIN_X+ibin*xbin_width;
		        	transbin(icam,xbin_det_lo,&xbin_width,ybin_det_lo,&ybin_width,thetax1,thetay1,&first,&ageom);
				y[counter]=first;
				y[counter]*=intensity;			
				
				error=((float)rand()/(float)RAND_MAX);
				error*=2;
				error-=1;
				error/=50;
				y[counter]+=error;
				counter++;
			}
					
		}
	
		}
	}
	else
	{
		
		for(k=0;k<3;k++)
		{
	
		for(i = 0; i <= NUM_WIRES-1; ++i)
		{
					
			ybin_det_lo = DET_ORIGIN_Y + i*ybin_width;
			for(ibin=0;ibin<=nbins;ibin++)
			{	
				xbin_det_lo = DET_ORIGIN_X+ibin*xbin_width;
		        	transbin(icam,xbin_det_lo,&xbin_width,ybin_det_lo,&ybin_width,thetax1,thetay1,&first,&ageom);
				first*=intensity;
				y[counter]+=first;

				error=((float)rand()/(float)RAND_MAX);
				error*=2;
				error-=1;
				error/=50;
				y[counter]+=error;
				counter++;
			}
					
		}
	
		}
	}					
}




void generateObservedShadow(float thetaX1,float thetaY1,float thetaX2,float thetaY2,float *y,int icam)
{
	
	int counter=0;
	float i=0;
	float xbin_det_lo=0,ybin_det_lo=0,thetax1,thetay1,thetax2,thetay2;
	float xbin_width=1,ybin_width=1,ageom;
	float deg2rad=0.01745329252;
	float ibin=0;
	int k=0;
	float first;
	transbininit();

	xbin_width = DET_WD_X/nbins;
	ybin_width = CELL_WD_Y;	
	
	thetax1=thetaX1*deg2rad;
	thetay1=thetaY1*deg2rad;

	thetax2=thetaX2*deg2rad;
	thetay2=thetaY2*deg2rad;
	for(k=0;k<3;k++)
	{
	
	for(i = 0; i <= NUM_WIRES-1; ++i)
	{
				
		ybin_det_lo = DET_ORIGIN_Y + i*ybin_width;
		for(ibin=0;ibin<=nbins;ibin++)
		{	
			xbin_det_lo = DET_ORIGIN_X+ibin*xbin_width;
		        transbin(icam,xbin_det_lo,&xbin_width,ybin_det_lo,&ybin_width,thetax1,thetay1,&first,&ageom);
			transbin(icam,xbin_det_lo,&xbin_width,ybin_det_lo,&ybin_width,thetax2,thetay2,&y[counter],&ageom);
			y[counter]+=first;
			//printf("%f \n",y[counter]);
			counter++;
		}
				
	}
	}					
}
void generateUpdatedShadow(float thetaX1,float thetaY1,float *y,int icam)
{
	
	int counter=0;
	float i=0;
	float xbin_det_lo=0,ybin_det_lo=0,thetax1,thetay1;
	float xbin_width=1,ybin_width=1,ageom;
	float deg2rad=0.01745329252;
	float ibin=0;
	int k=0;
	float first;
	transbininit();

	xbin_width = DET_WD_X/nbins;
	ybin_width = CELL_WD_Y;	
	
	thetax1=thetaX1*deg2rad;
	thetay1=thetaY1*deg2rad;

	for(k=0;k<3;k++)
	{
	
	for(i = 0; i <= NUM_WIRES-1; ++i)
	{
				
		ybin_det_lo = DET_ORIGIN_Y + i*ybin_width;
		for(ibin=0;ibin<=nbins;ibin++)
		{	
			xbin_det_lo = DET_ORIGIN_X+ibin*xbin_width;
		        transbin(icam,xbin_det_lo,&xbin_width,ybin_det_lo,&ybin_width,thetax1,thetay1,&first,&ageom);
			y[counter]=first;
			counter++;
		}
				
	}
	}					
}





