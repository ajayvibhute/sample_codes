/*
 * GetCorrelation.h
 *
 *  Created on: Apr 30, 2013
 *      Author: castro
 */

#ifndef GETCORRELATION_H_
#define GETCORRELATION_H_
#define NUMTHREADS 16
#define BUFSIZE 1024
#define ELEM 256
#define  NUMPIXPERQUAD 64
#define NUMPIXPERMODULE 16
#define ROWSPERMODULE 1950
#define COLSPERMODULE 1950
#define WINSIZE 3
#define NUMFITPARAM 3
#define FUNELEMENTS ((2*WINSIZE)+1)//*(1+(WINSIZE*2))
#define HEIGHT 477/(0.02)
#define THREASHOLD 3
int gensha(int moduleid,float*det,float xs,float ys,float zs,int *mask);
void calculateSigmaClippedMean(float* pixel_count,float *mean_out,float*rms_out);
void normalizeData(float *data,int numelements);
float getMaximum(float *data,int numelements);
int gensha_device(int moduleid,float*det,float xs,float ys,float zs,int*mask);
void readImage(char*filename,int *buffer,int hduNo);
void getminmax(int moduleno,int*x,int *y,float* x_center,float* y_center);
void do_svdfit(float * x,float *funx,float*a,int rows,int cols,float *chisq_out,void (*fpoly)(float,float*));
void getModule(int moduleNo,int *pixels);
int write_wcsaxis(fitsfile *imgfile, int axis, char *suffix,char *wcsname, char *wcstype, char *ctype, double crpix, double cdelt, double crval,char *cunit, int *status);
void writeFloatImage(char *filename,float*pixels,int rows,int cols,float resol);
void polyfunction(float x,float*funx);
void readFloatImage(char*filename,int hduno,float*data);
__global__ void CorrelationKernel(int *mask,float *dph,float *correlation,float *qemap,int numelem_x,int numelem_y,float x_start,float y_start,float resol,float zs,int moduleid);
__device__ void getminmax_device(int moduleno,int*x,int *y,float* x_center,float* y_center);
#endif /* GETCORRELATION_H_ */
