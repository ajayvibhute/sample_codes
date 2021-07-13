/*
 * statSubDefSer.h
 *
 *  Created on: Sep 16, 2013
 *  Author: Ajay Vibhute
 */

#ifndef STATSUBDEFSER_H_
#define STATSUBDEFSER_H_
/*#include "statSubDefSer.c"*/
#define DATASIZE 1000
#define BUFFSIZE 1024
struct lightcurveParam
{
	int numpoints;
	char id[256];
	//double mag[DATASIZE],mag_error[DATASIZE],mjd[DATASIZE],temp[DATASIZE],*buff;
	//double *mag,*mag_error,*mjd,*temp;
	double min_mag,max_mag,amplitude,fpr_mid20,fpr_mid35,fpr_mid50,fpr_mid65,fpr_mid80;
	double xm,mjdm,xsig,abbe,lsq[2],xsk,xkur,cusum,xmad,s,xbwmv,thiel_sen,xdw,xkt;
	double rchi,stetson_k,stetson_j,con,xvon,xdfa ;
};
struct lightcurve
{
	int numpoints;
	double *mag,*mag_error,*mjd,*temp;
};
int readCurve(char*infile,struct lightcurve *curves);
int parseFile(char*infile);
double getabs(double data);
int mean(double*data,double*mean,int n);
int sd(double*data,double*sd,int n,double mean);
int skew(double*x,double*xsk,int n,double xm);
int kurtosis(double*x,double*xkur,int n,double xm,double xs);
int median(double x[],double *xmed,int n);
int median(double *x,double *xmed,int n);
int  sort(double *x,int n);
int mad(double*x,double*xmad,int n,double*xtemp);
int sn(double*x,double *sn,int n);
int bwmv(double*x,double *xbwmv,int n,double*tmp);
int thiel_sen(double *x,double*y, double *xthiel_sen,double *med, int n);
int durbin_watson(double*x,double*xdw,int n,double*e,double xm);
int kendall_tau(double *x,double*y, double *xkt, int n);
int cusum(double*x,double*xcusum,int n,double xm,double xs);
int dfa(double*x,double*t,double *xdfa,double *xt,int n,double xm);
int von_neumann(double*x,double*xvon,int n,double xs);
int stetson_j(double*mjd,double*data,double*dataerr,double *stetson_j,int n,double datam);/*Arun*/
int stetson_k(double*data,double*dataerr,double *stetson_k,int n,double datam);/*Arun*/
int abbe(double*data, double *abbe,int n,double datam) ;/*Arun*/
int reducedchi(double*data,double *dataerr,double*rchi,int n,double datamean);/*Arun*/
double sum(double*a, int n);/*Arun*/
int lsq(double*data,double*mjd,double*lsq,int n,double datam,double mjdm);/*Arun*/
int subt(double *inarray,double*dataout,double subtr,int n);/*Arun*/
int con(double*data,double *con, int n,double datamean,double xsig);
int Rcode(double *mag,double*mag_err,double*timeo,int n,double*times,double*delta);
void processLightCurve(struct lightcurve *curves,struct lightcurveParam*param);
char *getFileName(char* mystr);
char *removeFileExt(char* mystr);
#endif /* STATSUBDEFSER_H_ */
