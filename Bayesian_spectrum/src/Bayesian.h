#define BUFFSIZE 1024
#define NUMQUADRANT 1
#define NUMELEM 4096
#define ELEMPERQUAD 64
#define PI 3.14159265
#define TORAD 0.0174532925

//#define LLDFILE "/Users/ajay/CALDB/data/as1/czti/bcf/AS1cztlld20160517v01.fits"
int getSourceList(char*,float*,float*,int*);
void printerror( int status);
void createChannelDPH(char*infile,float **count,int minpi,int maxpi,int binsize,int*LLD,int*pix_flag);
void generateMatrix(float*tx,float*ty,int numsources,float *x,float*sig,float **u, float**shadow,int rowcount);
int write_spectrum_file(char*specFilename, char*specTemplate,char* infile,int*pi,float* flux,float* error,int nbins,double tx,double ty); 
int create_empty_fitsfile(char* outFilename, char* outTemplate);
void copyheader(fitsfile *fptr,fitsfile *outfile);
void fitPowerLaw(float *energy,float *counts,int numbins,float *params);
void writeFitsDPH(char * outputfilename,int qid,float* count);
void getShadowSingleQuad(float *shadow_pixels,float tx,float ty);
void executeLikelihood(char*infile,char*outfile,float **params,float **shadow,int numsources,int minpi,int maxpi,float**,int*LLD,int*pix_flag,float**error);
void ReadShadow(char*filename,float *mask,int hduNo);
int read_lldfile(char* caldb_lld,int qid,int *lld);
void read_badpixfile(char*badpixfilename,int quadid, int*pix_flag);





/*void createDPH(char*,float **,float,float,int);
void writeFitsDPH(char * outputfilename,int qid,float* count);
float getmasktrans(float tx,float ty,float ekev);
double cztmasktrans(double ekev, double thetadeg);
void fun(float*fx,float minenergy,float binwidth,int numbins,float p);
void addIntensityRatio(char*infile,float *x,float **itensityratio,int numbis);
void polint(float *xa, float *ya, int n, float x, float *y, float *dy);
void executeLikelihood(char*infile,char*outfile,float **params,float **shadow,int numsources,float minenergy,float maxenergy);
void singlefit(float *fx,float*c,int n,float*a,float*chisq);
fitsfile *write_fptr;
void generateFitMatrix(float *x,float*sig,float **u,int cols,int rows );
void generateBackgroundMatrix(float*tx,float*ty,int numsources,float *x,float*sig,float **u );
*/
