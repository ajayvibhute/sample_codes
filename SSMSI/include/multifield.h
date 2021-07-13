
/*******************************************************************
            Mission: ASTROSAT
            Project: SSM

          File-Name: multifield.h ($Revision: 1.2 $)
#%%            Path: /data1/ravibt/ASTROSAT/SSM/multifield/src
         Created on: 2008 Sep 26 06:06:11 PM (IST: UTC+0530)
   Last Modified on: 2009 Feb 04 03:14:45 PM (IST: UTC+0530)
             Author: Ravishankar B.T. (ravibt@isac.gov.in)
        Supervisors: Dipankar Bhattacharya (dipankar@iucaa.ernet.in)
                     Seetha S (seetha@isac.gov.in)
*******************************************************************/




/*####################################################################

$Id: multifield.h,v 1.2 2009/04/06 12:33:04 ravibt Exp $
$Log: multifield.h,v $
Revision 1.2  2009/04/06 12:33:04  ravibt
Some additional definitions for writeshadow()

Revision 1.1  2008/10/27 10:05:46  ravibt
Initial revision

####################################################################*/




#ifndef ___MULTI_FIELD_H___
#define ___MULTI_FIELD_H___


#include "dim.h"
#include "pthread.h"
#include <stdio.h>		/* printf(3), sprintf(3), fopen(3), fclose(3)
						   fscanf(3)								*/
#include <stdlib.h>		/* exit(3), NULL, system(3), malloc(3),
						   free(3)									*/
#include <string.h>		/* strcmp(3)								*/
#include <unistd.h>		/* sleep(3) 								*/
#include <time.h>
#include <sys/time.h>
#include <stdio.h>		
#include <stdlib.h>
#include <math.h>	
#include <dirent.h>	
#include <sys/stat.h>	
#include <string.h>	

/* Misc */
#define PARSE_SUCCESS	0
#define PARSE_ERROR		1
#define SUCCESS			0
#define ERROR			1

/* unit = degree */
#define DEF_BOOM_RA		0.0
#define DEF_BOOM_DEC	0.0

/* unit = degree */
/*#define INCR_RA			0.1
#define INCR_DEC		0.1*/
#define INCR_RA			0.5
#define INCR_DEC		0.5

/*
  0.05deg=3arcmin;
  An increment of 0.05deg in both RA and dec =>
    nelem = (360/0.05) * (180/0.05) = 25920000
	=> 25920000* 8 = 207360000 bytes = 198MB!!
	For five different arrays = > 1GB!!!!

  May be 0.1deg (6arcmin) will be sufficient.
*/

#define NUM_RA			((int) (360.0/INCR_RA)+1)
#define NUM_DEC			((int) (180.0/INCR_DEC)+1)

/* unit = degree */
#define INCR_INCLI		60.0


#define RA_RANGE		261.0
#define DEC_RANGE		120.0

/* unit = degree;    0.4 degree = 24 arcmin */
#define INCR_RA_COARSE	0.4
#define INCR_DEC_COARSE	0.4


/*
  unit = degree; If dec of boom-pt or any of corners goes beyond this
  DEC_THRESH, the search-area will be extended to the complete sky
*/
#define DEC_THRESH		70.0


#define PATH_SHA_DIR	"./SHADOWS"
#define PATH_MFP		"../config/avg_mfp.dat"


/* used in writeshadow() */
#define GENSHA			0
#define INS_FLD_SEP		1
#define INS_REC_SEP		2
#define WRITE_ZEROES	3
#define GENSHA_CMD		"GenShadow.sh"
#define GENSHA_EXPTIME	"12500"

/* %s below will be replaced by one of the three MFP values */
#define NORM_SHA_NAME_FORMAT	"%s_Norm.sha"
#define NONORM_SHA_NAME_FORMAT	"%s_NoNorm.sha"

/*
 20160 = 63 * 4 * 8 * 10
 63 --> Number of bins per detector
  4 --> if we sample the detector bins four times smaller than now.
  8 --> number of wires
 10 --> Number of bytes per detector bin count (ascii representation)
*/
#define MAXLINELENGTH	20160


#define MFP_COLSIZE		50
#define KEYWORD_COLSIZE	100
#define PATH_COLSIZE	200


#define MFP_MIN			7.0
#define MFP_MAX			20.0
#define NUM_MFP			3

/*
  Number of results to be written: for each mfp an unnormalised and
  a normalised version of the shadow.
*/
#define NUM_RES_FILES	(NUM_MFP*2)


/* unit = mm */
#define DEF_SIGMA_X		0.8

/* Number of bins per anode */
#define DEF_NUM_BINS	63

/* unit = s */
#define DEF_EXP_TIME	1.0



#define MAX_THREAD_SETS      20
#define MAX_THREADS_PER_SET  NUM_MFP
#define MAX_CARDS_SUPPORTED	 2


#define IMAGEDATA_T float 

typedef struct {
float mask_thickness[NUM_CAM];
float elwidth[NUM_CAM];
float patwidth[NUM_CAM];
float patgap[NUM_CAM];
float xbin_width_max[NUM_CAM];
float ybin_width_max[NUM_CAM];
float zmasktop_0[NUM_CAM];
float zmasktopx_1[NUM_CAM];
float zmasktopy_1[NUM_CAM];
float xmin_mask_0[NUM_CAM];
float xmin_mask_1[NUM_CAM];
float ymin_mask_0[NUM_CAM];
float ymin_mask_1[NUM_CAM];
float ywinrod_0[NUM_CAM*NUM_WIN_SUP_RODS];
float ywinrod_1[NUM_CAM*NUM_WIN_SUP_RODS];
float winrod_dia[NUM_CAM*NUM_WIN_SUP_RODS];
float zcalwire_0[NUM_CAM];
float zcalwire_1[NUM_CAM];
float xcalwire_0[NUM_CAM];
float xcalwire_1[NUM_CAM];
float calwire_dia[NUM_CAM];
int npat[NUM_CAM];
int nelem[NUM_CAM];
int nwinrod[NUM_CAM];
int maskpat[NUM_PATTERNS*(NUM_MASK_ELEM+2)];
} multifield_consts_t;

typedef struct {
	IMAGEDATA_T lambda;
	IMAGEDATA_T beta;
	int n;
	int res;
	IMAGEDATA_T ThetaX;
	IMAGEDATA_T ThetaY;
	int ifverbose;
	IMAGEDATA_T xfrac;
	IMAGEDATA_T yfrac;
	IMAGEDATA_T CosThetaFactor;
	IMAGEDATA_T ht;
	IMAGEDATA_T MaskWdY;
	IMAGEDATA_T fovLimX;
	IMAGEDATA_T fovLimY;
} infov_struct_t;


/*  multifield: shadow generation context */
typedef struct {
	int inuse;
	void* 		     pctx;
	pthread_mutex_t th_consumer_mutex;
	pthread_cond_t th_consumer_sync;
	IMAGEDATA_T **bins;
 	IMAGEDATA_T **norm;
	IMAGEDATA_T **tbins;
	IMAGEDATA_T **backupbins; 
	infov_struct_t lambeta_infov;
/*device side*/
	IMAGEDATA_T *bins_d;
 	IMAGEDATA_T *norm_d;
	IMAGEDATA_T *tbins_d;
	IMAGEDATA_T *backupbins_d;
	IMAGEDATA_T *conv_d;
	IMAGEDATA_T *dummy;
	infov_struct_t *lambeta_infov_d;
	multifield_consts_t *consts_d;
/*...*/
	double		 **orig_backupbins; 
	double		 **orig_norm; 
	double		 **orig_tbins; 
	double		 **orig_bins; 
/*...*/
	pthread_t thread_id[MAX_THREADS_PER_SET];
	int iRA, idec, nRA, ndec;
	double boomRA_r, boomdec_r, incli;
	double RAval[NUM_RA], decval[NUM_DEC];
	//int nbins;
	int th_num;
	int state;
} mf_shagen_context_t;
enum th_state {
	TH_INIT= 0,
	TH_WAIT_WORK, 
	TH_WORK, 
};
/*
  Function Prototypes 
*/
	int parse(int, char *a[]);
	void SetDefaults();	
	void PrintInfo(char *);
	void PrintParam();
	void lambdabeta(double, double, double, double, double, double *, double *, double*, double *, double *);
	void infov(double *, double *, int, int *, double *, double *, double *, double *, double *, double, double, double, double);
	void RAdecval(double, double, double, double *, double *, int *, int *);
	void writeang(FILE *, double *, double *, int, int);
	void delnewline(char *);
	void freemem(double **mat);
	double** allocmem(double **mat, int nrow, int ncol);

typedef enum 
{
	NVIDIA_QUADROFX_3700=0,
}GPUcard_type;


/*  multifield: card object */
typedef struct  
{
	int use;
	/*device side*/
	GPUcard_type card_type;
	long 	cuda_cap_major;
	long 	cuda_cap_minor;
	long 	cuda_driver_version;
	long 	cuda_software_version;
	long 	card_ram;
	long 	max_threads;
	long 	max_threads_per_block;

	/* host side multithreading */
	int num_threadsets; 
	int num_threads_perset; 
	mf_shagen_context_t shagen_pool[MAX_THREAD_SETS] ; /* pool of shadow gen thread set context */

	/* host side card ops */
	int (*gpucard_init)(int cci);
	int (*gpucard_deinit)(int cci);
	int (*gpucard_get_shagen_numthreadsets)(int cci);
	int (*gpucard_get_numthread_persets)();
	int (*gpucard_create_shadowgen_pool)();
	int (*gpucard_join_shadowgen_pool)();
	int (*gpucard_shadow_gen_datacopy_fromhosttodevice)(mf_shagen_context_t *cnxt, int mvbins, int mvnorm, int mvtbins );
	int (*gpucard_shadow_gen_datacopy_fromdevicetohost)(mf_shagen_context_t *cnxt, int mvbins, int mvnorm, int mvtbins );
	int (*gpucard_shadow_gen_onetime_init)(mf_shagen_context_t *cnxt);
	int (*gpucard_shadow_gen_onetime_exit)(mf_shagen_context_t *cnxt);
	int (*gpucard_shadow_gen_perrun_init)(mf_shagen_context_t *cnxt);
	int (*gpucard_shadow_gen_perrun_exit)(mf_shagen_context_t *cnxt);
} 
gpu_card_base_t;

typedef struct {
	gpu_card_base_t cb;
}
MF_NVIDIA_card_base;

/*  Multifield: Card Manager */
typedef struct  
{
	int num_cards_detected; 
	int num_cards_active; 
	gpu_card_base_t cuda_cards[MAX_CARDS_SUPPORTED]; 
} 
GPUCOMPUTATION_card_manager;




/* INTERFACE: GPU COMPUTATION PRIMITIVES: GPUcard discovery, contract for init/deinit, shadowgen work partitioning etc.  */
int gpu_computation_init();
int gpu_computation_allocate_shadowgen_work(double incli,double RA_r,double dec_r,int nRA,int ndec,double*,double*);
int gpu_computation_wait_for_free_thread();
int gpu_computation_card_discovery();
int gpu_computation_cards_init();
int gpu_computation_create_worker_pool() ;
int gpu_computation_join_worker_pool();
/* CUDA specific impl: implements abstract contract for each card */
int generic_gpucard_create_shadowgen_pool(int cci);
/* NVIDIA_QUADROFX_3700*/
int NVIDIA_quadro_fx_3700_init(int cci);
int NVIDIA_quadro_fx_3700_deinit(int cci);

int cuda_create_context(mf_shagen_context_t *cnxt);
int cuda_attach_context(mf_shagen_context_t *cnxt);
int cuda_detach_context(mf_shagen_context_t *cnxt);
int cuda_destroy_context(mf_shagen_context_t *cnxt);
int cuda_init(mf_shagen_context_t *cnxt);
int cuda_test(mf_shagen_context_t *cnxt);

int multifield_cuda_host_alloc(mf_shagen_context_t *cnxt);
int multifield_cuda_device_alloc(mf_shagen_context_t *cnxt);
int multifield_cuda_device_to_host_xfer(mf_shagen_context_t *cnxt, int mv_bins, int mv_norm, int mv_tbins, int mv_infovparams);
int multifield_cuda_host_init(mf_shagen_context_t *cnxt);
int multifield_cuda_device_prepare(mf_shagen_context_t *cnxt);


/* Extern Globals */
	#define CUDA_CARDS g_mulfld_gpu_card_mgr.cuda_cards
	extern GPUCOMPUTATION_card_manager g_mulfld_gpu_card_mgr; 
	extern pthread_mutex_t g_producer_mutex;
	extern pthread_cond_t g_producer_sync;
	extern int g_shut;
	extern double g_boomRA, g_boomdec, g_SigmaX;
	extern double g_ExposureTime; /* for transbin - NEED TO BE REWORKED */
	extern int g_nbins, g_vFlag;
	extern char g_ShaDir[500];
	extern char g_mfp[NUM_MFP][MFP_COLSIZE];
	extern multifield_consts_t g_consts;

/* GPU COMPUTATION: Threads */
	void *image_shadow_processing_thread(void *context);
	void *image_shadow_processing_thread2(void *context);

/*Host side - CUDA card resource interacion interface */
	int cuda_shadow_gen_onetime_init(mf_shagen_context_t *cnxt);
	int cuda_shadow_gen_datacopy_fromdevicetohost(mf_shagen_context_t *, int mv_bn, int mv_nrm, int move_tbin,int mv_info);
	int cuda_shadow_gen_datacopy_fromhosttodevice(mf_shagen_context_t *, int mv_bn, int mv_nrm, int move_tbin,int mv_info);	
	int cuda_shadow_gen_perrun_exit(mf_shagen_context_t *cnxt);
	int cuda_shadow_gen_perrun_init(mf_shagen_context_t *cnxt);
	int cuda_shadow_gen_onetime_exit(mf_shagen_context_t *cnxt);

/*Host side - CUDA card init interface for transbin specific init*/
	int transbininit();
	void transbinexit();
	extern int writeshadow(mf_shagen_context_t *cnxt, int flag, FILE **fp, double ThetaX, double ThetaY, char mfp[][MFP_COLSIZE], 			double SigmaX, double ExpTime);
/*
  Multifield local - Function Prototypes 
*/
	void lambdabeta(double, double, double, double, double, double *, double *, double*, double *, double *);
	void infov(double *, double *, int, int *, double *, double *, double *, double *, double *, double, double, double, double);
	void RAdecval(double, double, double, double *, double *, int *, int *);
	void writeang(FILE *, double *, double *, int, int);

/*
  Util/Debugging Routines  - Function Prototypes 
*/
	void print_bins(int flag, float **bins, int nrow, int ncol);

#endif
