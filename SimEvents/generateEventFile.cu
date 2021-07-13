/*

	This program generates the DPH for given input angle thetaX and thetaY
	Here the total no of photons accepted is 2000000
	This program uses GPU computing, To execute it for all the four quadrant use the command

	mpirun -np 3 -host cn001 ./cudaevent 2.3215 2.3215 : -np 1 -host cn002 ./cudaevent
	This will use three gpu cards from first node and one gpu card from second node.
	Author: Ajay Vibhute
*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<time.h>
#include<sys/time.h>
#include<cuda.h>
#include <curand_kernel.h>
#include "fitsio.h"
#include <mpi.h>
#define PI 3.14159265
#define noBlock 10
#define noThread 256
#define TOTALACCEPTED 25600

//declaring the required functions
void printerror( int );
__global__ void kernel(float tx,float ty,float height,int *maskPattern,float * dphValues,int *,int *,int,int,int*,int*,int*,int*,int);
__device__ void generateEvent(float tx,float ty,float height,int *maskPattern,float * dphValues,int *,int*,int,int*,int*,int*,int*,int);
void executeKernel(float tx,float ty,int myrank,int gpuId,int*accepted,int*rejected,int);
__device__ void getDetectorIdPixelNo(int x,int y,int *moduleNo,int *pixelNo,int detectorId);
void getTimeEnergy(char * filename,int *time,int * energy);
void cudaInit(int gpuId);

//start of main
int main(int argc,char*argv[])
{

	int myrank=0,npes=0;
	float timeSpent=0.0;
	fitsfile *fptr=NULL;
	int status=0;
    	char output_filename[100] ="",temp[100],hostname[100];    
    	int bitpix   =  FLOAT_IMG; /* 16-bit unsigned short pixel values       */
 	float tx=0,ty=0;
	int gpuId=0,totalAccepted=0,totalRejected=0,totalGeneratedCount=0,totalGeneratedCountMean=0;
	int *x=NULL,*y=NULL,*detectorId=NULL,*pixelNo=NULL, *time=NULL,*energy=NULL;

	MPI_Status stat;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&npes);
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);	
	totalGeneratedCountMean++;
	totalGeneratedCount=TOTALACCEPTED;//2001920;
	//reading inputs at rank 0
	if(myrank==0)
	{
		int buf[TOTALACCEPTED];
		MPI_Buffer_attach( buf, TOTALACCEPTED+1000 );		
		timeSpent=MPI_Wtime();
				
		if(argv[1]==NULL||argv[2]==NULL)
		{
			printf("Enter value for ThetaX\n");
			scanf("%f",&tx);
			printf("Enter value for ThetaY\n");
			scanf("%f",&ty);
		}
		else
		{
			tx=(float)atof(argv[1]);
			ty=(float)atof(argv[2]);
		}
		
		/*
			Calculate the total generated Count by using poisson distribution.
		*/
		
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	//Broadcasting the inputs
	MPI_Bcast(&tx,1,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Bcast(&ty,1,MPI_FLOAT,0,MPI_COMM_WORLD);
	
	//Assigning GPU cards	
	if(myrank!=3)
	{
		gpuId=myrank;
	}
	else
	{
		//The last MPI process will run on another node, so first GPU card will be used.
		gpuId=0;
	}
	
	//calling kernal	
	executeKernel(tx,ty,myrank,gpuId,&totalAccepted,&totalRejected,totalGeneratedCount);

	gethostname(hostname,sizeof(hostname));
	if(myrank==0)
	{
		//writing event file

		time=(int*)malloc(sizeof(int)*totalGeneratedCount);
		energy=(int*)malloc(sizeof(int)*totalGeneratedCount);
		x=(int*)malloc(sizeof(int)*totalGeneratedCount);
		y=(int*)malloc(sizeof(int)*totalGeneratedCount);
		detectorId=(int*)malloc(sizeof(int)*totalGeneratedCount);
		pixelNo=(int*)malloc(sizeof(int)*totalGeneratedCount);
		int tfields   =6;      /* table will have 3 columns */
    		long nrows    = TOTALACCEPTED;       /* table will have 6 rows    */
		char extname[] = "EVENT";           /* extension name */
		char *ttype[] = { "TIME", "PHA", "DETID","PIXID","DETX","DETY" };
   		char *tform[] = { "1I",     "1I",       "1I"  ,"1I" ,"1I","1I" };
    		char *tunit[] = { "s",      "\0",       "\0"  , "\0" "\0","\0"};
		long firstrow=1, firstelem=1;
		status=0;	
		

		sprintf(temp,"%f",tx);
		strcat(output_filename,temp);
		strcat(output_filename,"_");
		bzero(temp,sizeof(temp));
		sprintf(temp,"%f",ty);
		strcat(output_filename,temp);
		strcat(output_filename,".event");
		bzero(temp,sizeof(temp));
		strcpy(temp,"rm ");
		strcat(temp,output_filename);
		system(temp);//to remove existing fits file
		bzero(temp,sizeof(temp));
			
 		if (fits_create_file(&fptr, output_filename, &status)) 
		{
        	 	printerror( status );           
			MPI_Finalize();
		}
		if ( fits_create_img(fptr,  bitpix, 0, 0, &status) )	
		{
        		printerror( status ); 	
			MPI_Finalize();
		}
		

		for(int i=0;i<npes;i++)
		{	getTimeEnergy("TimeEnergy",time,energy);
			MPI_Recv (x,totalGeneratedCount, MPI_INT, MPI_ANY_SOURCE,i+10, MPI_COMM_WORLD, &stat);
			MPI_Recv (y,totalGeneratedCount, MPI_INT, MPI_ANY_SOURCE,i+20, MPI_COMM_WORLD, &stat);
			MPI_Recv (detectorId,totalGeneratedCount, MPI_INT, MPI_ANY_SOURCE,i+30, MPI_COMM_WORLD, &stat);
			MPI_Recv (pixelNo,totalGeneratedCount, MPI_INT, MPI_ANY_SOURCE,i+40, MPI_COMM_WORLD, &stat);
			
			strcpy(extname,"Q");
			sprintf(temp,"%d",i);
			strcat(extname,temp);
			if ( fits_create_tbl( fptr, BINARY_TBL, nrows, tfields, ttype, tform,tunit, extname, &status) )
        	 		printerror( status );
			fits_write_col(fptr, TINT, 1, firstrow, firstelem, nrows, time,&status);
		        fits_write_col(fptr,TINT,2,firstrow, firstelem, nrows, energy,&status);				
			fits_write_col(fptr, TINT, 3, firstrow, firstelem, nrows, detectorId,&status);	
			fits_write_col(fptr, TINT, 4, firstrow, firstelem, nrows, pixelNo,&status);	
			fits_write_col(fptr, TINT, 5, firstrow, firstelem, nrows, x,&status);	
			fits_write_col(fptr, TINT, 6, firstrow, firstelem, nrows, y,&status);
		}
		if ( fits_close_file(fptr, &status) )                
        	printerror( status ); 
		timeSpent=MPI_Wtime()-timeSpent;
	}
	MPI_Finalize();
		
}

//Reading energy of event from input file
void getTimeEnergy(char * filename,int *time,int * energy)
{
	FILE *fp;
	int i=0;
	fp=fopen(filename,"r");
	if(fp==NULL)
	{
		printf("Error(%s:%d):Error while opening %s file\n",__FILE__,__LINE__,filename);
		exit(0);
	}
	for(i=0;i<TOTALACCEPTED;i++)
	{
		fscanf(fp,"%d",&time[i]);
		fscanf(fp,"%d",&energy[i]);
	}
	fclose(fp);
}

//Initilizing the cuda environment
void cudaInit(int gpuId)
{
	cudaSetDevice(gpuId);
	float *initmalloc;
	cudaMalloc(&initmalloc,sizeof(float));
}
//function to simulate event files
void executeKernel(float tx,float ty,int myrank,int gpuId,int *totalAccepted,int *totalRejected,int totalGeneratedCount)
{
	float height=481,*dphValues=NULL,*dphValues_device=NULL;
	int *maskPattern=NULL,*maskPattern_device=NULL,no_elements=64;
	int *acceptedCount=NULL,*rejectedCount=NULL,*acceptedCount_device=NULL,*rejectedCount_device=NULL;
	int i=0,j=0;
	char mask_fileName[100],temp[100];
	
	FILE *fp;
	cudaError_t cuerr;
	
	int noOfBlock=noBlock,noOfThread=noThread;

	int totalThreads=noOfBlock*noOfThread;	
	int totalElements=totalThreads*no_elements*no_elements;
	int countPerThread=totalGeneratedCount/totalThreads;
	int remainingCount=totalGeneratedCount%totalThreads;
	int *x,*y,*pixelNo,*detectorId;
	int *x_device,*y_device,*pixelNo_device,*detectorId_device;
	/*
		set and initilise the cuda device
	*/
	
	cudaInit(gpuId);
	
	/*
		Allocation of the memory for the host and device
	*/

	x=(int*)malloc(sizeof(int*)*(totalGeneratedCount));
	y=(int*)malloc(sizeof(int*)*(totalGeneratedCount));
	pixelNo=(int*)malloc(sizeof(int*)*(totalGeneratedCount));
	detectorId=(int*)malloc(sizeof(int*)*(totalGeneratedCount));

	maskPattern=(int*)malloc(sizeof(int)*no_elements*no_elements);
	dphValues=(float*)malloc(sizeof(float)*totalElements);
	acceptedCount=(int*)malloc(sizeof(int)*totalThreads);
	rejectedCount=(int*)malloc(sizeof(int)*totalThreads);

	/*
		Here add loop for to execute it for each quadrant 
	*/
		strcpy(mask_fileName,"maskpattern/Q");
		sprintf(temp,"%d",myrank);
		strcat(mask_fileName,temp);
		strcat(mask_fileName,"mask.dat");
		
		fp=fopen(mask_fileName,"r");
		if(fp==NULL)
		{
			printf("Error(%s:%d):%s file not exist\n",__FILE__,__LINE__,mask_fileName);	
			exit(0);
		}
		for(i=0;i<no_elements;i++)
		{
			for(j=0;j<no_elements;j++)
			{
				fscanf(fp,"%d",&maskPattern[((63-i)*no_elements)+j]);
			}
		}
		fclose(fp);	
					
		/*
			Memory Allocation for device
		*/

		if((cuerr = cudaGetLastError()) != cudaSuccess)
		{
			printf("\nError:(Pre Malloc) \"%s\"\n", cudaGetErrorString(cuerr));
		}

		cudaMalloc(&x_device,sizeof(int*)*(totalGeneratedCount));
		cudaMalloc(&y_device,sizeof(int*)*(totalGeneratedCount));
		cudaMalloc(&pixelNo_device,sizeof(int*)*(totalGeneratedCount));
		cudaMalloc(&detectorId_device,sizeof(int*)*(totalGeneratedCount));	
		cudaMalloc(&maskPattern_device,sizeof(int)*no_elements*no_elements);
		cudaMalloc(&dphValues_device,sizeof(float)*totalElements);
		cudaMalloc(&acceptedCount_device,sizeof(int)*totalThreads);
		cudaMalloc(&rejectedCount_device,sizeof(int)*totalThreads);
		if((cuerr = cudaGetLastError()) != cudaSuccess)
		{
			printf("\nError:(Pre CudaMemcpy) \"%s\"\n", cudaGetErrorString(cuerr));
		}		
		//coping mask pattern
		cudaMemcpy(maskPattern_device,maskPattern,sizeof(int)*no_elements*no_elements,cudaMemcpyHostToDevice);
		if((cuerr = cudaGetLastError()) != cudaSuccess)
		{
			printf("\nError:(Post memcpy) \"%s\"\n", cudaGetErrorString(cuerr));
		}
		//calling kernel
		kernel<<<noOfBlock,noOfThread>>>(tx,ty,height,maskPattern_device,dphValues_device,acceptedCount_device,rejectedCount_device,countPerThread,remainingCount,x_device,y_device,pixelNo_device,detectorId_device,myrank);

		//coping back the results
		cudaMemcpy(dphValues,dphValues_device,sizeof(float)*totalElements,cudaMemcpyDeviceToHost);
		cudaMemcpy(acceptedCount,acceptedCount_device,sizeof(int)*totalThreads,cudaMemcpyDeviceToHost);
		cudaMemcpy(rejectedCount,rejectedCount_device,sizeof(int)*totalThreads,cudaMemcpyDeviceToHost);

		cudaMemcpy(x,x_device,sizeof(int)*totalGeneratedCount,cudaMemcpyDeviceToHost);
		cudaMemcpy(y,y_device,sizeof(int)*totalGeneratedCount,cudaMemcpyDeviceToHost);
		cudaMemcpy(pixelNo,pixelNo_device,sizeof(int)*totalGeneratedCount,cudaMemcpyDeviceToHost);
		cudaMemcpy(detectorId,detectorId_device,sizeof(int)*totalGeneratedCount,cudaMemcpyDeviceToHost);

		if((cuerr = cudaGetLastError()) != cudaSuccess)
		{
			printf("\nError:(Post dph kernel) \"%s\"\n", cudaGetErrorString(cuerr));
		}

		MPI_Send (x,totalGeneratedCount, MPI_INT, 0,myrank+10, MPI_COMM_WORLD );
		MPI_Send (y,totalGeneratedCount, MPI_INT, 0, myrank+20, MPI_COMM_WORLD );
		MPI_Send (pixelNo,totalGeneratedCount, MPI_INT, 0, myrank+30, MPI_COMM_WORLD );
		MPI_Send (detectorId,totalGeneratedCount, MPI_INT, 0, myrank+40, MPI_COMM_WORLD );
		
		//releasing the memory
		cudaFree(&y_device);
		cudaFree(&x_device);
		cudaFree(&pixelNo_device);
		cudaFree(&detectorId_device);
		cudaFree(&maskPattern_device);
		cudaFree(&dphValues_device);	
		cudaFree(acceptedCount_device);	
		cudaFree(rejectedCount_device);	
		free(&maskPattern);
		free(&acceptedCount);
		free(&rejectedCount);

}//end of main
__global__ void kernel(float tx,float ty,float height,int *maskPattern,float * dphValues,int *acceptedCount,int *rejectedCount,int countPerThread,int remainingCount,int*x,int*y,int*pixelNo,int*detectorId,int quadrantId)
{	int index=blockIdx.x*blockDim.x+threadIdx.x;
	//Generating event using ray tracing	
	generateEvent(tx,ty,height,maskPattern,&dphValues[index*64*64],&acceptedCount[index],&rejectedCount[index],countPerThread,&x[index*countPerThread],&y[index*countPerThread],&pixelNo[index*countPerThread],&detectorId[index*countPerThread],quadrantId);
}

//function to convert length to detector module and pixel number
__device__ void getDetectorIdPixelNo(int x,int y,int *moduleNo,int *pixelNo,int quadrantId)
{
	
	int xtemp=0,ytemp=0;
	int moduleRows=16,moduleCols=16;

	*moduleNo=(x/16)+((y/16)*4);
	if(*moduleNo<4)
		*moduleNo+=12;
	else if(*moduleNo>=4 && *moduleNo<8)
		*moduleNo+=4;
	else if(*moduleNo>=8 && *moduleNo<12)			
		*moduleNo-=4;
	else
		*moduleNo-=12;
		
	xtemp=x;
	ytemp=y;
	while(xtemp>=16)
		xtemp-=16;
	while(ytemp>=16)
		ytemp-=16;
	*pixelNo=(moduleRows-1-ytemp)*moduleCols+xtemp;
	if(quadrantId==1||quadrantId==2)
	{
		*moduleNo=15-*moduleNo;
	}
}

//function to perform ray tracing
__device__ void generateEvent(float tx,float ty,float height,int *maskPattern,float * dphValues,int *accepted,int *rejected,int countPerThread,int*x_id,int*y_id,int *pixelNo,int *detectorId ,int quadrantId )
{	
	float binSize=2.5;
	float xmin=0,ymin=0,thetaX=0.0,thetaY=0.0;
	float mask_lower_left_x=0,mask_lower_left_y=0,x=0,y=0,x_dect=0,y_dect=0;;
	int acceptedCount=0,rejectedCount=0,totalCount=0;
	thetaX=ty*(PI/180);
	thetaY=tx*(PI/180);
	float x_mask=0,y_mask=0;
	//inilizing cuda random 
	curandState localState;
	curand_init(0,clock(), 0, &localState);

	for(int i=0;i<64;i++)
	{
		for(int j=0;j<64;j++)
			dphValues[i*64+j]=0;
	}
	while(acceptedCount<countPerThread)
	{
		totalCount++;
		xmin=0;
		ymin=0;
		//generating random lengths
		x=curand_uniform (&localState);
		x*=159.5;
		
		y=curand_uniform (&localState);
		y*=159.5;
		//getting pixel on the mask plate		
		x_mask=(x)*10;
		x_mask=(int)x_mask;
		x_mask/=10;
				
		y_mask=y*10;
		y_mask=(int)y_mask;
		y_mask/=10;	
				
		while(1)
		{
			if(fmod(x_mask,binSize)<0.01)
				break;
			x_mask-=0.01;
		}
		while(1)
		{
			if(fmod(y_mask,binSize)<0.01)
				break;
			y_mask-=0.01;
		}
		
		x_mask*=10;
		x_mask=(int)x_mask;
		x_mask/=10;

		y_mask*=10;
		y_mask=(int)y_mask;
		y_mask/=10;		
		
	//checking mask is open  or close	
	if(maskPattern[(int)(((x_mask/binSize)*64)+(y_mask/binSize))]==1)
	{
	
		x_dect=x-(height*(tan(thetaX)));
		y_dect=y-(height*(tan(thetaY)));
		while(x>=((xmin+16)*binSize))
		{
			if(xmin==48.0f)
				break;
			xmin+=16;
		}
		while(y>=((ymin+16)*binSize))
		{
			if(ymin==48.0f)
				break;
			ymin+=16;
		}
		
		if(x_dect<(xmin*binSize)|| (((xmin+16)*binSize)-x_dect)<0.1||y_dect<(ymin*binSize)||(((ymin+16)*binSize)-y_dect)<0.1)
		{
			rejectedCount++;
		}	
		else
		{
			mask_lower_left_x=(x_dect)*10;
			mask_lower_left_x=(int)mask_lower_left_x;
			mask_lower_left_x/=10;
				
			mask_lower_left_y=y_dect*10;
			mask_lower_left_y=(int)mask_lower_left_y;
			mask_lower_left_y/=10;	
				
			while(1)
			{
				if(fmod(mask_lower_left_x,binSize)<0.01)
					break;
				mask_lower_left_x-=0.01;
			}
			while(1)
			{
				if(fmod(mask_lower_left_y,binSize)<0.01)
					break;
				mask_lower_left_y-=0.01;
			}
			mask_lower_left_x*=10;
			mask_lower_left_x=(int)mask_lower_left_x;
			mask_lower_left_x/=10;

			mask_lower_left_y*=10;
			mask_lower_left_y=(int)mask_lower_left_y;
			mask_lower_left_y/=10;
			/*Get Pixel No and detector Id*/
			x_id[acceptedCount]=(int)(mask_lower_left_y/2.5);
			y_id[acceptedCount]=(int)(mask_lower_left_x/2.5);
			
			getDetectorIdPixelNo(x_id[acceptedCount],y_id[acceptedCount],&detectorId[acceptedCount],&pixelNo[acceptedCount],quadrantId);
			acceptedCount++;
			
		}//end of else i.e., pixel is outside of the detector
	}
	else
	{
		rejectedCount++;
	}
	}//end of the while
	accepted[0]=acceptedCount;
	rejected[0]=rejectedCount;
	
}

void printerror( int status)
{
    if (status)
    {
       fits_report_error(stderr, status); 
       exit( status );    
    }
    return;
}

