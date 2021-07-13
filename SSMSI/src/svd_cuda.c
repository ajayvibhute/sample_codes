//header file
#include"kernel.c"
#include "allocate_device.h"
/*
	set zero and null value for no of another function of x,and another function of x respectively and call SVD_CUDA
*/
float SVD_FIT(float *x,float *thetaX,float *thetaY,float* a,float *chisq,int *maskpat,int icam)
{
	int no_another=0;
	float x_another=0;
	float SVD_CUDA(float *x,float *thetaX,float *thetaY,float* a,float *chisq,int no_another,float *x_another,int *maskpat,int icam);
	return SVD_CUDA(x,thetaX,thetaY,a,chisq,no_another,&x_another,maskpat,icam);
}

/*
	This function allocate the memory on device and call kernel for the execution.
*/

float SVD_CUDA(float *x,float*thetaX,float*thetaY,float *a,float *chisq,int no_another,float *x_another,int *maskpat,int icam)
{
	cudaError_t cuerr;
	float *x_device=NULL,*y_device=NULL,*a_device=NULL,*u_device=NULL;
	float *chisq_device=NULL,*sig_device=NULL,*b=NULL,*func=NULL,*x_another_device=NULL;
	float *thetaX_device=NULL,*thetaY_device=NULL;
	int *maskpat_device=NULL;
	int noOfbins=rows;
	int col=colms;
	int noOfGrids=noOfDataSet;
	size_t size=(noOfbins+1)*noOfGrids*sizeof(float);
	cudaEvent_t start_event, stop_event;


	/* 
 	* 	Allocation of the memory on device
	*/


	cudaMalloc(&b,size);
	cudaMalloc(&func,(noOfGrids)*(col+1)*sizeof(float));
	cudaMalloc(&x_device,(noOfbins+2)*sizeof(float));
	cudaMalloc(&y_device,size);	
	cudaMalloc(&a_device,(noOfGrids+1)*(col+1)*sizeof(float));		;	
	cudaMalloc(&chisq_device,(noOfGrids+1)*sizeof(float));
	cudaMalloc(&sig_device,(noOfbins+2)*sizeof(float));	
	cudaMalloc(&u_device,size*(col+1));	
	cudaMalloc(&x_another_device,1);
	cudaMalloc(&thetaX_device,sizeof(float)*(noOfGrids+2));
	cudaMalloc(&thetaY_device,sizeof(float)*(noOfGrids+2));
	cudaMalloc(&maskpat_device,NUM_PATTERNS*(NUM_MASK_ELEM+2)*sizeof(int));

	if((cuerr = cudaGetLastError()) != cudaSuccess)
	{
		
		printf("\nError:(Pre SVD FIT) \"cudaMalloc %s\"\n", cudaGetErrorString(cuerr));
		return -1;
	}	

	/* 
		Cuda timer to calculate the time spent on GPU for execution
	*/

        cudaEventCreate(&start_event);
        cudaEventCreate(&stop_event);
	/*	
		Cuda error to check whether the execution is without error or with error
	*/	
	

	int noOfBlocks,threadsPerBlock,MaxThreads=512;	
	float timeSpent=0;
	
	/*copy observed shadow, x, y locations and another function of x from host to device*/
	cudaMemcpy(x_device,x,(noOfbins+1)*sizeof(float),cudaMemcpyHostToDevice);
	cudaMemcpy(thetaX_device,thetaX,sizeof(float)*(noOfGrids+2),cudaMemcpyHostToDevice);
	cudaMemcpy(thetaY_device,thetaY,sizeof(float)*(noOfGrids+2),cudaMemcpyHostToDevice);
	cudaMemcpy(maskpat_device,maskpat,NUM_PATTERNS*(NUM_MASK_ELEM+2)*sizeof(int),cudaMemcpyHostToDevice);
	/* Decides the no of Blocks and threads per block */
	if(noOfGrids<MaxThreads/2)
	{
		noOfBlocks=1;
		threadsPerBlock=noOfGrids;
	}
	else
	{
		threadsPerBlock=275;
		noOfBlocks=(int)ceil((float)noOfGrids/threadsPerBlock);	
	}
	
	if((cuerr = cudaGetLastError()) != cudaSuccess)
	{
		
		printf("\nError:(Pre SVD FIT) \"cudaMemcpy %s\"\n", cudaGetErrorString(cuerr));
		return -1;
	}
	
	//Start the cuda Timer
	cudaEventRecord(start_event,0);
	/*
		Call the kernel configuring with no of blocks and threads per block	
	*/
	Kernel<<<noOfBlocks,threadsPerBlock>>>(x_device,y_device,thetaX_device,thetaY_device,sig_device,a_device,chisq_device,u_device,b,func,no_another,x_another_device,maskpat_device,noOfGrids,icam);
	
	/*
		Cuda Error checking 
	*/
	if((cuerr = cudaGetLastError()) != cudaSuccess)
	{
		printf("\nError:(Post SVD_FIT) \"%s\" ", cudaGetErrorString(cuerr));
		
		return -1;
	}
		
	/*
		Stop the cuda timer and calculate the time spent on the GPU
	*/
	cudaEventRecord(stop_event, 0);
	cudaEventSynchronize(stop_event);   // block until the event is actually recorded
        cudaEventElapsedTime(&timeSpent, start_event, stop_event) ;
	cudaEventDestroy(start_event);
	cudaEventDestroy(stop_event);
         
	/* Copy calculated  Source strength,chisq from device to host memory */
	cudaMemcpy(a,a_device,(noOfGrids+1)*(colms+1)*sizeof(float),cudaMemcpyDeviceToHost);	
	cudaMemcpy(chisq,chisq_device,(noOfGrids+1)*sizeof(float),cudaMemcpyDeviceToHost);
	
	/*
		Free the allocated memory
	*/
	cudaFree(x_device);
	cudaFree(y_device);	
	cudaFree(a_device);	
	cudaFree(u_device);	
	cudaFree(chisq_device);
	cudaFree(sig_device);	
	cudaFree(b);
	cudaFree(func);
	cudaFree(thetaX_device);
	cudaFree(thetaY_device);
	cudaFree(maskpat_device);
	cudaFree(x_another_device);

	/* return the time spent on the execution of the program */
	return timeSpent;

}
