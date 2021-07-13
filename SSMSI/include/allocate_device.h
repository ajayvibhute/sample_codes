/* Declare Device side Variables*/

/*
float *x_device,*y_device,*a_device,*u_device,*v_device,*w_device;
float *chisq_device,*sig_device,*cvm_device,*b,*func,*x_another_device;
float *thetaX_device,*thetaY_device;
int *maskpat_device;
*/
void allocateDevice(int no_another)
{
	/* allocate memory on device*/
	/*
	int noOfbins=rows;
	int col=colms;
	int noOfGrids=noOfDataSet;
	size_t size=(noOfbins+1)*noOfGrids*sizeof(float);

	cudaMalloc(&b,size);
	cudaMalloc(&func,(noOfGrids)*(col+1)*sizeof(float));
	cudaMalloc(&x_device,(noOfbins+2)*sizeof(float));
	cudaMalloc(&y_device,size);	
	cudaMalloc(&a_device,(noOfGrids+1)*(col+1)*sizeof(float));		;	
	cudaMalloc(&chisq_device,(noOfGrids+1)*sizeof(float));
	cudaMalloc(&sig_device,(noOfbins+2)*sizeof(float));	
	cudaMalloc(&u_device,size*(col+1));	
	cudaMalloc(&v_device,noOfGrids*(col+1)*(col+1)*sizeof(float));	
	cudaMalloc(&w_device,noOfGrids*(col+1)*sizeof(float));
	
	cudaMalloc(&thetaX_device,sizeof(float)*(noOfGrids+2));
	cudaMalloc(&thetaY_device,sizeof(float)*(noOfGrids+2));
	cudaMalloc(&maskpat_device,NUM_PATTERNS*(NUM_MASK_ELEM+2)*sizeof(int));
	*/
}

void freeDevice()
{
	/* Free the device Memory*/
	/*
	cudaFree(x_device);
	cudaFree(y_device);	
	cudaFree(a_device);	
	cudaFree(u_device);	
	cudaFree(v_device);	
	cudaFree(w_device);		
	cudaFree(chisq_device);
	cudaFree(sig_device);	
	cudaFree(b);
	cudaFree(func);
	cudaFree(thetaX_device);
	cudaFree(thetaY_device);
	cudaFree(maskpat_device);
	*/
}
