#include<stdio.h>
#include<cuda.h>
#include<fitsio.h>
#include<pil.h>
#include<GetCorrelation.h>
#include"svdfit.c"
/*
 *  nvcc GetCorrelation.cu  $INCL $OTHERLIB -arch=sm_20
 *
 */
int main(int argc,char*argv[])
{
	//declaring the local variables
	char parfilename[BUFSIZE],dphfilename[BUFSIZE],correlationFilename[BUFSIZE],qeFilename[BUFSIZE],badpixfilename[BUFSIZE];
	int moddph_size=sizeof(float)*NUMPIXPERMODULE*NUMPIXPERMODULE;
	float *mod_dph,*qemap,*dphWithoutBadPix;
	int moduleid,llx,lly,numelements_x,numelements_y;
	int i=0,j=0,rv=0,ii=0,jj=0,numbadpix=0;
	int *mask,*badpix,*quad_dph,pixid=0;
	float x_center,y_center,*correlation,zs;// xs,ys*shadow,
	float xstart=-19.5,ystart=-19.5,xend=19.5,yend=19.5;
	float resol=1,mean,rms;
	time_t start;
	int tempsrc=0;
	float *dph_device,*correlation_device=NULL,*qemap_device=NULL;
	float fittedparam=0;
	int *mask_device=NULL,numblocks-0,threadsPerBlock=0;
	FILE *fp;
	start=clock();
	char  *CZTHOME;
	cudaError_t cuerr;
	float maxvalue=0;
	int max_x=0,max_y=0;
	void (*funx_pointer)(float,float*)=polyfunction;
	float *funy,*funxy,*coeff,chisq,*shadow;
	float xs,ys,minval=9999999999,min_x=0,min_y=0,xt_tmp=xstart,yt_tmp=ystart,xe_tmp=xend,ye_tmp=yend;
	int fitcounter=0;


	//getting path for CZTI environment
	CZTHOME = getenv ("CZTWORKSPACE");
	if(CZTHOME==NULL)
	{
		printf("CZTHOME Variable is not set\n");
		exit(0);
	}
	//Setting par file name
	strcpy(parfilename,CZTHOME);
	strcat(parfilename,"/GetCorrelationGPU/par/GetCorrelation");
	PILSetModuleName(parfilename);
	int r=PILInit(argc,argv);
	if(r<0)
	{
		printf("Error(%s:%d) : Error while loading par file\n",__FILE__,__LINE__);
		exit(0);
	}

	//reading inputs
	r=PILGetInt("moduleid",&moduleid);
	r=PILGetReal4("Zs",&zs);
	r=PILGetReal4("xstart",&xstart);
	r=PILGetReal4("xend",&xend);
	r=PILGetReal4("ystart",&ystart);
	r=PILGetReal4("yend",&yend);
	r=PILGetReal4("resolution",&resol);
	r=PILGetFname("dphfilename",dphfilename);
	r=PILGetFname("badpixfilename",badpixfilename);
	r=PILGetFname("qeFilename",qeFilename);
	r=PILGetFname("correlationfilename",correlationFilename);
	PILClose(r);

	//allocating the memory
	numelements_x=(int)(((xend-xstart)/resol)+1);
	numelements_y=(int)(((yend-ystart)/resol)+1);
	badpix=(int*)malloc(sizeof(int)*ROWSPERMODULE*COLSPERMODULE);
	mask=(int*)malloc(sizeof(int)*ROWSPERMODULE*COLSPERMODULE);
	correlation=(float*)malloc(sizeof(float)*numelements_y*numelements_x);
	quad_dph=(int*)malloc(sizeof(int)*NUMPIXPERQUAD*NUMPIXPERQUAD);
	mod_dph=(float*)malloc(moddph_size);
	qemap=(float*)malloc(sizeof(float)*NUMPIXPERMODULE*NUMPIXPERMODULE);
	dphWithoutBadPix=(float*)malloc(sizeof(float)*NUMPIXPERMODULE*NUMPIXPERMODULE);

	if(mask==NULL || badpix==NULL ||correlation==NULL ||quad_dph==NULL ||mod_dph==NULL ||qemap==NULL ||dphWithoutBadPix==NULL)
	{
		printf("Error while allocating memory\n");
		exit(0);
	}
	//Reading DPH, QE files
	readImage(dphfilename,quad_dph,2);
	readFloatImage(qeFilename,2,qemap);
	getminmax(moduleid,&llx,&lly,&x_center,&y_center);
	getModule(moduleid+1,mask);
	
	//initilizing the arrays
	for(i=0;i<numelements_x;i++)
	{
		for(j=0;j<numelements_y;j++)
			correlation[j*NUMPIXPERMODULE+i]=0;
	}
	for(i=llx,ii=0;i<llx+NUMPIXPERMODULE;i++,ii++)
	{
		for(j=lly,jj=0;j<lly+NUMPIXPERMODULE;j++,jj++)
			mod_dph[jj*NUMPIXPERMODULE+ii]=quad_dph[j*NUMPIXPERQUAD+i];
	}
	
	//reading bad pixel list
	fp=fopen(badpixfilename,"r");
	if(fp==NULL)
	{
		printf("%s file not found\n",badpixfilename);
		exit(0);
	}
	while(1)
	{
		rv=fscanf(fp,"%d",&badpix[numbadpix]);
		if(rv==-1)
			break;
		numbadpix++;
	}
	//here assign -1 for the bad pixels so they can be ignored from computations
	int isbad=0;
	for(i=0;i<NUMPIXPERMODULE;i++)
	{
		for(j=0,isbad=0;j<NUMPIXPERMODULE;j++)
		{
			pixid=(15-i)*16+j;
			for(ii=0;ii<numbadpix;ii++)
			{
					if(pixid==badpix[ii])
					{
						isbad=1;break;
					}
			}
			if(isbad)
			{
				dphWithoutBadPix[pixid]=-1;
			}
			else
			{
				dphWithoutBadPix[pixid]=mod_dph[pixid];
			}
		}
	}

	
	numblocks=(int)ceil(((numelements_x*numelements_y)/(float)NUMTHREADS));
	printf("NUM BLOCKS:%d\t",numblocks);
	printf("Threads per block:%d\n",threadsPerBlock);

	//allocating device memory
	cudaMalloc((void **) &dph_device,moddph_size );
	cudaMalloc((void **) &qemap_device,moddph_size );
	cudaMalloc((void **) &mask_device,sizeof(int)*ROWSPERMODULE*COLSPERMODULE );
	cudaMalloc((void **) &correlation_device, sizeof(float)*numelements_x*numelements_y);
	//checking error in memory allocation
	if((cuerr = cudaGetLastError()) != cudaSuccess)
	{
		printf("\nError: Cuda Malloc %s\n", cudaGetErrorString(cuerr));
		return -1;
	}
	//coping data from host to device memory
	cudaMemcpy(dph_device,dphWithoutBadPix,moddph_size,cudaMemcpyHostToDevice);
	cudaMemcpy(qemap_device,qemap,moddph_size,cudaMemcpyHostToDevice);
	cudaMemcpy(mask_device,mask,sizeof(int)*ROWSPERMODULE*COLSPERMODULE ,cudaMemcpyHostToDevice);
	if((cuerr = cudaGetLastError()) != cudaSuccess)
	{
		printf("\nError: Cuda Memcpy %s\n", cudaGetErrorString(cuerr));
		return -1;
	}
	
	//calling cross-correlation kernel
	CorrelationKernel<<<numblocks,threadsPerBlock>>>(mask_device,dph_device,qemap_device,correlation_device,numelements_x,numelements_y,xstart,ystart,resol,zs,moduleid);
	if((cuerr = cudaGetLastError()) != cudaSuccess)
	{
		printf("\nError: CUDA KERNEL ERROR: %s\n", cudaGetErrorString(cuerr));
		return -1;
	}
	//coping back the results from device to host memory
	cudaMemcpy(correlation,correlation_device,sizeof(float)*numelements_x*numelements_y,cudaMemcpyDeviceToHost);
	if((cuerr = cudaGetLastError()) != cudaSuccess)
	{
		printf("\nError: Cuda DEVICE TO HOST Memcpy: %s\n", cudaGetErrorString(cuerr));
		return -1;
	}
	//getting maximum cross-correlation value and preparing array to write FITS image 
	for(i=0;i<numelements_x;i++)
	{	for(j=0;j<numelements_y;j++)
		{
			if(correlation[i*numelements_y+j]>maxvalue){
				maxvalue=correlation[i*numelements_y+j];
				max_x=i;
				max_y=j;
			}
		}
	}
	printf("Max_x:%d\tMax_y:%d\n",max_x,max_y);
	printf("Peak:%f\tX:%f\tY:%f\n",maxvalue,xstart+(max_x*resol),ystart+(max_y*resol));
	writeFloatImage(correlationFilename,correlation,numelements_x,numelements_y,resol);
	
	//performing forward fitting to get source intensity
	shadow=(float*)malloc(sizeof(float)*(NUMPIXPERMODULE*NUMPIXPERMODULE)+2);
	funy=(float*)malloc(sizeof(float)*(NUMPIXPERMODULE*NUMPIXPERMODULE)+2);
	funxy=(float*)malloc(sizeof(float)*(NUMPIXPERMODULE*NUMPIXPERMODULE)+2);
	coeff=(float*)malloc(sizeof(float)*(NUMPIXPERMODULE*NUMPIXPERMODULE)+2);
	xstart=xstart+(max_x*resol);
	ystart=ystart+(max_y*resol);
	int windowsize=3;
	printf("Fitting source for iteration :%d\n",tempsrc);
	for(xs=xstart-windowsize;xs<=xstart+windowsize;xs+=0.25)
	{
		for(ys=ystart-windowsize;ys<=ystart+windowsize;ys+=0.25)
		{
			gensha(moduleid,shadow,xs,ys,zs,mask);
			fitcounter=0;
			for(i=0;i<NUMPIXPERMODULE*NUMPIXPERMODULE;i++){
				if(dphWithoutBadPix[i]!=-1)
				{
					funy[fitcounter+1]=shadow[i];
					funxy[fitcounter+1]=dphWithoutBadPix[i];
					fitcounter++;
				}
			}
			do_svdfit(funy,funxy,coeff,fitcounter,1,&chisq,funx_pointer);
			chisq/=fitcounter-1;
			if(chisq<minval)
			{
				minval=chisq;
				fittedparam=coeff[1];
				min_x=xs;
				min_y=ys;
			}
		}
	}
	char splot[BUFSIZ],matrix[BUFSIZE],xplot[BUFSIZE],yplot[BUFSIZE];
	FILE *fp1,*fp2,*fp3,*fp4;

	strcpy(splot,correlationFilename);
	strcpy(matrix,correlationFilename);
	strcpy(xplot,correlationFilename);
	strcpy(yplot,correlationFilename);

	strcat(splot,"_S.txt");
	strcat(matrix,"_Corr.txt");
	strcat(xplot,"_X.txt");
	strcat(yplot,"_Y.txt");

	fp1=fopen(splot,"w");
	fp2=fopen(matrix,"w");
	fp3=fopen(xplot,"w");
	fp4=fopen(yplot,"w");
	for(xs=xt_tmp,i=0;xs<=xe_tmp;xs+=0.25,i++)
	{
		for(ys=yt_tmp,j=0;ys<=ye_tmp;ys+=0.25,j++)
		{
			fprintf(fp1,"%f\t%f\t%f\n",xs,ys,correlation[i*numelements_y+j]);
			fprintf(fp2,"%f\t",correlation[i*numelements_y+j]);
			fprintf(fp3,"%f\t%f\n",xs,correlation[i*numelements_y+j]);
			fprintf(fp4,"%f\t%f\n",ys,correlation[i*numelements_y+j]);
		}
		fprintf(fp2,"\n");
	}
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	fclose(fp4);

	printf("Minimum Chisq is %f found at %f,%f\n",minval,min_x,min_y);
	printf("Fitted Param:%f\n",fittedparam);
	for(i=0;i<NUMPIXPERMODULE*NUMPIXPERMODULE;i++){
		if(dphWithoutBadPix[i]!=-1)
		{
			dphWithoutBadPix[i]-=shadow[i]*fittedparam;
		}
	}



	printf("TIME ELAPSED : %f seconds\n\n", ((double)clock() - start) / CLOCKS_PER_SEC);
}
void polyfunction(float x,float*funx)
{
	funx[1]=x;
}
__global__ void CorrelationKernel(int *mask,float *dph,float *qemap,float *correlation,int numelem_x,int numelem_y,float x_start,float y_start,float resol,float zs,int moduleid)
{
	int index=blockIdx.x*blockDim.x+threadIdx.x;
	int i=0;
	float xs,ys;
	float *shadow,sum=0;
	xs=(x_start)+((index/numelem_y)*resol);
	ys=(y_start)+((index%numelem_y)*resol);
	if(xs>19.5 || xs<-19.5 || ys>19.5 || ys<-19.5)
		return;
	//allocating memory
	shadow=(float*)malloc(sizeof(float)*ELEM);
	//computing shadow
	gensha_device(moduleid,shadow,xs,ys,zs,mask);
	
	//Computing cross-correlation value for a direction
	for(i=0;i<NUMPIXPERMODULE*NUMPIXPERMODULE;i++)
	{
		if(qemap[i]!=0 && dph[i]!=-1)
			sum+=((dph[i]*shadow[i])/qemap[i]);
	}
	correlation[index]=sum;
	free(shadow);
}

//function to generate shadow
int gensha(int moduleid,float*det,float xs,float ys,float zs,int *mask)
{
	int detx=0, dety=0;
	float x_center=0.0,y_center=0.0,tx=0.0,ty=0.0,xtemp=0,ytemp=0;
	int rowcounter=0,colcounter=0, x_index,y_index;
	int i=0,j=0,ii=0,jj=0;
	int llx=0,lly=0;


	getminmax(moduleid,&llx,&lly,&x_center,&y_center);
	dety=(llx/16)*ROWSPERMODULE+(llx/16)*100;
	detx=(lly/16)*COLSPERMODULE+(lly/16)*100;

	xtemp=19.5-xs;
	ytemp=19.5+ys;
	xs=ytemp;
	ys=xtemp;
	xs/=0.02;
	ys/=0.02;
	zs/=0.02;
	int tempx,tempy;
	for(i=0;i<256;i++)
			det[i]=0;
	for(i=detx,ii=0;i<detx+ROWSPERMODULE;i++,ii++)
	{
		tx=asin( (ys-ii) / (  sqrt( (zs*zs)+((ys-ii) * (ys-ii)) )  ) );
		tempx=(int)(ii+(HEIGHT*(tan(tx))));
		x_index=(ii/123);
		if(x_index>14)
			x_index=14+((ii-(14*123))/114);
		else
			x_index=x_index+((ii-(x_index*123))/114);
		for(j=dety,jj=0;j<dety+COLSPERMODULE;j++,jj++)
		{
			y_index=(jj/123);
			if(y_index>14)
				y_index=14+((jj-(14*123))/114);
			else
				y_index=y_index+((jj-(y_index*123))/114);
			ty=asin( (xs-jj) / (  sqrt( (zs*zs)+((xs-jj) * (xs-jj)) )  ) );
			tempy=(int)(jj+(HEIGHT*(tan(ty))));
			if(tempx>=0 && tempx<ROWSPERMODULE && tempy>=0 &&tempy<COLSPERMODULE)
			{
				if(mask[tempx*COLSPERMODULE+tempy]==1)
					det[x_index*16+y_index]+=1;
			}
		}
	}
	for(i=0;i<16;i++)
	{
		if(i%16==0 ||(i+1)%16==0)
			rowcounter=114;
		else
			rowcounter=123;
		for(j=0;j<16;j++)
		{
			if(j%16==0 || (j+1)%16==0)
				colcounter=114;
			else
				colcounter=123;
			det[i*16+j]=(float)det[i*16+j]/(float)(rowcounter*colcounter);
		}
	}

	return 0;
}

//Generatee shadow function for devicee
__device__ int gensha_device(int moduleid,float*det,float xs,float ys,float zs,int *mask)
{
	int detx=0, dety=0;
	float x_center=0.0,y_center=0.0,tx=0.0,ty=0.0,xtemp=0,ytemp=0;
	int rowcounter=0,colcounter=0, x_index,y_index;
	int i=0,j=0,ii=0,jj=0;
	int llx=0,lly=0;


	getminmax_device(moduleid,&llx,&lly,&x_center,&y_center);
	dety=(llx/16)*ROWSPERMODULE+(llx/16)*100;
	detx=(lly/16)*COLSPERMODULE+(lly/16)*100;

	xtemp=19.5-xs;
	ytemp=19.5+ys;
	xs=ytemp;
	ys=xtemp;
	xs/=0.02;
	ys/=0.02;
	zs/=0.02;
	int tempx,tempy;
	for(i=0;i<256;i++)
			det[i]=0;
	for(i=detx,ii=0;i<detx+ROWSPERMODULE;i++,ii++)
	{
		tx=asin( (ys-ii) / (  sqrt( (zs*zs)+((ys-ii) * (ys-ii)) )  ) );
		tempx=(int)(ii+(HEIGHT*(tan(tx))));
		x_index=(ii/123);
		if(x_index>14)
			x_index=14+((ii-(14*123))/114);
		else
			x_index=x_index+((ii-(x_index*123))/114);
		for(j=dety,jj=0;j<dety+COLSPERMODULE;j++,jj++)
		{
			y_index=(jj/123);
			if(y_index>14)
				y_index=14+((jj-(14*123))/114);
			else
				y_index=y_index+((jj-(y_index*123))/114);
			ty=asin( (xs-jj) / (  sqrt( (zs*zs)+((xs-jj) * (xs-jj)) )  ) );
			tempy=(int)(jj+(HEIGHT*(tan(ty))));
			if(tempx>=0 && tempx<ROWSPERMODULE && tempy>=0 &&tempy<COLSPERMODULE)
			{
				if(mask[tempx*COLSPERMODULE+tempy]==1)
					det[x_index*16+y_index]+=1;
			}
		}
	}
	for(i=0;i<16;i++)
	{
		if(i%16==0 ||(i+1)%16==0)
			rowcounter=114;
		else
			rowcounter=123;
		for(j=0;j<16;j++)
		{
			if(j%16==0 || (j+1)%16==0)
				colcounter=114;
			else
				colcounter=123;
			det[i*16+j]=(float)det[i*16+j]/(float)(rowcounter*colcounter);
		}
	}

	return 0;
}
//function to write fits images
void writeFloatImage(char *filename,float*pixels,int rows,int cols,float resol)
{
	int bitpix   =  FLOAT_IMG; /* 16-bit unsigned short pixel values       */
    long naxis    =   2;  /* 2-dimensional image                            */
    int fpixel = 1,status=0;
	long naxes[2] = { cols,rows };
	long nelements = naxes[0] * naxes[1];
	fitsfile *fptr;
	remove(filename);
	if (fits_create_file(&fptr, filename, &status))
	{
		printf("Error(%s:%d):Creating file\n",__FILE__,__LINE__);
	}

	if ( fits_create_img(fptr,  bitpix,0, naxes, &status) )
	{
		printf("Error(%s:%d):Creating image\n",__FILE__,__LINE__);
	}


	if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) )
	{
		printf("Error(%s:%d):Creating image\n",__FILE__,__LINE__);
	}
	write_wcsaxis(fptr,1,"","","","IMX",((double)(rows+1)/2.0),resol,0,"mm",&status);
	write_wcsaxis(fptr,2,"","","","IMY",((double)(cols+1)/2.0),resol,0,"mm",&status);
	if ( fits_write_img(fptr, TFLOAT, fpixel, nelements,pixels, &status) )
	{
			printf("Error(%s:%d):Wrting image\n",__FILE__,__LINE__);
	}
	if ( fits_close_file(fptr, &status) )
	{
			printf("Error(%s:%d):closing file\n",__FILE__,__LINE__);
	}
}

//function to read fits images
void readFloatImage(char*filename,int hduno,float*data)
{
		fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
	    int status,  nfound, anynull,hdutype;
	    long naxes[2], fpixel, nbuffer, npixels;
	    float  nullval;
	    int buffsize;
	    status = 0;
	    if ( fits_open_file(&fptr, filename, READONLY, &status) )
	    	printf("Error while reading fits file\n");

	   if ( fits_movabs_hdu(fptr, hduno, &hdutype, &status) )
	           printf("Error while moving HDU\n");
	    if ( fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) )
	    	  printf("Error while reading keys\n");
	    npixels  = naxes[0] * naxes[1];         /* number of pixels in the image */
	    buffsize=npixels;
	    fpixel   = 1;
	    nullval  = 0;                /* don't check for null values in the image */
	    while (npixels > 0)
	    {
	      nbuffer = npixels;
	      if (npixels > buffsize)
	        nbuffer = buffsize;
	      if ( fits_read_img(fptr, TFLOAT, fpixel, nbuffer, &nullval,data, &anynull, &status) )
	    	  printf("Error while reading fits image\n");
	      npixels -= nbuffer;
	      fpixel  += nbuffer;
	    }
	    if ( fits_close_file(fptr, &status) )
	    	  printf("Error while closing file\n");

}

//Writing WCS information to images
int write_wcsaxis(fitsfile *imgfile, int axis, char *suffix,char *wcsname, char *wcstype, char *ctype, double crpix, double cdelt, double crval,char *cunit, int *status)
{
	char key[20];

	if (status == 0) return NULL_INPUT_PTR;
	if (*status != 0) return (*status);
	if (imgfile == 0) return (*status = NULL_INPUT_PTR);

	if (wcsname && wcsname[0]) {
		sprintf(key, "WCSNAME%s", suffix);
		fits_update_key(imgfile, TSTRING, key, wcsname,"Coordinate system name", status);
	}
	if (wcstype && wcstype[0]){
		sprintf(key, "WCSTY%d%s", axis, suffix);
		fits_update_key(imgfile, TSTRING, key, wcstype,"Coordinate system axis", status);
	}
	sprintf(key, "CTYPE%d%s", axis, suffix);
	fits_update_key(imgfile, TSTRING, key, ctype,"Name of coordinate", status);

	if (cunit && cunit[0]) {
		sprintf(key, "CUNIT%d%s", axis, suffix);
		fits_update_key(imgfile, TSTRING, key, cunit,"Units of coordinate axis", status);

	}
	sprintf(key, "CRPIX%d%s", axis, suffix);
	fits_update_key(imgfile, TDOUBLE, key, &crpix,"Reference pixel position", status);
	sprintf(key, "CDELT%d%s", axis, suffix);
	fits_update_key(imgfile, TDOUBLE, key, &cdelt,"Pixel spacing in physical units", status);
	sprintf(key, "CRVAL%d%s", axis, suffix);
	fits_update_key(imgfile, TDOUBLE, key, &crval,"Coordinate value at reference pixel position", status);
	return (*status);
}

void normalizeData(float *data,int numelements)
{
	float peakval=getMaximum(data,numelements);
	int i=0;
	for(i=0;i<numelements;i++)
	{
		data[i]/=peakval;
	}
}
float getMaximum(float *data,int numelements)
{
	int i=0;
	float maxval=0;
	for(i=0;i<numelements;i++)
	{
		if(data[i]>maxval)
			maxval=data[i];
	}
	return maxval;
}

void readImage(char*filename,int *buffer,int hduNo)
{
	fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
	int status,  nfound, anynull,hdutype;
	long naxes[2], fpixel, nbuffer, npixels;
	float nullval;
	status = 0;

	if ( fits_open_file(&fptr, filename, READONLY, &status) )
		printf("Error while opening fits file\n");

	if ( fits_movabs_hdu(fptr, hduNo, &hdutype, &status) )
		printf("Error while moving module\n");

	if ( fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) )
		printf("Error while reading keys\n");

	npixels  = naxes[0] * naxes[1];         /* number of pixels in the image */
	fpixel   = 1;
	nullval  = 0;                /* don't check for null values in the image */
	while (npixels > 0)
	{
		nbuffer = npixels;
		if ( fits_read_img(fptr, TINT, fpixel, nbuffer, &nullval,buffer, &anynull, &status) )
		{
			printf("Error reading fits image\n");
		}
		npixels -= nbuffer;
		fpixel  += nbuffer;
	}
	if ( fits_close_file(fptr, &status) )
	{
		printf("Error while closing fits file\n");
	}
	return;
}
__device__ void getminmax_device(int moduleno,int*x,int *y,float* x_center,float* y_center)
{
	switch(moduleno)
	{
	case 12 :
		*x=0;
		*y=0;
		*x_center=19.5;
		*y_center=19.5;
		break;
	case 13 :
		*x=16;
		*y=0;
		*x_center=60.5;
		*y_center=19.5;
		break;
	case 14 :
		*x=32;
		*y=0;
		*x_center=101.5;
		*y_center=19.5;
		break;
	case 15 :
		*x=48;
		*y=0;
		*x_center=142.5;
		*y_center=19.5;
		break;
	case 8:
		*x=0;
		*y=16;
		*x_center=19.5;
		*y_center=60.5;
		break;
	case 9 :
		*x=16;
		*y=16;
		*x_center=60.5;
		*y_center=60.5;
		break;
	case 10 :
		*x=32;
		*y=16;
		*x_center=101.5;
		*y_center=60.5;
		break;
	case 11 :
		*x=48;
		*y=16;
		*x_center=142.5;
		*y_center=60.5;
		break;
	case 4 :
		*x=0;
		*y=32;
		*x_center=19.5;
		*y_center=101.5;
		break;
	case 5 :
		*x=16;
		*y=32;
		*x_center=60.5;
		*y_center=101.5;
		break;
	case 6 :
		*x=32;
		*y=32;
		*x_center=101.5;
		*y_center=101.5;
		break;
	case 7 :
		*x=48;
		*y=32;
		*x_center=142.5;
		*y_center=101.5;
		break;
	case 0 :
		*x=0;
		*y=48;
		*x_center=19.5;
		*y_center=142.5;
		break;
	case 1 :
		*x=16;
		*y=48;
		*x_center=60.5;
		*y_center=142.5;
		break;
	case 2 :
		*x=32;
		*y=48;
		*x_center=101.5;
		*y_center=142.5;
		break;
	case 3 :
		*x=48;
		*y=48;
		*x_center=142.5;
		*y_center=142.5;
		break;
	default:
		printf("Invalid module id\n");
		break;
	}
}
void getminmax(int moduleno,int*x,int *y,float* x_center,float* y_center)
{
	switch(moduleno)
	{
	case 12 :
		*x=0;
		*y=0;
		*x_center=19.5;
		*y_center=19.5;
		break;
	case 13 :
		*x=16;
		*y=0;
		*x_center=60.5;
		*y_center=19.5;
		break;
	case 14 :
		*x=32;
		*y=0;
		*x_center=101.5;
		*y_center=19.5;
		break;
	case 15 :
		*x=48;
		*y=0;
		*x_center=142.5;
		*y_center=19.5;
		break;
	case 8:
		*x=0;
		*y=16;
		*x_center=19.5;
		*y_center=60.5;
		break;
	case 9 :
		*x=16;
		*y=16;
		*x_center=60.5;
		*y_center=60.5;
		break;
	case 10 :
		*x=32;
		*y=16;
		*x_center=101.5;
		*y_center=60.5;
		break;
	case 11 :
		*x=48;
		*y=16;
		*x_center=142.5;
		*y_center=60.5;
		break;
	case 4 :
		*x=0;
		*y=32;
		*x_center=19.5;
		*y_center=101.5;
		break;
	case 5 :
		*x=16;
		*y=32;
		*x_center=60.5;
		*y_center=101.5;
		break;
	case 6 :
		*x=32;
		*y=32;
		*x_center=101.5;
		*y_center=101.5;
		break;
	case 7 :
		*x=48;
		*y=32;
		*x_center=142.5;
		*y_center=101.5;
		break;
	case 0 :
		*x=0;
		*y=48;
		*x_center=19.5;
		*y_center=142.5;
		break;
	case 1 :
		*x=16;
		*y=48;
		*x_center=60.5;
		*y_center=142.5;
		break;
	case 2 :
		*x=32;
		*y=48;
		*x_center=101.5;
		*y_center=142.5;
		break;
	case 3 :
		*x=48;
		*y=48;
		*x_center=142.5;
		*y_center=142.5;
		break;
	default:
		printf("Invalid module id\n");
		break;
	}
}
void getModule(int moduleNo,int *pixels)
{


	int i=0,j=0,ii=0,jj=0,temp=0;
	int cols=COLSPERMODULE;
	char moduleFileName[100]="module",blockageFileName[100]="bars_",tempBuff[100];
	FILE *fp,*blockage;
	int *data;
	char *CZTHOME = getenv ("CZTWORKSPACE");
	if(CZTHOME==NULL)
	{
		printf("CZTHOME Variable is not set\n");
		exit(0);
	}
	data=(int*)malloc(sizeof(int)*16*16);
	sprintf(tempBuff,"%d",moduleNo);
	strcpy(moduleFileName,CZTHOME);
	strcpy(blockageFileName,CZTHOME);
	strcat(moduleFileName,"/config/module");
	strcat(blockageFileName,"/config/bars_");
	strcat(moduleFileName,tempBuff);
	strcat(blockageFileName,tempBuff);
	fp=fopen(moduleFileName,"r");
	blockage=fopen(blockageFileName,"r");

	if(fp==NULL)
	{
		printf("Sorry Error while opening the module file\n%s\n",moduleFileName);
		exit(0);
	}
	if(blockage==NULL)
	{
		printf("Error while opening the blockage file\n%s\n",blockageFileName);
		exit(0);
	}
	int colPreVal=0,colInc=0,rowPreVal=0,rowInc=0;
	for(i=0;i<16;i++)
	{
		if(i==0||i==15)
		{
			rowInc=114;
		}
		else
		{
			rowInc=123;
		}

		for(j=0,colPreVal=0;j<16;j++)
		{
			if(j==0||j==15)
			{
				colInc=114;
			}
			else
			{
				colInc=123;
			}
			fscanf(fp,"%d",&data[i*16+j]);
			for(ii=rowPreVal;ii<rowPreVal+rowInc-10;ii++)
			{
				for(jj=colPreVal;jj<colPreVal+colInc;jj++)
				{

					pixels[ii*cols+jj]=data[i*16+j];
				}

			}
			temp=0;
			if(i!=15)
			{
				fscanf(blockage,"%d",&temp);
				for(ii=rowPreVal+rowInc-10;ii<rowPreVal+rowInc;ii++)
				{
					for(jj=colPreVal;jj<colPreVal+colInc;jj++)
					{
						pixels[ii*cols+jj]=temp;
					}
				}

			}
			else
			{
				for(ii=rowPreVal+rowInc-10;ii<rowPreVal+rowInc;ii++)
				{
					for(jj=colPreVal;jj<colPreVal+colInc;jj++)
					{
						pixels[ii*cols+jj]=data[i*16+j];
					}
				}


			}
			colPreVal+=colInc;
		}
		rowPreVal+=rowInc;
	}

}
void calculateSigmaClippedMean(float* pixel_count,float *mean_out,float*rms_out)
{
	int i=0;
	float mean_sum=0,rms_sum=0,mean,rms;
	int rows=NUMPIXPERMODULE,cols=NUMPIXPERMODULE;
	float temp_mean=0,temp_rms=0,temp_mean_sum=0,temp_rms_sum=0,mean_count=0,rms_count=0;
	mean_count=0;
	for(i=0;i<rows*cols;i++)
	{
		if(pixel_count[i]!=0)
		{
			mean_sum+=pixel_count[i];
			mean_count++;
		}
	}
	mean=mean_sum/(mean_count);
	rms_count=0;
	for(i=0;i<rows*cols;i++)
	{
		if(pixel_count[i]!=0)
		{
			rms_sum+=((pixel_count[i]-mean)*(pixel_count[i]-mean));
			rms_count++;
		}
	}
	rms=(rms_sum/rms_count);
	rms=sqrt(rms);

	while(1)
	{
		temp_mean_sum=0;
		temp_rms_sum=0;
		temp_mean=0;
		temp_rms=0;
		mean_count=0;
		rms_count=0;

		for(i=0;i<rows*cols;i++)
		{
			if((pixel_count[i])<((THREASHOLD*rms)+mean) || (pixel_count[i] > (mean-(THREASHOLD*rms))))
			{
				temp_mean_sum+=pixel_count[i];
				mean_count++;
			}
		}
		temp_mean=temp_mean_sum/mean_count;
		for(i=0;i<rows*cols;i++)
		{
			if((pixel_count[i])<((THREASHOLD*rms)+mean) || (pixel_count[i] > (mean-(THREASHOLD*rms))))
			{
				temp_rms_sum+=(pixel_count[i]-temp_mean)*(pixel_count[i]-temp_mean);
				rms_count++;
			}
		}
		if(mean_count!=rms_count)
			printf("Mean Rms count different\n");
		temp_rms=temp_rms_sum/rms_count;
		temp_rms=sqrt(temp_rms);

		float t1=0,t2=0;
		t1=((mean-temp_mean)/mean);
		t2=((rms-temp_rms)/rms);
		if(t1<0)
			t1*=-1;
		if(t2<0)
			t2*=-1;
		if(t1 <0.01 &&  t2<0.01 )
		{
			//mean=temp_mean;
			//rms=temp_rms;
			*mean_out=temp_mean;
			*rms_out=temp_rms;
			break;
		}
		else
		{

			mean=temp_mean;
			rms=temp_rms;
		}

	}
}
