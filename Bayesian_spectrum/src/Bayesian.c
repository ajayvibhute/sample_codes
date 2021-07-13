/*
 ============================================================================
 Name        : Bayesian.c
 Author      : Ajay Vibhute, Oct 23, 2017
 Version     :1.0
 Copyright   : 
 Description : This code is inherited from calcintensity ratio. Code is modifed to work on PI channels rather than
energy and generate output which can be later loaded into xspec. It bins the energy per channel
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <pil.h>
#include <fitsio.h>
#include<math.h>
#include <string.h>
#include "Bayesian.h"
#include "svdfit_old.c" /*Fits a function using svd fit*/
#include "common.c"

void getDetectorIdPixelNo(int x,int y,int *moduleNo,int *pixelNo);

void getShadow(float *shadow_pixels,float tx,float ty);
int main(int argc,char*argv[])
{

	char *basedir,parfilename[BUFFSIZE],infile[BUFFSIZE],sourcelistfile[BUFFSIZE],outfile[BUFFSIZE],backfile[BUFFSIZE],responsefile[BUFFSIZE],temp[BUFFSIZE],badpixfilename[BUFFSIZE],outspecfile[BUFFSIZE];
	int r,i=0,j=0,numsources,rows,cols,numbins,numpix=4*4096;
	float tx[10],ty[10];
	float *x=NULL,*x_good=NULL,**y=NULL,**u=NULL,**v=NULL,*w=NULL,*sig=NULL,chisq=0,*a=NULL,**intenRatio,*energy,*backinten;
	float*funy,**params,numparams=2,**background,**shadow,**error,**shadow_new;
	FILE *fp;
	float tmp,*channelf,*response,**channelDPH;
	int minchannel=0,maxchannel=512,binsize=1,*channel;
	int LLD[4096],pix_flag[4096];
	char specTemplate[BUFFSIZE],LLDFILE[BUFFSIZE];

	/*Get LLD file path from CALDB*/	
	basedir=getenv("CALDB");
	if(basedir==NULL)
	{
		printf("Error (%s:%d):CALDB variable is not set\n",__FILE__,__LINE__);
		exit(1);
	}
	strcpy(LLDFILE,basedir);
	strcat(LLDFILE,"/data/as1/czti/bcf/AS1cztlld20170324v02.fits");

	/*Get PATH for spectrum template*/

	basedir=getenv("as1czt");
	if(basedir==NULL)
	{
		printf("Error (%s:%d):as1czt variable is not set\n",__FILE__,__LINE__);
		exit(1);
	}
	strcpy(specTemplate,basedir);
	strcat(specTemplate,"/templates/spectrumTemplate");



	/*Get the working directory from environment */
	basedir=getenv("WORKSPACE");
	if(basedir==NULL)
	{
		printf("Error (%s:%d):WORKSPACE variable is not set\n",__FILE__,__LINE__);
		exit(1);
	}
	strcpy(parfilename,basedir);
	strcat(parfilename,"CalcIntensityRatio/par/CalcIntensityRatio");
	PILSetModuleName(parfilename);

	r=PILInit(argc,argv);//add error checking
	r=PILGetFname("infile",infile);
	r=PILGetFname("responsefile",responsefile);
	r=PILGetFname("sourcelistfile",sourcelistfile);
	r=PILGetFname("outfile",outfile);
	r=PILGetFname("badpixfile",badpixfilename);
	//r=PILGetFname("backfile",backfile);
	PILClose(r);
	
	numbins=(maxchannel-minchannel)/binsize;
	printf("Number of bins:%d\n",numbins);
	/* 
	Read list of sources in FOV. The list can be generated using any image reconstruction algorithm.
	Assuming there are n sources in FOV, format for sourcelist file will be,
	n
	ThetaX1 ThetaY1 
	ThetaX2 ThetaY2
	   .	   .
	ThetaXn ThetaYn 
	*/
	getSourceList(sourcelistfile,&tx[0],&ty[0],&numsources);
	printf("Total Number of sources:%d\n",numsources);
	printf("Theta x\tTheta y\n");
	for(i=0;i<numsources;i++)
	{
		printf("%0.2f\t%0.2f\n",tx[i],ty[i]);
		tx[i]*=-1.0;
		ty[i]*=-1.0;
	}
	/*
	We need to remove contribution of the background. Here we considering backaground as a source. Hence, increasing
	no of sources in FOV by 1.
	*/
	numsources+=1;
	cols=numsources;	
	rows=NUMELEM*NUMQUADRANT;
	/*Allocating memory required for the execution*/

	energy=(float*)malloc(sizeof(float)*(numbins+2));
	backinten=(float*)malloc(sizeof(float)*(numbins+2));
	intenRatio=(float**)malloc(sizeof(float*)*(cols+2));
	error=(float**)malloc(sizeof(float*)*(cols+2));
	x=(float*)malloc(sizeof(float)*rows+2);
	x_good=(float*)malloc(sizeof(float)*rows+2);
	u=(float**)malloc(sizeof(float*)*rows+2);
	v=(float**)malloc(sizeof(float*)*rows+2);
	y=(float**)malloc(sizeof(float*)*numbins+2);
	channelDPH=(float**)malloc(sizeof(float*)*numbins+2);
	channel=(int*)malloc(sizeof(int*)*numbins+2);
	channelf=(float*)malloc(sizeof(float*)*numbins+2);
	w=(float*)malloc(sizeof(float)*rows+2);
	sig=(float*)malloc(sizeof(float)*rows+2);
	a=(float*)malloc(sizeof(float)*cols+2);
	funy=(float*)malloc(sizeof(float)*rows+2);
	params=(float**)malloc(sizeof(float*)*numsources+2);
	shadow=(float**)malloc(sizeof(float*)*numsources+2);
	shadow_new=(float**)malloc(sizeof(float*)*numsources+2);
	response=(float*)malloc(sizeof(float)*300*512);
	if(energy==NULL ||backinten==NULL||intenRatio==NULL || x==NULL||u==NULL||v==NULL||w==NULL||a==NULL||funy==NULL||params==NULL||shadow==NULL)
	{
			printf("Error(%s:%d): Error while allocating the memory\n",__FILE__,__LINE__);
			exit(0);
	}
	for(i=0;i<numsources;i++)
	{
		params[i]=(float*)malloc(sizeof(float)*numparams);
		shadow[i]=(float*)malloc(sizeof(float)*numpix);
		shadow_new[i]=(float*)malloc(sizeof(float)*numpix);
	}

	for(i=0;i<numbins+1;i++)
	{
		y[i]=(float*)malloc(sizeof(float)*rows+1);
		channelDPH[i]=(float*)malloc(sizeof(float)*rows+1);

	}
	for(i=0;i<=cols;i++)
	{
		intenRatio[i]=(float*)malloc(sizeof(float)*(numbins+1));
		error[i]=(float*)malloc(sizeof(float)*numbins+1);
		
	}
	for(i=0;i<=rows;i++)
	{
		u[i]=(float*)malloc(sizeof(float)*(rows+1));
		v[i]=(float*)malloc(sizeof(float)*(cols+1));
		sig[i]=1;
	}

	for(i=0;i<=rows;i++)
		w[i]=0;
	for(i=0;i<=rows;i++)
		for(j=0;j<=cols;j++)
			v[i][j]=0;
	if(x==NULL || u==NULL || v==NULL||y==NULL||w==NULL||a==NULL||sig==NULL)
	{
		printf("\tError(%s:%d): Unable to allocate memory\n",__FILE__,__LINE__);
		exit(0);
	}

	for(i=0;i<numbins+1;i++)
	{
		for(j=0;j<rows+1;j++)
		{
				y[i][j]=1;
		}
	}


	for(i=0;i<numsources;i++)
	{
		for(j=0;j<numbins;j++)
		{
			error[i][j]=0;
		}
	}
	for(i=0;i<numsources;i++)
	{
		if(i==0)
		{

			//ReadShadow("backshadow.fits",shadow[0],2);
			for (j=0;j<4096;j++)
				shadow[i][j]=0.25;
			
		}
		else
		{
			printf("Getting shadow for %f\t%f\n",tx[i-1],ty[i-1]);
			getShadow(shadow[i],tx[i-1],ty[i-1]);
		}
		char shadowfile[1024];
		sprintf(shadowfile,"shadow_%d",i);
		
	}
	read_lldfile(LLDFILE,0,LLD);
	read_badpixfile(badpixfilename,0,pix_flag);
/******** Start Computing prior******************/


	printf("Computing the prior\n");

	createChannelDPH(infile,channelDPH,minchannel,maxchannel,binsize,LLD,pix_flag);

	int pixflagindex[4096],tmpdetid,tmppixid;
	for (i=0;i<64;i++)//stands for dety
	{
		for(j=0;j<64;j++)//stands for detx
		{
			getDetectorIdPixelNo(j,i,&tmpdetid,&tmppixid);
			pixflagindex[i*64+j]=(tmpdetid*256)+tmppixid;
		}
		
	}
	int ii=0,kk=0,jj=0,rowcount=0,rowcount1=0;

	for(jj=0;jj<64;jj++)//dety
	{

		for(kk=0;kk<64;kk++)//detx
		{
			if(pix_flag[pixflagindex[jj*64+kk]]==0)
			{
				for(ii=0;ii<numsources;ii++)
				{
					shadow_new[ii][rowcount]=shadow[ii][jj*64+kk];
				}
				rowcount++;
			}
		}
	}
 	for(i=0;i<numbins;i++)
	{
		rowcount1=1;
		for(jj=0;jj<64;jj++)//dety
		{

			for(kk=0;kk<64;kk++)//detx
			{

				if(pix_flag[pixflagindex[jj*64+kk]]==0)
				{
					y[i][rowcount1]=channelDPH[i][(jj*64+kk)+1];
					rowcount1++;
				}
			}	
		}

		generateMatrix(tx,ty,numsources,x,sig,u,shadow_new,rowcount);

		svdfit(x,y[i],sig, rowcount,a,cols,u,v,w,&chisq);
		channel[i]=minchannel+(i*binsize);	
		for(j=1;j<=cols;j++)
		{
			if(!isnan(a[j]) )
			{
				intenRatio[j][i+1]=a[j];
				error[j][i+1]=sqrt(a[j]);
				printf("%f\t",a[j]);
			}
			else
			{

				printf("nan\n");
				intenRatio[j][i+1]=0.0;
				error[j][i+1]=0.0;
			}
			printf("\n");
		}
	}
	printf("Writing Prior Spectrum\n");

	write_spectrum_file("prior2.pha",specTemplate,infile,channel,&intenRatio[2][1],&error[2][1],numbins,0,0);
	executeLikelihood(infile,"spectra.txt",params,shadow,numsources,minchannel,maxchannel,intenRatio,LLD,pix_flag,error);

	for(i=0;i<numsources;i++)
	{
		for(j=0;j<numbins;j++)
		{
			error[i][j]=sqrt(intenRatio[i][j]);
		}
	}
	for (i=0;i<numsources;i++)
	{
		bzero(outspecfile,sizeof(outfile));
		sprintf(outspecfile,"%s_%d.pha",outfile,i);
		if(i==1)
		{
			//writing background spectrum, hence tx, ty are set to 0.0
			write_spectrum_file(outspecfile,specTemplate,infile,channel,&intenRatio[i][0],&error[i][0],numbins,0.0,0.0);
		}
		else
		{
			//writing spectrum for source (i-1)
			write_spectrum_file(outspecfile,specTemplate,infile,channel,&intenRatio[i][0],&error[i][0],numbins,tx[i],ty[i]);
		}
		
	}
}
int readResponse(char*responsefile,float*response)
{
	fitsfile*fptr;
	int status,hdutype,hdunum;
	if(fits_open_file(&fptr,responsefile,READONLY,&status) )
	{
		printf("Error while opening file %s\n",responsefile);
		printf("Exiting...\n");
		exit(0);
	}

	if ( fits_movabs_hdu(fptr, hdunum+2, &hdutype, &status) )
		printf("Error while moving hdu\n");
}

void executeLikelihood(char*infile,char*outfile,float **params,float **shadow,int numsources,int minpi,int maxpi, float**outspectra, int * LLD,int*pix_flag,float**error)
{
	fitsfile *fptr;
	float binwidth=1,*area, masktrans=0;
	int status=0,hdunum=0,hdutype=0,nfound=0,anynull=0,toberead=0,numbins,k=0;
	int xcolno=0,ycolno=0,energycolno=0,i=0,chunksize=1024,totevt=0,j=0,iter=0;
	long naxes[2],rv=0, felem=1, nelem=0;
	char strnull[10],plottingCommand[2048],tempfile[1024];
	int *detx=NULL,*dety=NULL,binno=0,numcols=0,intencolno=0,ii=0;
	float*energy=NULL,**prob,**spectra,*sum,**oldparam,*compCount,**oldspectra,*midenergy;
	float maxdiff=0;
	FILE *fp,*fptmp;
	int picolno,*pi,detidcolno,pixidcolno,detpix_index;

	float normfrac[numsources];
	float width=0.0,height=0.0;
	int tmpdetid=0,tmppixid=0;
	int numpix[numsources];
	short *detid,*pixid;
	strcpy(strnull, " ");
	if ( fits_open_file(&fptr, infile, READONLY, &status) )
		printf("Error while opening fits file\n");
	numbins=(int)(maxpi-minpi)/binwidth;

	for(hdunum=0;hdunum<NUMQUADRANT;hdunum++)
	{

		if ( fits_movabs_hdu(fptr, hdunum+2, &hdutype, &status) )
			printf("Error while moving hdu\n");
		if (hdutype != BINARY_TBL)
		{
			printf("Error (%s:%d): %dth extension supposed to be binary table, but it is not. Exiting...\n",__FILE__,__LINE__,hdunum);
			exit(-1);
		}
		if ( fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) )
			printf("Error while getting naxis keyword");
		fits_get_num_cols(fptr,&numcols,&status);
		if(status)
		{
			fits_report_error(stderr,status);
			return ;
		}
		intencolno=numcols+1;
		fits_get_colnum(fptr,CASEINSEN,"DETX",&xcolno,&status);
		if(status)
		{
			fits_report_error(stderr,status);
			return ;
		}
		fits_get_colnum(fptr,CASEINSEN,"DETY",&ycolno,&status);
		if(status)
		{
			fits_report_error(stderr,status);
			return ;
		}
		fits_get_colnum(fptr,CASEINSEN,"PI",&picolno,&status);
		if(status)
		{
			fits_report_error(stderr,status);
			return ;
		}


		fits_get_colnum(fptr,CASEINSEN,"DETID",&detidcolno,&status);
		if(status)
		{
			fits_report_error(stderr,status);
			return ;
		}

		fits_get_colnum(fptr,CASEINSEN,"PIXID",&pixidcolno,&status);
		if(status)
		{
			fits_report_error(stderr,status);
			return ;
		}
		
		totevt=naxes[1];
		prob=(float**)malloc(sizeof(float*)*numsources);
		detx=(int*)malloc(sizeof(int)*totevt);
		dety=(int*)malloc(sizeof(int)*totevt);
		detid=(short*)malloc(sizeof(short)*totevt);
		pixid=(short*)malloc(sizeof(short)*totevt);
		energy=(float*)malloc(sizeof(float)*totevt);
		pi=(int*)malloc(sizeof(pi)*totevt);
		sum=(float*)malloc(sizeof(float)*totevt);
		spectra=(float**)malloc(sizeof(float*)*numsources+1);
		oldspectra=(float**)malloc(sizeof(float*)*numsources+1);
		area=(float*)malloc(sizeof(float)*numsources+1);
		oldparam=(float**)malloc(sizeof(float)*numsources+1);
		compCount=(float*)malloc(sizeof(float)*4096);
		chunksize=1024;
		for(i=0;i<numsources;i++)
		{
			prob[i]=(float*)malloc(sizeof(float)*totevt);
			spectra[i]=(float*)malloc(sizeof(float)*numbins+1);
			oldspectra[i]=(float*)malloc(sizeof(float)*numbins+1);
			area[i]=0;
			oldparam[i]=(float*)malloc(sizeof(float)*4);
		}
		midenergy=(float*)malloc(sizeof(float)*numbins+1);
		toberead=totevt;//(totevt>chunksize)?chunksize:totevt;
		fits_read_col(fptr, TINT, xcolno,1, felem, toberead, strnull, detx,&anynull, &status);
		if(status)
		{
			fits_report_error(stderr,status);
			return ;
		}
		fits_read_col(fptr, TINT, ycolno, 1, felem, toberead, strnull, dety,&anynull, &status);
		if(status)
		{
			fits_report_error(stderr,status);
			return ;
		}
		fits_read_col(fptr, TINT, picolno, 1, felem, toberead, strnull, pi,&anynull, &status);
		if(status)
		{
			fits_report_error(stderr,status);
			return ;
		}
		fits_read_col(fptr, TSHORT, detidcolno, 1, felem, toberead, strnull, detid,&anynull, &status);
		if(status)
		{
			fits_report_error(stderr,status);
			return ;
		}
		fits_read_col(fptr, TSHORT, pixidcolno, 1, felem, toberead, strnull, pixid,&anynull, &status);
		if(status)
		{
			fits_report_error(stderr,status);
			return ;
		}

		for(i=0;i<numsources;i++)
		{
			for(j=0;j<numbins;j++)
			{
				oldspectra[i][j]=0;
				error[i][j]=0.0;
			}
		}
		for(i=0;i<numbins;i++)
			midenergy[i]=minpi+((i+0.5)*binwidth);
		for(j=0;j<numbins;j++)
		{
			for(i=0;i<numsources;i++)
			{
				spectra[i][j]=outspectra[i+1][j+1];
			}
		}
//start of while
		while(iter<1)
		{
			//printf("\n\nProcessing iteration no %d\n",iter);
			for(j=0;j<totevt;j++)
				sum[j]=0;
			//compute the sum of probabilities as base
			for(j=0;j<totevt;j++)
				for(i=0;i<numsources;i++)
				{
					masktrans=0;//getmasktrans(dety[j],detx[j],energy[j]);

					binno=pi[j];//-minpi;

					detpix_index=detid[j]*256+pixid[j];			
					if(pi[j]>minpi && pi[j] <maxpi && pi[j]>LLD[detpix_index] && pix_flag[detpix_index]==0)
					{
						//sum[j]+=spectra[i][binno]*(shadow[i][((dety[j]*64)+detx[j])]);
						sum[j]+=spectra[i][binno]*(shadow[i][((dety[j]*64)+detx[j])]);
					}
				}



			masktrans=0;
			//compute the probabilities
			for(i=0;i<numsources;i++)
			{
				for(j=0;j<totevt;j++)
				{
					if(sum[j]!=0)
					{
						binno=pi[j];//-minpi;
						detpix_index=detid[j]*256+pixid[j];			
						if(pi[j]>minpi && pi[j] <maxpi && pi[j]>LLD[detpix_index] && pix_flag[detpix_index]==0)
						{

							if(sum[j]!=0)
							{
								prob[i][j]=(spectra[i][binno]*shadow[i][((dety[j]*64)+detx[j])])/sum[j];
							}
							else
							{
								prob[i][j]=0.0;			
							}	
						}
					}
				}
				
			}

			//code to bin the spectra
			for(i=0;i<numsources;i++)
				for(j=0;j<numbins;j++)
				{
					spectra[i][j]=0;
				}
			for(i=0;i<numsources;i++)
			{
				for(j=0;j<totevt;j++)
				{
					
					detpix_index=detid[j]*256+pixid[j];			
					
					if(pi[j]>minpi && pi[j] <maxpi && pi[j]>LLD[detpix_index] && pix_flag[detpix_index]==0)
					{
						binno=pi[j];//(int)((pi[j]-minpi)/binwidth);
						spectra[i][binno]+=prob[i][j];
						error[i][binno]+=(prob[i][j]*prob[i][j]);
					}
				}
				
			}


			for(i=0;i<numsources;i++)
			{
				normfrac[i]=0.0;
				numpix[i]=0.0;
					
				for(j=0;j<64;j++)//stands for dety
				{

						if(j%16==0|| j%16==15)
						{
							width=2.28;
						}
						else
						{
							width=2.46;
						}
						for(k=0;k<64;k++)//detx
						{
							if(k%16==0 || k%16==15)
							{
								height=2.28;
							}
							else
							{
								height=2.46;
							}

							getDetectorIdPixelNo(j,k,&tmpdetid,&tmppixid);
							detpix_index=(tmpdetid*256)+tmppixid;
							if(pix_flag[detpix_index]==0)
							{
								normfrac[i]+=(shadow[i][j*64+k]*width*height);
								numpix[i]++;
							}
						}	
				}



				normfrac[i]/=100.0;
				//normfrac[i]+=(normfrac[i]*0.4);

				printf("NormFrac %f\n",normfrac[i]);
				//normfrac[i]=115.553993;

			}	
			/*Area is not matching in common and bayesian*/
			//computing the flux by dividing the total area of czti i.e, 1024 cm^2
			for(j=0;j<numbins;j++)
				for(i=0;i<numsources;i++)
				{
					spectra[i][j]/=(normfrac[i]);
				}
			for(i=0;i<numsources;i++)
			{			

				printf("Sum of shadows for %d source is %f\t%d\n",i,normfrac[i],numpix[i]);
				for(j=0;j<numbins;j++)
				{
					error[i][j]=sqrt(error[i][j]/(normfrac[i]));///normfrac[i];
				}
			}
			iter++;
		}
		//end of while

		for(i=0;i<numsources;i++)
		{
			for(j=0;j<numbins;j++)
			{
				outspectra[i][j]=spectra[i][j];
			}
		}
	
	}

int computeError(float*error,float*spectra,int numbins)
{
	float norm=5.12763,index=-2.04863;
	int i=0;
	for(i=0;i<numbins;i++)
	{
		
	}
}

void ReadShadow(char*filename,float *mask,int hduNo)
{
	fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
    	int status,  nfound, anynull,hdutype,i,j,rsprows=0;
	int *buffer;
    	long naxes[2], fpixel, nbuffer, npixels, ii;
    	float datamin, datamax, nullval;
    	status = 0;
		
    	if ( fits_open_file(&fptr, filename, READONLY, &status) )
	{
    	
		printf("Unable to open background file\n");
		printerror( status );
	}
   	if ( fits_movabs_hdu(fptr, hduNo, &hdutype, &status) ) 
	{
		
		printf("Unable to move to background hdy\n");
   		printerror( status );
	}
    	if ( fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) )
	{
        
		printf("Unable to read keys\n");
		printerror( status );
	 }
    	npixels  = naxes[0] * naxes[1];         /* number of pixels in the image */
    	fpixel   = 1;
    	nullval  = 0;                /* don't check for null values in the image */
	while (npixels > 0)
    	{
      		nbuffer = npixels;
      		if ( fits_read_img(fptr, TFLOAT, fpixel, nbuffer, &nullval,mask, &anynull, &status) )
           		printerror( status ); 
      		npixels -= nbuffer;    
      		fpixel  += nbuffer;    
    	}
    	if ( fits_close_file(fptr, &status) )
	{
        
		printf("Unable to close background file\n");
		printerror( status );
	}
    	return;
}

void fitPowerLaw(float *energy,float *counts,int numbins,float *params)
{
		float *z,*y,sz=0,sy=0,szz=0,syz=0,delta=0;
		int i=0,j=0,n=1,isbreak=0;
		z=(float*)malloc(sizeof(float)*numbins+1);
		y=(float*)malloc(sizeof(float)*numbins+1);
		for(j=0;j<numbins;j++)
		{
			if(counts[j]>0)
			{

				z[n]=log(energy[j]);
				y[n]=log(counts[j]);
				n++;
			}
		}

		n-=1;
		sz=0;
		sy=0;
		szz=0;
		syz=0;
		for(j=1;j<n+1;j++)
		{
			sz+=z[j];
			sy+=y[j];

			szz+=(z[j]*z[j]);
			syz+=(y[j]*z[j]);
		}
		delta=(n*szz)-(sz*sz);
		params[0]=exp(((szz*sy)-(sz*syz))/delta);
		params[1]=-((n*syz)-(sy*sz))/delta;

		free(z);
		free(y);
}


int write_spectrum_file(char*specFilename, char*specTemplate,char* infile,int*pi,float* flux,float* error,int nbins,double tx,double ty) 
{
    int status = 0;
    int colnum = 0;
    long nrows = 0;
    fitsfile *fspec;

    int hdunum=0;
    int hdutype=0;
    fitsfile *fevt;

    hdunum=2;//spec.quadrantid+2;
    fits_open_file(&fevt,infile, READONLY, &status);
   
    if (create_empty_fitsfile(specFilename, specTemplate)) {
        printf("Error in creating spectrum file from corresponding template");
    }

    fits_open_file(&fspec, (char*) specFilename, READWRITE, &status);
    if (status) {
        fits_report_error(stderr, status);
    }


   if ( fits_movabs_hdu(fspec, 1, &hdutype, &status) )
   {
        fits_report_error(stderr, status);
    }

    copyheader(fevt ,fspec );


   if ( fits_movabs_hdu(fspec, 2, &hdutype, &status) )
   {
        fits_report_error(stderr, status);
    }

    copyheader(fevt ,fspec );
   

    //Writing PI channel
    if ( fits_write_col(fspec, TINT, 1, 1, 1, nbins, pi,&status) ) 
    {
        printf("Error in writing column CHANNEL of spectrum file %s" ,specFilename);
    }
    //Writing flux
    if ( fits_write_col(fspec, TFLOAT, 2, 1, 1, nbins, flux,&status) ) 
    {

        printf("Error in writing column COUNTS of spectrum file %s" ,specFilename);
    }
    //Writing Error
    if ( fits_write_col(fspec, TFLOAT, 3, 1, 1, nbins, error,&status) ) 
    {

        printf("Error in writing column ERR of spectrum file %s" ,specFilename);
    }
    fits_close_file(fspec, &status);
    if (status) {
        fits_report_error(stderr, status);
    }

return 0;
}

void copyheader(fitsfile *fptr,fitsfile *outfile)
{
	int status=0;
	long tstarti=0,tstopi=0;
	double tstart=0,tstop=0,exposure=0;
	char gtitype[100];
/*	
	fits_read_key(fptr,TLONG,"STARTI",&tstarti,NULL, &status);
        if (status) 
	{
		printf("Error while readiing tstarti key\n");
		status=0;
	}
	fits_read_key(fptr,TLONG,"STOPI",&tstopi,NULL, &status);
        if (status) 
	{
		printf("Error while readiing tstopi key\n");
		status=0;
	}
*/
	fits_read_key(fptr,TDOUBLE,"TSTART",&tstart,NULL, &status);

        if (status) 
	{
		printf("Error while readiing tstart key\n");

		status=0;
	}
	fits_read_key(fptr,TDOUBLE,"TSTOP",&tstop,NULL, &status);
        if (status) 
	{
		printf("Error while readiing tstop key\n");
		status=0;
	}

	fits_read_key(fptr,TDOUBLE,"EXPOSURE",&exposure,NULL, &status);

        if (status) 
	{
		printf("Error while readiing exposure key\n");
		status=0;
	}
	fits_read_key(fptr,TSTRING,"GTITYPE",&gtitype,NULL, &status);

        if (status) 
	{
		printf("Error while readiing gtitype key\n");
		status=0;
	}
	fits_update_key(outfile, TDOUBLE, "TSTART",&tstart, "", &status);

        if (status) 
	{
		printf("error while writing tstart key\n");
		status=0;
	}
	fits_update_key(outfile, TDOUBLE, "TSTOP",&tstop, "", &status);

        if (status) 
	{
		printf("error while writing tstop key\n");
		status=0;
	}
	fits_update_key(outfile, TDOUBLE, "EXPOSURE",&exposure, "", &status);

        if (status) 
	{
		printf("error while writing exposure key\n");
		status=0;
	}
	fits_update_key(outfile, TSTRING, "GTITYPE",gtitype, "", &status);

        if (status) 
	{
		printf("error while writing gtitype key\n");
		status=0;
	}

	
}
/*This function calls generate shadow function and creates a matrix u which is used for svd fit*/
void generateMatrix(float*tx,float*ty,int numsources,float *x,float*sig,float **u ,float **shadow,int numrows)
{
	float tmp=0;
	int j=0,i=0;
	int counter=0;
	int rows=NUMELEM*NUMQUADRANT;
	float **srcshadow,*backshadow;
	srcshadow=(float**)malloc(sizeof(float**)*numsources+1);
	backshadow=(float*)malloc(sizeof(float*)*numsources+1);

	for(i=0;i<numsources;i++)
	{
		srcshadow[i]=(float*)malloc(sizeof(float)*NUMELEM*NUMQUADRANT);
	}



	for (j=1;j<=numsources;j++)
	{

			for (i=1;i<=numrows;i++)
			{
				u[i][j]=shadow[counter][i-1];
			}
			counter++;
	}
}


/*
This functions creates DPH (histogram of counts in each pixel) for every channel. 
*/
void createChannelDPH(char*infile,float **count,int minpi,int maxpi, int binsize,int*LLD,int*pixflag)
{
	fitsfile *fptr,*badfptr;
	int status=0,hdunum=0,hdutype=0,nfound=0,anynull=0,toberead=0;
	int detidcolno=0,pixidcolno=0;
	int xcolno=0,ycolno=0,picolno=0,i=0,chunksize=1024,totevt=0,j=0,detpix_index=0;
	long naxes[2],rv=0, felem=1, nelem=0;
	char strnull[10];
	int *detx,*dety,tmp=0,numbins=0,binno=0;
	int*pi;
	short *detid,*pixid;
//	int maxpi=512,minpi=0,binsize=1;
	numbins=(maxpi-minpi)/binsize;
	strcpy(strnull, " ");

	if ( fits_open_file(&fptr, infile, READONLY, &status) )
		printerror( status );



	detx=(int*)malloc(sizeof(int)*chunksize);
	dety=(int*)malloc(sizeof(int)*chunksize);
	pi=(int*)malloc(sizeof(float)*chunksize);
	detid=(short*)malloc(sizeof(short)*chunksize);
	pixid=(short*)malloc(sizeof(short)*chunksize);

	for(hdunum=0;hdunum<NUMQUADRANT;hdunum++)
	{

		if ( fits_movabs_hdu(fptr, hdunum+2, &hdutype, &status) )
				printerror( status );
		if (hdutype != BINARY_TBL)
		{
			printf("Error (%s:%d): %dth extension supposed to be binary table, but it is not. Exiting...\n",__FILE__,__LINE__,hdunum);
			exit(-1);
		}
		if ( fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) )
				printerror( status );
		/*Get column no for detx*/
		fits_get_colnum(fptr,CASEINSEN,"DETX",&xcolno,&status);
		if(status)
		{
			fits_report_error(stderr,status);
			return ;
		}

		/*Get column no for dety*/
		fits_get_colnum(fptr,CASEINSEN,"DETY",&ycolno,&status);
		if(status)
		{
			fits_report_error(stderr,status);
			return ;
		}

		/*Get column no for PI*/
		fits_get_colnum(fptr,CASEINSEN,"PI",&picolno,&status);
		if(status)
		{
			printf("Inside get energy column\n");
			fits_report_error(stderr,status);
			return ;
		}

		/*Get column no for detid*/
		fits_get_colnum(fptr,CASEINSEN,"DETID",&detidcolno,&status);
		if(status)
		{
			printf("Inside get energy column\n");
			fits_report_error(stderr,status);
			return ;
		}



		/*Get column no for pixid*/
		fits_get_colnum(fptr,CASEINSEN,"PIXID",&pixidcolno,&status);
		if(status)
		{
			printf("Inside get energy column\n");
			fits_report_error(stderr,status);
			return ;
		}


		/*
		naxes[1] stores the number of rows in a fits file. In this case total no of rows
		in a quadrant corrosponds to total no of recorded events.	
		*/
		totevt=naxes[1];
		/*
		Inoder to avoid memory related issue, complete file has not been read in memory. Instead we are
		reading few rows (bufisize) and processing them. One having machine with higher ram can change this value which 
		will help to improve performance of the routine.
		*/
		chunksize=1024;
		tmp=0;
		/*
		count array stores the channel wise dph; initilizing counts to zero.
		*/	

		for(i=0;i<numbins;i++)
				for(j=0;j<4096*NUMQUADRANT;j++)
					count[i][j]=0;




		for(i=0;i<naxes[1];i+=chunksize)
		{
			/*
			calculting how many rows to read during the current pass. For all passes except last pass it will be bufsize. But
			for the last pass it will read only remaining events.
			*/
			toberead=(totevt>chunksize)?chunksize:totevt;
			/*
			Reading detx column and storing it in a variable detx
			*/
			fits_read_col(fptr, TINT, xcolno,i+1, felem, toberead, strnull, detx,&anynull, &status);
			if(status)
			{
				fits_report_error(stderr,status);
				return ;
			}

			/*
			Reading dety column and storing it in a variable dety
			*/
			fits_read_col(fptr, TINT, ycolno, (i+1), felem, toberead, strnull, dety,&anynull, &status);
			if(status)
			{
				fits_report_error(stderr,status);
				return ;
			}

			/*
			Reading pi column and storing it in a variable pi
			*/
			fits_read_col(fptr, TINT, picolno, (i+1), felem, toberead, strnull, pi,&anynull, &status);
			if(status)
			{
				fits_report_error(stderr,status);
				return ;
			}

			/*
			Reading detid column and storing it in a variable detid
			*/
			fits_read_col(fptr, TSHORT, detidcolno, (i+1), felem, toberead, strnull, detid,&anynull, &status);
			if(status)
			{
				printf("Error while readind detid colno\n");
				fits_report_error(stderr,status);
				return ;
			}


			/*
			Reading pixid column and storing it in a variable pixid
			*/
			fits_read_col(fptr, TSHORT, pixidcolno, (i+1), felem, toberead, strnull, pixid,&anynull, &status);
			if(status)
			{
				printf("Error while readind pixid colno\n");
				
				fits_report_error(stderr,status);
				
				exit(0) ;
			}





			/* Keeping track of number of events yet to be read. After reading complete file, totevt will be zero
			*/
			totevt-=toberead;

//printf("Okay %d\n",__LINE__);
			for(j=0;j<toberead;j++)
			{
				detpix_index=detid[j]*256+pixid[j];			
				/* Chacking range of the PI value*/
				if(pi[j]<maxpi && pi[j]>=minpi && pi[j]>LLD[detpix_index] && pixflag[detpix_index]==0 )
				{
					/*
					Incrementing count of pixel corrosponding to the current pi channel
					*/
					/*if (pi[j]>=100 && pi[j] <=120)
					{
					}
					else
					{*/
						binno=(pi[j]-minpi)/binsize;	
						count[binno][(hdunum*NUMELEM)+(dety[j]*ELEMPERQUAD)+detx[j]+1]++;
					//}
				}
				else
				{
					//printf("Error(%s:%d) Unexpected pi value: %d\n. Terminating the executio...",__FILE__,__LINE__,pi[j]);		
					//exit(0);
				}

			}
		}
	}
	if ( fits_close_file(fptr, &status) )
		printerror( status );
	free(detx);
	free(dety);
	free(pi);
}


int create_empty_fitsfile(char* outFilename, char* outTemplate)
{
    int status=0; //status variable
    fitsfile *fptr;
   remove(outFilename); 
    fits_create_template(&fptr, (char*) outFilename, (char*) outTemplate, &status);
    if(status){
        fits_report_error(stderr, status);
    }
    
    fits_close_file(fptr, &status);
    if (status) {
        fits_report_error(stderr, status);
    }
return 0;
}


void writeFitsDPH(char * outputfilename,int qid,float* count)
{
	long naxes[2] = { 64,64};
	long nelements = naxes[0] * naxes[1];
	int bitpix   =  FLOAT_IMG;
	long naxis    =   2;
	int fpixel = 1,status=0;
	float energy=10;
	char ext[100]="Q",temp[100];
	sprintf(temp,"%d",qid);
	strcat(ext,temp);
	int i=0,j=0;
	fitsfile*write_fptr;
//	if(qid==0)
	{
		remove(outputfilename);
		if (fits_create_file(&write_fptr,outputfilename, &status))
		{
			printf("\tError(%s:%d):\n",__FILE__,__LINE__);
			printerror( status );
		}
		if ( fits_create_img(write_fptr,  bitpix, 0, naxes, &status) )
		{
			printf("\tError(%s:%d):\n",__FILE__,__LINE__);
			printerror( status );
		}
	}

	if ( fits_create_img(write_fptr,  bitpix, naxis, naxes, &status) )
	{
		printf("\tError(%s:%d):\n",__FILE__,__LINE__);
		printerror( status );
	}

	if ( fits_write_img(write_fptr, TFLOAT, fpixel, nelements,count, &status) )
	{
		printf("\tError(%s:%d):\n",__FILE__,__LINE__);
		printerror( status );
	}
	fits_update_key(write_fptr, TSTRING, "EXTNAME", ext,"", &status) ;
//	if(qid==3)
	{
		if ( fits_close_file(write_fptr, &status) )
		{
			printf("\tError(%s:%d):\n",__FILE__,__LINE__);
			printerror( status );
		}
	}

}




int read_lldfile(char* caldb_lld,int qid,int *lld)
{
    fitsfile *fptr;
    int status=0,hdutype=0;
    int lld_col,i;
    long nrows;

    fits_open_file(&fptr, (char *)caldb_lld, READONLY, &status);
    if (status)
    {
    	printf("Error (%s:%d): Error while opening %s file\n",__FILE__,__LINE__,caldb_lld);
    	fits_report_error(stderr, status);
    	return (EXIT_FAILURE);
    }

    fits_movabs_hdu(fptr, qid+2, &hdutype, &status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    fits_get_num_rows(fptr, &nrows, &status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    fits_get_colnum(fptr,CASEINSEN,"LLD",&lld_col,&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    for(i=0;i<4096;i++)
    {

        fits_read_col(fptr, TINT, lld_col, i+1, 1, 1, NULL, &lld[i],NULL, &status);
        if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    }

    fits_close_file(fptr,&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    return(EXIT_SUCCESS);

}

void read_badpixfile(char*badpixfilename,int quadid, int*pixflag)
{
	fitsfile*badfptr;
	int pixflagcolno;
	   int status=0,hdutype=0;

	if ( fits_open_file(&badfptr, badpixfilename, READONLY, &status) )
	printerror( status );
	
	if ( fits_movabs_hdu(badfptr, quadid+2, &hdutype, &status) )
		printerror( status );

		

		printf("Pixel flag read\n");
	/*Get column no for PIX_FLAG from badpixel file*/
	fits_get_colnum(badfptr,CASEINSEN,"PIX_FLAG",&pixflagcolno,&status);
	if(status)
	{
		printf("Unable to read PIX_FLAG column from the badpixel file\n");
		fits_report_error(stderr,status);
		return ;
	}


	/*
	Reading pixidflag column from badpix file and storing it in a variable pixflag
	*/
	fits_read_col(badfptr, TINT, pixflagcolno, 1, 1, 4096, NULL, pixflag,NULL, &status);
	if(status)
	{
		printf("Error (%s:%d) while reading pix flag column from badpixel file\n",__FILE__,__LINE__);
		fits_report_error(stderr,status);
			exit(0) ;
	}
    	fits_close_file(badfptr,&status);
	int i=0,tmppix=0;
	for(i=0;i<4096;i++)
		if(pixflag[i]==0)
			tmppix++;

}
void getDetectorIdPixelNo(int x,int y,int *moduleNo,int *pixelNo)
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


}

