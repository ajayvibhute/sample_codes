/************************************************************************************************************************
 Simultaneous  Linear Least Square Fit via Singular Value Decomposition for Multiple Data sets

Technologies:C,CUDA.
Author:Ajay Vibhute
Project Guide:Prof.Dipankar Bhattacharya.
Project For: IUCAA
*************************************************************************************************************************/
//Header Files.
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<cuda.h>
#include<math.h>
#include<pthread.h>
#include<time.h>
#include<unistd.h>
#include "PreProcessors.h"
#include"generate_xy.h"
#include"svd_cuda.c"
extern "C"
void startExecution(char *observationId); 
#define NUM_THREADS 3

#define print(format,...) fprintf(format,__VA_ARGS__)

//function to initilise the cuda environment and map GPU to a MPI process
void cudaInit(int myrank)
{

	float *initmalloc=NULL;

	cudaSetDevice(myrank);
	cudaMalloc(&initmalloc,sizeof(float));
}

//This function reads the data and performs the forward fitting for a given sky range
void*svdfit(void * data)
{
	//creating object to store the source information 
	InfoData * info=(InfoData*)data;

	//creating local variable and initilising them

	char filename[BUFSIZE],command[BUFSIZE],observation_id[BUFSIZE];	
	char string[100];
	char fileName[BUFSIZE],xtemp[BUFSIZE],ytemp[BUFSIZE],stemp[BUFSIZE];

	int myrank=info->myrank;
	int xstart1=0,xend1=0,ystart1=0,yend1=0;
	int intensityCounter=0;
	int rv=0,jj;
	int trueCondition=1;
	int index=0,updateCounter=0;
	int threadId=myrank;
	int counter=0, no_another=0,*maskpat=NULL;
	int noOfGrids=noOfDataSet;
	int noOfbins=rows;
	int ielem=0, ipat=0;	
	int icam=myrank;//threadId;

	float *x=NULL,*a=NULL,*chisq=NULL,*x_another=NULL,*x_updated=NULL,*x_passing=NULL;
	float intensity=1,thetaX_Updated=0,thetaY_Updated=0.0;
	float *thetaX=NULL,*thetaY=NULL;
	float xdiff1=0.0,ydiff1=0.0;
	float tolerance=0.001;
	float xCheck=0.0,yCheck=0.0,tempChi=0.0;
	float timeSpent=0;

	FILE  *print_pointer=NULL;
	FILE *fp=NULL,*read_observed_shadow=NULL;
	FILE *pipe=NULL; 


	strcpy(observation_id,info->observationId);
	if(myrank<3)
	{
		cudaInit(myrank);

		//checking mode of the excution
		#ifdef __SIM__
	
			print_pointer=stdout;
	
		#else
			//populating the directory structure
			strcpy(filename,OBSERVED_SHADOW_PATH);
			strcat(filename,observation_id);
			strcat(filename,"/output");
			strcat(filename,"/");

			//populating the log file name. Log of individual camera will be stored in the different file
			if(myrank==0)
				strcat(filename,"BOOM_Execution_Log");
	
			else if(myrank==1)
				strcat(filename,"SLANT1_Execution_Log");
			else
				strcat(filename,"SLANT2_Execution_Log");

			//creating log file
			print_pointer=fopen(filename,"w");
	
			if(print_pointer==NULL)
			{
				print(stderr,"%s:%d:Error while creating output file %s \n",__FILE__,__LINE__,filename);
				printf("%s:%d:Error while creating output file %s \n",__FILE__,__LINE__,filename);	
				return (void*)-1;
			}
		#endif

		//initilizing the grid information
		xstart1 =info->xs;             //starting of thetaX
		xend1 =info->xe;              //end of the thetax
		ystart1=info->ys;            //starting of thetaY
		yend1 =info->ye;  
		xdiff1=xdiff;	
		ydiff1=ydiff;

		//getting size of the sky grid
		xgridSize=(int)((xend1-xstart1)/xdiff1);
		ygridSize=(int) ((yend1-ystart1)/ydiff1);
		ygridSize++;
		if(xdiff1>0.1)
			xgridSize++;
		noOfDataSet=xgridSize*ygridSize;
		thetaX=(float*)malloc(sizeof(float)*(xgridSize+1)*(ygridSize+1));
		thetaY=(float*)malloc(sizeof(float)*(xgridSize+1)*(ygridSize+1));
		//Running a loop over the grid and initilizing the pointing information
		for(float ty=ystart1;ty<=yend1;ty+=ydiff1)
		{	
			for(float tx=xstart1;tx<=xend1;tx+=xdiff1)
			{	
				thetaX[counter]=tx;
				thetaY[counter++]=ty;
			}
		}
		counter=0;

		no_another=0;	

		size_t size=sizeof(float)*(noOfbins+2)*noOfGrids;

		/* Allocaing host (CPU) memory  */
		x=(float*)malloc(size);
		x_updated=(float*)malloc(size);
		x_passing=(float*)malloc(size);
		a=(float*)malloc(size);
		chisq=(float*)malloc(size);
		x_another=(float*)malloc(size*no_another);
		maskpat=(int*)malloc(NUM_PATTERNS*(NUM_MASK_ELEM+2)*sizeof(int));

		//Reading Maskpattern from the file ./config/MASKPATTERNS/....
		for(ipat = 0; ipat <= NUM_PATTERNS-1; ++ipat)
		{
			sprintf(string, MASK_PATH_FORMAT, ipat+1);
		
			if((fp = fopen(string, "r")) == NULL)
			{
				print(stderr,"\nError: (%s:%d) Mask Pattern file \"%s\" could not be opened for reading\n", __FILE__, 	__LINE__, string);
			}
	
			//setting up the mask pattern
			maskpat[ipat*(NUM_MASK_ELEM+2)+0] = 0; /* pad closed elem in the beginning */
			maskpat[ipat*(NUM_MASK_ELEM+2)+(NUM_MASK_ELEM+1)] = 0;     /* and at the end */
		
			//reading the mask elements	
			for(ielem = 1; ielem <= NUM_MASK_ELEM; ++ielem)
				fscanf(fp, "%1d", &maskpat[ipat*(NUM_MASK_ELEM+2)+ielem]);
	
			//closing the file
			fclose(fp);
		}
		/*
		Generating the shadow for the given location with the given intensity
		when flag is set to 0 then new shadow is generated and return (void*)-1;
		when flag is set to 1 then new shadow is genereated and added to the previous shadow which is 
		in the buffer.
		
		If 5 or 6 info->sources out of 10 are at the edge then it will give nearly correct location of the all info->sources 
		if number of info->sources at the edge is increased beyound that then it will not work properly i.e., it will not 
                give the correct strength and correct location of the info->sources
		If info->sources are not at the edge of the grid and maximum source strength is 2000 and minimum source strength
 		is 10 then also it's able to detect all the info->sources with it's relative strength
		*/


	
		sprintf(fileName,OBSERVED_SHADOW_PATH);
		strcat(fileName,observation_id);
		//One telescope is assigned to one node and each node will run the processing on 3 GPUs
		if(myrank==0)
			strcat(fileName,"/BOOM");
		else if(myrank==1)
			strcat(fileName,"/SLANT1");
		else 
			strcat(fileName,"/SLANT2");

		//reading observation data/shadows
		read_observed_shadow=fopen(fileName,"r");
		if(read_observed_shadow==NULL)
		{
			print(print_pointer,"%s (%d):%s file not found\n",__FILE__ ,__LINE__,fileName);
			fclose(print_pointer);
			return (void*)-1;;
		}
		for( int k=1;k<=noOfbins;k++)
		{
			rv=fscanf(read_observed_shadow,"%f",&x[k]);
			if(rv==-1)
			{
				print(print_pointer,"%s (%d):%s file have insufficient data to proceed\n",__FILE__ ,__LINE__,fileName);
				fclose(print_pointer);
				return (void*)-1;;
			
			}
		}
		fclose(read_observed_shadow);
		

	
/************************************************************************************************************************
	Pass User defined Function of X
	x-another contains the values for another function of x noOfbins wise not column wise.
	If input for the x_another is read from the file then read all values of first function for all column then all 	values of the second function and so on... if you pass values column wise then it will give error or the output 	values would not be correct
************************************************************************************************************************/
	/*
	for(k=0;k<noOfGrids;k++)
	{
		for(i=1;i<=noOfbins;i++)
		{
			temp=no_another;	
			for(int l=0;l<no_another;l++)
			{
				x_another[counter++]=powf(x[k*noOfbins+i],colms-temp);	//change the counter++
				temp--;
			}
			
		}
		
		
	}
	*/
/************************************************************************************************************************

	Call SVD_FIT to find coordinates,chisquare and cvm matrix
	It provide two overloaded function one with user defined function and one without user defined function.
	1)timeSpent=SVD_FIT(x,y,noOfbins,colms,sig,a,chisq,cvm,noOfGrids,no_another,x_another);
	2)timeSpent=SVD_FIT(x,y,noOfbins,colms,sig,a,chisq,cvm,noOfGrids);
	Function will return (void*)-1; the time spent on GPU for exection 
*************************************************************************************************************************/

	info->numSources=0;

	//running a loop till we discover a source in a field of view. Sometimes a source may take multiple
	//iterations to get detected. In such cases windowing is performed and source intensity is populated accordingly.
	while(trueCondition)
	{
	
		updateCounter++;
		intensityCounter=0;
		for(int kk=1;kk<noOfbins;kk++)
		{
			x_passing[kk]=0;
			x_passing[kk]=x[kk]-(x_updated[kk]);	
		}
		//performing fitting on thee GPU
		timeSpent=SVD_FIT(x_passing,thetaX,thetaY,a,chisq,maskpat,icam);

		bzero(fileName,sizeof(fileName));
		strcpy(fileName,"../output/");
		sprintf(string,"%ld",threadId+1);
		strcat(fileName,string);
		strcat(fileName,"_output.txt");
		index=0;
		tempChi=chisq[0];	
		counter=1;
	
		//checking chisquare value of each grid point and computing lowest chisquare value
		for(int i=0;i<noOfGrids;i++)
		{
			tempChi=chisq[i];
			index=i;
			for(int j=0;j<noOfGrids;j++)
			{
				if(tempChi>chisq[j])
				{
					tempChi=chisq[j];
					index=j;	
					
				}
			
			}
			if(intensityCounter==0)
			{
				intensity=a[index*colms+2];
				thetaX_Updated=thetaX[index];
				thetaY_Updated=thetaY[index];
				intensityCounter++;
			
			}
		chisq[index]=10000;
		
		}
	
				
		//updating the source information
		if(info->numSources==0)
		{
			info->sources[0].x_position=thetaX_Updated;
			info->sources[0].y_position=thetaY_Updated;
			info->sources[0].sourceStrength=intensity;
			info->numSources++;
		}	
		else
		{
		

			
			for(jj=0;jj<info->numSources;jj++)
			{
				xCheck=info->sources[jj].x_position-thetaX_Updated;
				yCheck=info->sources[jj].y_position-thetaY_Updated;
				if(xCheck<0)
					xCheck*=-1;
				if(yCheck<0)
					yCheck*=-1;
				//updating source intensity if source is already listed
				if(xCheck<=(0.5*xdiff1) && yCheck<=(2.5) )	
				{
					info->sources[jj].sourceStrength+=intensity;	
					break;
				}
			}
			if(jj==info->numSources)
			{
				if(intensity>tolerance)
				{	
					if(intensity>0.1)
					{
						info->sources[info->numSources].x_position=thetaX_Updated;
						info->sources[info->numSources].y_position=thetaY_Updated;
						info->sources[info->numSources].sourceStrength=intensity;
						info->numSources++;
					}
					else	
					{
						trueCondition=0;
					}
				
				}
				else
				{	
					trueCondition=0;
				}
		
			}
		}
		//generating shadow for updated location	
		generateShadow(thetaX_Updated,thetaY_Updated,intensity,1,x_updated,icam);
		print(print_pointer,"...................%d th pass rank %d ..\n",updateCounter,myrank);
		print(print_pointer,"\tThetaX=%f\tThetaY=%f\n",thetaX_Updated,thetaY_Updated);
		print(print_pointer,"\tSource Strength=%f\n",intensity);

		}//end of updateCounter loop
		print(print_pointer,"\n-----------------------Processing Complete by thread %d -------------------\n",myrank);
	
		//Printing the results
		print(print_pointer,"NoOfGrids processed=%d\n",noOfGrids);
		print(print_pointer,"X Resolution=%f\tY Resolution=%f\n",xdiff1,ydiff1);
		print(print_pointer,"No of sources Found =%d\n",info->numSources);
		print(print_pointer,"sources with their position and strength\n\n");
		print(print_pointer,"\tThetaX\t\tThetaY\t\tsourcestrength\n\n");
		for(int kk=0;kk<info->numSources;kk++)
		{
			print(print_pointer,"\t%f \t %f \t %f\n",info->sources[kk].x_position,info->sources[kk].y_position,info->sources[kk].sourceStrength);
		}
	
		strcpy(filename,OBSERVED_SHADOW_PATH);
		strcat(filename,observation_id);
		strcat(filename,"/output");
		strcat(filename,"/");	
		if(myrank==0)
			strcat(filename,"BOOM_SOURCE_LIST");
		else if(myrank==1)
			strcat(filename,"SLANT1_SOURCE_LIST");
		else
		strcat(filename,"SLANT2_SOURCE_LIST");
			
		//creating a file and writing source information
		fp=fopen(filename,"w");
		if(fp==NULL)
		{
			print(print_pointer,"%s (%d):Error while creating file %s \n",__FILE__ ,__LINE__,fileName);
			return (void*)-1;;
		}
		for(int i=0;i<info->numSources;i++)
		{
			fprintf(fp,"\t%.2f",info->sources[i].x_position);
			fprintf(fp,"\t\t%.2f",info->sources[i].y_position);
			fprintf(fp,"\t\t%.2f\n",info->sources[i].sourceStrength);
		}
		fclose(fp);
	
		//Generating a plot indicating the source location in the field of view using GNU plot	

		pipe= popen("gnuplot -persist","w");
		strcpy(command,"set title 'Located sources By ");
		if(myrank==0)
			strcat(command," Boom Camera '\n");
		else if(myrank==1)
			strcat(command," SLANT1 Camera '\n");
		else
			strcat(command," SLANT2 Camera '\n");

		fprintf(pipe,command);
		fprintf(pipe, "set xlabel 'X'\n");
		fprintf(pipe, "set ylabel 'Y'\n");
		for(int i=0;i<info->numSources;i++)
		{
			strcpy(command,"set obj ");
			sprintf(xtemp,"%d",i+1);
			strcat(command,xtemp);
			strcat(command," circle  center ");
			sprintf(xtemp,"%.2f",info->sources[i].x_position);
			strcat(command,xtemp);
			strcat(command,",");
			sprintf(ytemp,"%.2f",info->sources[i].y_position);
			strcat(command,ytemp);
			strcat(command," radius 0.1 fill solid fc rgbcolor 'red'\n");
			fprintf(pipe,command);
			sprintf(stemp,"%0.2f",info->sources[i].sourceStrength);
			bzero(command,sizeof(command));
			strcpy(command,"set label ' ( ");
			strcat(command,xtemp);
			strcat(command,",");
			strcat(command,ytemp);
			strcat(command,",");
			strcat(command,stemp);
			strcat(command,")' at ");
			strcat(command,xtemp);
			strcat(command,",");
			strcat(command,ytemp);
			fprintf(pipe,"%s\n",command);
		}
		fprintf(pipe, "plot -40\n");
		fclose(pipe);

#ifdef __SIM__
#else
	fclose(print_pointer);
#endif

	//Free allocated memory
	free(x);
	free(chisq);	
	free(x_another);
	free(thetaX);
	free(thetaY);
	free(maskpat);
	freeDevice();	
		
}	
	return (void*) 0;
}//main

void startExecution(char *observationId)
{
	char hostname[100];
	gethostname(hostname,sizeof(hostname));
	printf("Executing %s on  node %s\n",observationId,hostname);	
	InfoData info[NUM_THREADS+1];	
	
	time_t t1,t2;	
	(void) time(&t1);
	pthread_t  thread[NUM_THREADS];
	for(int i=0;i<NUM_THREADS;i++)
	{
		info[i].myrank=i;
		info[i].numSources=i;
		if(i==0)
		{
			info[i].xs=-10;
			info[i].xe=10;
		}
		else
		{
			info[i].xs=-13;
			info[i].xe=13;
		}
		info[i].ys=-45;
		info[i].ye=45;
		strcpy(info[i].observationId,observationId);
		//printf("ystart=%d\tyend=%d\n",info[i].ys,info[i].ye);
		pthread_create(&thread[i],NULL,svdfit,&info[i]);
		sleep(5);
		
	}
	for(int i=0;i<NUM_THREADS;i++)
	{
		pthread_join(thread[i],NULL);
	}
	(void) time(&t2);
	
	//printf("Total time spent by all threads is=%f\n",(float)(t2-t1));
	printf("Execution of %s is complete on node %s\n",observationId,hostname);		
	return ;
	//pthread_exit(NULL);
	//printf("Exit");
	
}
