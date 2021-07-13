/*
 * mcmc_simulation.c
 *
 *  Created on: Sept 29, 2018
 *  Author: Ajay Vibhute (IUCAA)
 *  Ver1.0 
 *  A client server model to run multiple  mcmc simulation to get error function for bayesian spectrum reconstruction algorithm.
 *  Ver1.1, Ajay Dec 20, 2018
 *  Modified code to work with multiple sources given list of sources in an
 *  input file. Modified the argument list
 */
#include<mpi.h>
#include <unistd.h>
#include<fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include<time.h>
#include<unistd.h>
#include<stdlib.h>
#define BUFFSIZE 1024
#define FILEPATHSIZE 4096
int main(int argc,char*argv[])
{




	int rv=0,i=0,isfirst=1,freeThreadId=0,*jobmap=NULL;
	char command[FILEPATHSIZE],temp[BUFFSIZE];
	FILE *fp;
	int rank,size,counter=0,sampleid=0;
	int numsamples=0;

	MPI_Status stat;
	MPI_Init(&argc,&argv);
	MPI_Comm_size( MPI_COMM_WORLD, &size );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	
/*	Simulation configuration: 
 *	number of sources in the fov
 *	source pointing for each source: thetax thetay
 *	source properties: powerlaw_index powerlaw_normalization
 *
 * */	


	double exposure=80000;
	char responsefile[BUFFSIZE],badpixfile[BUFFSIZE],outfile[BUFFSIZE],tmpoutfile[BUFFSIZE],backgroundfile[FILEPATHSIZE],sourcelistfile[BUFFSIZE];
	char outphafile[BUFFSIZE],outputfile[BUFFSIZE],inlistfile[BUFFSIZE];
	strcpy(outfile,"simulation");

	if(argc<5)
	{
		printf("Please enter commandline arguments\n1)sourcelist file 2) num samples 3) background event file 4)badpixfile 5)outputfile\n");
		exit(0);
	}
	


	strcpy(inlistfile,argv[1]);
	numsamples=atoi(argv[2]);
	strcpy(backgroundfile,argv[3]);
	strcpy(badpixfile,argv[4]);
	strcpy(outputfile,argv[5]);
	


	/*	
 *	Reading the input file containing the information about the sources to
 *	be simulated. The input file should be in following format
 *	Number_of_sources
 *	ThetaX	ThetaY	PowerLawIndex	Normalization	SourceName RA Dec
 *	ThetaX	ThetaY	PowerLawIndex	Normalization	SourceName RA Dec
 *	.
 *	.
 * */
	fp=fopen(inlistfile,"r");
	if(fp==NULL)
	{
		printf("Error(%s:%d): %s file not found\nExiting...\n",__FILE__,__LINE__,inlistfile);
	}	

	int numsources=0;
	fscanf(fp,"%d",&numsources);

	printf("Number of sources %d\n",numsources);
	double tx[numsources+1],ty[numsources+1],photindex[numsources+1],norm[numsources+1],ra[numsources+1],dec[numsources+1];
	char sourcename[numsources+1];
	for(i=0;i<numsources;i++)
	{

		fscanf(fp,"%lf",&tx[i]);
		fscanf(fp,"%lf",&ty[i]);
		fscanf(fp,"%lf",&photindex[i]);
		fscanf(fp,"%lf",&norm[i]);
		//fscanf(fp,"%s",&sourcename[i]);
		fscanf(fp,"%lf",&ra[i]);
		fscanf(fp,"%lf",&dec[i]);
		printf("%f\t%f\t%f\t%f\t%s\t%f\t%f\n",tx[i],ty[i],photindex[i],norm[i],sourcename[i],ra[i],dec[i]);
	}
	fclose(fp);

/*	Creating sourcelist file for all processes
 * */
	bzero(sourcelistfile,sizeof(sourcelistfile));
	sprintf(sourcelistfile,"sourcelist_%d.txt",rank);
	fp=fopen(sourcelistfile,"w");
	if(fp==NULL)
	{
		printf("Error(%s:%d): %d rank is unable to create sourcelist file",__FILE__,__LINE__,rank);
		exit(0);
	}
	fprintf(fp,"%d\n",numsources);
	for(i=0;i<numsources;i++)
	{
		fprintf(fp,"%f\t%f\n",tx[i],ty[i]);	
	}
	fclose(fp);	
/*	Sourcelist file created
 * */	
	jobmap=(int*)malloc(sizeof(int)*size);

	char allmasterid[size][256];
	if(rank==0)
	{
/*		// Simulating the response file from rank 0. Continue the executation once response file is ready.
		bzero(command,sizeof(command));
		sprintf(command,"cztrspgen thetax=%f thetay=%f rspfile=%s badpixfile=%s history=y clobber=y",tx[0],ty[0],responsefile,badpixfile);
		printf("%s\n",command);
		system(command);
*/
	}
	printf("%d rank waiting barrier\n",rank);
	srand (time(NULL));
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0)
	{

			for(i=0;i<size;i++)
				jobmap[i]=-1;
			counter=100;

				
			for(i=0;i<numsamples;i++)
			{

				/*
 				* Get free processor ID in a pool		
 				* */
				MPI_Recv ( &freeThreadId,1, MPI_INT, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &stat);
				if(jobmap[freeThreadId]!=-1)
				{

					jobmap[freeThreadId]=-1;

				}
				/* Assign configuration file to free processor 
 				* */
				MPI_Send ( &i,1, MPI_INT, freeThreadId, 1, MPI_COMM_WORLD );
				jobmap[freeThreadId]=counter;
			}
			
			sampleid=-1;

			counter=-1;
			for(i=1;i<size;i++)
			{
				if(jobmap[i]!=-1)
				{
					/*
 					*Send free processors id
 					*  		*/
					MPI_Recv ( &freeThreadId,1, MPI_INT, i, 2, MPI_COMM_WORLD, &stat);

					jobmap[i]=-1;
				}
				/*
 				*	Send terminate flag
 				* */
				MPI_Send ( &sampleid,1, MPI_INT, i, 1, MPI_COMM_WORLD );
			}
		}//end of rank 0 if
		else
		{
			isfirst=1;
			while(1)
			{
				MPI_Send ( &rank, 1, MPI_INT, 0, 2, MPI_COMM_WORLD );
				/*
				if(!isfirst)
				{
				//send result
				}*/
				isfirst=0;

				/*
 				* Receive configuration file name to process
 				* 		*/
				MPI_Recv ( &sampleid,1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &stat);
				printf("Rank %d\tReceived Job:%d\n",rank,sampleid);
				/* Processed all configuration file and terminating
 				* */
				if(sampleid==-1)
				{
					break;
				}
				else
				{
					bzero(temp,sizeof(temp));
					//strcpy(temp,"addevt");
					sprintf(temp,"addevt %f",exposure);
					/*	
 *					Simulating event files
 *
 * */
					for(i=0;i<numsources;i++)
					{

						// Simulating the response file from rank 0. Continue the executation once response file is ready.

						sprintf(responsefile,"simulation_%0.2f_%0.2f.rsp",tx[i],ty[i]);
						bzero(command,sizeof(command));
						sprintf(command,"cztrspgen thetax=%f thetay=%f rspfile=%s badpixfile=%s history=y clobber=y",tx[i],ty[i],responsefile,badpixfile);
						printf("%s\n",command);
						system(command);

						bzero(command,sizeof(command));
						sprintf(tmpoutfile,"%s_%d_%d.evt",outfile,i,rank);	
						sprintf(command,"simevt thetax=%f thetay=%f p=%f a=%f exposure=%f responsefile=%s badpixfile=%s outfile=%s",tx[i],ty[i],photindex[i],norm[i],exposure,responsefile,badpixfile,tmpoutfile);
						printf("%s\n",command);
						system(command);
						sprintf(temp,"%s %s",temp,tmpoutfile);
						bzero(command,sizeof(command));
						bzero(tmpoutfile,sizeof(tmpoutfile));
					}
					/*Event file simulation complete*/
					/*Adding event files*/
					bzero(tmpoutfile,sizeof(tmpoutfile));
					sprintf(tmpoutfile,"%s_%d_bg.evt",outfile,rank);	
					sprintf(command,"%s %s %s",temp,backgroundfile,tmpoutfile);
					printf("%s\n",command);
					system(command);


					bzero(command,sizeof(command));
					bzero(outphafile,sizeof(outphafile));
					sprintf(outphafile,"simulation_%d",rank);	
					sprintf(command,"bayesian infile=%s responsefile=%s sourcelistfile=%s outfile=%s badpixfile=%s",tmpoutfile,responsefile,sourcelistfile,outphafile,badpixfile);
					printf("%s\n",command);	
					system(command);	
					sprintf(tmpoutfile,"%s_%d_bg.evt",outfile,rank);	

					for(i=1;i<numsources+1;i++)
					{
						bzero(temp,sizeof(temp));

						sprintf(responsefile,"simulation_%0.2f_%0.2f.rsp",tx[i],ty[i]);
						sprintf(temp,"%s_%d.pha",outphafile,i);
						//sprintf(command,"fit_spectrum.py  %s  %s --outfile %s_%d.txt",temp,responsefile,outputfile,i);	
						printf("%s\n",command);
						system(command);
						remove(temp);
					
					}
				}
			}

		}
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0)
	{
		printf("Execution complete\n");
	}
	free(jobmap);

	for(i=0;i<numsources;i++)
	{
		bzero(temp,sizeof(temp));
		sprintf(temp,"%s_%d.evt",outfile,rank);	
		//removing simulated event files
		remove(temp);
	}
	//removing background added event file
	sprintf(tmpoutfile,"%s_%d_bg.evt",outfile,rank);	
	remove(tmpoutfile);
	remove (sourcelistfile);
	MPI_Finalize();
}


