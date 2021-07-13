/*
 * lightCurveParam_mpi.c
 *
 *  Created on: Oct 21, 2013
 *  Author: Ajay Vibhute (IUCAA)
 *  Version: 1.0
 *  Description: This is the mpi version of the code which calculate various parameters of the
 *  light curve. This program will divide the job in multiple processors available to use.
 */
#include<mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include<time.h>
#include "statSubDef_mpi.c"
int rank,size,counter=0,num_curves=0,cpp=0;

int main(int argc,char*argv[])
{
	int data_size,rv=0,i=0,isfirst=1,index=0,totpoints=0,freeThreadId,*jobmap,numpoints=0;
	double start=0;
	double temp_mag,temp_magerr,temp_mjd;
	char last_id[256],id[256],infile[BUFFSIZE],outfile[BUFFSIZE],**filenames;
	FILE *fp;
	struct lightcurve*cur_curve=NULL;
	struct lightcurveParam *paramin=NULL;
	size_t structsize=sizeof(struct lightcurveParam);
	void    *buffer=NULL;
	
	MPI_Status stat;
	MPI_Init(&argc,&argv);
	MPI_Comm_size( MPI_COMM_WORLD, &size );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	jobmap=(int*)malloc(sizeof(int)*size);
	paramin=(struct lightcurveParam*)malloc(sizeof(struct lightcurveParam));
	cur_curve=(struct lightcurve*) malloc(sizeof(struct lightcurve));
	cur_curve->mag=(double*)malloc(sizeof(double));
	cur_curve->mag_error=(double*)malloc(sizeof(double));
	cur_curve->mjd=(double*)malloc(sizeof(double));
	cur_curve->temp=(double*)malloc(sizeof(double));
	if(rank==0)
	{
		if(argc<4)
		{
			printf("Error: insufficient inputs. Please provide all command line inputs\n");
			printf("Inputs:\n\t1: Number of light curves\n\t2: Input file\n\t3: Output file\n");
			printf("Exiting....\n");
			exit(0);//check alternative for normal termination

		}
		strcpy(infile,argv[2]);
		strcpy(outfile,argv[3]);
		start = MPI_Wtime();
		num_curves=parseFile(infile);
		if(num_curves==-1)
		{
			//printf("Error (%s:%d): in estimating number of light curves\n",__FILE__,__LINE__);
			exit(0);
		}
		printf("Number of light curves:%d\n",num_curves);
		fp=fopen(infile,"r");
		if(fp==NULL)
		{
			printf("Error (%s:%d): %s file not found\n",__FILE__,__LINE__,infile);
			exit(0);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Bcast(&infile[0],sizeof(infile),MPI_BYTE,0,MPI_COMM_WORLD);
	MPI_Bcast(&num_curves,1,MPI_INT,0,MPI_COMM_WORLD);
	filenames=(char**)malloc(sizeof(char*)*num_curves);
	if(filenames==NULL)
	{
		printf("Error (%s:%d) : %d process out of memory\n",__FILE__,__LINE__,rank);
		printf("Exit...\n");
		exit(0);
	}
	for(i=0;i<num_curves;i++)
		filenames[i]=(char*)malloc(sizeof(char)*BUFFSIZE);
	fp=fopen(infile,"r");
	if(fp==NULL)
	{
		printf("Error: %s file not found\n",infile);
		printf("Exiting....\n");
		exit(0);
	}
	for(i=0;i<num_curves;i++)
	{
		fscanf(fp,"%s",&filenames[i][0]);
	}
	fclose(fp);

	if(rank==0)
	{
			fp=fopen(outfile,"w");
			if(fp==NULL)
			{
				printf("Error (%s:%d) Error while creating %s file\n",__FILE__,__LINE__,outfile);
				printf("Writing output to standard output\n");
				fp=stdout;
			}
			for(i=0;i<size;i++)
				jobmap[i]=-1;
			counter=0;
			fprintf(fp,"MasterId,count,mean,sd,skew,kurtosis,mad,bwmv,thiel_sen,durbin_watson,stetson_j,stetson_k,kendall_tau,cusum,con,abbe,reducedchi\n");
			while(counter!=num_curves)
			{
				freeThreadId=2;

				MPI_Recv ( &freeThreadId,1, MPI_INT, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &stat);
				if(jobmap[freeThreadId]!=-1)
				{
					MPI_Recv(paramin,structsize , MPI_BYTE, freeThreadId, 3, MPI_COMM_WORLD, &stat);
					fprintf(fp,"%s,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",removeFileExt( basename(filenames[jobmap[freeThreadId]])),paramin->numpoints,paramin->xm,paramin->xsig,paramin->xsk,paramin->xkur,paramin->xmad,paramin->xbwmv,paramin->thiel_sen,paramin->xdw,paramin->stetson_j,paramin->stetson_k,paramin->xkt,paramin->cusum,paramin->con,paramin->abbe,paramin->rchi);
					jobmap[freeThreadId]=-1;
				}
				if(counter<num_curves)
				{
					MPI_Send ( &counter,1, MPI_INT, freeThreadId, 1, MPI_COMM_WORLD );
					jobmap[freeThreadId]=counter++;
					//Putting counter value as current job for requesting thread and incrementing counter value
				}
			}
			counter=-1;
			for(i=1;i<size;i++)
			{
				if(jobmap[i]!=-1)
				{
					MPI_Recv ( &freeThreadId,1, MPI_INT, i, 2, MPI_COMM_WORLD, &stat);
					MPI_Recv(paramin,structsize , MPI_BYTE, i, 3, MPI_COMM_WORLD, &stat);
					fprintf(fp,"%s,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",removeFileExt( basename(filenames[jobmap[freeThreadId]])),paramin->numpoints,paramin->xm,paramin->xsig,paramin->xsk,paramin->xkur,paramin->xmad,paramin->xbwmv,paramin->thiel_sen,paramin->xdw,paramin->stetson_j,paramin->stetson_k,paramin->xkt,paramin->cusum,paramin->con,paramin->abbe,paramin->rchi);
					jobmap[i]=-1;
				}
				MPI_Send ( &counter, 1, MPI_INT,i, 1, MPI_COMM_WORLD );
			}
			fclose(fp);
		}
		else
		{
			isfirst=1;
			while(1)
			{
				MPI_Send ( &rank, 1, MPI_INT, 0, 2, MPI_COMM_WORLD );
				if(!isfirst)
				{
					rv = MPI_Send(paramin, structsize, MPI_BYTE, 0, 3, MPI_COMM_WORLD);
					if (rv != MPI_SUCCESS)
					{
						printf("Buffer attach failed. Return code= %d Terminating\n", rv);
						MPI_Finalize();
					}
				}
				isfirst=0;
				MPI_Recv ( &counter,1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &stat);
				if(counter==-1 || counter>num_curves)
				{
					break;
				}
				else
				{
					numpoints=parseFile(filenames[counter]);//add error bit in param and if num point is -1 then set error flag for that bit in param
					cur_curve->numpoints=numpoints;
					cur_curve->mag=(double*)realloc(cur_curve->mag,sizeof(double)*numpoints);
					cur_curve->mag_error=(double*)realloc(cur_curve->mag_error,sizeof(double)*numpoints);
					cur_curve->mjd=(double*)realloc(cur_curve->mjd,sizeof(double)*numpoints);
					cur_curve->temp=(double*)realloc(cur_curve->temp,sizeof(double)*numpoints*numpoints);

					if(readCurve(filenames[counter],cur_curve)!=-1)
					{
						paramin->numpoints=numpoints;
						processLightCurve(cur_curve,paramin);
						//printf("Processor %d processing light curve named %s with index %d\n",rank,filenames[counter],counter);
					}
					//printf("Processor %d processing light curve named %s at %d with %d datapoints\tmean:%f \n",rank,filenames[counter],counter,numpoints,paramin->xm);
					cpp++;

				}
			}

		}
		MPI_Barrier(MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0)
		printf("Execution time:%f\n",MPI_Wtime()-start);

	//printf("%d processor work done waiting for finalize\n",rank);
	MPI_Finalize();
}
int readCurve(char*infile,struct lightcurve *curves)
{
	FILE *fp;
	int i=0;
	double tmpmjd,tmpmag,tmpmagerr;
	fp=fopen(infile,"r");
	if(fp==NULL)
	{
		printf("Error (%s:%d) %s file not found\n",__FILE__,__LINE__,infile);
		return -1;
	}
	/*Modify to accept it in any format*/
	for(i=0;i<curves->numpoints;i++)
	{
		fscanf(fp,"%lf",&tmpmjd);
		curves->mjd[i]=tmpmjd;
		fscanf(fp,"%lf",&tmpmag);
		curves->mag[i]=tmpmag;
		fscanf(fp,"%lf",&tmpmagerr);
		curves->mag_error[i]=tmpmagerr;
	}
	fclose(fp);
	return 0;
}
int parseFile(char*infile)
{
	FILE *fp;
	int lccounter=0,rv=0;
	int c;
	fp=fopen(infile,"r");
	if(fp==NULL)
	{
		printf("Error (%s:%d) %s file not found\n",__FILE__,__LINE__,infile);
		return -1;
	}
	while ( (c=fgetc(fp)) != EOF ) {
	        if ( c == '\n' )
	            lccounter++;
	    }
	fclose(fp);
	return lccounter;

}
