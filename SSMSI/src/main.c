/* Tool to perform imaging on the SSM data. This code generates simulated data for SLANT1, SLANT2, and BOOM camera. The simulated data is stored in the shadows directory.
 * After generating the simulated data, code perfroms forward fitting and tries to reconstruct the position of the simulated sources. Logs of the reconstruction is dumped in the log directory and output is stroed in the output directory.
 *
 * This code uses all available resources, i.e., 4 compute nodes, 12 Tesla M-2050 GPUs  on hyades cluster. Code uses MPI, CUDA and pthreads for job distribution.
 * 
 * 
 * Author: Ajay Vibhute, ajay@iucaa.in
 * Guide: Prof. Dipankar Bhattacharya
 */

#include<stdio.h>
#include<string.h>
#include<mpi.h>
#include<time.h> 
#include<stdlib.h>
//defined scripts to manage the queues
#define updateJobQueue "../scripts/updateQueue.csh"
#define jobQueuePath "../shadows/log/jobQueue"
#define moveJob "../scripts/moveJob.csh"
extern void startExecution(char *observationId);

int main(int argc,char*argv[])
{
	//defining required local variables

	char workQueue[10][50],removeCommand[100];
	char currentWork[10][50];
	char hostname[100],freeNodeList[100]; 
	char job[100],command[100],currentJob[100];
	char directoryQueue[10][100];	

	int myrank=0,npes=0;
	int i=0,counter=0,freeThreadId=0,rv=0;
	int directoryCount=0;
	float timeSpent=0.0;

	FILE *getCurrentJob=NULL;
	time_t now=NULL;

	//initializing MPI environment 
	MPI_Status stat;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&npes);
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

	// Rank 0 acts as a master and distributes job to other avaiable processes.
	if(myrank==0)
	{	

			
		gethostname(hostname,sizeof(hostname));
		printf("Master is in node %s\n",hostname);

		//updating the data queue
		system(updateJobQueue);				

		//starting the counter
		timeSpent=MPI_Wtime();
 		time(&now);
		counter=0;
		rv=0;
		//Running infinite loop, loop will stop once all data is processed.
		while(1)	
		{
			//Getting id of the free thread and assinging job to it
			MPI_Recv ( &freeThreadId,1, MPI_INT, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &stat);
			if(freeThreadId>0)
			{
				//getting job id of the first available job
				getCurrentJob=fopen(jobQueuePath,"r");
				if(getCurrentJob==NULL)
				{
					printf("Error(%s:%d): No jobs found\n",__FILE__,__LINE__);
					break;
				}
				rv=fscanf(getCurrentJob,"%s",currentJob);
				fclose(getCurrentJob);	
				if(rv==-1)
				{

					printf("Error(%s:%d): No jobs found\n",__FILE__,__LINE__);
					break;
				}
				else
				{
					//performing queue management
					strcpy(command,moveJob);
					strcat(command," process ");
					strcat(command,currentJob);
					system(command);
				}
					
			        //Sending job id to the avaiable processor for further processing					
				MPI_Send ( currentJob,sizeof(workQueue[counter]), MPI_CHAR, freeThreadId, 1, MPI_COMM_WORLD );
				//printf("sending %s job to %d node\n",currentJob,freeThreadId);
			}
			counter++;
		}
			
		//Once processing is complete send message to avaiable processes and cleanup the environment
		for(i=1;i<npes;i++)
		{
			char done[100];
			strcpy(done,"done");
			MPI_Send ( &done, sizeof(done), MPI_CHAR,i, 1, MPI_COMM_WORLD );
		}
		
	}
	else
	{

		//Run in a loop till all observation ids are processed
		while(1)
		{	
			//Request for a job from Master process	
			MPI_Send ( &myrank, 1, MPI_INT, 0, 2, MPI_COMM_WORLD );

			//Receive a job from master process
			MPI_Recv ( &job,sizeof(job), MPI_CHAR, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &stat);

			//come out of loop if all ids are processed
			if(strcmp(job,"done")==0)
			{
				break;
			}		
			else
			{	
				//start processing the observation id and update queues once processing is done.
				//In the processing code processes data for all three camera, and perform reconstuction,  
				//find streangth of the source using forward fitting via SVD fit.

				startExecution(job);
				//startExecution("20110720:0123456789"); 

				//performing the queue mangement
				strcpy(command,moveJob);
				strcat(command," done ");
				strcat(command,job);
				
				system(command);
				strcpy(command,moveJob);
				strcat(command," repository ");
				strcat(command,job);
				system(command);
				
		
				
			}
		}
		
	}
	//wait till all processes come to this stage	
	MPI_Barrier(MPI_COMM_WORLD);
	if(myrank==0)
	{
		/*
			Stop the timer and calculate the total time spent in the execution 
		*/
		timeSpent=MPI_Wtime()-timeSpent;
		printf("Time spent=%f\n",timeSpent);
	
	}
	//clean up MPI environment
	MPI_Finalize();	
	
}
