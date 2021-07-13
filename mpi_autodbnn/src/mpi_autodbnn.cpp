/*
 * Name			: mpi_autodbnn.cpp
 * Author		: Ajay
 * Created on	: Jan 20, 2013
 * Description 	:This is the MPI version of the the autobdnn originally written by Prof. Sajith Philip
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <mpi.h>
using namespace std;
#include <stdlib.h>
#include<sys/times.h> // times() fun. is here.
#include <time.h>
#include <vector>
using std::vector;
#define max_resol 1600
#define features 100
#define classes 500  
#define master_rank 0
static double bgain,gain,dmyclass[classes+2],classval[classes+2],cmax,c2max,c3max,c4max,tmp2_wts,totprob,oldj;
static double LoC=0.65;
static double nLoC=0.0;
//NSP_Added new variable jx=0;
static int jx=0,resol=100,nresol=0,nerror=0,nLoCcnt=1,skpchk=0,MissingDat=-9999;
static double omax,omin,rslt,rslt2,orslt,orslt2,prslt,nrslt,fst_gain;
clock_t start,stop;
static int argfnd, oneround=100,kmax,k2max,k3max,k4max,ans1,tcnt,rnn,rnd,i,j,k,l,n,m,p,c1cnt,c2cnt,pcnt,pocnt,invcnt,innodes=100,outnodes=100,send_status=0;
char fln[256],fltmp[256],urchoice,urchoicex,bcchoice,savedpar,datfilename[256];
FILE *fl1,*fl2,*fl3,*fl4,*fl5,*fl6,*fl7,*fl8,*fl9,*fl10;
double *arr_anti_wts_temp=(double*) malloc(512*sizeof(double));
MPI_Status status;
int main(int argv, char *argp[256])
{
	int rank, size;
	double start_t,end_t;
	int testsize=0,rc=0;
	int jc=0,jpt,jobsperthread=0,exjpt=0;;
	int ii=0,index=0,remaining=0,tobesent=0,tobereceived=0;
	MPI_Init (&argv, &argp);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	MPI_Comm_size (MPI_COMM_WORLD, &size);
	strcpy(fln,argp[1]);
	rnd=rank;
	if(rank==master_rank)
	{
		start_t=MPI_Wtime();
		if(argv > 3)
		{
			argfnd=1;
			cout << "The selected option is " << *argp[3] <<"\n";
			switch(*argp[3])
			{
				case '0':
					ans1=0;
					if((fl2=fopen("0.par","r"))!=NULL)
					{
						fscanf(fl2,"%c\n",&bcchoice);  //Handle missing or out of range values? Y if yes. NEW in Ver 7
						fscanf(fl2,"%c\n",&urchoice);
						fscanf(fl2,"%c\n",&savedpar);
						fscanf(fl2,"%c\n",&urchoicex);
						if(bcchoice == 'Y'||bcchoice =='y')
						{
							fscanf(fl2,"%d\n",&skpchk);
							if(skpchk <0) MissingDat=skpchk;
							cout << "System  is configured for handling missing data with missing data indicator" << MissingDat <<"\n";
						}
						fclose(fl2);
					}
				else
				{
					cout << "No Parameter File... existing..";
					exit(1);
				}
				break;
				case '1':
					ans1=1;
					if((fl2=fopen("0.par","r"))!=NULL)
					{
						fscanf(fl2,"%c\n",&bcchoice);  //Handle missing or out of range values? Y if yes. NEW in Ver 7
						fscanf(fl2,"%c\n",&urchoice);
						fscanf(fl2,"%c\n",&savedpar);
						fscanf(fl2,"%c\n",&urchoicex);
						if(bcchoice == 'Y'||bcchoice =='y')
						{
							fscanf(fl2,"%d\n",&skpchk);
							if(skpchk <0) MissingDat=skpchk;
							cout << "System  is configured for handling missing data with missing data indicator" << MissingDat <<"\n";
						}
						fclose(fl2);
					}
					else
					{
						cout << "No Parameter File... existing..";
						exit(1);
					}
					if((fl2=fopen("1.par","r"))!=NULL)
					{
						fscanf(fl2,"%lf",&gain);
						fscanf(fl2,"%d",&oneround);
						fclose(fl2);
					}
					else
					{
						cout << "No Parameter File... existing..";
						exit(1);
					}
				break;
				case '2':
					ans1=2;
					if((fl2=fopen("0.par","r"))!=NULL)
					{
						fscanf(fl2,"%c\n",&bcchoice);  //Handle missing or out of range values? Y if yes. NEW in Ver 7
						fscanf(fl2,"%c\n",&urchoice);
						fscanf(fl2,"%c\n",&savedpar);
						fscanf(fl2,"%c\n",&urchoicex);
						if(bcchoice == 'Y'||bcchoice =='y')
						{
							fscanf(fl2,"%d\n",&skpchk);
							if(skpchk <0) MissingDat=skpchk;
							cout << "System  is configured for handling missing data with missing data indicator" << MissingDat <<"\n";
						}
						fclose(fl2);
					}
					else
					{
						cout << "No Parameter File... existing..";
						exit(1);
					}
				break;
				case '3':
					ans1=3;
					if((fl2=fopen("0.par","r"))!=NULL)
					{
						fscanf(fl2,"%c\n",&bcchoice);  //Handle missing or out of range values? Y if yes. NEW in Ver 7
						fscanf(fl2,"%c\n",&urchoice);
						fscanf(fl2,"%c\n",&savedpar);
						fscanf(fl2,"%c\n",&urchoicex);
						if(bcchoice == 'Y'||bcchoice =='y')
						{
							fscanf(fl2,"%d\n",&skpchk);
							if(skpchk <0) MissingDat=skpchk;
							cout << "System  is configured for handling missing data with missing data indicator" << MissingDat <<"\n";
						}
						fclose(fl2);
					}
					else
					{
						cout << "No Parameter File... existing..";
						exit(1);
					}
				break;
					default:
						cout << "Create the APF file(0) or Create the Weights file (1) or Classify Data(2,3) ?";
						cin >> ans1;
				break;
			}
		}
		else
		{
			argfnd=0;
			cout << "Create the APF file(0) or Create the Weights file (1) or Classify Data(2,3) ?";
			cin >> ans1;
		}
		if(ans1 == 2)
		{
			if(argfnd==1)
			bgain=0.0;
			else
			{
				cout << "Allowed relaxation on the boundary (in % use 0 for default from training data) :";
				cin >> bgain;
				bgain=bgain*1.0;
			}
		}
		else
		bgain= 0;  // During training we are strict on boundary constraints.
		if(argv < 3)
		{
			cout << "Enter the name of the input file without extension (dat) :";
			cin >> fln;
		}
		else
		{
			strcpy(fln,argp[1]);
		}
		if(argp[4]!=NULL)
			strcpy(datfilename,argp[4]);
		else
		{
			strcpy(datfilename,argp[1]);
			strcat(datfilename,".dat");
		}
		strcpy(fltmp,datfilename);

		if((fl1=fopen(fltmp,"r"))!=NULL)
		{
			strcpy(fltmp,fln);
			strcat(fltmp,".inf");
			if((fl2=fopen(fltmp,"r"))!=NULL)
			{
				i=0;
				fscanf(fl2,"%d",&innodes);
				fscanf(fl2,"%d",&outnodes);
				for (i=0;i<=outnodes;i++) // dmyclass[0] contains margin others are expected values.
				fscanf(fl2,"%lf\n",&dmyclass[i]);
				fscanf(fl2,"%lf",&LoC);   // New parameter to specify the Line Of Control
				fscanf(fl2,"%d",&nresol);
				fscanf(fl2,"%d",&nerror);
				cout <<"You have "<< innodes << " input nodes and " << outnodes <<" Output nodes with " << "margin set to " << LoC << "\n";
				cout << "The target outputs are\n";
				for (i=0;i<=outnodes;i++) cout << dmyclass[i] <<"\n";
				if(nresol >0)
				{
					resol=nresol;cout << "The maximum binsize is: " << resol <<"\n";
				}
				else
				{
					cout << "The maximum binsize is: " << resol<<"\n";
				}
				fst_gain=1.0/outnodes;
			}
			else
			{
				cout << "Unable to find the Info file. Exiting !!";
				exit(1);
			}

		} // program ends.
		else   // data file read error.
		{
			cout << "Unable to open the data file";
			exit(1);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(datfilename,256,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&innodes,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&outnodes,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&resol,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&argfnd,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&ans1,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&MissingDat,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&dmyclass,outnodes+2,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&gain,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&LoC,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&oneround,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nerror,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&savedpar,1,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

 /**************** Let us Define the Network Structure *********************************/
	double vects[innodes+outnodes+2],vectso[innodes+outnodes+2],tmpv,max[innodes+2],min[innodes+2];
	double err1vects[innodes+2], err2vects[innodes+2];
	//vector<vector<vector<vector<vector<int> > > > > anti_net(innodes+2,vector<vector<vector<vector<int> > > >(resol+2, vector<vector<vector<int> > >(innodes+2,vector<vector<int> >(resol+4,vector<int>(outnodes+2,0)))));

	int totsize=(innodes+2)*(resol+2)*(innodes+2)*(resol+4)*(outnodes+2);
	int totsendreceivesize=(innodes+1)*(resol+2)*(innodes+1)*(resol+1)*(outnodes+1);
	double *arr_anti_wts=(double*) malloc(totsize*sizeof(double));
	int *arr_anti_net=(int*)malloc(sizeof(int)*totsize);
	int *anti_net_temp=(int*) malloc(totsize*sizeof(int));
	int ik=innodes,jk=resol,lk=innodes,mk=resol,kk=outnodes;
	// Max Threshold
	int resolution[innodes+8];
	double classtot[innodes+2][resol+2];           // Total Prob. computed
	if(classtot==NULL){cout << "Out of Memory to Run Code at classtot.. Exiting\n";exit(1);}
	double binloc[innodes+4][resol+8];
	if(binloc==NULL){cout << "Out of Memory to Run Code at binloc.. Exiting\n";exit(1);}
  /***************************Let us put up the Network***********************************/
//    Start the counter for case 2 here.................
	//if(rank==0)
	{
		strcpy(fltmp,datfilename);

		fl1=fopen(fltmp,"r");
		start = times(NULL);

		if (ans1==0)
		{
			n=0;
			omax=-400;
			omin=400;
			while (!feof(fl1))
			{
				skpchk=0;
				for(i=1;i<=innodes;i++)
				if (n==0)
				{
				fscanf(fl1,"%lf",&vects[i]);
				if(nerror ==2){fscanf(fl1,"%lf",&err1vects[i]);fscanf(fl1,"%lf",&err2vects[i]);}else
				if(nerror ==1){fscanf(fl1,"%lf",&err1vects[i]); err2vects[i]=err1vects[i];}
				if(vects[i] != MissingDat)
				{
					min[i]=vects[i];
					max[i]=vects[i];
				}
					else max[i]=MissingDat;
				}
				else
				{
					fscanf(fl1,"%lf",&vects[i]);
					if(vects[i] != MissingDat)
					{
						if( vects[i]> max[i]) max[i]=vects[i];
						if (min[i] > vects[i]) min[i]=vects[i];
					}
				}
				fscanf(fl1,"%lf\n",&tmpv);
				if(tmpv>omax) omax = tmpv;
				if(tmpv<omin) omin =tmpv;
				k=1;
				j=0;
				n++;
			}
			if(rank==master_rank)
			cout << "No of vectors =" << n <<" and i/n is= " << 1.0/n << "\n";
			for(i=1;i<=innodes;i++)
			{
				if(min[i]==max[i])if(min[i]!=0){min[i]= -1.0*max[i];}else{min[i]=0.0; max[i]=1.0;}
			}
			if(argfnd==0)
			{
				cout <<"Do you want to use the saved parameters (Y/N)? ";
				cin >>savedpar;
			}
			if (savedpar == 'y') savedpar='Y';
			else
			if(savedpar == 'n') savedpar='N';
			if((savedpar == 'Y') || (savedpar=='y'))
			{
				strcpy(fltmp,fln);
				strcat(fltmp,".apf");
				fl2=NULL;
				if((fl2=fopen(fltmp,"r"))!=NULL)
				{
					if(rank==master_rank)
					cout << "Reading from the saved information\n";
					for (i=1;i<=innodes;i++)
					{
						fscanf(fl2,"%d",&resolution[i]);
						for(j=0;j<=resolution[i];j++) binloc[i][j+1]=j*1.0;
					}
					if(rank==master_rank)
					cout << innodes << " items read from " << fltmp <<"\n";
				}
				else
				{
					cout << "ERROR: File " << fltmp << " not found" << "\n";
					exit(1);
				}
			}
			else
			for(i=1;i<=innodes;i++)
			{
				if(min[i]==max[i])if(min[i]!=0){min[i]= -1.0*max[i];}else{min[i]=0.0; max[i]=1.0;}
				cin >> resolution[i];
				for(j=0;j<=resolution[i];j++) binloc[i][j+1]=j*1.0;
			}
			for(k=1;k<=outnodes;k++)
			for(i=1;i<=innodes;i++)
			for(j=0;j<=resolution[i];j++)
			for(l=1;l<=innodes;l++)
			for(m=0;m<=resolution[l];m++)
			{
				arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+k)]=1;
				arr_anti_wts[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+k)]=(double)(1.0);
				anti_net_temp[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+k)]=1;
			}
			  // Start the counter now...............
			start = times(NULL);
			rewind(fl1);
			tcnt=0;



			while (!feof(fl1))
			{
				tcnt++;
				for (i=1;i<=innodes;i++)
				{
					fscanf(fl1,"%lf",&vects[i]);
					if(nerror ==2){fscanf(fl1,"%lf",&err1vects[i]);fscanf(fl1,"%lf",&err2vects[i]);}else
					if(nerror ==1){fscanf(fl1,"%lf",&err1vects[i]);err2vects[i]=err1vects[i];}
				}
				fscanf(fl1,"%lf\n",&tmpv);
				for(i=1;i<=innodes;i++)
				{
					if((vects[i] != MissingDat)&&(max[i] !=MissingDat))
					{
						vectso[i]=vects[i];
						vects[i]=round((vects[i]-min[i])/(max[i]-min[i])*resolution[i]);
						err1vects[i]=round((err1vects[i])/(max[i]-min[i])*resolution[i]);
						err2vects[i]=round((err2vects[i])/(max[i]-min[i])*resolution[i]);
					}
				}

				jobsperthread=innodes/size;
				jpt=innodes/size;
				exjpt=0;
				if(innodes%size!=0)
				{
					if(rank<innodes%size)
					{
						jobsperthread+=1;
						exjpt=rank;
					}
					else
						exjpt=innodes%size;
				}

				for (jc=0;jc<jobsperthread;jc++)
				{
					i=(jc*size)+rank+1;

				/*for (i=1;i<=innodes;i++)
				{*/
					j=0;
					if(vects[i] != MissingDat)
					{
						while ((fabs(vects[i]-binloc[i][j+1]) >=1.0)&& (j<= resolution[i]))
						{
							j++;
						}
						for (l=1;l<=innodes;l++)
						{
							m=0;
							if(i!=l)
							{
								while ((fabs(vects[l]-binloc[l][m+1]) >=1.0)&& (m<= resolution[l]))
								{
									m++;
								}
								k=1;
								while ((k<=outnodes)&&(fabs(tmpv - dmyclass[k])) > dmyclass[0]) k++;
								arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+k)]++;
								arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+0)]++;


							}
						}
					}
				}

			}

			if(rank!=0)
			{

				MPI_Pack_size (totsendreceivesize, MPI_DOUBLE, MPI_COMM_WORLD, &testsize);
				testsize = testsize +  MPI_BSEND_OVERHEAD;
				void* buffer = (void*)malloc(testsize);
				rc = MPI_Buffer_attach(buffer, testsize);
				if (rc != MPI_SUCCESS) {
					printf("Buffer attach failed. Return code= %d Terminating\n", rc);
					MPI_Finalize();
				}
				MPI_Bsend((void*)arr_anti_net,totsendreceivesize,MPI_INT,0,11,MPI_COMM_WORLD);
				/*add mpi dettach here*/
				 MPI_Buffer_detach(&buffer, &testsize);
			}
			else
			{
				ii=1;
				index=0;
				for(ii=1;ii<size;ii++)
				{
					jc=0,jpt;
					jobsperthread=innodes/size;
					jpt=innodes/size;
					exjpt=0;
					if(innodes%size!=0)
					{
						if(ii<innodes%size)
						{
							jobsperthread+=1;
							exjpt=ii;
						}
						else
							exjpt=innodes%size;
					}
					remaining=totsendreceivesize;
					tobereceived=128;
					k=0;

					MPI_Recv(anti_net_temp,totsendreceivesize,MPI_INT,ii,11,MPI_COMM_WORLD,&status);

					for (jc=0;jc<jobsperthread;jc++)
					{
						i=(jc*size)+ii+1;
						for(j=0;j<=jk;j++)
							for(m=0;m<=mk;m++)
								for(l=0;l<=lk;l++)
									for(k=0;k<=kk;k++)
									{

										arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+k)]=anti_net_temp[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+k)];
										anti_net_temp[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+k)]=0;
									}
					}
				}

			}
			free(anti_net_temp);
			fclose(fl1);
			fclose(fl2);
			MPI_Barrier(MPI_COMM_WORLD);

			 /*
				The conditional Probability,
			P(A|B) = P(A intersection B)/P(B) is the
			probability for the occurance of A(k) if B(ij) has happened =
			Share of B(ij) that is held by A(k) / Probability of total B(ij)
			in that particular feature i with resolution j.

						  */
			if(rank==master_rank)
			{
				strcpy(fltmp,fln);
				strcat(fltmp,".awf");      // This file holds the weights
				fl6=fopen(fltmp,"w+");
				strcpy(fltmp,fln);
				strcat(fltmp,".apf");     // This file holds the estimated probability
				if((fl1=fopen(fltmp,"w+"))!=NULL)
				{
					for(i=1;i<=innodes;i++) fprintf(fl1,"%d ",resolution[i]);
					fprintf(fl1,"\n%lf %lf \n",omax,omin);
					for(i=1;i<=innodes;i++) fprintf(fl1,"%lf ",max[i]);
					fprintf(fl1,"\n");
					for(i=1;i<=innodes;i++) fprintf(fl1,"%lf ",min[i]);
					fprintf(fl1,"\n");
					for(k=1;k<=outnodes;k++)
					{
						for(i=1;i<=innodes;i++)
						for(j=0;j<=resolution[i];j++)
						{
							for(l=1;l<=innodes;l++)
							if(i!=l)
							{
								for(m=0;m<=resolution[l];m++)
								{
									fprintf(fl1,"%d ",arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+k)]);
									fprintf(fl6,"%lf ",(double)arr_anti_wts[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+k)]);

								}
								fprintf(fl6,"\n");
								fprintf(fl1,"\n");
							}
						}
						fprintf(fl6,"\n");
						fprintf(fl1,"\n");
					}
					fprintf(fl6,"\n");
					fprintf(fl1,"\n");
				}
				else
				{
					cout << "Unable to create file for output\n";
					exit(1);
				}
				for(i=1;i<=innodes;i++)
				for(j=1;j<=resolution[i];j++)
				fprintf(fl6,"%lf\n", (double)binloc[i][j]);                 /// Let us print the bins.
				fclose(fl1);
				fclose(fl6);
				fflush(NULL);
				cout << "Creating the Anticipated Weights data file\n";
			}
		}
	}
/**********************************End of Case 0 ******************************/


	if(ans1==1)
	{
		if(argfnd==0)
		{
			cout << "Please enter the gain:";
			cin >> gain;
			cout << "Please enter the number of training epochs:";
			cin >> oneround;
		}
		// Start the counter in this round here...................
		strcpy(fltmp,fln);
		strcat(fltmp,".awf");
		if((fl6=fopen(fltmp,"r"))!=NULL)
		{
			strcpy(fltmp,fln);
			strcat(fltmp,".apf");
			fl2=NULL;
			if((fl2=fopen(fltmp,"r"))!=NULL)
			{
				for (i=1;i<=innodes;i++)
				{
					fscanf(fl2,"%d",&resolution[i]);
					for(j=0;j<=resolution[i];j++) binloc[i][j+1]=j*1.0;
				}
				fscanf(fl2,"\n%lf",&omax);
				fscanf(fl2,"%lf",&omin);
				fscanf(fl2,"\n");
				for(i=1;i<=innodes;i++) fscanf(fl2,"%lf",&max[i]);
				fscanf(fl2,"\n");
				for(i=1;i<=innodes;i++) fscanf(fl2,"%lf",&min[i]);
				fscanf(fl2,"\n");
				for(i=1;i<=innodes;i++)for(j=0;j<=resolution[i];j++)
				for(l=1;l<=innodes;l++)for(m=0;m<=resolution[l];m++) arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+0)] =0;
				for(k=1;k<=outnodes;k++)
				{
					for(i=1;i<=innodes;i++)
					for(j=0;j<=resolution[i];j++)
					{
						for(l=1;l<=innodes;l++)
						if(i!=l)
						{
							for(m=0;m<=resolution[l];m++)
							{
								fscanf(fl2,"%d",&arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+k)]);
								arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+0)]+=arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+k)];
								if(rank==master_rank)
								fscanf(fl6,"%lf",&arr_anti_wts[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+k)]);
								else
								{
									arr_anti_wts[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+k)]=0;
									double temp=0;
									fscanf(fl6,"%lf",&temp);
								}
							}
							fscanf(fl2,"\n");
							fscanf(fl6,"\n");
						}
					}
					fscanf(fl2,"\n");
					fscanf(fl6,"\n");
				}
				for(i=1;i<=innodes;i++)
				for(j=1;j<=resolution[i];j++)
				fscanf(fl6,"%lf\n", &binloc[i][j]);                 /// Let us print the bins.
			}
			else
			{
				cout << "Unable to Open the APF information file\n";
				exit(1);
			}
			fclose(fl2);
		}
		else
		{
			cout << "Unable to Open the AWF information file\n";
			exit(1);
		}
		fclose(fl6);
		// Training round starts here....
		//for(rnd=0;rnd<=oneround;rnd++)
		jobsperthread=(oneround+1)/size;
		jpt=(oneround+1)/size;
		exjpt=0;
		if((oneround+1)%size!=0)
		{
			if(rank<(oneround+1)%size)
			{
				jobsperthread+=1;
				exjpt=rank;
			}
			else
				exjpt=(oneround+1)%size;
		}

		for (jc=0;jc<jobsperthread;jc++)
		{
			rnd=(jc*size)+rank;
/*			printf("%d\t%d\n",rnd,rank);*/
			strcpy(fltmp,datfilename);


			fl1=fopen(fltmp,"r");
			n=0;
			rslt=0.0;
			rslt2=0.0;
			pcnt=0;
			pocnt=0;
			orslt=rslt;
			orslt2=rslt2;

			ii=0;

			//if(rank!=master_rank)
			if(rnd!=0)
			{
				if(rank==master_rank)
					ii=size-1;
				else
					ii=rank-1;
				k=0;
				remaining=totsendreceivesize;
				tobereceived=256;
				while(remaining!=0)
				{
					if(remaining<256)
						tobereceived=remaining;
					MPI_Recv(&arr_anti_wts_temp[0],tobereceived,MPI_DOUBLE,ii,1,MPI_COMM_WORLD,&status);
					for(i=0,l=k;i<tobereceived;i++,l++)
					{
						arr_anti_wts[l]+=arr_anti_wts_temp[i];
						arr_anti_wts_temp[i]=0;
					}
					k+=tobereceived;
					remaining-=tobereceived;
				}
			}
			while (!feof(fl1))
			{
				for(k=1;k<=outnodes;k++) classval[k]=1.0;
				n++;
				if(ans1==3)
				{
					for (i=1;i<=innodes;i++)
					{
						fscanf(fl1,"%lf",&vects[i]);
						if(nerror ==2){fscanf(fl1,"%lf",&err1vects[i]);fscanf(fl1,"%lf",&err2vects[i]);}else
						if(nerror ==1){fscanf(fl1,"%lf",&err1vects[i]);err2vects[i]=err1vects[i];}
					}
					fscanf(fl1,"\n");
				}
				else
				{
					for (i=1;i<=innodes;i++)
					{
						fscanf(fl1,"%lf",&vects[i]);
						if(nerror ==2){fscanf(fl1,"%lf",&err1vects[i]);fscanf(fl1,"%lf",&err2vects[i]);}else
						if(nerror ==1){fscanf(fl1,"%lf",&err1vects[i]);err2vects[i]=err1vects[i];}
					}
					fscanf(fl1,"%lf\n",&tmpv);
				}
				for(i=1;i<=innodes;i++)
				{
					if((vects[i] != MissingDat)&&(max[i]!=MissingDat))
					{
						vectso[i]=vects[i];
						vects[i]=round((vects[i]-min[i])/(max[i]-min[i])*resolution[i]);
						err1vects[i]=round((err1vects[i])/(max[i]-min[i])*resolution[i]);
						err2vects[i]=round((err2vects[i])/(max[i]-min[i])*resolution[i]);
					}
				}
				for(k=1;k<=outnodes;k++) classval[k]=1.0;
				for (i=1;i<=innodes;i++)
				{
					j=0;
					k=1;
					if(vects[i] != MissingDat)
					{
						//NSP_edited 2013 oldj is no longer needed
						//oldj=(double)2*resolution[i];
						// NSP_Edited resolution[i]+1 replaced with resolution[i] removed
						while ((fabs(vects[i]-binloc[i][j+1]) >=1.0)&& (j<= resolution[i]))
						{
							//oldj=fabs(vects[i]-binloc[i][j+1]);
							j++;
						}
						//if(j>0)j--;
						// NSP_added 	if(fabs(vects[i]-binloc[i][j+1]) <=1){jx=0;} else{jx=-1;}
						if(fabs(vects[i]-binloc[i][j+1]) <= 1){jx=0;} else{jx=-1;}
						for (l=1;l<=innodes;l++)
						{
						  if(i!=l)
						  {
							m=0;
							k=1;
							if(vects[l] != MissingDat)
							{
								//NSP_edited 2013 oldj is no longer needed
								//oldj=(double)2*resolution[l];
								// NSP_Edited resolution[l]+1 replaced with resolution[l] removed
								//NSP_---> forgot to replace oldj in next line (Feb 2013)
								while ((fabs(vects[l]-binloc[l][m+1]) >=1.0)&& (m<= resolution[l]))
								{
									//oldj=fabs(vects[l]-binloc[l][m+1]);
									m++;
								}
								//if(m>0)m--;
								for (k=1;k<=outnodes;k++)
								{
									//NSP_added conditional tmp2_wts
									if(jx==0)tmp2_wts=(double)arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+k)];
									else tmp2_wts=1.0/outnodes;
									if(nerror ==2)
									{
										for(p=(m-(int)err1vects[l]);p<=(m+(int)err2vects[l]);p++)
										{
											if(p<0) p=0; if(p>resolution[l]) break;
											//NSP replaced tmp2_wts by (double)[i][j][l][m][k]
											if ((double)arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+p*(kk+1)+k)] > (double)arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+k)])
											m=p;
									    }
									}
									if(nerror ==1)
									{
										for(p=(m-(int)err1vects[l]);p<=(m+(int)err1vects[l]);p++)
										{
											if(p<0) p=0; if(p>resolution[l]) break;
											//NSP replaced tmp2_wts by (double)arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+k)]
											if ((double)arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+p*(kk+1)+k)] > (double)arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+k)])
											m=p;
									    }
									}
									if(arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+0)] > 0)
									{
										//NSP_added conditional tmp2_wts
										if(jx==0)
										tmp2_wts=(double)arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+k)]*arr_anti_wts[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+k)]/(arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+0)]);
										else
										tmp2_wts=1.0/outnodes;
									}
									else
										tmp2_wts= 1.0/outnodes;
									classval[k]*=(double)tmp2_wts;
								}
								totprob=0;
								for(k=1;k<=outnodes;k++) totprob+=classval[k];
								//NSP Replaced totprob=1; with totprob=innodes*outnodes;
								if (totprob==0) {totprob=innodes*outnodes; cout <<"Caution\n";}
								for(k=1;k<=outnodes;k++) classval[k]=classval[k]/totprob;
							 }
						  }
						}
					}
				}
				kmax=1;
				cmax=0;
				for (k=1;k<=outnodes;k++)
				{
					if (classval[k] > cmax)
					{
						cmax=classval[k];
						kmax=k;
					}
				}
				if ((fabs(dmyclass[kmax]-tmpv) >= dmyclass[0]) && (rnd >0))
				{
					for (i=1;i<=innodes;i++)
					{
						j=0;
						k=1;
						if(vects[i] != MissingDat)
						{
							//NSP_edited 2013 oldj is no longer needed
							//oldj=(double)2*resolution[i];
							// NSP_Edited resolution[i]+1 replaced with resolution[i] removed
							while ((fabs(vects[i]-binloc[i][j+1]) >=1.0)&& (j<= resolution[i]))
							{
								//oldj=fabs(vects[i]-binloc[i][j+1]);
								j++;
							}
							//if(j>0)j--;
							for (l=1;l<=innodes;l++)
							{
							if(i!=l)
							{
								m=0;
								k=1;
								if(vects[l] != MissingDat)
								{
									//NSP_edited 2013 oldj is no longer needed
									//oldj=(double)2*resolution[l];
									// NSP_Edited resolution[l]+1 replaced with resolution[l] removed
									while ((fabs(vects[l]-binloc[l][m+1]) >=1.0)&& (m<= resolution[l]))
									{
										//oldj=fabs(vects[l]-binloc[l][m+1]);
										m++;
									}
									//if(m>0)m--;
									while ((k<=outnodes)&&fabs(dmyclass[k]-tmpv) > dmyclass[0]) k++;
									if((classval[(int)kmax] >0)&&(classval[k]<classval[(int)kmax]))
									{
										arr_anti_wts[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+k)]+=(double)gain*(1.0-(classval[k]/classval[(int)kmax]));
									}
									if(arr_anti_wts[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+k)]<= 0.0)
										cout << k << " "<< tmpv << "[" << dmyclass[1] << "]" << dmyclass[outnodes] << "\n";
								}
							}
							}
						}
					}
				} // kmax che
			} // while not eof check
		 // Now save the wieights
			fclose(fl1);

			//if(rank!=size-1)
			if(rnd!=(oneround))
			{
				remaining=totsendreceivesize;
				tobesent=256;
				k=0;
				if(rank==size-1)
					ii=0;
				else
					ii=rank+1;

				while(remaining!=0)
				{
					if(remaining<256)
						tobesent=remaining;
					MPI_Send(&arr_anti_wts[k],tobesent,MPI_DOUBLE,ii,1,MPI_COMM_WORLD);
					k+=tobesent;
					remaining-=tobesent;
				}
			}

			strcpy(fltmp,datfilename);

			fl1=fopen(fltmp,"r");
			m=n;
			n=0;
			rslt=0.0;
			rslt2=0.0;
			pcnt=0;
			while (!feof(fl1))                    // Test round...
			{
				n++;
				kmax=1;
				cmax=0;
				if(ans1==3)
				{
					for (i=1;i<=innodes;i++)
					{
						fscanf(fl1,"%lf",&vects[i]);
						if(nerror ==2){fscanf(fl1,"%lf",&err1vects[i]);fscanf(fl1,"%lf",&err2vects[i]);}else
						if(nerror ==1){fscanf(fl1,"%lf",&err1vects[i]);err2vects[i]=err1vects[i];}
					}
					fscanf(fl1,"\n");
				}
				else
				{
					for (i=1;i<=innodes;i++)
					{
						fscanf(fl1,"%lf",&vects[i]);
						if(nerror ==2){fscanf(fl1,"%lf",&err1vects[i]);fscanf(fl1,"%lf",&err2vects[i]);}else
						if(nerror ==1){fscanf(fl1,"%lf",&err1vects[i]);err2vects[i]=err1vects[i];}
					}
					fscanf(fl1,"%lf\n",&tmpv);
				}
				for(i=1;i<=innodes;i++)
				{
					if((vects[i] != MissingDat)&&(max[i]!=MissingDat))
					{
						vects[i]=round((vects[i]-min[i])/(max[i]-min[i])*resolution[i]);
						err1vects[i]=round((err1vects[i])/(max[i]-min[i])*resolution[i]);
						err2vects[i]=round((err2vects[i])/(max[i]-min[i])*resolution[i]);
						if (vects[i] < 0) vects[i]=0;             // let us be bounded. #Oct 2001.
					}
				}
				for(k=1;k<=outnodes;k++) classval[k]=1.0;
				for (i=1;i<=innodes;i++)
				{
					j=0;
					k=1;
					if(vects[i] != MissingDat)
					{
						//NSP_edited 2013 oldj is no longer needed
						//oldj=(double)2*resolution[i];
						// NSP_Edited resolution[i]+1 replaced with resolution[i] removed
						while ((fabs(vects[i]-binloc[i][j+1]) >=1.0)&& (j<= resolution[i]))
						{
							//oldj=fabs(vects[i]-binloc[i][j+1]);
							j++;
						}
						//if(j>0)j--;
						//NSP_added if(fabs(vects[i]-binloc[i][j+1])<=1){jx=0;} else{jx=-1;}
						if(fabs(vects[i]-binloc[i][j]) <=1){jx=0;} else{jx=-1;}
						for (l=1;l<=innodes;l++)
						{
						 if (i !=l)
						 {
							m=0;
							k=1;
							if(vects[l] != MissingDat)
							{
								//NSP_edited 2013 oldj is no longer needed
								//oldj=(double)2*resolution[l];
								// NSP_Edited resolution[l]+1 replaced with resolution[l] removed
								while ((fabs(vects[l]-binloc[l][m+1]) >=1.0)&& (m<= resolution[l]))
								{
									//oldj=fabs(vects[l]-binloc[l][m+1]);
									m++;
								}
								//if(m>0)m--;
								for (k=1;k<=outnodes;k++)
								{
									//NSP_added conditional tmp2_wts
									if(jx==0)tmp2_wts=(double)arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+k)];
									else tmp2_wts=1.0/outnodes;
									if(nerror ==2)
									{
										for(p=(m-(int)err1vects[l]);p<=(m+(int)err2vects[l]);p++)
										{
											if(p<0) p=0; if(p>resolution[l]) break;
											//NSP replaced tmp2_wts by (double)[i][j][l][m][k]
											if ((double)arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+p*(kk+1)+k)]> (double)arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+k)])
											m=p;
										}
									}
									if(nerror ==1)
									{
										for(p=(m-(int)err1vects[l]);p<=(m+(int)err1vects[l]);p++)
										{
											if(p<0) p=0; if(p>resolution[l]) break;
											//NSP replaced tmp2_wts by (double)[i][j][l][m][k]
											if ((double)arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+p*(kk+1)+k)]> (double)arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+k)])
											m=p;
										}
									}
									if(arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+0)]> 0)
									{
										if(jx==0)
										tmp2_wts=(double)arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+k)]*arr_anti_wts[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+k)]/(arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+0)]);
										else
										tmp2_wts=1.0/outnodes;
									}
									else
										tmp2_wts=(double)1.0/outnodes;
									classval[k]*=(double)tmp2_wts;
								}
								totprob=0;
								for(k=1;k<=outnodes;k++) totprob+=classval[k];
								//NSP Replaced totprob=1; with totprob=innodes*outnodes;
								if (totprob==0) {totprob=innodes*outnodes; cout << "Caution!!\n";}
								for(k=1;k<=outnodes;k++) classval[k]=classval[k]/totprob;
							}
						 }
						}
					}
				}
				for (k=1;k<=outnodes;k++)
				{
					if (classval[k] > cmax)
					{
						cmax=classval[k];
						kmax=k;
					}
				}
				if (fabs(dmyclass[kmax]-tmpv) < dmyclass[0])
				{
					rslt2+=cmax;
					pcnt++;
				}
				else
				{
					k=1;
					while ((k<=outnodes)&&fabs(dmyclass[k]-tmpv) > dmyclass[0]) k++;
					rslt+=cmax-classval[k];
				}
			} // while not eof check
		// Now save the wieights
			fclose(fl1);
			kmax=1;
			//if(rank!=0)
			if(rnd!=0)
			{
				if(rank==master_rank)
					ii=size-1;
				else
					ii=rank-1;
				MPI_Recv(&orslt,1,MPI_DOUBLE,rank-1,1,MPI_COMM_WORLD,&status);
				MPI_Recv(&orslt2,1,MPI_DOUBLE,rank-1,2,MPI_COMM_WORLD,&status);
				MPI_Recv(&pocnt,1,MPI_INT,rank-1,3,MPI_COMM_WORLD,&status);
			}
			if(orslt2==0) orslt2=rslt2;
			if(orslt==0) orslt=rslt;
			prslt=(rslt2-orslt2);
			if(rslt > 0)
			nrslt=(orslt/rslt);
			send_status=0;
			if(pcnt>pocnt)
			{
				rnn=rnd;
				pocnt=pcnt;   // The best result is now saved in pocnt
				strcpy(fltmp,fln);
				strcat(fltmp,".awf");
				fl6=fopen(fltmp,"w+");
				kmax=1;
				for(k=1;k<=outnodes;k++)
				{
					for(i=1;i<=innodes;i++)
					for(j=0;j<=resolution[i];j++)
					{
						for(l=1;l<=innodes;l++)
						if(i!=l)
						{
							for(m=0;m<=resolution[l];m++)
							{
								fprintf(fl6,"%lf ",arr_anti_wts[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+k)]);
							}
							fprintf(fl6,"\n");
						}
					}
					fprintf(fl6,"\n");
				}
				fprintf(fl6,"\n");
				for(i=1;i<=innodes;i++)
				for(j=1;j<=resolution[i];j++)
				fprintf(fl6,"%lf\n", binloc[i][j]);                 /// Let us print the bins.
				fflush(fl6);
				fclose(fl6);
				cout << "Round:" << rnn << "| TProb["<<prslt<<"," <<nrslt<<"] | Passed count:" << pocnt << endl;
				if(orslt2 <rslt2) orslt2=rslt2;
				if(rslt < orslt) orslt=rslt;
			}
			n=m;
			/*if(rank!=size-1)*/
			if(rnd!=oneround)
			{
				if(rank!=size-1)
					ii=rank+1;
				else
					ii=0;
				MPI_Send(&orslt,1,MPI_DOUBLE,ii,1,MPI_COMM_WORLD);
				MPI_Send(&orslt2,1,MPI_DOUBLE,ii,2,MPI_COMM_WORLD);
				MPI_Send(&pocnt,1,MPI_INT,ii,3,MPI_COMM_WORLD);
			}
			for(i=0;i<totsize;i++)
				arr_anti_wts[i]=0;
		}  //rnd inc.
		MPI_Barrier(MPI_COMM_WORLD);
	}  // ans <> 1


/***********************************End of Case 1*******************************/
	MPI_Barrier(MPI_COMM_WORLD);
	//if(rank==master_rank)
	{
		strcpy(fltmp,datfilename);

		fl1=fopen(fltmp,"r");
		strcpy(fltmp,fln);
		strcat(fltmp,".awf");
		fl6=NULL;
		fl6=fopen(fltmp,"r");
		strcpy(fltmp,fln);
		strcat(fltmp,".apf");
		fl2=NULL;
		if((fl2=fopen(fltmp,"r"))!=NULL)
		{
			//cout << "Creating the Anticipated Network outputs\n";
			for (i=1;i<=innodes;i++)
			{
				fscanf(fl2,"%d",&resolution[i]);
				for(j=0;j<=resolution[i];j++) binloc[i][j+1]=j*1.0;
			}
			fscanf(fl2,"%lf",&omax);
			fscanf(fl2,"%lf",&omin);
			fscanf(fl2,"\n");
			for(i=1;i<=innodes;i++) fscanf(fl2,"%lf",&max[i]);
			fscanf(fl2,"\n");
			for(i=1;i<=innodes;i++) fscanf(fl2,"%lf",&min[i]);
			fscanf(fl2,"\n");
			for(i=1;i<=innodes;i++)for(j=0;j<=resolution[i];j++)
			for(l=1;l<=innodes;l++)for(m=0;m<=resolution[l];m++) arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+0)]=0;
			for(k=1;k<=outnodes;k++)
			{
				for(i=1;i<=innodes;i++)
				for(j=0;j<=resolution[i];j++)
				{
					for(l=1;l<=innodes;l++)
					if(i!=l)
					{
						for(m=0;m<=resolution[l];m++)
						{
							fscanf(fl2,"%d",&arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+k)]);
							fscanf(fl6,"%lf",&arr_anti_wts[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+k)]);
							arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+0)]+=(double)(arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+k)]);
						}
						fscanf(fl2,"\n");
						fscanf(fl6,"\n");
					}
				}
				fscanf(fl2,"\n");
				fscanf(fl6,"\n");
			}
		}
		else
		{
			cout << "Unable to Open the APF information file";
			exit(1);
		}
		for(i=1;i<=innodes;i++)
		for(j=1;j<=resolution[i];j++)
		{
			fscanf(fl6,"%lf\n",&binloc[i][j]);                 /// Let us print the bins.
		}
		fclose(fl6);

		//cout << "Read all input parameters\n";
	// *********** case 3 ***********************************************
		if(rank==master_rank)
		{
			fl4=fopen("output.dat","w+");  // Network Output values
			if (ans1 !=3)
			{
				fl5=fopen("actual.dat","w+");  // Expected Output Values
				strcpy(fltmp,fln);
				strcat(fltmp,argp[2]);
				strcpy(fltmp,fln);
				strcat(fltmp,argp[2]);
				strcat(fltmp,".cmp");         // Lets see how well the classification went.
				fl7=fopen(fltmp,"w+");
				fprintf(fl7,"Sample         Predicted     Actual            Prediction \n");
				fprintf(fl7," No.       Ist 2nd  3rd  4th  item             Confidence\n");

			}
		}
		c1cnt=0;
		c2cnt=0;
		invcnt=0;
		n=0;
	 // Create classtot values ***********************
		while (!feof(fl1))
		{
			n++;
			cmax= 0.0;
			c2max=0.0;
			c3max=0.0;
			c4max=0.0;
			kmax=0;
			k2max=0;
			k3max=0;
			k4max=0;
			classval[0]=0.0;
			if(ans1==3)
			{
				for (i=1;i<=innodes;i++)
				{
					fscanf(fl1,"%lf",&vects[i]);
					if(nerror ==2){fscanf(fl1,"%lf",&err1vects[i]);fscanf(fl1,"%lf",&err2vects[i]);}else
					if(nerror ==1){fscanf(fl1,"%lf",&err1vects[i]);err2vects[i]=err1vects[i];}
				}
				fscanf(fl1,"\n");
			}
			else
			{
				for (i=1;i<=innodes;i++)
				{
					fscanf(fl1,"%lf",&vects[i]);
					if(nerror ==2){fscanf(fl1,"%lf",&err1vects[i]);fscanf(fl1,"%lf",&err2vects[i]);}else
					if(nerror ==1){fscanf(fl1,"%lf",&err1vects[i]);err2vects[i]=err1vects[i];}
				}
				fscanf(fl1,"%lf\n",&tmpv);
			}
			skpchk=0;
			for(i=1;i<=innodes;i++)
			{
				vectso[i]=vects[i];
				if((((max[i]-min[i]) >0)&& (vects[i] !=MissingDat))&&(max[i] !=MissingDat))
				{
					vects[i]=round(((vects[i]-min[i])/(max[i]-min[i]))*resolution[i]);
					err1vects[i]=round((err1vects[i])/(max[i]-min[i])*resolution[i]);
					err2vects[i]=round((err2vects[i])/(max[i]-min[i])*resolution[i]);
					skpchk=0;
				}
				else
				skpchk=1;
			}
			for(k=1;k<=outnodes;k++) classval[k]=1.0;

			jc=0,jpt=0;
			jobsperthread=innodes/size;
			jpt=innodes/size;
			exjpt=0;
			if(innodes%size!=0)
			{
				if(rank<innodes%size)
				{
					jobsperthread+=1;
					exjpt=rank;
				}
				else
					exjpt=innodes%size;
			}
			for (jc=0;jc<jobsperthread;jc++)
			{
				i=(jc*size)+rank+1;
				if(i!=1)
				{

					if(rank==0)
					{
						MPI_Recv(&classval,outnodes+2,MPI_DOUBLE,size-1,1,MPI_COMM_WORLD,&status);

					}
					else
					{
						MPI_Recv(&classval,outnodes+2,MPI_DOUBLE,rank-1,1,MPI_COMM_WORLD,&status);
					}
				}
			//for (i=1;i<=innodes;i++)                     // NSP_Feb13_2013 : Let each node be a node on the cluster.
			//{
				j=0;
				k=1;
				if(vects[i]==MissingDat)
					skpchk=1;
				else
					skpchk=0;
				if ((resolution[i] >= vects[i]) &&(skpchk==0))
				{
					while ((fabs(vects[i]-binloc[i][j+1]) >=1.0)&& (j<= resolution[i]))
					{
						j++;
					}
					jx=0;
				}
				  else
				  {
				  jx=-1;
				  }

				for (l=1;l<=innodes;l++)
				{
					 if(i!=l)
					 {
						m=0;
						k=1;
						if((vects[l]==MissingDat)||(vects[i]==MissingDat))
							skpchk=1;
						else
							skpchk=0;
						if ((resolution[l] >= vects[l]) &&(skpchk==0))
						{
							while ((fabs(vects[l]-binloc[l][m+1]) >=1.0)&& (m<= resolution[l]))
							{
								m++;
							}
						}
						for (k=1;k<=outnodes;k++)
						{
							if(jx==0) {tmp2_wts=(double)arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+k)];} else{tmp2_wts=1.0/outnodes;}
							if(nerror ==2)
							{
								for(p=(m-(int)err1vects[l]);p<=(m+(int)err2vects[l]);p++)
								{
									if(p<0) p=0; if(p>resolution[l]) break;
									if ((double)arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+p*(kk+1)+k)] > (double)arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+k)])
									m=p;
								}
							}
							if(nerror ==1)
							{
								for(p=(m-(int)err1vects[l]);p<=(m+(int)err1vects[l]);p++)
								{
									if(p<0) p=0; if(p>resolution[l]) break;
									if ((double)arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+p*(kk+1)+k)]> (double)arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+k)])
									m=p;
								}
							}
							if((arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+0)] > 0) && (resolution[i]>= vects[i])&& (resolution[l]>= vects[l])&&(skpchk==0))
							{
							 if(jx==0){tmp2_wts=(double)arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+k)]*arr_anti_wts[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+k)]*1.0/(arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+0)]);}
							 else{tmp2_wts=1.0/outnodes;}
							}
							else
							if(skpchk == 1)   // || bcchoice == 'y')
							{
								tmp2_wts= 1.0; //(double)1.0/outnodes; //1.0; //
							}
							else
							{
								tmp2_wts=(double)1.0/outnodes;
							}
							if((resolution[i] >= vects[i])&& (resolution[l]>= vects[l])&&(skpchk==0))
							{
								classval[k]*=(double)tmp2_wts;
							}
						}
						totprob=0;
						for(k=1;k<=outnodes;k++) totprob+=classval[k];
						if (totprob==0) {totprob=innodes*outnodes; cout<<"["<< i<<j<<l<<m<<"]"<<arr_anti_net[(i*(jk+1)*(lk+1)*(mk+1)*(kk+1)+j*(lk+1)*(mk+1)*(kk+1)+l*(mk+1)*(kk+1)+m*(kk+1)+k)]<<"Caution!!\n";}
						for(k=1;k<=outnodes;k++) classval[k]=classval[k]/totprob;
					 }
				}
				if(i!=innodes)
				{
					if(rank==size-1)
					{
						MPI_Send(classval,outnodes+2,MPI_DOUBLE,0,1,MPI_COMM_WORLD);

					}
					else
					{
						MPI_Send(classval,outnodes+2,MPI_DOUBLE,rank+1,1,MPI_COMM_WORLD);
					}
				}
				else
				{
					MPI_Send(classval,outnodes+2,MPI_DOUBLE,0,2,MPI_COMM_WORLD);

				}
			}                                // NSP_Feb13_2013 : Each node returns classval[k] for one feature.
			MPI_Barrier(MPI_COMM_WORLD);
			if(rank==master_rank)
			{


				MPI_Recv(&classval,outnodes+2,MPI_DOUBLE,MPI_ANY_SOURCE,2,MPI_COMM_WORLD,&status);
				cmax=0.0;
				c2max=0.0;
				c3max=0.0;
				k3max=0.0;
				kmax=0.0;
				k2max=0.0;
				totprob=0.0;
				for (k=1;k<=outnodes;k++)
				{
					if (classval[k] > cmax)
					{
						c4max=c3max;
						k4max=k3max;
						c3max=c2max;
						k3max=k2max;
						c2max=classval[kmax];
						k2max=kmax;
						cmax=classval[k];
						kmax=k;
					}
					else
					if (classval[k]>c2max)
					{
						c4max=c3max;
						k4max=k3max;
						c3max=c2max;
						k3max=k2max;
						c2max=classval[k];
						k2max=k;
					}
				   else
					if (classval[k]>c3max)
					{
						c4max=c3max;
						k4max=k3max;
						c3max=classval[k];
						k4max=k;
					}
				   else
					if (classval[k]>c4max)
					{
						c4max=classval[k];
						k4max=k;
					}
					totprob += (double)classval[k];
				}
				if(totprob <=0.0) totprob=innodes*outnodes;
				if(ans1 ==3)
				{
					if (dmyclass[(int)kmax]- (int)dmyclass[(int)kmax] ==0.0)
					{
						fprintf(fl4,"%d  %d %-5.2f %d %-5.2f %d %-5.2f %d %-5.2f",n, (int)dmyclass[(int)kmax],100.0*((classval[kmax])/totprob),(int)dmyclass[(int)k2max],100.0*((classval[k2max])/totprob),(int)dmyclass[(int)k3max],100.0*((classval[k3max])/totprob),(int)dmyclass[(int)k4max],100.0*((classval[k4max])/totprob));
					}
					else
					{
						fprintf(fl4,"%d  %lf %-5.2f %lf %-5.2f %lf %-5.2f %lf %-5.2f",n, dmyclass[(int)kmax],100.0*((classval[kmax])/totprob),dmyclass[(int)k2max],100.0*((classval[k2max])/totprob),dmyclass[(int)k3max],100.0*((classval[k3max])/totprob),dmyclass[(int)k4max],100.0*((classval[k4max])/totprob));
					}
					if((fabs(classval[kmax]-classval[k2max]))<0.01*classval[kmax]) //classval[kmax])
					{
						nLoC+=classval[kmax]/totprob;
						nLoCcnt++;
						if(classval[kmax]>totprob*LoC)    //LoC)
						{
							fprintf(fl4, " <-- Either of it");
						}
						else
						{
							fprintf(fl4, " <-- Rejected");
						}
					}
					else
					{
						if(classval[kmax]>totprob*LoC)    //LoC)
						{
							fprintf(fl4, " <-- confident");
						}
						else
						{
							fprintf(fl4, " <-- Rejected");
						}
					}
					fprintf(fl4,"\n");
				}
				if(ans1 !=3)
				{
					if (dmyclass[(int)kmax]- (int)dmyclass[(int)kmax] ==0.0)
					{
						fprintf(fl4,"%d  %d\n",n, (int)dmyclass[(int)kmax]);
						fprintf(fl7, "%-8d    %d   %d     %d  %d     %d   ",n,(int)dmyclass[(int)kmax],(int)dmyclass[(int)k2max],(int)dmyclass[(int)k3max],(int)dmyclass[(int)k4max],(int)tmpv);
					}
					else
					{
						fprintf(fl4,"%d  %lf\n",n, dmyclass[(int)kmax]);
						fprintf(fl7, "%-8d    %lf   %lf     %lf    %lf     %lf    ",n,dmyclass[(int)kmax],dmyclass[(int)k2max],dmyclass[(int)k3max],dmyclass[(int)k4max],tmpv);
					}
					if(fabs(dmyclass[kmax]-tmpv) >= dmyclass[0])
					{
						if (classval[kmax]==0.0)
						{
							invcnt++;
							fprintf(fl7, "%-5.2f %% <-Out of range %-5.2f %% \n",100.0*((classval[kmax])/totprob),100.0*((classval[k2max])/totprob));
						}
						else
						{
							if (fabs(dmyclass[k2max]-tmpv) < dmyclass[0])
							{
								// The Next line defines the margin level required to pick up less confident examples!
								if((fabs(classval[kmax]-classval[k2max]))<0.01*classval[k2max]) //classval[kmax])
								{
									nLoC+=classval[kmax]/totprob;
									nLoCcnt++;
									if (classval[kmax]>totprob*LoC) // LoC)
									{
										c2cnt++;  // No more differences. NSP (OCT 2001)
										fprintf(fl7, "%-5.2f %% <-F(1)P(2) %-5.2f %%  %-5.2f %% %-5.2f %% \n",100.0*((classval[kmax])/totprob),100.0*((classval[k2max])/totprob),100.0*((classval[k3max])/totprob),100.0*((classval[k4max])/totprob));
									}
									else
									{
										fprintf(fl7, "%-5.2f %%  <-FMC %-5.2f %% %-5.2f %% %-5.2f %% \n",100.0*((classval[kmax])/totprob),100.0*((classval[k2max])/totprob),100.0*((classval[k3max])/totprob),100.0*((classval[k4max])/totprob));
										invcnt++;
									}
								}
								else
								{
									if (classval[kmax]>totprob*LoC) // LoC)
									{
										fprintf(fl7, "%-5.2f %% <-Failed %-5.2f %%  %-5.2f %%  %-5.2f %% \n",100.0*((classval[kmax])/totprob),100.0*((classval[k2max])/totprob),100.0*((classval[k3max])/totprob),100.0*((classval[k4max])/totprob));
									}
									else
									{
										fprintf(fl7, "%-5.2f %% <-FMC %-5.2f %%  %-5.2f %%  %-5.2f %% \n",100.0*((classval[kmax])/totprob),100.0*((classval[k2max])/totprob),100.0*((classval[k3max])/totprob),100.0*((classval[k4max])/totprob));
										invcnt++;
									}
								}
							}
							else
							{
								if (classval[kmax]>totprob*LoC) // LoC)
								{
									fprintf(fl7, "%-5.2f %% <-Failed %-5.2f %%  %-5.2f %% %-5.2f %% \n",100.0*((classval[kmax])/totprob),100.0*((classval[k2max])/totprob),100.0*((classval[k3max])/totprob),100.0*((classval[k4max])/totprob));
								}
								else
								{
									fprintf(fl7, "%-5.2f %% <-FMC %-5.2f %% %-5.2f %% %-5.2f %% \n",100.0*((classval[kmax])/totprob),100.0*((classval[k2max])/totprob),100.0*((classval[k3max])/totprob),100.0*((classval[k4max])/totprob));
									invcnt++;
								}
							}
						}
					}
					else
					{
						if((fabs(classval[kmax]-classval[k2max]))<0.01*classval[kmax])
						{
							nLoC+=classval[kmax]/totprob;
							nLoCcnt++;
							if (classval[kmax]>totprob*LoC) // LoC)
							{
								fprintf(fl7, "%-5.2f %% <-P(1)F(2) %-5.2f %% %-5.2f %% %-5.2f %% \n",100.0*((classval[kmax])/totprob),100.0*((classval[k2max])/totprob),100.0*((classval[k3max])/totprob),100.0*((classval[k4max])/totprob));
								c1cnt++;
							}
							else
							{
								invcnt++;
								fprintf(fl7, "%-5.2f %% <-PMC %-5.2f %% %-5.2f %% %-5.2f %% \n",100.0*((classval[kmax])/totprob),100.0*((classval[k2max])/totprob),100.0*((classval[k3max])/totprob),100.0*((classval[k4max])/totprob));
							}
						}
						else
						{
							if (classval[kmax]>totprob*LoC) // LoC)
							{
								fprintf(fl7, "%-5.2f %% <-Passed %-5.2f %% %-5.2f %% %-5.2f %% \n",100.0*((classval[kmax])/totprob),100.0*((classval[k2max])/totprob),100.0*((classval[k3max])/totprob),100.0*((classval[k4max])/totprob));
								c1cnt++;
							}
							else
							{
								invcnt++;
								fprintf(fl7, "%-5.2f %% <-PMC %-5.2f %% %-5.2f %% %-5.2f %% \n",100.0*((classval[kmax])/totprob),100.0*((classval[k2max])/totprob),100.0*((classval[k3max])/totprob),100.0*((classval[k4max])/totprob));
							}
						}
					}
					fprintf(fl5,"%d %e \n",n,(double) tmpv);
				} // ans1 != 3 ends here ******************
			}


		}
		if(rank==master_rank)
		{
			cout << "The suggested LoC is " << nLoC/nLoCcnt << "\n";
			fclose(fl1);
			fclose(fl2);
			fclose(fl4);
			if(ans1 < 3)
			{
				strcpy(fltmp,fln);
				fclose(fl5);
				fprintf(fl7,"*________________________________________________________________________\n");
				fprintf(fl7,"*Total    Success in   Success in   Non classified   Real success in    \n");
				cout << "*________________________________________________________________________\n";
				cout << "*Total    Success in   Success in   Non classified   Real success in    \n";
				if (outnodes > 3)
				{
					fprintf(fl7,"* No.    Ist Choice  2nd Choice     items           two chances    \n");
					fprintf(fl7,"* %d       %d          %d           %d             %-5.2f %% \n",n,c1cnt,c2cnt,invcnt,(double)100.0*(c1cnt+c2cnt)/(n-invcnt));
					cout << "* No.    Ist Choice  2nd Choice     items           two chances    \n";
					printf("* %d       %d          %d           %d             %-5.2f %% \n",n,c1cnt,c2cnt,invcnt,(double)100.0*(c1cnt+c2cnt)/(n-invcnt));
				}
				else
				{
					fprintf(fl7,"* No.    Ist Choice  2nd Choice     items           First chance    \n");
					fprintf(fl7,"* %d       %d          %d           %d             %-5.2f %% \n",n,c1cnt,c2cnt,invcnt,(double)100.0*(c1cnt)/(n-invcnt));
					cout << "* No.    Ist Choice  2nd Choice     items           First chance    \n";
					printf("* %d       %d          %d           %d             %-5.2f %% \n",n,c1cnt,c2cnt,invcnt,(double)100.0*(c1cnt)/(n-invcnt));
				}
				fprintf(fl7,"*________________________________________________________________________\n");
				printf("*________________________________________________________________________\n");
				fclose(fl7);
			} // ******** ans1!=3 ends here *************
			cout << "Done.\n";
			stop = times(NULL);
			//printf("Anti net Time = %f\n", (double)(MPI_Wtime()-anit_net));
		}
	}//end of master rank if
	if(rank==master_rank)
	{
		end_t = MPI_Wtime();
		printf("Execution Time = %f\n", (double)(end_t-start_t));
	}
	/*free(anti_net_temp);*/
	free(arr_anti_net);
	free(arr_anti_wts);
	free(arr_anti_wts_temp);
	MPI_Finalize();
} //end main



