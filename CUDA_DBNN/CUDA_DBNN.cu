/* CUDA version of the DBNN code for classification of stars, galaxies. The code was originally written by Prof. Sajeeth

Author: Ajay Vibhute

*/


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>
using namespace std;
#include <stdlib.h>
#include<sys/times.h> // times() fun. is here.
#include <time.h>
#include <vector>
#define classes 500
#define max_resol 1600
#define features 100
#include"kernel.cu"
using std::vector;
/**************************
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

#include "vsipc.h"
***************************/

//#define oneround 100       //  Memory size.
//#define fst_gain 1.0 Moved to be a floating variable
static float bgain,gain,dmyclass[classes+2],classval[classes+2],cmax,c2max,c3max,c4max,tmp2_wts,totprob,oldj;
static float LoC=0.65;
static float nLoC=0.0;
static int jx=0, resol=100,nresol=0,nerror=0,nLoCcnt=1,skpchk=0,MissingDat=-9999;
static float omax,omin,rslt,rslt2,orslt,orslt2,prslt,nrslt,fst_gain;
clock_t start,stop;
static int argfnd, oneround=100,kmax,k2max,k3max,k4max,ans1,tcnt,rnn,rnd,i,j,k,l,n,m,p,c1cnt,c2cnt,pcnt,pocnt,invcnt,innodes=100,outnodes=100;
char fln[256],fltmp[256],urchoice,urchoicex,bcchoice,savedpar;
FILE *fl1,*fl2,*fl3,*fl4,*fl5,*fl6,*fl7,*fl8,*fl9,*fl10;
int main(int argv, char *argp[256])
/*
 Important note for revision in ver. 4.02
 Compute the bin center of gravities and the Kernel fit that holds the probability to
 find them. The slope of it should be used as bgain for each bin.
*/
/*
 You can now run dbnn in automated mode by specifying the parameters in 0.par and 1.par
 files. Also dbnn can now use the bin values from the saved apf file.
*/
{
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
					fscanf(fl2,"%f",&gain);
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
    strcpy(fltmp,fln);
    strcat(fltmp,".dat");
/*
  The structure of the data file is:
  Feature1 Feature2 Feature3 ....(etc upto innodes) ActualClass
  Feature1 Feature2 Feature3 ....(etc upto innodes) ActualClass
  Feature1 Feature2 Feature3 ....(etc upto innodes) ActualClass
  The delimiters are spaces and not tabs!!
  ActualClass should be a numeric > 0
*/
    if((fl1=fopen(fltmp,"r"))!=NULL)
    {
		strcpy(fltmp,fln);
	    strcat(fltmp,".inf");
/*
  The format of the info file is: (in each line enter)
  innodes
  outnodes
  margin   <- This addition is required for regression problems.
  1.0       <- You can give any real positive value here. It is just a label.
  2.0
  ... (etc. upto no of classes)
  0.65 <- The Margin or Line of Control for marginal values.
  100 <- By default, the maximum bin size is set to 100. You can change this if required.
  0,1,2 <- no error bars, uniform error bar, upper lower separate error values per entry.
*/	
     	if((fl2=fopen(fltmp,"r"))!=NULL)
	    {
			i=0;
			fscanf(fl2,"%d",&innodes);
			fscanf(fl2,"%d",&outnodes);
			for (i=0;i<=outnodes;i++) // dmyclass[0] contains margin others are expected values.
			fscanf(fl2,"%f",&dmyclass[i]);
			fscanf(fl2,"%f",&LoC);   // New parameter to specify the Line Of Control
			fscanf(fl2,"%d",&nresol);
			fscanf(fl2,"%d",&nerror);
			cout <<"You have "<< innodes << " input nodes and " << outnodes <<" Output nodes with " << "margin set to " << LoC << " and error levels set to "<< nerror <<"\n";
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
			fst_gain*=1.0/outnodes;
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
    cout << "Going to initialise the arrays\n";
 /**************** Let us Define the Network Structure *********************************/
//float mask_disp_maxres; // Space to save max resol for normalisation of mask_dist
    strcpy(fltmp,fln);
    strcat(fltmp,".dat");

	int numlines=getNumlines(fltmp);
	printf("NUMLINES:%d\n",numlines);

	float vectso[innodes+outnodes+2],tmpv,max[innodes+2],min[innodes+2],vects[innodes+outnodes+2];
	float err1vects[innodes+2], err2vects[innodes+2];
	//float arr_vects[numlines][innodes+outnodes+2];
	float *arr_tmpv=(float*)malloc(sizeof(float)*(numlines+2));
	float *arr_vects=(float*)malloc(sizeof(float)*(numlines+2)*(innodes+outnodes+2));
	float *arr_err1vects=(float*)malloc(sizeof(float)*(numlines+2)*(innodes+2));
	float *arr_err2vects=(float*)malloc(sizeof(float)*(numlines+2)*(innodes+2));
	int totsize=(innodes+2)*(resol+2)*(innodes+2)*(resol+4)*(outnodes+2);
	int totsendreceivesize=(innodes+1)*(resol+2)*(innodes+1)*(resol+1)*(outnodes+1);
	float *arr_anti_wts=(float*) malloc(totsize*sizeof(float));
	int *arr_anti_net=(int*)malloc(sizeof(int)*totsize);
	int ik=innodes+1,jk=resol+1,lk=innodes+1,mk=resol+1,kk=outnodes+1;
	int resolution[innodes+8];
	float classtot[innodes+2][resol+2];           // Total Prob. computed
	if(classtot==NULL){cout << "Out of Memory to Run Code at classtot.. Exiting\n";exit(1);}
	//float binloc[innodes+4][resol+8];
	float *arr_binloc=(float*)malloc(sizeof(float)*(innodes+4)*(resol+8));
	int rn=resol+1;
	int iin=innodes+4;





  /***************************Let us put up the Network***********************************/
//    Start the counter for case 2 here.................
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
			fscanf(fl1,"%f",&vects[i]); 
			if(nerror ==2){fscanf(fl1,"%f",&err1vects[i]);fscanf(fl1,"%f",&err2vects[i]);}else
			if(nerror ==1){fscanf(fl1,"%f",&err1vects[i]); err2vects[i]=err1vects[i];} 
			if(vects[i] != MissingDat)
				{
					min[i]=vects[i];
					max[i]=vects[i];
				}
				else max[i]=MissingDat;
			}
			else
			{
				fscanf(fl1,"%f",&vects[i]);
				if(vects[i] != MissingDat)
				{
					if( vects[i]> max[i]) max[i]=vects[i];
					if (min[i] > vects[i]) min[i]=vects[i];
				}
			}
			fscanf(fl1,"%f\n",&tmpv);
			if(tmpv>omax) omax = tmpv;
			if(tmpv<omin) omin =tmpv;
			k=1;
			j=1;
			n++;
		}
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
				cout << "Reading from the saved information\n";
				for (i=1;i<=innodes;i++)
				{
					fscanf(fl2,"%d",&resolution[i]);
					for(j=0;j<=resolution[i];j++) arr_binloc[(i*rn)+(j+1)]=j*1.0;
				}
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
			for(j=0;j<=resolution[i];j++) arr_binloc[(i*rn)+(j+1)]=j*1.0;

		}
		for(k=1;k<=outnodes;k++)
		for(i=1;i<=innodes;i++)
		for(j=0;j<=resolution[i];j++)
		for(l=1;l<=innodes;l++)
		for(m=0;m<=resolution[l];m++)
		{
			//anti_net[i][j][l][m][k]=1;
			//anti_wts[i][j][l][m][k]=(float)(1.0);
			arr_anti_wts[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(m*kk)+k)]=(double)(1.0);
			arr_anti_net[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(m*kk)+k)]=1;
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
				fscanf(fl1,"%f",&vects[i]); 
				if(nerror ==2){fscanf(fl1,"%f",&err1vects[i]);fscanf(fl1,"%f",&err2vects[i]);}else
				if(nerror ==1){fscanf(fl1,"%f",&err1vects[i]);err2vects[i]=err1vects[i];} 
			}
			fscanf(fl1,"%f\n",&tmpv);
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
			for (i=1;i<=innodes;i++)
			{
				j=0;
				if(vects[i] != MissingDat)
				{
//					oldj=(float)2*resolution[i];

					while ((fabs(vects[i]-arr_binloc[(i*rn)+(j+1)]) >=1.0 )&& (j<= resolution[i]))
					{
	//					oldj=fabs(vects[i]-binloc[i][j+1]);
						j++;
					}
					for (l=1;l<=innodes;l++)
					{
						m=0;
						if(i!=l)
						{
//							oldj=(float)2*resolution[l];
							while ((fabs(vects[l]-arr_binloc[(l*rn)+(m+1)]) >=1.0)&& (m<= resolution[l]))
							{
	//							oldj=fabs(vects[l]-binloc[l][m+1]);
								m++;
							}
							k=1;
							while ((k<=outnodes)&&(fabs(tmpv - dmyclass[k])) > dmyclass[0]) k++;
							//(anti_net[i][j][l][m][k])++;
							//(anti_net[i][j][l][m][0])++;
							(arr_anti_net[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(m*kk)+k)])++;
							(arr_anti_net[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(m*kk)+0)])++;


						}

					}
				}
			}
		}//end of while
		fclose(fl1);
		fclose(fl2);
		stop = times(NULL);
		cout << "The computation took " << fabs(start - stop)*10000/(CLOCKS_PER_SEC) << " Secs.\n";
         /*
            The conditional Probability,
	    P(A|B) = P(A intersection B)/P(B) is the
	    probability for the occurance of A(k) if B(ij) has happened =
	    Share of B(ij) that is held by A(k) / Probability of total B(ij)
	    in that particular feature i with resolution j.

                      */
		strcpy(fltmp,fln);
		strcat(fltmp,".awf");      // This file holds the weights
		fl6=fopen(fltmp,"w+");
		strcpy(fltmp,fln);
		strcat(fltmp,".apf");     // This file holds the estimated probability
		if((fl1=fopen(fltmp,"w+"))!=NULL)
		{
			for(i=1;i<=innodes;i++) fprintf(fl1,"%d ",resolution[i]);
			fprintf(fl1,"\n%f %f \n",omax,omin);
			for(i=1;i<=innodes;i++) fprintf(fl1,"%f ",max[i]);
			fprintf(fl1,"\n");
			for(i=1;i<=innodes;i++) fprintf(fl1,"%f ",min[i]);
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
							//fprintf(fl1,"%d ",anti_net[i][j][l][m][k]);
							//fprintf(fl6,"%f ",(float)anti_wts[i][j][l][m][k]);

							fprintf(fl1,"%d ",arr_anti_net[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(m*kk)+k)]);
							fprintf(fl6,"%f ",(float)arr_anti_wts[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(m*kk)+k)]);
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
		fprintf(fl6,"%f\n", (float)arr_binloc[(i*rn)+(j)]);                 /// Let us print the bins.
		fclose(fl1);
		fclose(fl6);
		fflush(NULL);
		cout << "Creating the Anticipated Weights data file\n";
	}
/**********************************End of Case 0 ******************************/












	if(ans1==1)
	{
		start = times(NULL);
		pcnt=0;
		pocnt=0;
		rslt=0.0;
		rslt2=0.0;
		orslt=rslt;
		orslt2=rslt2;

		for(i=0;i<totsize;i++)
			arr_anti_wts[i]=0;
		cout << "The programe will now modify the compensatory weights\n";
		if(argfnd==0)
		{
			cout << "Please enter the gain:";
			cin >> gain;
			cout << "Please enter the number of training epochs:";
			cin >> oneround;
		}
		// Start the counter in this round here...................
		start = times(NULL);
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
					for(j=0;j<=resolution[i];j++) arr_binloc[(i*rn)+(j+1)]=j*1.0;
				}

				fscanf(fl2,"\n%f",&omax);
				fscanf(fl2,"%f",&omin);
				fscanf(fl2,"\n");
				for(i=1;i<=innodes;i++) fscanf(fl2,"%f",&max[i]);
				fscanf(fl2,"\n");
				for(i=1;i<=innodes;i++) fscanf(fl2,"%f",&min[i]);
				fscanf(fl2,"\n");
				for(i=1;i<=innodes;i++)for(j=0;j<=resolution[i];j++)
				//for(l=1;l<=innodes;l++)for(m=0;m<=resolution[l];m++) anti_net[i][j][l][m][0] =0;
				for(l=1;l<=innodes;l++)for(m=0;m<=resolution[l];m++) arr_anti_net[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(m*kk)+0)] =0;
				int ijk=0;
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
								ijk++;
								//fscanf(fl2,"%d",&anti_net[i][j][l][m][k]);
								//anti_net[i][j][l][m][0]+=anti_net[i][j][l][m][k];
								//fscanf(fl6,"%f",&anti_wts[i][j][l][m][k]);
								fscanf(fl2,"%d",&arr_anti_net[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(m*kk)+k)]);
								arr_anti_net[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(m*kk)+0)]+=arr_anti_net[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(m*kk)+k)];
								fscanf(fl6,"%f",&arr_anti_wts[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(m*kk)+k)]);
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
					fscanf(fl6,"%f\n", &arr_binloc[(i*rn)+(j)]);                 /// Let us print the bins.

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
		/*GPU Memory allocation*/

		int *d_arr_anti_net,*d_resolution;
		float *d_arr_anti_wts,*d_arr_tmpv,*d_arr_vects,*d_arr_err1vects,*d_arr_err2vects,*d_min,*d_max,*d_arr_binloc,*d_dmyclass;
		float *d_rslt,*d_rslt2,*tmp_rslt,*tmp_rslt2;
		int *d_pcnt,*tmp_pcnt,chunksize=0;
		cudaError_t status ;
		tmp_pcnt=(int*)malloc(sizeof(int)*numlines);
		tmp_rslt=(float*)malloc(sizeof(float)*numlines);
		tmp_rslt2=(float*)malloc(sizeof(float)*numlines);
		for(i=0;i<numlines;i++)
			tmp_pcnt[i]=0;
		
		//allocate memory on GPU
		status=cudaMalloc((void **)&d_pcnt,sizeof(int)*(numlines+10));
		status=cudaMemset(d_pcnt, 0, sizeof(int)*(numlines+10));
		status=cudaMalloc((void **)&d_rslt,sizeof(float)*(numlines+10));
		status=cudaMemset(d_rslt, 0, sizeof(float)*(numlines+10));
		status=cudaMalloc((void **)&d_rslt2,sizeof(float)*(numlines+10));
		status=cudaMemset(d_rslt2, 0, sizeof(float)*(numlines+10));

		status=cudaMalloc((void **)&d_arr_anti_net,sizeof(int)*totsize);
		status=cudaMemset(d_arr_anti_net, 0, sizeof(int)*totsize);
		status=cudaMalloc((void **)&d_resolution,sizeof(int)*(innodes+8));
		status=cudaMemset(d_resolution, 0, sizeof(int)*(innodes+8));

		status=cudaMalloc((void **)&d_arr_anti_wts,sizeof(float)*totsize);
		status=cudaMemset(d_arr_anti_wts, 0, sizeof(int)*totsize);

		status=cudaMalloc((void **)&d_arr_tmpv,sizeof(float)*(numlines+2));
		status=cudaMemset(d_arr_tmpv, 0, sizeof(float)*(numlines+2));

		status=cudaMalloc((void **)&d_arr_vects,sizeof(float)*(numlines+2)*(innodes+outnodes+2));
		status=cudaMemset(d_arr_vects, 0, sizeof(float)*(numlines+2)*(innodes+outnodes+2));

		status=cudaMalloc((void **)&d_arr_err1vects,sizeof(float)*(numlines+2)*(innodes+2));
		status=cudaMemset(d_arr_err1vects, 0, sizeof(float)*(numlines+2)*(innodes+2));

		status=cudaMalloc((void **)&d_arr_err2vects,sizeof(float)*(numlines+2)*(innodes+2));
		status=cudaMemset(d_arr_err2vects, 0, sizeof(float)*(numlines+2)*(innodes+2));

		status=cudaMalloc((void **)&d_min,sizeof(float)*(innodes+2));
		status=cudaMemset(d_min, 0, sizeof(float)*(innodes+2));

		status=cudaMalloc((void **)&d_max,sizeof(float)*(innodes+2));
		status=cudaMemset(d_max, 0,sizeof(float)*(innodes+2));

		status=cudaMalloc((void **)&d_arr_binloc,sizeof(float)*(innodes+4)*(resol+8));
		status=cudaMemset(d_arr_binloc, 0, sizeof(float)*(innodes+4)*(resol+8));

		status=cudaMalloc((void **)&d_dmyclass,sizeof(float)*(classes+2));
		status=cudaMemset(d_dmyclass, 0, sizeof(float)*(classes+2));

		if (status != cudaSuccess)
			printf("Error in cuda memory allocation\n");


		for(rnd=0;rnd<=oneround;rnd++)     // Training round starts here....
		{
			if((n==pocnt)&& (n>0)){ printf("breaking\n"); break;}
			strcpy(fltmp,fln);
			strcat(fltmp,".dat");
			fl1=fopen(fltmp,"r");
			n=0;
			rslt=0.0;
			rslt2=0.0;
			pcnt=0;
			int cindex=0;
			for(cindex=0;cindex<numlines;cindex++)
			{
				for(k=1;k<=outnodes;k++) classval[k]=1.0;
				n++;
				if(ans1==3)
				{
					for (i=1;i<=innodes;i++)
					{
						fscanf(fl1,"%f",&arr_vects[(cindex*innodes)+i]);
						if(nerror ==2){fscanf(fl1,"%f",&arr_err1vects[(cindex*innodes)+i]);fscanf(fl1,"%f",&err2vects[(cindex*innodes)+i]);}else
						if(nerror ==1){fscanf(fl1,"%f",&arr_err1vects[(cindex*innodes)+i]);err2vects[(cindex*innodes)+i]=arr_err1vects[(cindex*innodes)+i];}
					}
					fscanf(fl1,"\n");
				}
				else
				{
					for (i=1;i<=innodes;i++)
					{
						fscanf(fl1,"%f",&arr_vects[(cindex*innodes)+i]);
						if(nerror ==2){fscanf(fl1,"%f",&arr_err1vects[(cindex*innodes)+i]);fscanf(fl1,"%f",&arr_err2vects[(cindex*innodes)+i]);}else
						if(nerror ==1){fscanf(fl1,"%f",&arr_err1vects[(cindex*innodes)+i]);arr_err2vects[(cindex*innodes)+i]=arr_err1vects[(cindex*innodes)+i];}
					}
					fscanf(fl1,"%f\n",&arr_tmpv[cindex]);
				}
			}
			fclose(fl1);

			cudaMemcpy(d_arr_anti_net,arr_anti_net,sizeof(int)*totsize,cudaMemcpyHostToDevice);
			cudaMemcpy(d_resolution,resolution,sizeof(int)*(innodes+8),cudaMemcpyHostToDevice);
			cudaMemcpy(d_arr_anti_wts,arr_anti_wts,sizeof(float)*totsize,cudaMemcpyHostToDevice);
			cudaMemcpy(d_arr_tmpv,arr_tmpv,sizeof(float)*(numlines+2),cudaMemcpyHostToDevice);
			cudaMemcpy(d_arr_vects,arr_vects,sizeof(float)*(numlines+2)*(innodes+outnodes+2),cudaMemcpyHostToDevice);
			cudaMemcpy(d_arr_err1vects,arr_err1vects,sizeof(float)*(numlines+2)*(innodes+2),cudaMemcpyHostToDevice);
			cudaMemcpy(d_arr_err2vects,arr_err2vects,sizeof(float)*(numlines+2)*(innodes+2),cudaMemcpyHostToDevice);
			cudaMemcpy(d_min,min,sizeof(float)*(innodes+2),cudaMemcpyHostToDevice);
			cudaMemcpy(d_max,max,sizeof(float)*(innodes+2),cudaMemcpyHostToDevice);
			cudaMemcpy(d_arr_binloc,arr_binloc,sizeof(float)*(innodes+4)*(resol+8),cudaMemcpyHostToDevice);
			cudaMemcpy(d_dmyclass,dmyclass,sizeof(float)*(classes+2),cudaMemcpyHostToDevice);


			int numblocks=ceil(numlines/512.0);
			kernel1(numlines,arr_tmpv,min,max,resolution,arr_vects,arr_err1vects,arr_err2vects,arr_binloc,rn,arr_anti_net,arr_anti_wts,dmyclass,gain,innodes,resol,outnodes,nerror,rnd);

			strcpy(fltmp,fln);
			strcat(fltmp,".dat");
			fl1=fopen(fltmp,"r");

			m=n;
			n=0;
			rslt=0.0;
			rslt2=0.0;
			pcnt=0;


//			while (!feof(fl1))                    // Test round...
			for(cindex=0;cindex<numlines;cindex++)
			{
				n++;

				if(ans1==3)
				{
					for (i=1;i<=innodes;i++)
					{
						fscanf(fl1,"%f",&arr_vects[(cindex*innodes)+i]);
						if(nerror ==2){fscanf(fl1,"%f",&arr_err1vects[(cindex*innodes)+i]);fscanf(fl1,"%f",&arr_err2vects[(cindex*innodes)+i]);}else
						if(nerror ==1){fscanf(fl1,"%f",&arr_err1vects[(cindex*innodes)+i]);arr_err2vects[(cindex*innodes)+i]=arr_err1vects[(cindex*innodes)+i];}
					}
					fscanf(fl1,"\n");
				}
				else
				{
					for (i=1;i<=innodes;i++)
					{
						fscanf(fl1,"%f",&arr_vects[(cindex*innodes)+i]);
						if(nerror ==2){fscanf(fl1,"%f",&arr_err1vects[(cindex*innodes)+i]);fscanf(fl1,"%f",&arr_err2vects[(cindex*innodes)+i]);}else
						if(nerror ==1){fscanf(fl1,"%f",&arr_err1vects[(cindex*innodes)+i]);arr_err2vects[(cindex*innodes)+i]=arr_err1vects[(cindex*innodes)+i];}
					}
					fscanf(fl1,"%f\n",&arr_tmpv[cindex]);
				}

			}
			fclose(fl1);
			i=0;

			status=cudaMemset(d_pcnt, 0, sizeof(int)*(numlines+10));
			status=cudaMemset(d_rslt, 0, sizeof(float)*(numlines+10));
			status=cudaMemset(d_rslt2, 0, sizeof(float)*(numlines+10));
			//copy the results
			cudaMemcpy(d_resolution,resolution,sizeof(int)*(innodes+8),cudaMemcpyHostToDevice);
			cudaMemcpy(d_arr_tmpv,arr_tmpv,sizeof(float)*(numlines+2),cudaMemcpyHostToDevice);
			cudaMemcpy(d_arr_vects,arr_vects,sizeof(float)*(numlines+2)*(innodes+outnodes+2),cudaMemcpyHostToDevice);
			cudaMemcpy(d_arr_err1vects,arr_err1vects,sizeof(float)*(numlines+2)*(innodes+2),cudaMemcpyHostToDevice);
			cudaMemcpy(d_arr_err2vects,arr_err2vects,sizeof(float)*(numlines+2)*(innodes+2),cudaMemcpyHostToDevice);
			cudaMemcpy(d_min,min,sizeof(float)*(innodes+2),cudaMemcpyHostToDevice);
			cudaMemcpy(d_max,max,sizeof(float)*(innodes+2),cudaMemcpyHostToDevice);
			cudaMemcpy(d_arr_binloc,arr_binloc,sizeof(float)*(innodes+4)*(resol+8),cudaMemcpyHostToDevice);
			cudaMemcpy(d_dmyclass,dmyclass,sizeof(float)*(classes+2),cudaMemcpyHostToDevice);

			if((status = cudaGetLastError()) != cudaSuccess)
			{
				printf("Error(%s:%d) %s\n",__FILE__,__LINE__,cudaGetErrorString(status));
			}

			numblocks=ceil(numlines/512.0);
			kernel2<<<numblocks,512>>>(numlines,d_arr_tmpv,d_min,d_max,d_resolution,d_arr_vects,d_arr_err1vects,d_arr_err2vects,d_arr_binloc,rn,d_arr_anti_net,d_arr_anti_wts,d_dmyclass,gain,innodes,resol,outnodes,nerror,rnd,d_pcnt,d_rslt,d_rslt2);


			if((status = cudaGetLastError()) != cudaSuccess)
			{
				printf("%s\n",cudaGetErrorString(status));
			}
			cudaMemcpy(tmp_pcnt,d_pcnt,sizeof(int)*numlines,cudaMemcpyDeviceToHost);
			if((status = cudaGetLastError()) != cudaSuccess)
			{
				printf("%s\n",cudaGetErrorString(status));
			}
			cudaMemcpy(tmp_rslt,d_rslt,sizeof(float)*numlines,cudaMemcpyDeviceToHost);
			if((status = cudaGetLastError()) != cudaSuccess)
			{
				printf("%s\n",cudaGetErrorString(status));
			}
			cudaMemcpy(tmp_rslt2,d_rslt2,sizeof(float)*numlines,cudaMemcpyDeviceToHost);
			if((status = cudaGetLastError()) != cudaSuccess)
			{
				printf("%s\n",cudaGetErrorString(status));
			}

			pcnt=0;
			rslt=0;
			rslt2=0;
			for(i=0;i<numlines;i++)
			{
				pcnt+=tmp_pcnt[i];
				rslt+=tmp_rslt[i];
				rslt2+=tmp_rslt2[i];
			}



			printf("rnd:%d\trslt:%f\tRslt2:%f\tOrslt2:%f\tpcnt:%d\n",rnd,rslt,rslt2,orslt2,pcnt);

			kmax=1;
			if(orslt2==0) orslt2=rslt2;
			if(orslt==0) orslt=rslt;

			prslt=(rslt2-orslt2);
			if(rslt > 0)
			nrslt=(orslt/rslt);
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
								//fprintf(fl6,"%f ",anti_wts[i][j][l][m][k]);
								fprintf(fl6,"%f ",arr_anti_wts[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(m*kk)+k)]);

							}
							fprintf(fl6,"\n");
						}
					}
					fprintf(fl6,"\n");
				}
				fprintf(fl6,"\n");
				for(i=1;i<=innodes;i++)
				for(j=1;j<=resolution[i];j++)
				fprintf(fl6,"%f\n", arr_binloc[(i*rn)+(j)]);                 /// Let us print the bins.
				fflush(fl6);
				fclose(fl6);
				cout << "Round:" << rnn << "| TProb["<<prslt<<"," <<nrslt<<"] | Passed count:" << pocnt << endl;
				if(orslt2 <rslt2) orslt2=rslt2;
				if(rslt < orslt) orslt=rslt;
			}
			n=m;
		}  //rnd inc.
		fl6=NULL;
		cout << "Best result at round " << rnn<< endl;
	}  // ans <> 1
/***********************************End of Case 1*******************************/







    strcpy(fltmp,fln);
    strcat(fltmp,".dat");
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
		cout << "Creating the Anticipated Network outputs\n";
		for (i=1;i<=innodes;i++)
		{ 
			fscanf(fl2,"%d",&resolution[i]);
			for(j=0;j<=resolution[i];j++) arr_binloc[(i*rn)+(j+1)]=j*1.0;
		}
		fscanf(fl2,"%f",&omax);
		fscanf(fl2,"%f",&omin);
		fscanf(fl2,"\n");
        for(i=1;i<=innodes;i++) fscanf(fl2,"%f",&max[i]);
        fscanf(fl2,"\n");
    	for(i=1;i<=innodes;i++) fscanf(fl2,"%f",&min[i]);
		fscanf(fl2,"\n");
		for(i=1;i<=innodes;i++)for(j=0;j<=resolution[i];j++) 
		for(l=1;l<=innodes;l++)for(m=0;m<=resolution[l];m++) arr_anti_net[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(m*kk)+0)] =0;
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

						fscanf(fl2,"%d",&arr_anti_net[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(m*kk)+k)]);
						fscanf(fl6,"%f",&arr_anti_wts[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(m*kk)+k)]);
						arr_anti_net[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(m*kk)+0)]+=(float)(arr_anti_net[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(m*kk)+k)]);
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
		fscanf(fl6,"%f\n",&arr_binloc[(i*rn)+(j)]);                 /// Let us print the bins.
    }
    fclose(fl6);
    fl4=fopen("output.dat","w+");  // Network Output values
    cout << "Read all input parameters\n";
// *********** case 3 ***********************************************
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
		c1cnt=0;
		c2cnt=0;
		invcnt=0;
		n=0;
    }
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
				fscanf(fl1,"%f",&vects[i]);
				if(nerror ==2){fscanf(fl1,"%f",&err1vects[i]);fscanf(fl1,"%f",&err2vects[i]);}else
				if(nerror ==1){fscanf(fl1,"%f",&err1vects[i]);err2vects[i]=err1vects[i];} 
			}
			fscanf(fl1,"\n");
	    }
	    else
	    {
			for (i=1;i<=innodes;i++) 
			{
				fscanf(fl1,"%f",&vects[i]);
				if(nerror ==2){fscanf(fl1,"%f",&err1vects[i]);fscanf(fl1,"%f",&err2vects[i]);}else
				if(nerror ==1){fscanf(fl1,"%f",&err1vects[i]);err2vects[i]=err1vects[i];} 
			}
			fscanf(fl1,"%f\n",&tmpv);
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
		for(k=1;k<=outnodes;k++) classval[k]=1.0; tmp2_wts=1.0;
 		for (i=1;i<=innodes;i++)
		{
			j=0;
            if(vects[i]==MissingDat)
				skpchk=1;
            else
				skpchk=0;
            if ((resolution[i] >= vects[i]) &&(skpchk==0))
            {
                while ((fabs(vects[i]-arr_binloc[(i*rn)+(j+1)]) >=1.0)&& (j<= resolution[i]))
                {
                    j++;
                }
                jx=0;
            }
			else
			{
//NSP_added jx=-1;				  
			  jx=1;
			}
	        for (l=1;l<=innodes;l++)
	        {
				if((i!=l) && (jx==0))
				{
					m=0;
					if((vects[l]==MissingDat)||(vects[i]==MissingDat))
						skpchk=1;
					else
						skpchk=0;
					if ((resolution[l] >= vects[l]) &&(skpchk==0))
					{
						while ((fabs(vects[l]-arr_binloc[(l*rn)+(m+1)]) >=1.0)&& (m<= resolution[l]))
						{       
							m++;
						}
					}
					for (k=1;k<=outnodes;k++)
					{

						if(jx==0){tmp2_wts=(float)arr_anti_net[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(m*kk)+k)];}else{tmp2_wts=1.0/outnodes;}
						if(nerror ==2) 
						{
							for(p=(m-(int)err1vects[l]);p<=(m+(int)err2vects[l]);p++)
							{

								if(p<0) p=0; if(p>resolution[l]) break;
								if ((float)arr_anti_net[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(p*kk)+k)] > tmp2_wts)
								m=p;
							}
						}
						if(nerror ==1) 
						{
							for(p=(m-(int)err1vects[l]);p<=(m+(int)err1vects[l]);p++)
							{

								if(p<0) p=0; if(p>resolution[l]) break;
								if ((float)arr_anti_net[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(p*kk)+k)] > tmp2_wts)
								m=p;
							}
						}

						if((arr_anti_net[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(m*kk)+0)] > 0) && (resolution[i]>= vects[i])&& (resolution[l]>= vects[l])&&(skpchk==0))
						{

							if(jx==0){tmp2_wts=(float)arr_anti_net[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(m*kk)+k)]*arr_anti_wts[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(m*kk)+k)]*1.0/(arr_anti_net[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(m*kk)+0)]);}
							else{tmp2_wts=1.0/outnodes;}
						}
						else
						if(skpchk == 1)   // || bcchoice == 'y')
						{
							tmp2_wts= 1.0; //(float)1.0/outnodes; //1.0; //
						}
						else
						{
							tmp2_wts=(float)1.0/outnodes;
						}
						if((resolution[i] >= vects[i])&& (resolution[l]>= vects[l])&&(skpchk==0))
						{
							classval[k]*=(float)tmp2_wts;
						}
					}
					totprob=0;
					for(k=1;k<=outnodes;k++) totprob+=classval[k];
					if (totprob==0) {totprob=innodes*outnodes; cout <<"Caution!! Item did not have known types\n";}
					for(k=1;k<=outnodes;k++) classval[k]=classval[k]/totprob;
 				
				}
			}
	    }
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
	        totprob += (float)classval[k];
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
				fprintf(fl4,"%d  %f %-5.2f %f %-5.2f %f %-5.2f %f %-5.2f",n, dmyclass[(int)kmax],100.0*((classval[kmax])/totprob),dmyclass[(int)k2max],100.0*((classval[k2max])/totprob),dmyclass[(int)k3max],100.0*((classval[k3max])/totprob),dmyclass[(int)k4max],100.0*((classval[k4max])/totprob));
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
				fprintf(fl4,"%d  %f\n",n, dmyclass[(int)kmax]);
				fprintf(fl7, "%-8d    %f   %f     %f    %f     %f    ",n,dmyclass[(int)kmax],dmyclass[(int)k2max],dmyclass[(int)k3max],dmyclass[(int)k4max],tmpv);
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
			fprintf(fl5,"%d %e \n",n,(float) tmpv);
		} // ans1 != 3 ends here ******************
	}
	cout << "The suggested LoC is " << nLoC/nLoCcnt << "\n";
	fclose(fl1);
	fclose(fl2);
	fclose(fl4);
	if(ans1 < 3)
	{
		strcpy(fltmp,fln);
//	   tmp2_wts=0.0;
		fclose(fl5);
		fprintf(fl7,"*________________________________________________________________________\n");
		fprintf(fl7,"*Total    Success in   Success in   Non classified   Real success in    \n");
		cout << "*________________________________________________________________________\n";
		cout << "*Total    Success in   Success in   Non classified   Real success in    \n";
		if (outnodes > 3)
	    {
			fprintf(fl7,"* No.    Ist Choice  2nd Choice     items           two chances    \n");
     		fprintf(fl7,"* %d       %d          %d           %d             %-5.2f %% \n",n,c1cnt,c2cnt,invcnt,(float)100.0*(c1cnt+c2cnt)/(n-invcnt));
     		cout << "* No.    Ist Choice  2nd Choice     items           two chances    \n";
     		printf("* %d       %d          %d           %d             %-5.2f %% \n",n,c1cnt,c2cnt,invcnt,(float)100.0*(c1cnt+c2cnt)/(n-invcnt));
     	}
     	else
     	{
			fprintf(fl7,"* No.    Ist Choice  2nd Choice     items           First chance    \n");
     		fprintf(fl7,"* %d       %d          %d           %d             %-5.2f %% \n",n,c1cnt,c2cnt,invcnt,(float)100.0*(c1cnt)/(n-invcnt));
     		cout << "* No.    Ist Choice  2nd Choice     items           First chance    \n";
     		printf("* %d       %d          %d           %d             %-5.2f %% \n",n,c1cnt,c2cnt,invcnt,(float)100.0*(c1cnt)/(n-invcnt));
     	}
	  	fprintf(fl7,"*________________________________________________________________________\n");
		printf("*________________________________________________________________________\n");
		fclose(fl7);
	} // ******** ans1!=3 ends here *************
	cout << "Done.\n";
	stop = times(NULL);
	cout << "The computation took " << fabs(start - stop)*10000/(CLOCKS_PER_SEC) << " Secs.\n";
} //end main



