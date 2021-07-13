/*
 * kernel.cu
 *
 *  Created on: May 22, 2015
 *  Author: Ajay Vibhute
 */




void kernel1(int numlines,float *arr_tmpv,float *min,float *max,int* resolution,float *arr_vects,float *arr_err1vects,float *arr_err2vects,float*arr_binloc,int rn,int *arr_anti_net,float *arr_anti_wts,float*dmyclass, float gain,int innodes,int resol,int outnodes,int nerror,int rnd)
{
	int cindex=0,i,j,k,l,m,p,jx=0,MissingDat=-9999,kmax=0;
	float tmp2_wts=0,tmpv=0,totprob=0,cmax=0;
	int jk=resol+1,lk=innodes+1,mk=resol+1,kk=outnodes+1;
	float classval[classes+2];
	for(cindex=0;cindex<numlines;cindex++)
	{
		tmpv=arr_tmpv[cindex];
		for(i=1;i<=innodes;i++)
		{
			if((arr_vects[(cindex*innodes)+i] != MissingDat)&&(max[i]!=MissingDat))
			{
				arr_vects[(cindex*innodes)+i]=round((arr_vects[(cindex*innodes)+i]-min[i])/(max[i]-min[i])*resolution[i]);
				arr_err1vects[(cindex*innodes)+i]=round((arr_err1vects[(cindex*innodes)+i])/(max[i]-min[i])*resolution[i]);
				arr_err2vects[(cindex*innodes)+i]=round((arr_err2vects[(cindex*innodes)+i])/(max[i]-min[i])*resolution[i]);
			}
		}
		for(k=1;k<=outnodes;k++) classval[k]=1.0;
		for (i=1;i<=innodes;i++)
		{
			j=0;
			k=1;
			if(arr_vects[(cindex*innodes)+i] != MissingDat)
			{

				while ((fabs(arr_vects[(cindex*innodes)+i]-arr_binloc[(i*rn)+(j+1)]) >=1.0)&& (j<= resolution[i]))
				{
					j++;
				}
				//NSP_added if(fabs(vects[i]-binloc[i][j+1]) <= binloc[i][0]){jx=0;} else{jx=-1;}
				if(fabs(arr_vects[(cindex*innodes)+i]-arr_binloc[(i*rn)+(j+1)]) <= 1.0){jx=0;} else{jx=-1;}

				for (l=1;l<=innodes;l++)
				{
				  if(i!=l)
				  {
					m=0;
					k=1;
					if(arr_vects[(cindex*innodes)+l] != MissingDat)
					{
						while ((fabs(arr_vects[(cindex*innodes)+l]-arr_binloc[(l*rn)+(m+1)]) >=1.0)&& (m<= resolution[l]))
						{
							m++;
						}
						for (k=1;k<=outnodes;k++)
						{

							if(jx==0)tmp2_wts=(float)arr_anti_net[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(m*kk)+k)];
							else tmp2_wts=1.0/outnodes;
							if(nerror ==2)
							{
								for(p=(m-(int)arr_err1vects[(cindex*innodes)+l]);p<=(m+(int)arr_err1vects[(cindex*innodes)+l]);p++)
								{
									if(p<0) p=0; if(p>resolution[l]) break;

									if ((float)arr_anti_net[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(p*kk)+k)] > tmp2_wts)
									m=p;
								}
							}
							if(nerror ==1)
							{
								for(p=(m-(int)arr_err1vects[(cindex*innodes)+l]);p<=(m+(int)arr_err1vects[(cindex*innodes)+l]);p++)
								{

									if(p<0) p=0; if(p>resolution[l]) break;
									if ((float)arr_anti_net[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(p*kk)+k)] > tmp2_wts)
									m=p;
								}
							}

							if(arr_anti_net[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(m*kk)+0)] > 0)
							{
								if(jx==0)
									tmp2_wts=(float)arr_anti_net[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(m*kk)+k)]*arr_anti_wts[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(m*kk)+k)]/(arr_anti_net[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(m*kk)+0)]);
								else
								tmp2_wts=1.0/outnodes;
							}
							else
								tmp2_wts= 1.0/outnodes;
							classval[k]*=(float)tmp2_wts;
						}
						totprob=0;
						for(k=1;k<=outnodes;k++) totprob+=classval[k];
						if (totprob==0) {totprob=innodes*outnodes; /*cout <<"Caution : Item has no representation type\n";*/}
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
				if(arr_vects[(cindex*innodes)+i] != MissingDat)
				{
					while ((fabs(arr_vects[(cindex*innodes)+i]-arr_binloc[(i*rn)+(j+1)]) >=1.0)&& (j<= resolution[i]))
					{
						j++;
					}
					for (l=1;l<=innodes;l++)
					{
					if(i!=l)
					{
						m=0;
						k=1;
						if(arr_vects[(cindex*innodes)+l] != MissingDat)
						{
							while ((fabs(arr_vects[(cindex*innodes)+l]-arr_binloc[(l*rn)+(m+1)]) >=1.0)&& (m<= resolution[l]))
							{
								m++;
							}
							while ((k<=outnodes)&&fabs(dmyclass[k]-tmpv) > dmyclass[0]) k++;
							if((classval[(int)kmax] >0)&&(classval[k]<classval[(int)kmax]))
							{

								arr_anti_wts[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(m*kk)+k)]+=(float)gain*(1.0-(classval[k]/classval[(int)kmax]));
							}
							/*if(arr_anti_wts[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(m*kk)+k)] <= 0.0)
								cout << k << " "<< tmpv << "[" << dmyclass[1] << "]" << dmyclass[outnodes] << "\n";*/
						}
					}
					}
				}
			}
		} // kmax che
	} // while not eof check
}









__global__ void kernel2(int numlines,float *arr_tmpv,float *min,float *max,int* resolution,float *arr_vects,float *arr_err1vects,float *arr_err2vects,float*arr_binloc,int rn,int *arr_anti_net,float *arr_anti_wts,float*dmyclass,
		float gain,int innodes,int resol,int outnodes,int nerror,int rnd,int*pcnt,float*rslt,float*rslt2)
{

	int cindex=0,i,j,k,l,m,p,jx=0,MissingDat=-9999,kmax=0;
	float tmp2_wts=1,tmpv=0,totprob=0,cmax=0,oldj;
	int jk=resol+1,lk=innodes+1,mk=resol+1,kk=outnodes+1;
	float classval[classes];
	cindex=threadIdx.x +( blockIdx.x * blockDim.x);
	if(cindex<numlines)
	{

		kmax=1;
		cmax=0;
		tmpv=arr_tmpv[cindex];
		for(i=1;i<=innodes;i++)
		{
			if((arr_vects[(cindex*innodes)+i] != MissingDat)&&(max[i]!=MissingDat))
			{
				arr_vects[(cindex*innodes)+i]=round((arr_vects[(cindex*innodes)+i]-min[i])/(max[i]-min[i])*resolution[i]);
				arr_err1vects[(cindex*innodes)+i]=round((arr_err1vects[(cindex*innodes)+i])/(max[i]-min[i])*resolution[i]);
				arr_err2vects[(cindex*innodes)+i]=round((arr_err2vects[(cindex*innodes)+i])/(max[i]-min[i])*resolution[i]);
				if (arr_vects[(cindex*innodes)+i] < 0) arr_vects[(cindex*innodes)+i]=0;             // let us be bounded. #Oct 2001.
			}
		}
		for(k=1;k<=outnodes;k++) classval[k]=1.0;
		for (i=1;i<=innodes;i++)
		{
			j=0;
			k=1;
			if(arr_vects[(cindex*innodes)+i] != MissingDat)
			{
				while ((fabs(arr_vects[(cindex*innodes)+i]-arr_binloc[(i*rn)+(j+1)]) >=1.0)&& (j<= resolution[i]))
				{
					j++;
				}
				//NSP_added if(fabs(vects[i]-binloc[i][j+1]) <= binloc[i][0]){jx=0;} else{jx=-1;}
				if(fabs(arr_vects[(cindex*innodes)+i]-arr_binloc[(i*rn)+(j+1)]) < 1.0){jx=0;} else{jx=-1;}
				for (l=1;l<=innodes;l++)
				{
				 if (i !=l)
				 {
					m=0;
					k=1;
					if(arr_vects[(cindex*innodes)+l] != MissingDat)
					{
						oldj=(float)2*resolution[l];
						while ((fabs(arr_vects[(cindex*innodes)+l]-arr_binloc[(l*rn)+(m+1)]) >=1.0)&& (m<= resolution[l]))
						{
							oldj=fabs(arr_vects[(cindex*innodes)+l]-arr_binloc[(l*rn)+(m+1)]);
							m++;
						}
						for (k=1;k<=outnodes;k++)
						{

							if(jx==0)
							tmp2_wts=(float)arr_anti_net[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(m*kk)+k)];
							else
							tmp2_wts=1.0/outnodes;
							if(nerror ==2)
							{
								for(p=(m-(int)arr_err1vects[(cindex*innodes)+l]);p<=(m+(int)arr_err2vects[(cindex*innodes)+l]);p++)
								{
									if(p<0) p=0; if(p>resolution[l]) break;

									if ((float)arr_anti_net[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(p*kk)+k)] > tmp2_wts)
									m=p;
								}
							}
							if(nerror ==1)
							{
								for(p=(m-(int)arr_err1vects[(cindex*innodes)+l]);p<=(m+(int)arr_err1vects[(cindex*innodes)+l]);p++)
								{

									if(p<0) p=0; if(p>resolution[l]) break;
									if ((float)arr_anti_net[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(0*kk)+k)] > tmp2_wts)
									m=p;
								}
							}

								if(arr_anti_net[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(m*kk)+0)] > 0)
								{
									if(jx==0)
									tmp2_wts=(float)arr_anti_net[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(m*kk)+k)]*arr_anti_wts[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(m*kk)+k)]/(arr_anti_net[((i*jk*lk*mk*kk)+(j*lk*mk*kk)+(l*mk*kk)+(m*kk)+0)]);
									else
									tmp2_wts=1.0/outnodes;
								}
							else
								tmp2_wts=(float)1.0/outnodes;
							classval[k]*=(float)tmp2_wts;
						}

						totprob=0;
						for(k=1;k<=outnodes;k++) totprob+=classval[k];
						if (totprob==0) {totprob=innodes*outnodes; }
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
			rslt2[cindex]+=cmax;
			pcnt[cindex]+=1;
		}
		else
		{
			k=1;
			while ((k<=outnodes)&&fabs(dmyclass[k]-tmpv) > dmyclass[0]) k++;
			rslt[cindex]+=(cmax-classval[k]);
		}

		} // while not eof check

}

int getNumlines(char*filename)
{
	FILE *fp=fopen(filename,"r");
	char ch;
	int lines=0;
	while(!feof(fp))
	{
	  ch = fgetc(fp);
	  if(ch == '\n')
	  {
	    lines++;
	  }
	}
	fclose(fp);
	return lines;


}
