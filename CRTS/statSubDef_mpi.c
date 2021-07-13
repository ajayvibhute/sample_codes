/*
 * statSubDefSer.c
 *
 *  Created on: Sep 16, 2013
 *  Author: Ajay Vibhute
 */

/*Function to calculate mean.
 * */
#include "statSubDef_mpi.h"
#include<stdio.h>
void processLightCurve(struct lightcurve *curves,struct lightcurveParam*param)
{
		int index=0;
		if(mean(&curves[index].mag[0],&param[index].xm,curves[index].numpoints)==-1)
		{
			printf("Error: error in mean calculation\n");
			printf("Exiting......\n");
			//exit(0);
		}
		if(mean(&curves[index].mjd[0],&param[index].mjdm,curves[index].numpoints)==-1)
		{
			printf("Error: error in mjd mean calculation\n");
			printf("Exiting......\n");
			//exit(0);
		}
		if(sd(&curves[index].mag[0],&param[index].xsig,curves[index].numpoints, param[index].xm)==-1)
		{
			printf("Error: error in standard deviation calculation\n");
			printf("Exiting......\n");
			//exit(0);
		}

		if(abbe(&curves[index].mag[0],&param[index].abbe,curves[index].numpoints,param[index].xm)==-1)
		{
			printf("Error: error in abbe calculation\n");
			printf("Exiting......\n");
			//exit(0);
		}
		if(skew(&curves[index].mag[0],&param[index].xsk,curves[index].numpoints, param[index].xm)==-1)
		{
			printf("Error: error in skew calculation\n");
			printf("Exiting......\n");
			//exit(0);
		}
		if(kurtosis(&curves[index].mag[0],&param[index].xkur,curves[index].numpoints, param[index].xm,param[index].xsig)==-1)
		{
			printf("Error: error in kurtosis calculation\n");
			printf("Exiting......\n");
			//exit(0);
		}
		if(mad(&curves[index].mag[0],&param[index].xmad,curves[index].numpoints,&curves[index].temp[0]) ==-1)
		{
			printf("Error: error in Median absolute deviation calculation\n");
			printf("Exiting......\n");
			//exit(0);
		}
		if(bwmv(&curves[index].mag[0],&param[index].xbwmv,curves[index].numpoints,&curves[index].temp[0]) ==-1)
		{
			printf("Error: error in Biweight midvariance calculation\n");
			printf("Exiting......\n");
			//exit(0);
		}
		if(thiel_sen(&curves[index].mjd[0],&curves[index].mag[0],&param[index].thiel_sen,&curves[index].temp[0], curves[index].numpoints) ==-1)
		{
			printf("Error: error in Thiel-Sen estimator calculation\n");
			printf("Exiting......\n");
			//exit(0);
		}
		if(durbin_watson(&curves[index].mag[0],&param[index].xdw,curves[index].numpoints,&curves[index].temp[0], param[index].xm)==-1)
		{
			printf("Error: error in Durbin-Watson calculation\n");
			printf("Exiting......\n");
			//exit(0);
		}
		if(kendall_tau(&curves[index].mjd[0],&curves[index].mag[0],&param[index].xkt,curves[index].numpoints) ==-1)
		{
			printf("Error: error in Kendall tau rank statistic calculation\n");
			printf("Exiting......\n");
			//exit(0);
		}
		if(cusum(&curves[index].mag[0],&param[index].cusum,curves[index].numpoints,param[index].xm,param[index].xsig) ==-1)
		{
			printf("Error: error in Kendall tau rank statistic calculation\n");
			printf("Exiting......\n");
			//exit(0);
		}
		//not calling due to some memory allocation problem
		/*if(dfa(&curves[index].mag[0],&curves[index].mjd[0],&curves[index].xdfa,&curves[index].temp[0], curves[index].numpoints,curves[index].xm) ==-1)
		{
			printf("Error: error in Kendall tau rank statistic calculation\n");
			printf("Exiting......\n");
			//exit(0);
		}*/

		if(von_neumann(&curves[index].mag[0],&param[index].xvon,curves[index].numpoints,param[index].xsig) ==-1)
		{
			printf("Error: error in cusum calculation\n");
			printf("Exiting......\n");
			//exit(0);
		}
		if(lsq(&curves[index].mag[0],&curves[index].mjd[0],&param[index].lsq[0], curves[index].numpoints,param[index].xm,  param[index].mjdm) ==-1)

		{
			printf("Error: error in lsq calculation\n");
			printf("Exiting......\n");
			//exit(0);
		}

		if(stetson_j(&curves->mjd[0],&curves->mag[0],&curves->mag_error[0],&param->stetson_j, curves->numpoints,param->xm) ==-1)
		{
			printf("Error: error in stetson_j calculation\n");
			printf("Exiting......\n");
			//exit(0);
		}
		if(stetson_k(&curves[index].mag[0],&curves[index].mag_error[0],&param[index].stetson_k, curves[index].numpoints,param[index].xm) ==-1)
		{
			printf("Error: error in stetson_k calculation\n");
			printf("Exiting......\n");
			//exit(0);
		}
		if(con(&curves[index].mag[0],&param[index].con, curves[index].numpoints,param[index].xm,param[index].xsig) ==-1)
		{
			printf("Error: error in con calculation\n");
			printf("Exiting......\n");
			//exit(0);
		}

		if(reducedchi(&curves[index].mag[0],&curves[index].mag_error[0],&param[index].rchi, curves[index].numpoints,param[index].xm) ==-1)
		{
			printf("Error: error in reducedchi calculation\n");
			printf("Exiting......\n");
			//exit(0);
		}

}
int getPercentileRank(double percentile,int n )
{
	return((int)round(((percentile/100.0)*n)+0.5))-1;
}

int mean(double*data,double*mean,int n)
{
	int i=0;
	*mean=0;
	if(n<0)
		return -1;

	for(i=0;i<n;i++)
	{
		*mean+=data[i];
	}
	*mean/=n;
	return 0;
}

/*Function to calculate sd.
 * */
int sd(double*data,double*sd,int n,double xm=0)
{
	int i=0;

	*sd=0;
	if(n<=1)/*check because denom is n-1*/
		return 0;

	if(xm==0)
		mean(data,&xm,n);
	for(i=0;i<n;i++)
		*sd+=((data[i]-xm)*(data[i]-xm));
	*sd=sqrt(*sd/(n-1));
	return 0;
}

int skew(double*x,double*xsk,int n,double xm=0)
{
	int i=0;
	double m3=0,m2=0,tmp=0;
	*xsk=0;
	if(n<=0)
		return -1;
	if(xm==0)
		mean(x,&xm,n);
	for(i=0;i<n;i++)
	{
		tmp=x[i]-xm;
		m3+=pow(tmp,3);
		m2+=pow(tmp,2);
	}
	m3/=n;
	m2/=n;
	m2=pow(m2,1.5);
	*xsk=(sqrt(n*(n-1))*(m3/m2))/(n-2); /*n-2 may cause be zero discuss what to do in that situation to avoid nan*/

	return 0;
}

int kurtosis(double*x,double*xkur,int n,double xm=0,double xs=0)
{
	int i=0;
	double m3=0,m2=0,tmp=0;
	*xkur=0;
	if(n<=0)
		return -1;
	if(n<4)
		return 0;
	if(xm==0)
		mean(x,&xm,n);
	if(xs==0)
		sd(x,&xs,n,xm);

	for(i=0;i<n;i++)
	{
		tmp+= pow ((x[i]-xm)/xs,4);
	}
	*xkur=( (n*(n+1)*tmp)/((n-1)* (n-2)*(n-3)) )  - 3*( ((n-1)*(n-1)) / ((n-2)*(n-3)) );

	return 0;

}
int  sort(double *x,int n)
{
	int i=0,j=0;
	double tempval=0;
	for (i = 0; i < n; ++i)
	{
		for (j = i + 1; j < n; ++j)
		{
			if (x[i] > x[j])
			{
				tempval =  x[i];
				x[i] = x[j];
				x[j] = tempval;
			}

		}

	}
	return 0;
}


int median(double *x,double *xmed,int n)/*check fortran is zero based or not change code accordingly*/
{
	int i=0;
	double *tmp;
	*xmed=0;
	if(n<=1)
		return -1;
	tmp=(double*)malloc(sizeof(double)*n);
	if(tmp==NULL)
	{
		printf("Error (%s:%d):Out of memory\n",__FILE__,__LINE__);
		return -1;
	}
	for(i=0;i<n;i++)
		tmp[i]=x[i];
	sort(tmp,n);
	if(n%2==0)
		*xmed=( (tmp[(n/2)+1]+tmp[n/2])/2.0 );
	else
		*xmed=(tmp[(n/2)]);
	free(tmp);
	return 0;
}


/* Median absolute deviation*/
int mad(double*x,double*xmad,int n,double*xtemp)
{
	double xmed=0,tmp=0;
	int i=0;

	if(n<=0)
		return -1;

	if(xtemp==NULL)
	{
		printf("Error : Out of memory from mad\n");
		printf("Exiting.....\n");
		exit(0);
	}
	median(x,&xmed,n);
	for(i=0;i<n;i++)
	{
		tmp=(x[i]-xmed);
		xtemp[i]=(tmp>0)?tmp:tmp*-1;
	}
	tmp=0;
	median(xtemp,&tmp,n);
	*xmad = 1.4826 * tmp;

	return 0;
}

/*Absolute pairwise difference statistics*/
int sn(double*x,double *sn,int n)
{
	int i=0,j=0,k=0,count=0;
	double *dx,*imed,*diff,dn;
	if(n<=0)
		return -1;
	dx=(double*)malloc(sizeof(double)*n);
	imed=(double*)malloc(sizeof(double)*n);
	diff=(double*)malloc(sizeof(double)*n*(n-1)/2);
	count=0;
	for(i=0;i<n;i++)
	{
		for(k=0;k<n;k++)
			dx[j]=x[i]-x[j];
		median(dx,&imed[i],n);
		for(j=i+1;j<n;j++)
			diff[count++]=dx[j];
	}

	median(imed,&sn[0],n);
	sn[0]*=1.1926;
	if(n%2==1)
		dn = n / (n + 1.4);
	else
	    dn = n / (n + 3.8);

	 //call sortx_h(diff, idx) need to implement this method. discuss with arun
	  /*sn(2) = dn * diff(idx(nint(0.25 * n)))*/
	//free(dx);
	//free(imed);
	//free(diff);
	return 0;
}
/* Biweight midvariance*/
int bwmv(double*x,double *xbwmv,int n,double*buffer)
{
	int i=0;
	double u,  xmad, q, num, denom, usq,tmp;
	*xbwmv=0;
	if(n<=0)
		return -1;

	median(x,&q,n);
	mad(x,&xmad,n,buffer);
	num = 0;
	denom = 0;
	for(i=0;i<n;i++)
	{
		u=(x[i]-q)/(9*xmad);
		tmp=u<0?u*-1:u;
		if(tmp<1)
		{
			usq=u*u;
			num+=pow( (x[i] - q), 2) * pow((1 - usq), 4);
			denom = denom + (1 - usq) * (1 - 5 * usq);
		}
	}
	*xbwmv =(n * num) / (denom * denom);
	return 0;
}
/*Thiel-Sen estimator of median slope*/
int thiel_sen(double *x,double*y, double *xthiel_sen,double *med, int n)
{

	int count=0,i=0,j=0;
	*xthiel_sen=0;

	if(n<=0)
		return -1;
	med=(double*)malloc(sizeof(double)* (n * (n - 1) / 2) );
	if(med==NULL)
	{
		printf("Error : Out of memory from thiel_sen\n");
		printf("Exiting.....\n");
		exit(0);
	}
	count=0;
	for(i=0;i<n;i++)
	{
		for(j=i+1;j<n;j++)
		{
			med[count++]=(y[j] - y[i]) / (x[j] - x[i]);
		}
	}
	median(med,xthiel_sen,count);
	return 0;
}

/*Durbin-Watson statistic*/
int durbin_watson(double*x,double*xdw,int n,double*e,double xm=0)
{
	int i=0;
	double tmp=0,esqsum=0;
	*xdw=0;
	if(n<=0)
		return -1;
	if(xm==0)
		mean(x,&xm,n);
	for(i=0;i<n;i++)
		esqsum+=pow((x[i] - xm),2);
	*xdw=0;
	for(i=1;i<n;i++)
		tmp +=pow((double)((x[i] - xm) - (x[i-1] - xm)),2);
	tmp/=esqsum;
	*xdw=tmp;
	return 0;
}
/*Kendall tau rank statistic*/
int kendall_tau(double *x,double*y, double *xkt, int n)
{
	int   i=0, j=0, con=0, dis=0, tj=0, uj=0;
	*xkt=0;
	if(n<=0)
		return -1;
	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
		{
			if (i !=j){
				if (((x[i] > x[j]) && (y[i] > y[j])) || ((x[i] < x[j]) &&(y[i] < y[j])))
					con++;
				else if (x[i] == x[j])
					tj++;
				else if (y[i] == y[j])
					uj++;
				else
					dis++;
			}
		}
	/*Wikipedia tau-a
		kt = (con - dis) / (0.5 * n * (n - 1))
	Wikipedia tau-b
		n0 = 0.5 * n * (n - 1)
		kt = (con - dis) / sqrt(double((n0 - tj) * (n0 - uj)))
	 Python version of tau-b (takes into account ties)*/
	 *xkt = (double)(con - dis) / (double)sqrt(((con + dis + tj) * (con + dis + uj)));
	 return 0;
}
/*Cumulative sum range*/
int cusum(double*x,double*xcusum,int n,double xm=0,double xs=0)
{
	int i=0,j=0;
	double min=99999999,max=0,tmp=0;
	*xcusum=0;
	if(n<=0)
		return -1;
	if(n==1)
		return 1;
	if(xm==0)
		mean(x,&xm,n);
	if(xs==0)
		sd(x,&xs,n,xm);
	for(i=0;i<n;i++)
	{
		tmp=0;
		for(j=0;j<=i;j++)
			tmp+=(x[j] - xm);
		tmp/= (n * xs);
		if(tmp<min)
			min=tmp;
		if(tmp>max)
			max=tmp;
	}
	*xcusum=max-min;
	return 0;


}
/*von Neumann variability index*/
int von_neumann(double*x,double*xvon,int n,double xs=0)
{
	int i=0,j=0;
	double eta=0;
	*xvon=0;
	if(n<=0)
		return -1;
	if(n==1)
		return 0;
	if(xs==0)
		sd(x,&xs,n);
	for(i=0;i<n-1;i++)
		eta+=((x[i+1]-x[i])*(x[i+1]-x[i]));
	*xvon= eta / ((n - 1) * xs * xs);
	return 0;
}


int stetson_j(double*mjd,double*data,double*dataerr,double *stetson_j,int n,double datam=0)/*Arun*/
{
	int i;
	double dt=0,sum=0,w_sm=0,wk,pk;
	double tmp;
	*stetson_j=0;
	/* Return zero if number of data points is null  */
	if(n<=0)
		return -1;
	/* Check for mean*/
	if(datam==0)
	   mean(data,&datam,n);

	for (i=0;i<n-1;i++)
		dt += (double)mjd[i+1] - mjd[i];
	dt = dt/(n-1);
	double tval=(double)sqrt( ((double)n/(double)(n-1)) );
	for (i=0;i<(n-1);i++)
	{
		wk =(double) exp(-(mjd[i+1] - mjd[i])/ dt);
		pk = (tval * ((double)(data[i] -datam) / dataerr[i]))  * (tval * ((double)(data[i+1] -datam) / dataerr[i+1])) ;
		tmp=pk;
		pk<0?pk*=-1:pk;
		sum += (double)(wk * (tmp/pk) * sqrt(pk));/*Check this*/
		w_sm = (double)(w_sm + wk);
	}
	*stetson_j =(double) sum /(double)w_sm;
	return 0;
}

/* Function to calculate Stetson -K Statistics*/

int stetson_k(double*data,double*dataerr,double *stetson_k,int n,double datam=0)/*Arun*/
{
	int i;
	double delta=0,deltasum=0,deltasq=0,deltasqsum=0,datatemp=0,denom=0,tmp;
	/* Return zero if number of data points is null  */
	if(n<=0)
	{
		*stetson_k=0;
		return -1;
	}
	/* Check for mean*/
	if(datam==0)
	   mean(data,&datam,n);
	tmp=sqrt(((double)n / (double)(n-1)));
	for (i=0;i<n;i++)
	{
		delta = ((data[i] -datam)/dataerr[i])*tmp;
		deltasum+=delta<0?delta*-1:delta;
		deltasqsum+=(delta*delta);
	}
	deltasqsum=sqrt((double)deltasqsum)*sqrt((double)n);
	*stetson_k = deltasum/deltasqsum;
	return 0;

}

int con(double*data,double *con, int n,double datamean=0,double xsig=0)
 {
	int i,total=0;
	double tol;
	*con=0;
	/* Return zero if number of data points is null  */
	if(n<=1)
		return 0;
		/* Check for mean*/
     if(datamean==0)
	   mean(data,&datamean,n);
     if(xsig==0)
    	 sd(data,&xsig,n,datamean);
    tol = 2.0* xsig;
	for (i=0;i<(n-2);i++)
	{
		//check
		if(((getabs(data[i])-datamean) > tol) && (( getabs(data[i+1])-datamean) > tol ) && ((getabs(data[i+2])-datamean ) >tol) )
		{
			total = total+1;
		}
	}
	*con = total/(n-1);
	return 0;

}



/*Reduced chi-square*/
int reducedchi(double*data,double *dataerr,double*rchi,int n,double datamean=0)/*Arun*/
{
	int i;
	double datasum=0;
	*rchi=0;
	/* Return zero if number of data points is null  */
	if(n<=0)
		return -1;
	if(n==1)
		return 1;
	/* Check for mean*/
	if(datamean==0)
	   mean(data,&datamean,n);

	for(i=0;i<n;i++)
	    datasum += ((data[i] -datamean) *(data[i] -datamean) )/dataerr[i];

	*rchi =	 datasum / (n-1);

	return 0;

}
/*Detrended fluctuation analysis*/
int dfa(double*x,double*t,double *xdfa,double *xt,int n,double xm=0)
{
	int  k=0, i=0, j=0, m=0, nbins=0,si=0,xlengh=0;
		*xdfa=0;
		double *f,*l,err=0,*lfit;
		if(n<=0)
			return -1;

		/* Check for mean*/
		if(xm==0)
		   mean(x,&xm,n);
		nbins = n / 11 ;/* No. of bins (min. bin occupancy = 11)*/
		f=(double*)malloc(sizeof(double)*nbins);
		l=(double*)malloc(sizeof(double)*nbins);
		lfit=(double*)malloc(sizeof(double)*2);
		for(i=0;i<n;i++)
		{
			xt[i]=0;
			for(j=0;j<i+1;j++)
			{
				xt[j]+=x[j] - xm;
			}
		}
		for(k = 1;k<=nbins;k++)
		{
			m=n/k;
		    err = 0.0;
		    for(j = 1;j<= k;j++)
		    {
		    	si=(j-1)*m;
		    	xlengh=(j*m)-si+1;
		    	lsq(&t[si],&xt[si],lfit,xlengh,0,0 );
		        for(i=0;i<n;i++)
		        	err+=pow(    ( xt[i]-lfit[0]-(lfit[1]*t[i]))    ,2);//check this expression not sure about this
		    }
		    if (k >0)
		    	f[k] = log10( (f[k - 1] +sqrt(err / k)) );
		    else
		    	f[0] = log10(sqrt(err));

		    l[k] = log10((double)k);
		}
		lsq(f,l,lfit,nbins,0,0 );
		*xdfa= lfit[1];

		free(f);
		free(l);
		free(lfit);
		return 0;
}

int lsq(double*x,double*y,double*lsq,int n,double xm=0,double ym=0)
{
	int i=0;
	double alpha,beta;
	double nsum=0,dsum=0,tmp=0;
	/* Return zero if number of data points is null  */
	lsq[0]=0;
	lsq[1]=0;
	if(n<=0)
		return -1;
	/* Check for mean*/
	if(xm==0)
	   mean(x,&xm,n);

	if(ym==0)
	   mean(y,&ym,n);
	/* Fit coefficients alpha and beta */
	for(i=0;i<n;i++)
	{
		tmp=(x[i]-xm);
		nsum+=(tmp*(y[i]-ym));
		dsum+=tmp*tmp;
	}
	beta=(dsum==0)?0:(nsum/dsum);
	alpha = ym-(beta *xm);
	lsq[0] = alpha;
	lsq[1] = beta;

	return 0;
}



double getabs(double data)
{
	return (data<0?data*-1:data);
}

int subt(double *inarray,double*dataout,double subtr,int n)
{
	int i=0;
	for (i=0;i<n;i++)
		dataout[i]=inarray[i]- subtr;
	return 0;

}



int abbe(double*data, double *abbe,int n,double datam =0) /*Arun*/
{
	int i;
	double  var=0, denom=0 ;
	/* Return zero if number of data points is null  */
	if(n<=0)
	{
		*abbe=0;
		return -1;
	}
	if(datam==0)
		mean(data,&datam,n);
	for (i=0;i<n-1;i++)
	{
		var += ( pow((double)(data[i+1] - data[i]),2));
		denom += pow((double)(data[i+1] - datam) ,2);
	}
	*abbe = var/(2*denom);
	return 0;
}

double sum(double*a, int n)/*Arun*/
{

   int i;
   double sum=0;
   for (i=0; i<n; i++)
	 sum += a[i];
   return(sum);
}








int Rcode(double *mag,double*mag_err,double*timeo,int n,double*times,double*delta)
{
	int i=0,j=0,cons=1,maxScale=10,*smj,*slj;
	double meanDelta=0,mu=0,*dpss,*lamplus;


	smj=(int*)malloc(sizeof(int)*maxScale);
	slj=(int*)malloc(sizeof(int)*maxScale);
	dpss=(double*)malloc(sizeof(double)*maxScale);
	lamplus=(double*)malloc(sizeof(double)*maxScale);
	if(smj==NULL || slj==NULL || dpss==NULL || lamplus==NULL)
	{
		printf("Error: Out of memory\n");
		exit(0);
	}

	for(i=0;i<n;i++)
		times[i]=timeo[i]-timeo[0];

	for(i=0;i<n-1;i++)
		delta[i]=times[i+1]-times[i];
	meanDelta=times[n-1]/(n-1);
	for(j=0;j<maxScale;j++)
	{
		slj[j]=(int)(pow(2,j+1)*cons);
		 smj[j] = n-slj[j]-1;//this is because R codes modifies n value to n-2
	}
	mu=meanDelta;
	for(i=0;i<maxScale;i++)
	{

	}
	return 0;
}
char *removeFileExt(char* mystr) {
    char *retstr;
  //  char *lastdot;
    if (mystr == NULL)
         return NULL;
    if ((retstr = (char*)malloc (strlen (mystr) + 1)) == NULL)
        return NULL;
    strcpy (retstr, mystr);
    char *lastdot = strrchr (retstr, '.');
    if (lastdot != NULL)
        *lastdot = '\0';
    return retstr;
}
