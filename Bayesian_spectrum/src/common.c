/*
 * common.c
 *
 *  Created on: Nov 18, 2014
 *      Author: ajayvibhute
 */




/*This function read the input file and initilize the source list for further processing*/
int getSourceList(char*sourcelistfile,float *tx,float*ty,int*numsources)
{
	FILE *fp;
	int i=0;
	fp=fopen(sourcelistfile,"r");
	if(fp==NULL)
	{
		printf("Error (%s:%d): %s file not found\n",__FILE__,__LINE__,sourcelistfile);
		exit(-1);
	}
	fscanf(fp,"%d",numsources);
	for(i=0;i<*numsources;i++)
		fscanf(fp,"%f %f",&tx[i],&ty[i]);
	fclose(fp);
	return 0;
}
