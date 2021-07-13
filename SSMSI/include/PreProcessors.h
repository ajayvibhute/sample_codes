/*
 This  file contains all the # defines which common two or more files.
 
*/
#define BUFSIZE 1024            //size of the buffer
#define rows 241*NUM_WIRES
#define colms 2
#define xstart -13              //starting of thetaX
#define xend 13	              //end of the thetax
#define ystart -40           //starting of thetaY
#define yend 40             //end of the thetaY




float xdiff=0.1;         //step size for thetaX while iterating from xstart to xend
float   ydiff=1 ;     //step size for thetaY while iterating from ystart to yend  
int xgridSize;     	//grid size for x ((xend-xstart)/xdiff)+1
int ygridSize;	       //grid size for y ((yend-ystart)/ydiff)+1
int noOfDataSet;      // Total no of images to be fitted
struct sourceInfo
{
	float x_position,y_position,sourceStrength;
};
struct InfoData
{
	int myrank;
	int xs,xe,ys,ye,numSources;
	int icam;
	struct sourceInfo sources[100];
	char observationId[100];
};


