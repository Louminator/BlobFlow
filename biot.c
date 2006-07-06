#include <math.h>
#include "biot_global.h"

#define npts_max 1

/*
 ------------------------------------------------------------------------
  Evaluates the Biot-Savart integral (or its dervatives)
 ------------------------------------------------------------------------
 double a -> parameter
 int npts -> number of points for the evaluation
 double x[] -> x[0] to x[npts-1] contains the x values for the evaluationg
 double y[] -> y[0] to y[npts-1] contains the y values for the evaluation 
 int nderx -> number of derivatives to be taken in the x direction
 int ndery -> number of derivatives to be taken in the y direction
 double phi[]-> output:
             phi[0] to phi[npts-1] values of the solution (or its derivatives)


 Note: Current implementation is for
  ************************************************************
                       1/10 < a <10
     |x(k)| < 15*max(a,1/a) and  |y(k)| < 15*max(a,1/a) 
  ************************************************************

 by Rodrigo Platte 06/24/2006
*/

void biot(double a, int npts, double x[], double y[], int nderx, int ndery, double phi[])
{
 
  int k,idom1[npts_max],idom2[npts_max],idom3[npts_max]; 
  
  if (a<1)   /* a=1/a and interchange x and y */
     biot(1/a, npts, y, x, ndery, nderx, phi);
  else
  {
    /* Find domains ---------------- */
    int j1=0., j2=0., j3=0.;
    for (k=0; k<npts; ++k)
	{
	if (fabs(x[k])<=5*a && fabs(y[k])<=5)
		{idom1[j1]=k; j1+=1;}
	else
		{
	    	if (fabs(x[k])>5*a && fabs(y[k])<=5)
			{idom2[j2]=k; j2+=1;}
		else
			{idom3[j3]=k; j3+=1;}
		}
	}

	if (j1>0)
		biot_dom1(a, j1, x, y,nderx , ndery ,phi ,idom1);
	
	if (j2>0)
		biot_dom2(a, j2, x, y,nderx , ndery ,phi ,idom2);
	
	if (j3>0)
		biot_dom3(a, j3, x, y,nderx , ndery ,phi ,idom3);
  }

  }


