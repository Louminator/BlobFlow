#include <math.h>

/*
 ------------------------------------------------------------------------
  Evaluates the derivatives the Biot-Savart integral 
 ------------------------------------------------------------------------
 double a -> parameter
 double x -> the x value for the evaluation
 double y -> the y value for the evaluation 
 double phi[9]-> output:
			phi[0]=phi_x
			phi[1]=phi_y
			phi[2]=phi_xy
			phi[3]=phi_xx
			phi[4]=phi_yy
			phi[5]=phi_xxx
			phi[6]=phi_yyy
			phi[7]=phi_xxy
			phi[8]=phi_xyy

 Note: Current implementation is for
  ************************************************************
                       1/10 < a <10
     |x(k)| < 15*max(a,1/a) and  |y(k)| < 15*max(a,1/a) 
  ************************************************************

 by Rodrigo Platte 09/20/2006
*/

void eval_biot(double a, double x, double y, double phi[9])
{
 
  int k; 
  double temp;

  if (a<1-1e-11)   /* a=1/a and interchange x and y */
  {
	 eval_biot(1/a, y, x, phi);
	 /* Reorganize derivatives */
	 temp=phi[0];
	 phi[0]=phi[1];
	 phi[1]=temp;
	 temp=phi[3];
	 phi[3]=phi[4];
	 phi[4]=temp;
	 temp=phi[5];
	 phi[5]=phi[6];
	 phi[6]=temp;
	 temp=phi[7];
	 phi[7]=phi[8];
	 phi[8]=temp;
  }
  else 
  {
	 eval_doms(a, fabs(x), fabs(y), phi);
	 if (x<0)
	{
		phi[0]=-phi[0];
		phi[2]=-phi[2];
		phi[5]=-phi[5];
		phi[8]=-phi[8];
	}
	if (y<0)
	{
		phi[1]=-phi[1];
		phi[2]=-phi[2];
		phi[6]=-phi[6];
		phi[7]=-phi[7];
	}	
  }
}
