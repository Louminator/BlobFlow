
// ------------------------------------------------------------------------
//  Evaluates the solution (or its dervatives) of the Biot-Savart integral
// 			DOMAIN 2
// ------------------------------------------------------------------------
//
//
//
// by Rodrigo Platte 06/24/2006

#include <math.h>
#include "biot_global.h"


void biot_dom2 (double a, int npts, double x[], double y[], int nderx, int ndery, double phi[],
                int idom[])
{
  
  double Cheb(double z, double acosz, double sqrt1z, int d, int nder);

  double c[1000], Ak[1000];
  static double a1_2_old=0.0, a2_2_old=0.0;
  double sx1, sx2, sy1, sy2, xa, xb, ya, yb, coefder, z, acosz, sqrt1z;
  int n,m,k,i,j;
  
	if (a<=1.5)
		{  
		if ((a1_2_old>a)||(a2_2_old<a)) data_1a1p5_dom2();
		a1_2_old=1.; a2_2_old=1.5;
		xa=5*a1_2_old-0.05*a2_2_old; xb=15.1*a2_2_old; ya=-0.1; yb=5.1;
		}
	else 
		if (a<=2)
		{
		if ((a1_2_old>a)||(a2_2_old<a)) data_1p5a2_dom2();
		a1_2_old=1.5; a2_2_old=2.;
		xa=5*a1_2_old-0.05*a2_2_old; xb=15.1*a2_2_old; ya=-0.1; yb=5.1;
		}
		else
			if (a<=3)
			{
			if ((a1_2_old>a)||(a2_2_old<a)) data_2a3_dom2();
			a1_2_old=2; a2_2_old=3;
			xa=5*a1_2_old-0.05*a2_2_old; xb=15.1*a2_2_old; ya=-0.1; yb=5.1;
			}
			else
				if (a<=4)
				{
				if ((a1_2_old>a)||(a2_2_old<a)) data_3a4_dom2();
				a1_2_old=3; a2_2_old=4;
				xa=5*a1_2_old-0.05*a2_2_old; xb=15.1*a2_2_old; ya=-0.1; yb=5.1;
				}
				else
					if (a<=5)
					{
					if ((a1_2_old>a)||(a2_2_old<a)) data_4a5_dom2();
					a1_2_old=4; a2_2_old=5;
					xa=5*a1_2_old-0.05*a2_2_old; xb=15.1*a2_2_old; ya=-0.1; yb=5.1;
					}
					else
						if (a<=6)
						{
						if ((a1_2_old>a)||(a2_2_old<a)) data_5a6_dom2();
						a1_2_old=5; a2_2_old=6;
						xa=5*a1_2_old-0.05*a2_2_old; xb=15.1*a2_2_old; ya=-0.1; yb=5.1;
						}
						else
							if (a<=8)
							{
							if ((a1_2_old>a)||(a2_2_old<a)) data_6a8_dom2();
							a1_2_old=6; a2_2_old=8;
							xa=5*a1_2_old-0.05*a2_2_old; xb=15.1*a2_2_old; ya=-0.1; yb=5.1;
							}
							else
							//	if (a<=10)
								{
								if ((a1_2_old>a)||(a2_2_old<a)) data_8a10_dom2();
								a1_2_old=8; a2_2_old=10;
								xa=5*a1_2_old-0.05*a2_2_old; xb=15.1*a2_2_old; ya=-0.1; yb=5.1;
								}
      

  n=size2[0];
  m=size2[1];
  
  for (k=0; k<=m; ++k)
    {
      c[k]=Coefs2[0][k];
      for (j=1; j<=n; ++j)
			c[k]+=cos(j*acos(sa2[0]*a-sa2[1]))*Coefs2[j][k];
    }

sx1=2/(xb-xa);
sx2=(xb+xa)/(xb-xa);
sy1=2/(yb-ya);
sy2=(yb+ya)/(yb-ya);
  

for (i=0; i<npts; ++i)
	{

  // x-part
	if ( nderx==0 )	
		{
		acosz=acos(sx1*fabs(x[idom[i]])-sx2);
		for (k=0; k<=m; ++k) Ak[k]=cos(Deg2[k][0]*acosz);
		}
	else
		{
		if (x[idom[i]]>0)
			coefder=pow(sx1,nderx);
		else
			coefder=pow(-sx1,nderx);
		z=(sx1*fabs(x[idom[i]])-sx2);
		acosz=acos(z);
		sqrt1z=sqrt(1-z*z);
		for (k=0; k<=m; ++k)
			Ak[k]=coefder*Cheb(z,acosz,sqrt1z,Deg2[k][0],nderx);
		}

 // y-part
	if ( ndery==0 )
		{
		acosz=acos(sy1*fabs(y[idom[i]])-sy2);
		for (k=0; k<=m; ++k) Ak[k]*=cos(Deg2[k][1]*acosz);
		}
	else
		{
		if (y[idom[i]]>0)
			coefder=pow(sy1,ndery);
		else
			coefder=pow(-sy1,ndery);
		z=(sy1*fabs(y[idom[i]])-sy2);
		acosz=acos(z);
		sqrt1z=sqrt(1-z*z);
		for (k=0; k<=m; ++k)
			Ak[k]*=coefder*Cheb(z,acosz,sqrt1z,Deg2[k][1],ndery);
		}		

	phi[idom[i]]=Ak[0]*c[0];
	for (k=1;k<=m; ++k) 
		phi[idom[i]]+=Ak[k]*c[k];    

    	}
}

//------------------------------------------------------------------------
