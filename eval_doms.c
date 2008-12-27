#include <stdio.h>
#include <math.h>
#include "global_matrices.h"

extern void find_coeffs (double, double[], double[][700]);

void eval_doms(double a, double xp, double yp, double phi[])
{
 
	int k, nx, ny; 
	double coeffs[289], xmin, xmax, ymin, ymax, xr, yr, xpow, ypow;
	double sx1, sx2, sy1, sy2, r2;
	double c_far=0.15915494309190;

	r2=(xp*xp+yp*yp);

	if (sqrt(r2)>40*a) /* Use far field solution */
	{	
		phi[0]=c_far*xp/(r2);
		phi[1]=c_far*yp/(r2);
		phi[2]=c_far*(-2.0*xp*yp)/(r2*r2);
		phi[3]=c_far*(yp*yp-xp*xp)/(r2*r2);
		phi[4]=c_far*(xp*xp-yp*yp)/(r2*r2);
		phi[5]=c_far*2.0*xp*(xp*xp-3*yp*yp)/(r2*r2*r2);
		phi[6]=c_far*2.0*yp*(yp*yp-3*xp*xp)/(r2*r2*r2);
		phi[7]=c_far*2.0*yp*(3*xp*xp-yp*yp)/(r2*r2*r2);
		phi[8]=c_far*2.0*xp*(3*yp*yp-xp*xp)/(r2*r2*r2);
	}
	else /* Use near field solution */
	{

	/* Find domains ---------------- */
	if (xp<5*a)
	{
		xmin=0.0; xmax=5.0*a;
		if (yp<0.5)
		{
			find_coeffs (a, coeffs, Md1);
			ymin=0.0; ymax=0.5;
		}
			else if (yp<1.)
			{
				find_coeffs (a, coeffs, Md2);
				ymin=0.5; ymax=1.0;
			}
				else if (yp<3.)
				{
					find_coeffs (a, coeffs, Md3);
					ymin=1.0; ymax=3.0;
				}
					else if (yp<7*sqrt(a))
					{
						find_coeffs (a, coeffs, Md4);
						ymin=3.0; ymax=7.0*sqrt(a);
					}
						else if (yp<11.5*a)
						{
							find_coeffs (a, coeffs, Md5);
							ymin=7.0*sqrt(a); ymax=11.5*a;
						}
							else 
							{
								find_coeffs (a, coeffs, Md12);
								ymin=11.5*a; ymax=40.5*a;
								xmin=0.0; xmax=11.5*a;

							}
	}
		else if (xp<11.5*a)
		{
			xmin=5.0*a; xmax=11.5*a;
			if (yp<0.5)
			{
				find_coeffs (a, coeffs, Md6);
				ymin=0.0; ymax=0.5;
			}
				else if (yp<1)
				{
					find_coeffs (a, coeffs, Md7);
					ymin=0.5; ymax=1.0;
				}
					else if (yp<3.)
					{
						find_coeffs (a, coeffs, Md8);
						ymin=1.0; ymax=3.0;
					}
						else if (yp<7*sqrt(a))
						{
							find_coeffs (a, coeffs, Md9);
							ymin=3.0; ymax=7.0*sqrt(a);
						}
							else if (yp<11.5*a)
							{
								find_coeffs (a, coeffs, Md10);
								ymin=7.0*sqrt(a); ymax=11.5*a;
							}
								else
								{
									find_coeffs (a, coeffs, Md12);
									ymin=11.5*a; ymax=40.5*a;
									xmin=0.0; xmax=11.5*a;
								}
		}
			else if (yp<11.5*a)
			{
				find_coeffs (a, coeffs, Md11);
				ymin=0.0; ymax=11.5*a;
				xmin=11.5*a; xmax=40.5*a;				
			}
				else
				{
					find_coeffs (a, coeffs, Md13);
					ymin=11.5*a; ymax=40.5*a;
					xmin=11.5*a; xmax=40.5*a;				
				}

	
	/******************* Compute derivatives!!! ***************/
	
	/* Rescale */
	sx1=(2.0/(xmax-xmin)); sx2=(xmax+xmin)/(xmax-xmin);
	sy1=(2.0/(ymax-ymin)); sy2=(ymax+ymin)/(ymax-ymin);
	xr=sx1*xp-sx2; yr=sy1*yp-sy2;
	
	/* Phi_x **************************************************/
	xpow=1.0;  phi[0]=0.0; k=17;
	for (nx=1; nx<17; nx++)
	{
		ypow=1.0;
		for (ny=0; ny<17; ny++)
		{
			phi[0]+=coeffs[k]*nx*xpow*ypow;
			ypow*=yr;
			k+=1;
		}
		xpow*=xr;
	}
	phi[0]*=sx1;
	
	/* Phi_y *************************************************/
	xpow=1.0;  phi[1]=0.0; k=0;
	for (nx=0; nx<17; nx++)
	{
		ypow=1.0; k+=1;
		for (ny=1; ny<17; ny++)
		{
			phi[1]+=coeffs[k]*ny*xpow*ypow;
			ypow*=yr;
			k+=1;
		}
		xpow*=xr;
	}
	phi[1]*=sy1;
	
	/* Phi_xy *************************************************/
	xpow=1.0;  phi[2]=0.0; k=17;
	for (nx=1; nx<17; nx++)
	{
		ypow=1.0; k+=1;
		for (ny=1; ny<17; ny++)
		{
			phi[2]+=coeffs[k]*nx*ny*xpow*ypow;
			ypow*=yr;
			k+=1;
		}
		xpow*=xr;
	}
	phi[2]*=sy1*sx1;
	
	/* Phi_xx *************************************************/
	xpow=1.0;  phi[3]=0.0; k=34;
	for (nx=2; nx<17; nx++)
	{
		ypow=1.0;
		for (ny=0; ny<17; ny++)
		{
			phi[3]+=coeffs[k]*nx*(nx-1)*xpow*ypow;
			ypow*=yr;
			k+=1;
		}
		xpow*=xr;
	}
	phi[3]*=sx1*sx1;
	
	/* Phi_yy *************************************************/
	xpow=1.0;  phi[4]=0.0; k=0;
	for (nx=0; nx<17; nx++)
	{
		ypow=1.0; k+=2;
		for (ny=2; ny<17; ny++)
		{
			phi[4]+=coeffs[k]*ny*(ny-1)*xpow*ypow;
			ypow*=yr;
			k+=1;
		}
		xpow*=xr;
	}
	phi[4]*=sy1*sy1;

	/* Phi_xxx *************************************************/
	xpow=1.0;  phi[5]=0.0; k=51;
	for (nx=3; nx<17; nx++)
	{
		ypow=1.0; 
		for (ny=0; ny<17; ny++)
		{
			phi[5]+=coeffs[k]*nx*(nx-1)*(nx-2)*xpow*ypow;
			ypow*=yr;
			k+=1;
		}
		xpow*=xr;
	}
	phi[5]*=sx1*sx1*sx1;


	/* Phi_yyy *************************************************/
	xpow=1.0;  phi[6]=0.0; k=0;
	for (nx=0; nx<17; nx++)
	{
		ypow=1.0; k+=3;
		for (ny=3; ny<17; ny++)
		{
			phi[6]+=coeffs[k]*ny*(ny-1)*(ny-2)*xpow*ypow;
			ypow*=yr;
			k+=1;
		}
		xpow*=xr;
	}
	phi[6]*=sy1*sy1*sy1;

	/* Phi_xxy *************************************************/
	xpow=1.0;  phi[7]=0.0; k=34;
	for (nx=2; nx<17; nx++)
	{
		ypow=1.0; k+=1;
		for (ny=1; ny<17; ny++)
		{
			phi[7]+=coeffs[k]*nx*(nx-1)*(ny)*xpow*ypow;
			ypow*=yr;
			k+=1;
		}
		xpow*=xr;
	}
	phi[7]*=sx1*sx1*sy1;

	/* Phi_xyy *************************************************/
	xpow=1.0;  phi[8]=0.0; k=17;
	for (nx=1; nx<17; nx++)
	{
		ypow=1.0; k+=2;
		for (ny=2; ny<17; ny++)
		{
			phi[8]+=coeffs[k]*ny*(ny-1)*(nx)*xpow*ypow;
			ypow*=yr;
			k+=1;
		}
		xpow*=xr;
	}
	phi[8]*=sx1*sy1*sy1;

	} /* End if r>40a **********************/

}
