 /* STREAM_FN.C */
 /* Copyright (c) 2004 Louis F. Rossi                                    *
  * Elliptical Corrected Core Spreading Vortex Method (ECCSVM) for the   *
  * simulation/approximation of incompressible flows.                    *

  * This program is free software; you can redistribute it and/or modify *
  * it under the terms of the GNU General Public License as published by *
  * the Free Software Foundation; either version 2 of the License, or    *
  * (at your option) any later version.					*

  * This program is distributed in the hope that it will be useful,      *
  * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        *
  * GNU General Public License for more details.                         *

  * You should have received a copy of the GNU General Public License    *
  * along with this program; if not, write to the Free Software          *
  * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 *
  * USA                                                                  *
  *                                                                      *
  * Louis Rossi                                                          *
  * Department of Mathematical Sciences                                  *
  * University of Delaware                                               *
  * Newark, DE 19716                                                     */

 #include <stdio.h>
 #include <math.h>
 #include "global.h"

/* This is (3.47) of UDel Math Tech Report 2001-4 */
double expint(c,s2,r2)
     double c[maxpolyn],s2,r2;
{
  return(-0.5/s2*
	 (r2*((c[0]+s2*(8.0*c[1]+s2*(96.0*c[2]+s2*(1536.0*c[3]+30720.0*c[4]*s2))))+
	      r2*((c[1]+s2*(12.0*c[2]+s2*(192.0*c[3]+3840.0*c[4]*s2)))+
		  r2*((c[2]+s2*(16.0*c[3]+320.0*c[4]*s2))+
		      r2*((c[3]+20.0*c[4]*s2)+
			  c[4]*r2))))*exp(-r2/(4.0*s2))+
	  (-s2)*(4.0*c[0]+s2*(32.0*c[1]+
			    s2*(384.0*c[2]+s2*(6144.0*c[3]+s2*122880.0*c[4]))))*
	  (1.0-exp(-r2/(4.0*s2)))));
}

void build_rt(dx,dy,eps,r,t)
     double dx,dy,eps,r[maxexp],t[maxexp];
{
  int i;

  r[0] = SQR(dx*(1-eps)/(1+eps))+SQR(dy*(1+eps)/(1-eps));
  t[0] = (SQR(dx*(1-eps)/(1+eps))-SQR(dy*(1+eps)/(1-eps)))/r[0];
   
  r[1] = SQR(r[0]);
  t[1] = SQR(t[0]);

  for (i=2; i<maxexp; ++i)
    {
      r[i] = r[i-1]*r[0];
      t[i] = t[i-1]*t[0];
    }
}     

void build_psiRT(dx,dy,eps,r,t,RT)
  double dx,dy,eps,r[maxexp],t[maxexp],RT[maxpolyn][maxexp];
{
  int p,m;

  build_rt(dx,dy,eps,r,t);

  for (p=0; p<maxpolyn; ++p)
    for (m=0; m<maxexp; ++m)
      RT[p][m] = 0.0;

  RT[0][0] = 0.25 - 4.0*SQR(SQR(eps))+SQR(eps);
  RT[0][1] = eps - 5.0*SQR(eps)*eps;
  RT[0][2] = 20.0*SQR(SQR(eps)) - 2.0*SQR(eps);
  RT[0][3] = 16.0/3.0*SQR(eps)*eps;
  RT[0][4] = -16.0*SQR(SQR(eps));

  RT[1][0] = -2.0*SQR(eps) + 20.0*SQR(SQR(eps));
  RT[1][1] = 37.0/2.0*SQR(eps)*eps - eps/2.0;
  RT[1][2] = 4.0*SQR(eps) - 136.0*SQR(SQR(eps));
  RT[1][3] = -24.0*SQR(eps)*eps;
  RT[1][4] = 128.0*SQR(SQR(eps));

  RT[2][0] = SQR(eps) - 42.0*SQR(SQR(eps));
  RT[2][1] = -24.0*SQR(eps)*eps;
  RT[2][2] = -2.0*SQR(eps) + 324.0*SQR(SQR(eps));
  RT[2][3] = 32.0*SQR(eps)*eps;
  RT[2][4] = -320.0*SQR(SQR(eps));

  RT[3][0] = 40.0*SQR(SQR(eps));
  RT[3][1] = 10.0*SQR(eps)*eps;
  RT[3][2] = -320.0*SQR(SQR(eps));
  RT[3][3] = -40.0/3.0*SQR(eps)*eps;
  RT[3][4] = 320.0*SQR(SQR(eps));

  RT[4][0] = -14.0*SQR(SQR(eps));
  RT[4][2] = 112.0*SQR(SQR(eps));
  RT[4][4] = -112.0*SQR(SQR(eps));

  /* Fifth order corrections */

#if PSI_ORDER > 4
  RT[0][1] +=  29.0*eps*SQR(SQR(eps));
  RT[0][3] += -80.0*eps*SQR(SQR(eps));
  RT[0][5] +=  256.0/5.0*eps*SQR(SQR(eps));

  RT[1][1] += -509.0/2.0*eps*SQR(SQR(eps));
  RT[1][3] +=  872.0*eps*SQR(SQR(eps));
  RT[1][5] += -640*eps*SQR(SQR(eps));

  RT[2][1] +=  872.0*eps*SQR(SQR(eps));
  RT[2][3] += -3296.0*eps*SQR(SQR(eps));
  RT[2][5] +=  2560.0*eps*SQR(SQR(eps));

  RT[3][1] += -1430.0*eps*SQR(SQR(eps));
  RT[3][3] +=  5640.0*eps*SQR(SQR(eps));
  RT[3][5] += -4480.0*eps*SQR(SQR(eps));

  RT[4][1] +=  1120.0*eps*SQR(SQR(eps));
  RT[4][3] += -4480.0*eps*SQR(SQR(eps));
  RT[4][5] +=  3584.0*eps*SQR(SQR(eps));

  RT[5][1]  = -336.0*eps*SQR(SQR(eps));
  RT[5][3]  =  1344.0*eps*SQR(SQR(eps));
  RT[5][5]  = -5376.0/5.0*eps*SQR(SQR(eps));
#endif

  /* Sixth order corrections */

#if PSI_ORDER > 5
  RT[0][0] += 49.0/3.0*SQR(SQR(eps))*SQR(eps);
  RT[0][2] += -166.0*SQR(SQR(eps))*SQR(eps);
  RT[0][4] += 320.0*SQR(SQR(eps))*SQR(eps);
  RT[0][6] += -512.0/3.0*SQR(SQR(eps))*SQR(eps);

  RT[1][0] += -166.0*SQR(SQR(eps))*SQR(eps);
  RT[1][2] += 2252.0*SQR(SQR(eps))*SQR(eps);
  RT[1][4] += -5120.0*SQR(SQR(eps))*SQR(eps);
  RT[1][6] += 3072.0*SQR(SQR(eps))*SQR(eps);

  RT[2][0] += 723.0*SQR(SQR(eps))*SQR(eps);
  RT[2][2] += -11366.0*SQR(SQR(eps))*SQR(eps);
  RT[2][4] += 28160.0*SQR(SQR(eps))*SQR(eps);
  RT[2][6] += -17920.0*SQR(SQR(eps))*SQR(eps);

  RT[3][0] += -4960/3.0*SQR(SQR(eps))*SQR(eps);
  RT[3][2] += 28160.0*SQR(SQR(eps))*SQR(eps);
  RT[3][4] += -72960.0*SQR(SQR(eps))*SQR(eps);
  RT[3][6] += 143360.0/3.0*SQR(SQR(eps))*SQR(eps);

  RT[4][0] +=  2072.0*SQR(SQR(eps))*SQR(eps);
  RT[4][2] += -36736.0*SQR(SQR(eps))*SQR(eps);
  RT[4][4] +=  97216.0*SQR(SQR(eps))*SQR(eps);
  RT[4][6] += -64512.0*SQR(SQR(eps))*SQR(eps);

  RT[5][0] += -1344.0*SQR(SQR(eps))*SQR(eps);
  RT[5][2] +=  24192.0*SQR(SQR(eps))*SQR(eps);
  RT[5][4] += -64512.0*SQR(SQR(eps))*SQR(eps);
  RT[5][6] +=  43008.0*SQR(SQR(eps))*SQR(eps);

  RT[6][0] =  352.0*SQR(SQR(eps))*SQR(eps);
  RT[6][2] = -6336.0*SQR(SQR(eps))*SQR(eps);
  RT[6][4] =  16896.0*SQR(SQR(eps))*SQR(eps);
  RT[6][6] = -11264.0*SQR(SQR(eps))*SQR(eps);
#endif
}

void build_psi(dx,dy,eps,r,t,psi_c,psi_RT)
  double dx,dy,eps,r[maxexp],t[maxexp];
  double psi_c[maxpolyn],psi_RT[maxpolyn][maxexp];
{
  int p,m;

  /*** Build psi polynomial ***/

  build_psiRT(dx,dy,eps,r,t,psi_RT);

  for (p=0; p<maxpolyn; ++p)
    {
      psi_c[p] = 0.0;
      for (m=0; m<maxexp; ++m)
	{
	  if (p==0)
	    if (m == 0)
	      psi_c[p] += psi_RT[p][m];
	    else
	      psi_c[p] += psi_RT[p][m]*t[m-1];
	  else
	    if (m == 0)
	      psi_c[p] += psi_RT[p][m]/r[p-1];
	    else
	      psi_c[p] += psi_RT[p][m]*t[m-1]/r[p-1];
	}
    }
  
  /*** End psi calculation ***/

}



/* These results must be multiplied by dx */
void dx_coeff(oldRT,RT,offset_exp,a2)
     double oldRT[maxexp],RT[maxexp],a2;
     int    offset_exp;
{
  int i,m;

  for (i=0; i<maxexp; ++i)
    RT[i] = 0.0;

  for (m=0; m<maxexp; ++m)
    {
      if (m>0)
	RT[m-1] += (2.0/a2)*m*oldRT[m];
      RT[m] -= (2.0/a2)*(m+offset_exp)*oldRT[m];
    }
}

/* These results must be multiplied by dy */
void dy_coeff(oldRT,RT,offset_exp,a2)
     double oldRT[maxexp],RT[maxexp],a2;
     int    offset_exp;
{
  int i,m;

  for (i=0; i<maxexp; ++i)
    RT[i] = 0.0;

  for (m=0; m<maxexp; ++m)
    {
      if (m>0)
	RT[m-1] -= (2.0*a2)*m*oldRT[m];
      RT[m] -= (2.0*a2)*(m+offset_exp)*oldRT[m];
    }
}

void dxx_coeff(oldRT,RT,offset_exp,a2)
     double oldRT[maxexp],RT[maxexp],a2;
     int    offset_exp;
{
  int i,m;

  for (i=0; i<maxexp; ++i)
    RT[i] = 0.0;

  for (m=0; m<maxexp; ++m)
    {
      if (m>1)
	RT[m-2] += (2.0/a2)*m*(m-1)*oldRT[m];
      if (m>0)
	RT[m-1] += (2.0/a2)*m*(m-2*(m+offset_exp))*oldRT[m];
      RT[m] += (2.0/a2)*(m+offset_exp)*(m+offset_exp-2*m)*oldRT[m];
      if (m < maxexp-1)
	RT[m+1] += (2.0/a2)*(m+offset_exp)*(m+offset_exp+1)*oldRT[m];
    }
}

void dyy_coeff(oldRT,RT,offset_exp,a2)
     double oldRT[maxexp],RT[maxexp],a2;
     int    offset_exp;
{
  int i,m;

  for (i=0; i<maxexp; ++i)
    RT[i] = 0.0;

  for (m=0; m<maxexp; ++m)
    {
      if (m>1)
	RT[m-2] += 2.0*a2*m*(m-1)*oldRT[m];
      if (m>0)
	RT[m-1] += 2.0*a2*m*(2*(m+offset_exp)-m)*oldRT[m];
      RT[m] += 2.0*a2*(m+offset_exp)*(m+offset_exp-2*m)*oldRT[m];
      if (m < maxexp-1)
	RT[m+1] -= 2.0*a2*(m+offset_exp)*(m+offset_exp+1)*oldRT[m];
    }
}

/* These results must be multiplied by dx*dy */
void dxy_coeff(oldRT,RT,offset_exp,a2)
     double oldRT[maxexp],RT[maxexp],a2;
     int    offset_exp;
{
  int i,m;

  for (i=0; i<maxexp; ++i)
    RT[i] = 0.0;

  for (m=0; m<maxexp; ++m)
    {
      if (m>1)
	RT[m-2] -= 4.0*m*(m-1)*oldRT[m];
      RT[m] += 4.0*(m+offset_exp)*(m+offset_exp+1)*oldRT[m];
    }
}

double build_psi_x(dx,dy,str,s2,s4,a2,eps,r,t,psi_c,psi_x_c,psi_RT,psi_x_RT)
  double dx,dy,str,s2,s4,a2,eps,r[maxexp],t[maxexp];
  double psi_c[maxpolyn],psi_x_c[maxpolyn];
  double psi_RT[maxpolyn][maxexp],psi_x_RT[maxpolyn][maxexp];
{
  int p,m;
  double psi_x;

  /*** Build psi_x polynomial ***/

  for (p=0; p<maxpolyn; ++p)
    dx_coeff(psi_RT[p],psi_x_RT[p],p,a2);

  /* Don't forget the log term. */
  psi_x_RT[0][0] += 0.5/a2;

  for (p=0; p<maxpolyn; ++p)
    {
      psi_x_c[p] = 0.0;
      for (m=0; m<maxexp; ++m)
	{
	  if (m == 0)
	    psi_x_c[p] += psi_x_RT[p][m]/r[p];
	  else
	    psi_x_c[p] += psi_x_RT[p][m]*t[m-1]/r[p];
	}
    }
  
  for (p=0; p<maxpolyn; ++p)
    psi_x_c[p] *= dx;
  
  psi_x = str*expint(psi_x_c,s2,r[0]);

  /* Plus the psi_2 component */

  psi_x += str/2.0/s2*exp(-r[0]/4.0/s2)*
    dx*psi2coeffA;

  /*** End psi_x calculation ***/
  return(psi_x);

}

double build_psi_y(dx,dy,str,s2,s4,a2,eps,r,t,psi_c,psi_y_c,psi_RT,psi_y_RT)
  double dx,dy,str,s2,s4,a2,eps,r[maxexp],t[maxexp];
  double psi_c[maxpolyn],psi_y_c[maxpolyn];
  double psi_RT[maxpolyn][maxexp],psi_y_RT[maxpolyn][maxexp];
{
  int p,m;
  double psi_y;

  /*** Build psi_y polynomial ***/

  for (p=0; p<maxpolyn; ++p)
    dy_coeff(psi_RT[p],psi_y_RT[p],p,a2);

  /* Don't forget the log term. */
  psi_y_RT[0][0] += 0.5*a2;

  for (p=0; p<maxpolyn; ++p)
    {
      psi_y_c[p] = 0.0;
      for (m=0; m<maxexp; ++m)
	{
	  if (m == 0)
	    psi_y_c[p] += psi_y_RT[p][m]/r[p];
	  else
	    psi_y_c[p] += psi_y_RT[p][m]*t[m-1]/r[p];
	}
    }
  
  for (p=0; p<maxpolyn; ++p)
    psi_y_c[p] *= dy;

  psi_y = str*expint(psi_y_c,s2,r[0]);

  /* Plus the psi_2 component */

  psi_y += str/2.0/s2*exp(-r[0]/4.0/s2)*
    dy*psi2coeffB;

  /*** End psi_y calculation ***/

  return(psi_y);
}

double build_psi_xx(dx,dy,str,s2,s4,a2,eps,r,t,
		    psi_x_c,psi_xx_c,psi_RT,psi_xx_RT)
  double dx,dy,str,s2,s4,a2,eps,r[maxexp],t[maxexp];
  double psi_x_c[maxpolyn],psi_xx_c[maxpolyn];
  double psi_RT[maxpolyn][maxexp],psi_xx_RT[maxpolyn][maxexp];
{
  int p,m;
  double psi_xx;

  /*** Build psi_xx polynomial ***/

  for (p=0; p<maxpolyn; ++p)
    dxx_coeff(psi_RT[p],psi_xx_RT[p],p,a2);

  /* Don't forget the log term. */
  psi_xx_RT[0][1] -= 0.5/a2;

  for (p=0; p<maxpolyn; ++p)
    {
      psi_xx_c[p] = 0.0;
      for (m=0; m<maxexp; ++m)
	{
	  if (m == 0)
	    psi_xx_c[p] += psi_xx_RT[p][m]/r[p];
	  else
	    psi_xx_c[p] += psi_xx_RT[p][m]*t[m-1]/r[p];
	}
    }
  
  psi_xx = str*expint(psi_xx_c,s2,r[0]);

  /* Plus the psi_2 component */

  psi_xx += str/2.0/s2*exp(-r[0]/4.0/s2)*psi2coeffA;

  /*** End psi_xx calculation ***/

  return(psi_xx);
}

double build_psi_yy(dx,dy,str,s2,s4,a2,eps,r,t,
		    psi_y_c,psi_yy_c,psi_RT,psi_yy_RT)
  double dx,dy,str,s2,s4,a2,eps,r[maxexp],t[maxexp];
  double psi_y_c[maxpolyn],psi_yy_c[maxpolyn];
  double psi_RT[maxpolyn][maxexp],psi_yy_RT[maxpolyn][maxexp];
{
  int p,m;
  double psi_yy;

  /*** Build psi_yy polynomial ***/

  for (p=0; p<maxpolyn; ++p)
    dyy_coeff(psi_RT[p],psi_yy_RT[p],p,a2);

  /* Don't forget the log term. */
  psi_yy_RT[0][1] += 0.5*a2;

  for (p=0; p<maxpolyn; ++p)
    {
      psi_yy_c[p] = 0.0;
      for (m=0; m<maxexp; ++m)
	{
	  if (m == 0)
	    psi_yy_c[p] += psi_yy_RT[p][m]/r[p];
	  else
	    psi_yy_c[p] += psi_yy_RT[p][m]*t[m-1]/r[p];
	}
    }
  
  psi_yy = str*expint(psi_yy_c,s2,r[0]);

  /* Plus the psi_2 component */

  psi_yy += str/2.0/s2*exp(-r[0]/4.0/s2)*psi2coeffB;

  /*** End psi_yy calculation ***/

  return(psi_yy);
}

double build_psi_xy(dx,dy,str,s2,s4,a2,eps,r,t,
		    psi_y_c,psi_xy_c,psi_RT,psi_xy_RT)
  double dx,dy,str,s2,s4,a2,eps,r[maxexp],t[maxexp];
  double psi_y_c[maxpolyn],psi_xy_c[maxpolyn];
  double psi_RT[maxpolyn][maxexp],psi_xy_RT[maxpolyn][maxexp];
{
  int p,m;
  double psi_xy;

  /*** Build psi_xy polynomial ***/

  for (p=0; p<maxpolyn; ++p)
    dxy_coeff(psi_RT[p],psi_xy_RT[p],p,a2);

  /* Don't forget the log term. */
  psi_xy_RT[0][0] -= 1.0;

  for (p=0; p<maxpolyn; ++p)
    {
      psi_xy_c[p] = 0.0;
      for (m=0; m<maxexp; ++m)
	{
	  if (m == 0)
	    psi_xy_c[p] += psi_xy_RT[p][m]/r[p]/r[0];
	  else
	    psi_xy_c[p] += psi_xy_RT[p][m]*t[m-1]/r[p]/r[0];
	}
    }

  for (p=0; p<maxpolyn; ++p)
    psi_xy_c[p] *= dx*dy;

  psi_xy = str*expint(psi_xy_c,s2,r[0]);

  /*** End psi_xy calculation ***/

  
  return(psi_xy);
}

double build_psi_xxx(dx,dy,str,s2,s4,a2,eps,r,t,
		    psi_xx_c,psi_xxx_c,psi_xx_RT,psi_xxx_RT)
  double dx,dy,str,s2,s4,a2,eps,r[maxexp],t[maxexp];
  double psi_xx_c[maxpolyn],psi_xxx_c[maxpolyn];
  double psi_xx_RT[maxpolyn][maxexp],psi_xxx_RT[maxpolyn][maxexp];
{
  int p,m;
  double psi_xxx;

  /*** Build psi_xxx polynomial ***/

  for (p=0; p<maxpolyn; ++p)
    dx_coeff(psi_xx_RT[p],psi_xxx_RT[p],p+1,a2);

  for (p=0; p<maxpolyn; ++p)
    {
      psi_xxx_c[p] = 0.0;
      for (m=0; m<maxexp; ++m)
	{
	  if (m == 0)
	    psi_xxx_c[p] += psi_xxx_RT[p][m]/r[p]/r[0];
	  else
	    psi_xxx_c[p] += psi_xxx_RT[p][m]*t[m-1]/r[p]/r[0];
	}
    }
  
  for (p=0; p<maxpolyn; ++p)
    psi_xxx_c[p] *= dx;
  
  psi_xxx = str*expint(psi_xxx_c,s2,r[0]);

  /* Leibnitz rule correction for axisymmetric solution.*/

  psi_xxx += str/4.0/s4*exp(-r[0]/4.0/s2)*
    dx/a2*
    (psi_xx_c[0]*r[0]+psi_xx_c[1]*r[1]+psi_xx_c[2]*r[2]+
     psi_xx_c[3]*r[3]+psi_xx_c[4]*r[4]-
     psi2coeffA);

  /*** End psi_xxx calculation ***/

  return(psi_xxx);
}

double build_psi_xxy(dx,dy,str,s2,s4,a2,eps,r,t,
		    psi_xx_c,psi_xxy_c,psi_xx_RT,psi_xxy_RT)
  double dx,dy,str,s2,s4,a2,eps,r[maxexp],t[maxexp];
  double psi_xx_c[maxpolyn],psi_xxy_c[maxpolyn];
  double psi_xx_RT[maxpolyn][maxexp],psi_xxy_RT[maxpolyn][maxexp];
{
  int p,m;
  double psi_xxy;

 /*** Build psi_xxy polynomial ***/

  for (p=0; p<maxpolyn; ++p)
    dy_coeff(psi_xx_RT[p],psi_xxy_RT[p],p+1,a2);

  for (p=0; p<maxpolyn; ++p)
    {
      psi_xxy_c[p] = 0.0;
      for (m=0; m<maxexp; ++m)
	{
	  if (m == 0)
	    psi_xxy_c[p] += psi_xxy_RT[p][m]/r[p]/r[0];
	  else
	    psi_xxy_c[p] += psi_xxy_RT[p][m]*t[m-1]/r[p]/r[0];
	}
    }
  
  for (p=0; p<maxpolyn; ++p)
    psi_xxy_c[p] *= dy;
  
  psi_xxy = str*expint(psi_xxy_c,s2,r[0]);

  /* Leibnitz rule correction for axisymmetric solution.*/

  psi_xxy += str/4.0/s4*exp(-r[0]/4.0/s2)*
    dy*a2*
    (psi_xx_c[0]*r[0]+psi_xx_c[1]*r[1]+psi_xx_c[2]*r[2]+
     psi_xx_c[3]*r[3]+psi_xx_c[4]*r[4]-
     psi2coeffA);

  /*** End psi_xxy calculation ***/

  return(psi_xxy);
}

double build_psi_xyy(dx,dy,str,s2,s4,a2,eps,r,t,
		    psi_yy_c,psi_xyy_c,psi_yy_RT,psi_xyy_RT)
  double dx,dy,str,s2,s4,a2,eps,r[maxexp],t[maxexp];
  double psi_yy_c[maxpolyn],psi_xyy_c[maxpolyn];
  double psi_yy_RT[maxpolyn][maxexp],psi_xyy_RT[maxpolyn][maxexp];
{
  int p,m;
  double psi_xyy;

  /*** Build psi_xyy polynomial ***/

  for (p=0; p<maxpolyn; ++p)
    dx_coeff(psi_yy_RT[p],psi_xyy_RT[p],p+1,a2);

  for (p=0; p<maxpolyn; ++p)
    {
      psi_xyy_c[p] = 0.0;
      for (m=0; m<maxexp; ++m)
	{
	  if (m == 0)
	    psi_xyy_c[p] += psi_xyy_RT[p][m]/r[p]/r[0];
	  else
	    psi_xyy_c[p] += psi_xyy_RT[p][m]*t[m-1]/r[p]/r[0];
	}
    }
  
  for (p=0; p<maxpolyn; ++p)
    psi_xyy_c[p] *= dx;
  
  psi_xyy = str*expint(psi_xyy_c,s2,r[0]);

  /* Leibnitz rule correction for axisymmetric solution.*/

  psi_xyy += str/4.0/s4*exp(-r[0]/4.0/s2)*
    dx/a2*
    (psi_yy_c[0]*r[0]+psi_yy_c[1]*r[1]+psi_yy_c[2]*r[2]+
     psi_yy_c[3]*r[3]+psi_yy_c[4]*r[4]-
     psi2coeffB);

  /*** End psi_xyy calculation ***/

  return(psi_xyy);
}

double build_psi_yyy(dx,dy,str,s2,s4,a2,eps,r,t,
		    psi_yy_c,psi_yyy_c,psi_yy_RT,psi_yyy_RT)
  double dx,dy,str,s2,s4,a2,eps,r[maxexp],t[maxexp];
  double psi_yy_c[maxpolyn],psi_yyy_c[maxpolyn];
  double psi_yy_RT[maxpolyn][maxexp],psi_yyy_RT[maxpolyn][maxexp];
{
  int p,m;
  double psi_yyy;

  /*** Build psi_yyy polynomial ***/

  for (p=0; p<maxpolyn; ++p)
    dy_coeff(psi_yy_RT[p],psi_yyy_RT[p],p+1,a2);

  for (p=0; p<maxpolyn; ++p)
    {
      psi_yyy_c[p] = 0.0;
      for (m=0; m<maxexp; ++m)
	{
	  if (m == 0)
	    psi_yyy_c[p] += psi_yyy_RT[p][m]/r[p]/r[0];
	  else
	    psi_yyy_c[p] += psi_yyy_RT[p][m]*t[m-1]/r[p]/r[0];
	}
    }
  
  for (p=0; p<maxpolyn; ++p)
    psi_yyy_c[p] *= dy;
  
  psi_yyy = str*expint(psi_yyy_c,s2,r[0]);

  /* Leibnitz rule correction for axisymmetric solution.*/

  psi_yyy += str/4.0/s4*exp(-r[0]/4.0/s2)*
    dy*a2*
    (psi_yy_c[0]*r[0]+psi_yy_c[1]*r[1]+psi_yy_c[2]*r[2]+
     psi_yy_c[3]*r[3]+psi_yy_c[4]*r[4]-
     psi2coeffB);

  /*** End psi_yyy calculation ***/

  return(psi_yyy);
}

