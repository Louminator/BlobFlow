/* SOLVE.C */
/* Copyright (c) 2000 Louis F. Rossi                                    *
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

#include "global.h"

#define axisymmtol 1.0e-3
#define stifftol 1.0e-6
#define rk4itertol 14
#define rk4l2errtol 1.0e-14
#define rkckreltol 1.0e-7
#define rkck_err_dereg 1.0e-12

void set_blob(the_blobguts,parms)
Blob_internal *the_blobguts;
Blob_parms *parms;
{
    (*parms).costh = cos((*the_blobguts).th);
    (*parms).sinth = sin((*the_blobguts).th);
    (*parms).cos2 = SQR((*parms).costh);
    (*parms).sin2 = SQR((*parms).sinth);
    (*parms).sincos = (*parms).sinth*(*parms).costh;
}

void th_slave(blobguts,parms)
Blob_internal *blobguts;
Blob_parms     *parms;
{
  double eps,th1;

  if ( ((*blobguts).a2-1/(*blobguts).a2)<0 )
    (*blobguts).th = 
      atan2(-(*parms).du11+
	    sqrt(SQR((*parms).du11)+
		 SQR((*parms).du12+(*parms).du21)/4),
	    ((*parms).du12+(*parms).du21)/2);
  else
    (*blobguts).th = 
      atan2(-(*parms).du11-
	    sqrt(SQR((*parms).du11)+
		 SQR((*parms).du12+(*parms).du21)/4),
	    ((*parms).du12+(*parms).du21)/2);

  eps = (1/(*blobguts).a2-(*blobguts).a2)/(1/(*blobguts).a2+(*blobguts).a2);

  th1 = -((*parms).du21-(*parms).du12)/2/
    (4*sin((*blobguts).th)*cos((*blobguts).th)*
     ((*parms).du12+(*parms).du21)/2*
     (1+SQR((*parms).du11*2/(((*parms).du12+(*parms).du21)))));

  /* Check for special case where correction is second order. */

  if ( (isinf(th1) == 0) && (isnan(th1) == 0) )
    (*blobguts).th = (*blobguts).th + eps*th1 + 2/3*eps*SQR(eps)*th1*SQR(th1);

  set_blob(blobguts,parms);
}

void dy_slave(blobguts,parms,ds2,da2)
Blob_internal *blobguts;
Blob_parms     *parms;
double        *ds2,*da2;
{
  *ds2 = dts2(blobguts);
  *da2 = dta2(blobguts,parms);
}

void dy(blobguts,parms,ds2,da2,dth,prefstep)
Blob_internal *blobguts;
Blob_parms     *parms;
double        *ds2,*da2,*dth,prefstep;
{
  *dth = dtth(blobguts,parms,prefstep,axisymmtol);
  *ds2 = dts2(blobguts);
  *da2 = dta2(blobguts,parms);
}

void rk4(double prefstep)
{
   Blob_internal  tempguts[NMAX];
   Blob_parms     tempparms[NMAX];
   double         ds2[4][NMAX],da2[4][NMAX],dth[4][NMAX],
                  dxdt[4][NMAX],dydt[4][NMAX];
   int            j;

  /* Take a step of RK4 */

  /* Note: blob1, tempparms and tempguts stores the initial velocity field at t=0. */

  printf("Stage 1a\n");


  printf("Stage 1b\n");

  write_vorts(9990);

  vel_field();

  for (j=0; j<N; ++j)
    if (fabs(blobguts[j].a2-1/blobguts[j].a2)/prefstep < axisymmtol)
      th_slave(blobguts+j,tmpparms+j);

  printf("Stage 1c\n");

  for (j=0; j<N; ++j)
    {
      mblob[j].blob1 = mblob[j].blob0;
      tempguts[j]    = blobguts[j];
      tempparms[j]   = tmpparms[j];
    }

  /* Half step of FCE predictor. */

  for (j=0; j<N; ++j)
    {
      dxdt[0][j] = mblob[j].blob0.dx;
      dydt[0][j] = mblob[j].blob0.dy;
      set_blob(blobguts+j,tmpparms+j);
      dy(blobguts+j,tmpparms+j,ds2[0]+j,da2[0]+j,dth[0]+j,prefstep);
    }

  printf("Stage 2a\n");

  for (j=0; j<N; ++j)
    {
      mblob[j].blob0.x = mblob[j].blob1.x + 0.5*prefstep*dxdt[0][j];
      mblob[j].blob0.y = mblob[j].blob1.y + 0.5*prefstep*dydt[0][j];
      blobguts[j].th   = tempguts[j].th   + 0.5*prefstep*dth[0][j];
      blobguts[j].s2   = tempguts[j].s2   + 0.5*prefstep*ds2[0][j];
      blobguts[j].a2   = tempguts[j].a2   + 0.5*prefstep*da2[0][j];
      set_blob(blobguts+j,tmpparms+j);
    }

  printf("Stage 2b\n");

  write_vorts(9991);

  vel_field();

  for (j=0; j<N; ++j)
    if (fabs(blobguts[j].a2-1/blobguts[j].a2)/prefstep < axisymmtol)
      th_slave(blobguts+j,tmpparms+j);

  printf("Stage 2c\n");

  /* Half step of BCE predictor. */

  for (j=0; j<N; ++j)
    {
      dxdt[1][j] = mblob[j].blob0.dx;
      dydt[1][j] = mblob[j].blob0.dy;
      dy(blobguts+j,tmpparms+j,ds2[1]+j,da2[1]+j,dth[1]+j,prefstep);
    }

  printf("Stage 2d\n");

  for (j=0; j<N; ++j)
    {
      mblob[j].blob0.x = mblob[j].blob1.x + 0.5*prefstep*dxdt[1][j];
      mblob[j].blob0.y = mblob[j].blob1.y + 0.5*prefstep*dydt[1][j];
      blobguts[j].th   = tempguts[j].th   + 0.5*prefstep*dth[1][j];
      blobguts[j].s2   = tempguts[j].s2   + 0.5*prefstep*ds2[1][j];
      blobguts[j].a2   = tempguts[j].a2   + 0.5*prefstep*da2[1][j];
      set_blob(blobguts+j,tmpparms+j);
    }

  printf("Stage 3\n");

  write_vorts(9992);

  vel_field();

  for (j=0; j<N; ++j)
    if (fabs(blobguts[j].a2-1/blobguts[j].a2)/prefstep < axisymmtol)
      th_slave(blobguts+j,tmpparms+j);

  /* Full step of midpoint predictor. */

  for (j=0; j<N; ++j)
    {
      dxdt[2][j] = mblob[j].blob0.dx;
      dydt[2][j] = mblob[j].blob0.dy;
      dy(blobguts+j,tmpparms+j,ds2[2]+j,da2[2]+j,dth[2]+j,prefstep);
    }

  for (j=0; j<N; ++j)
    {
      mblob[j].blob0.x = mblob[j].blob1.x + prefstep*dxdt[2][j];
      mblob[j].blob0.y = mblob[j].blob1.y + prefstep*dydt[2][j];
      blobguts[j].th   = tempguts[j].th   + prefstep*dth[2][j];
      blobguts[j].s2   = tempguts[j].s2   + prefstep*ds2[2][j];
      blobguts[j].a2   = tempguts[j].a2   + prefstep*da2[2][j];
      set_blob(blobguts+j,tmpparms+j);
    }

  write_vorts(9993);

  vel_field();

  for (j=0; j<N; ++j)
    if (fabs(blobguts[j].a2-1/blobguts[j].a2)/prefstep < axisymmtol)
      th_slave(blobguts+j,tmpparms+j);

  /* Simpson's rule corrector. */

  for (j=0; j<N; ++j)
    {
      dxdt[3][j] = mblob[j].blob0.dx;
      dydt[3][j] = mblob[j].blob0.dy;
      dy(blobguts+j,tmpparms+j,ds2[3]+j,da2[3]+j,dth[3]+j,prefstep);
    }

  for (j=0; j<N; ++j)
    {
      mblob[j].blob0.x = mblob[j].blob1.x + 
	prefstep*(dxdt[0][j]+2.0*dxdt[1][j]+2.0*dxdt[2][j]+dxdt[3][j])/6.0;
      mblob[j].blob0.y = mblob[j].blob1.y + 
	prefstep*(dydt[0][j]+2.0*dydt[1][j]+2.0*dydt[2][j]+dydt[3][j])/6.0;
      blobguts[j].th   = tempguts[j].th   + 
	prefstep*(dth[0][j]+2.0*dth[1][j]+2.0*dth[2][j]+dth[3][j])/6.0;
      blobguts[j].s2   = tempguts[j].s2   + 
	prefstep*(ds2[0][j]+2.0*ds2[1][j]+2.0*ds2[2][j]+ds2[3][j])/6.0;
      blobguts[j].a2   = tempguts[j].a2   + 
	prefstep*(da2[0][j]+2.0*da2[1][j]+2.0*da2[2][j]+da2[3][j])/6.0;
      set_blob(blobguts+j,tmpparms+j);
    }

  printf("Stage 4\n");

  write_vorts(9994);
}
