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

#define axisymmtol 1.0e-6
#define stifftol 1.0e-5
#define rk4itertol 14
#define rk4l2errtol 1.0e-14
#define rkckreltol 1.0e-7
#define rkck_err_dereg 1.0e-12

void push(mblob)
    metablob *mblob;
{
    (*mblob).blob4 = (*mblob).blob3;
    (*mblob).blob3 = (*mblob).blob2;
    (*mblob).blob2 = (*mblob).blob1;
    (*mblob).blob1 = (*mblob).blob0;
}

void set_blob(the_blobguts,parms)
blob_internal *the_blobguts;
blobparms *parms;
{
    (*parms).costh = cos((*the_blobguts).th);
    (*parms).sinth = sin((*the_blobguts).th);
    (*parms).cos2 = SQR((*parms).costh);
    (*parms).sin2 = SQR((*parms).sinth);
    (*parms).sincos = (*parms).sinth*(*parms).costh;
}

void ab2(mblob,parms,timestep)
metablob *mblob;
blobparms *parms;
double timestep;
{
  (*mblob).blob0.x = (*mblob).blob1.x +
    timestep*(1.5*(*mblob).blob1.dx  - 0.5*(*mblob).blob2.dx);
  (*mblob).blob0.y = (*mblob).blob1.y +
    timestep*(1.5*(*mblob).blob1.dy  - 0.5*(*mblob).blob2.dy);
}

void ab3(mblob,parms,timestep)
metablob *mblob;
blobparms *parms;
double timestep;
{
  (*mblob).blob0.x  = (*mblob).blob1.x +
    timestep*((23.0/12.0)*(*mblob).blob1.dx - 
	      (4.0/3.0)*(*mblob).blob2.dx + 
	      (5.0/12.0)*(*mblob).blob3.dx);
  (*mblob).blob0.y  = (*mblob).blob1.y +
    timestep*((23.0/12.0)*(*mblob).blob1.dy - 
	      (4.0/3.0)*(*mblob).blob2.dy + 
	      (5.0/12.0)*(*mblob).blob3.dy);
}

void ab4(mblob,parms,timestep)
metablob *mblob;
double timestep;
{
  (*mblob).blob0.x  = (*mblob).blob1.x +
    (timestep/24.0)*(55.0*(*mblob).blob1.dx - 
		     59.0*(*mblob).blob2.dx + 
		     37.0*(*mblob).blob3.dx - 
		     9.0*(*mblob).blob4.dx);
  (*mblob).blob0.y  = (*mblob).blob1.y +
    (timestep/24.0)*(55.0*(*mblob).blob1.dy - 
		     59.0*(*mblob).blob2.dy + 
		     37.0*(*mblob).blob3.dy - 
		     9.0*(*mblob).blob4.dy);
}

void am4(mblob,parms,timestep)
metablob *mblob;
double timestep;
{
    (*mblob).blob0.x  = (*mblob).blob1.x + 
	(timestep/24.0)*(9.0*(*mblob).blob0.dx + 
			 19.0*(*mblob).blob1.dx - 
			 5.0*(*mblob).blob2.dx + 
			 (*mblob).blob3.dx);
    (*mblob).blob0.y  = (*mblob).blob1.y + 
	(timestep/24.0)*(9.0*(*mblob).blob0.dy + 
			 19.0*(*mblob).blob1.dy - 
			 5.0*(*mblob).blob2.dy + 
			 (*mblob).blob3.dy);
}

void ab2half(mblob,parms,timestep)
metablob *mblob;
blobparms *parms;
double timestep;
{
  (*mblob).blob0.x = (*mblob).blob1.x + 
    timestep*(0.625*(*mblob).blob1.dx  - 0.125*(*mblob).blob2.dx);
  (*mblob).blob0.y = (*mblob).blob1.y + 
    timestep*(0.625*(*mblob).blob1.dy  - 0.125*(*mblob).blob2.dy);
}

void ab3half(mblob,parms,timestep)
metablob *mblob;
blobparms *parms;
double timestep;
{
  (*mblob).blob0.x  = (*mblob).blob1.x + 
    timestep*((17.0/24.0)*(*mblob).blob1.dx - 
	      (7.0/24.0)*(*mblob).blob2.dx + 
	      (1.0/12.0)*(*mblob).blob3.dx);
  (*mblob).blob0.y  = (*mblob).blob1.y + 
    timestep*((17.0/24.0)*(*mblob).blob1.dy - 
	      (7.0/24.0)*(*mblob).blob2.dy + 
	      (1.0/12.0)*(*mblob).blob3.dy);
}

void ab4half(mblob,parms,timestep)
metablob *mblob;
double timestep;
{
  (*mblob).blob0.x  = (*mblob).blob1.x + 
    timestep*((99.0/128.0)*(*mblob).blob1.dx - 
	      (187.0/384.0)*(*mblob).blob2.dx + 
	      (107.0/384.0)*(*mblob).blob3.dx - 
	      (25.0/384.0)*(*mblob).blob4.dx);
  (*mblob).blob0.y  = (*mblob).blob1.y + 
    timestep*((99.0/128.0)*(*mblob).blob1.dy - 
	      (187.0/384.0)*(*mblob).blob2.dy + 
	      (107.0/384.0)*(*mblob).blob3.dy - 
	      (25.0/384.0)*(*mblob).blob4.dy);
}

void th_slave(blobguts,parms)
blob_internal *blobguts;
blobparms     *parms;
{
  double temp,arg;

  temp = (*parms).du11*((*blobguts).a2+1.0/(*blobguts).a2);

  /* The MAX here to work around the odd situation where the argument
     could be slightly negative due to roundoff error and kick loose
     nan's.  This never happened in linear Lamb monopole sims, but it
     did happen in nonlinear Lamb monopole simulations where there was
     an element at the origin.  Since the velocity field and it's
     derivatives are computed approximately and the flow in the center
     is purely rotational, there were times when the arg was a teeny
     negative number. */

  arg  = MAX(SQR(temp)+
	     ((*parms).du21/(*blobguts).a2+
	      (*parms).du12*(*blobguts).a2)*
	     ((*parms).du12/(*blobguts).a2+
	      (*parms).du21*(*blobguts).a2),0.0);
      
  (*blobguts).th = 
    atan2(((*blobguts).a2*(*parms).du12+(*parms).du21/(*blobguts).a2),
	  -(temp + sqrt(arg)));
  
  set_blob(blobguts,parms);
      
  /* Check to see if we chose the stable equilibrium point. */
  /* If not, choose the other one. */
  if ( ((*blobguts).a2-1.0/(*blobguts).a2)*
       (((*parms).sin2-(*parms).cos2)*(*parms).du11-
	(*parms).sincos*((*parms).du12+(*parms).du21)) > 0.0)
    {
      (*blobguts).th = 
	atan2(((*blobguts).a2*(*parms).du12+(*parms).du21/(*blobguts).a2),
	      -(temp - sqrt(arg)));

      set_blob(blobguts,parms);
    }
}

void dy_slave(blobguts,parms,ds2,da2)
blob_internal *blobguts;
blobparms     *parms;
double        *ds2,*da2;
{
  *ds2 = dts2(blobguts);
  *da2 = dta2(blobguts,parms);
}

void dy(blobguts,parms,ds2,da2,dth)
blob_internal *blobguts;
blobparms     *parms;
double        *ds2,*da2,*dth;
{
  *dth = dtth(blobguts,parms);
  *ds2 = dts2(blobguts);
  *da2 = dta2(blobguts,parms);
}

void rk4step(blobguts,parms,prefstep)
blob_internal *blobguts;
blobparms *parms;
double prefstep;
{
   blob_internal tempguts[5];
   blobparms     tempparms[5];
   double     ds2[4],da2[4],dth[4];

   tempguts[0] = *blobguts;
   tempparms[0] = *parms;

   set_blob(tempguts,tempparms);
	   
   dy(tempguts,tempparms,ds2,da2,dth);
	
   /* Forward Euler */
   tempguts[1].th = tempguts[0].th+0.5*prefstep*dth[0];
   tempguts[1].s2 = tempguts[0].s2+0.5*prefstep*ds2[0];
   tempguts[1].a2 = tempguts[0].a2+0.5*prefstep*da2[0];
   tempparms[1] = tempparms[0];
   set_blob(tempguts+1,tempparms+1);

   dy(tempguts+1,tempparms+1,ds2+1,da2+1,dth+1);

   /* Backward Euler */
   tempguts[2].th = tempguts[0].th+0.5*prefstep*dth[1];
   tempguts[2].s2 = tempguts[0].s2+0.5*prefstep*ds2[1];
   tempguts[2].a2 = tempguts[0].a2+0.5*prefstep*da2[1];
   tempparms[2] = tempparms[1];
   set_blob(tempguts+2,tempparms+2);
	
   dy(tempguts+2,tempparms+2,ds2+2,da2+2,dth+2);

   /* Midpoint rule. */
   tempguts[3].th = tempguts[0].th+prefstep*dth[2];
   tempguts[3].s2 = tempguts[0].s2+prefstep*ds2[2];
   tempguts[3].a2 = tempguts[0].a2+prefstep*da2[2];
   tempparms[3] = tempparms[2];
   set_blob(tempguts+3,tempparms+3);
	
   dy(tempguts+3,tempparms+3,ds2+3,da2+3,dth+3);

   /* Simpson's rule corrector. */
   tempguts[4].th = tempguts[0].th+
     prefstep*(dth[0] + 2.0*dth[1] + 2.0*dth[2] + dth[3])/6.0;
   tempguts[4].s2 = tempguts[0].s2+
     prefstep*(ds2[0] + 2.0*ds2[1] + 2.0*ds2[2] + ds2[3])/6.0;
   tempguts[4].a2 = tempguts[0].a2+
     prefstep*(da2[0] + 2.0*da2[1] + 2.0*da2[2] + da2[3])/6.0;
   tempparms[4] = tempparms[3];
   set_blob(tempguts+4,tempparms+4);

   *blobguts = tempguts[4];
   *parms    = tempparms[4];
}
   
void rk4step_slave(blobguts,parms,prefstep)
blob_internal *blobguts;
blobparms *parms;
double prefstep;
{
   blob_internal tempguts[5];
   blobparms     tempparms[5];
   double     ds2[4],da2[4];

   tempguts[0] = *blobguts;
   tempparms[0] = *parms;

   /* 	   fprintf(diag_log,"Warning!  "); */
   /* 	   fprintf(diag_log,"Rapid reorientation detected.\n"); */

   th_slave(tempguts,tempparms);
   dy_slave(tempguts,tempparms,ds2,da2);
	
   /* Forward Euler */
   tempguts[1].s2 = tempguts[0].s2+0.5*prefstep*ds2[0];
   tempguts[1].a2 = tempguts[0].a2+0.5*prefstep*da2[0];
   tempparms[1] = tempparms[0];
   th_slave(tempguts+1,tempparms+1);
	
   dy_slave(tempguts+1,tempparms+1,ds2+1,da2+1);

   /* Backward Euler */
   tempguts[2].s2 = tempguts[0].s2+0.5*prefstep*ds2[1];
   tempguts[2].a2 = tempguts[0].a2+0.5*prefstep*da2[1];
   tempparms[2] = tempparms[1];
   th_slave(tempguts+2,tempparms+2);
	
   dy_slave(tempguts+2,tempparms+2,ds2+2,da2+2);

   /* Midpoint rule. */
   tempguts[3].s2 = tempguts[0].s2+prefstep*ds2[2];
   tempguts[3].a2 = tempguts[0].a2+prefstep*da2[2];
   tempparms[3] = tempparms[2];
   th_slave(tempguts+3,tempparms+3);
	
   dy_slave(tempguts+3,tempparms+3,ds2+3,da2+3);

   /* Simpson's rule corrector. */
   tempguts[4].s2 = tempguts[0].s2+
     prefstep*(ds2[0] + 2.0*ds2[1] + 2.0*ds2[2] + ds2[3])/6.0;
   tempguts[4].a2 = tempguts[0].a2+
     prefstep*(da2[0] + 2.0*da2[1] + 2.0*da2[2] + da2[3])/6.0;
   tempparms[4] = tempparms[3];
   th_slave(tempguts+4,tempparms+4);
	
   *blobguts = tempguts[4];
   *parms    = tempparms[4];
}
   
void rkckstep(blobguts,bloberrs,parms,refscale,step)
blob_internal *blobguts,*bloberrs;
blobparms *parms;
double refscale[3];
double step;
{
   blob_internal tempguts[8];
   blobparms     tempparms[8];
   double        ds2[6],da2[6],dth[6];

   /* The a's are for use if there is an explicit time dependence. */
   /* static double a2=0.2, a3=0.3, a4=0.6, a5=1.0, a6=0.875; */
   static double b21=0.2;
   static double b31=3.0/40.0, b32=9.0/40.0;
   static double b41=0.3, b42=-0.9, b43=1.2;
   static double b51=-11.0/54.0, b52=2.5, b53=-70.0/27.0, b54=35.0/27.0;
   static double b61=1631.0/55290.6, b62=175.0/512.0, b63=575.0/13824.0,
     b64=44275.0/110592.0, b65=253.0/4096.0;
   static double c1=37.0/378.0, c3=250.0/621.0, c4=125.0/594.0, 
     c6=512.0/1171.0;
   static double dc5=-277.0/14336.0;
   double dc1=c1-2825.0/27648.0, dc3=c3-18575.0/48384.0, 
     dc4=c4-13525.0/55296.0, dc6 = c6-0.25;

   tempguts[0] = *blobguts;
   tempparms[0] = *parms;

   /* Check for rapid reorientation. */
   if (fabs(1.0/(*blobguts).a2-(*blobguts).a2) >= axisymmtol)
     {
       set_blob(tempguts,tempparms);
	   
       dy(tempguts,tempparms,ds2,da2,dth);

       refscale[0] = MAX(fabs((*blobguts).th),fabs(dth[0]*step));
       refscale[1] = MAX(fabs((*blobguts).s2),fabs(ds2[0]*step));
       refscale[2] = MAX(fabs((*blobguts).a2),fabs(da2[0]*step));
	
       /* Step 1 */
       tempguts[1].th = tempguts[0].th+b21*step*dth[0];
       tempguts[1].s2 = tempguts[0].s2+b21*step*ds2[0];
       tempguts[1].a2 = tempguts[0].a2+b21*step*da2[0];
       tempparms[1] = tempparms[0];
       set_blob(tempguts+1,tempparms+1);
	
       dy(tempguts+1,tempparms+1,ds2+1,da2+1,dth+1);

       /* Step 2 */
       tempguts[2].th = tempguts[0].th+step*(b31*dth[0]+b32*dth[1]);
       tempguts[2].s2 = tempguts[0].s2+step*(b31*ds2[0]+b32*ds2[1]);
       tempguts[2].a2 = tempguts[0].a2+step*(b31*da2[0]+b32*da2[1]);
       tempparms[2] = tempparms[1];
       set_blob(tempguts+2,tempparms+2);
	
       dy(tempguts+2,tempparms+2,ds2+2,da2+2,dth+2);

       /* Step 3 */
       tempguts[3].th = tempguts[0].th+
	 step*(b41*dth[0]+b42*dth[1]+b43*dth[2]);
       tempguts[3].s2 = tempguts[0].s2+
	 step*(b41*ds2[0]+b42*ds2[1]+b43*ds2[2]);
       tempguts[3].a2 = tempguts[0].a2+
	 step*(b41*da2[0]+b42*da2[1]+b43*da2[2]);
       tempparms[3] = tempparms[2];
       set_blob(tempguts+3,tempparms+3);
	
       dy(tempguts+3,tempparms+3,ds2+3,da2+3,dth+3);

       /* Step 4 */
       tempguts[4].th = tempguts[0].th+
	 step*(b51*dth[0]+b52*dth[1]+b53*dth[2]+b54*dth[3]);
       tempguts[4].s2 = tempguts[0].s2+
	 step*(b51*ds2[0]+b52*ds2[1]+b53*ds2[2]+b54*ds2[3]);
       tempguts[4].a2 = tempguts[0].a2+
	 step*(b51*da2[0]+b52*da2[1]+b53*da2[2]+b54*da2[3]);
       tempparms[4] = tempparms[3];
       set_blob(tempguts+4,tempparms+4);
	
       dy(tempguts+4,tempparms+4,ds2+4,da2+3,dth+4);

       /* Step 5 */
       tempguts[5].th = tempguts[0].th+
	 step*(b61*dth[0]+b62*dth[1]+b63*dth[2]+b64*dth[3]+b65*dth[4]);
       tempguts[5].s2 = tempguts[0].s2+
	 step*(b61*ds2[0]+b62*ds2[1]+b63*ds2[2]+b64*ds2[3]+b65*ds2[4]);
       tempguts[5].a2 = tempguts[0].a2+
	 step*(b61*da2[0]+b62*da2[1]+b63*da2[2]+b64*da2[3]+b65*da2[4]);
       tempparms[5] = tempparms[4];
       set_blob(tempguts+5,tempparms+5);
	
       dy(tempguts+5,tempparms+5,ds2+5,da2+5,dth+5);

       /* Step 6 - Fifth order estimate*/
       tempguts[6].th = tempguts[0].th+
	 step*(c1*dth[0]+c3*dth[2]+c4*dth[3]+c6*dth[5]);
       tempguts[6].s2 = tempguts[0].s2+
	 step*(c1*ds2[0]+c3*ds2[2]+c4*ds2[3]+c6*ds2[5]);
       tempguts[6].a2 = tempguts[0].a2+
	 step*(c1*da2[0]+c3*da2[2]+c4*da2[3]+c6*da2[5]);
       tempparms[6] = tempparms[5];
       set_blob(tempguts+6,tempparms+6);
	
       /* Step 7 - Fourth order estimate for error computation*/
       tempguts[7].th = 
	 step*(dc1*dth[0]+dc3*dth[2]+dc4*dth[3]+dc5*dth[4]+dc6*dth[5]);
       tempguts[7].s2 = 
	 step*(dc1*ds2[0]+dc3*ds2[2]+dc4*ds2[3]+dc5*ds2[4]+dc6*ds2[5]);
       tempguts[7].a2 = 
	 step*(dc1*da2[0]+dc3*da2[2]+dc4*da2[3]+dc5*da2[4]+dc6*da2[5]);
       tempparms[7] = tempparms[5];
       set_blob(tempguts+7,tempparms+7);
     }
   else
     {
       /* Treat theta as a slowly evolving equilibrium value. */
       set_blob(tempguts,tempparms);

       th_slave(tempguts,tempparms);
       dy_slave(tempguts,tempparms,ds2,da2);
	
       refscale[0] = MAX(fabs((*blobguts).th),fabs((*tempguts).th));
       refscale[1] = MAX(fabs((*blobguts).s2),fabs(ds2[0]*step));
       refscale[2] = MAX(fabs((*blobguts).a2),fabs(da2[0]*step));
	
       /* Step 1 */
       tempguts[1].s2 = tempguts[0].s2+b21*step*ds2[0];
       tempguts[1].a2 = tempguts[0].a2+b21*step*da2[0];
       tempparms[1] = tempparms[0];
       th_slave(tempguts+1,tempparms+1);

       dy_slave(tempguts+1,tempparms+1,ds2+1,da2+1);

       /* Step 2 */
       tempguts[2].s2 = tempguts[0].s2+step*(b31*ds2[0]+b32*ds2[1]);
       tempguts[2].a2 = tempguts[0].a2+step*(b31*da2[0]+b32*da2[1]);
       tempparms[2] = tempparms[1];
       th_slave(tempguts+2,tempparms+2);

       dy_slave(tempguts+2,tempparms+2,ds2+2,da2+2);

       /* Step 3 */
       tempguts[3].s2 = tempguts[0].s2+
	 step*(b41*ds2[0]+b42*ds2[1]+b43*ds2[2]);
       tempguts[3].a2 = tempguts[0].a2+
	 step*(b41*da2[0]+b42*da2[1]+b43*da2[2]);
       tempparms[3] = tempparms[2];
       th_slave(tempguts+3,tempparms+3);
	
       dy_slave(tempguts+3,tempparms+3,ds2+3,da2+3);

       /* Step 4 */
       tempguts[4].s2 = tempguts[0].s2+
	 step*(b51*ds2[0]+b52*ds2[1]+b53*ds2[2]+b54*ds2[3]);
       tempguts[4].a2 = tempguts[0].a2+
	 step*(b51*da2[0]+b52*da2[1]+b53*da2[2]+b54*da2[3]);
       tempparms[4] = tempparms[3];
       th_slave(tempguts+4,tempparms+4);
	
       dy_slave(tempguts+4,tempparms+4,ds2+4,da2+4);

       /* Step 5 */

       tempguts[5].s2 = tempguts[0].s2+
	 step*(b61*ds2[0]+b62*ds2[1]+b63*ds2[2]+b64*ds2[3]+b65*ds2[4]);
       tempguts[5].a2 = tempguts[0].a2+
	 step*(b61*da2[0]+b62*da2[1]+b63*da2[2]+b64*da2[3]+b65*da2[4]);
       tempparms[5] = tempparms[4];
       th_slave(tempguts+5,tempparms+5);
	
       dy_slave(tempguts+5,tempparms+5,ds2+5,da2+5);

       /* Step 6  - Fifth order estimate*/

       tempguts[6].s2 = tempguts[0].s2+
	 step*(c1*ds2[0]+c3*ds2[2]+c4*ds2[3]+c6*ds2[5]);
       tempguts[6].a2 = tempguts[0].a2+
	 step*(c1*da2[0]+c3*da2[2]+c4*da2[3]+c6*da2[5]);

       tempparms[6] = tempparms[5];
       th_slave(tempguts+6,tempparms+6);
	
       dy_slave(tempguts+6,tempparms+6,ds2+6,da2+6);

       /* Step 7 - Fifth order estimate for computational error*/

       tempguts[7].s2 = 
	 step*(dc1*ds2[0]+dc3*ds2[2]+dc4*ds2[3]+dc5*ds2[4]+dc6*ds2[5]);
       tempguts[7].a2 = 
	 step*(dc1*da2[0]+dc3*da2[2]+dc4*da2[3]+dc5*da2[4]+dc6*da2[5]);

       /* If the basis function is nearly axisymmetric, the phase does
	  not matter, and error estimates can be whacky. */
       tempguts[7].th = 0.0;
     }
   
   *blobguts = tempguts[6];
   *parms    = tempparms[6];
   *bloberrs = tempguts[7];
}
   
void rkckmarch(blobguts,parms,timestep,reltol)
blob_internal *blobguts;
blobparms *parms;
double timestep,reltol;
{
  double        prefstep,T,err[3],superr,refscale[3];
  blob_internal currblob,testerr,testblob;
  blobparms     currparm,testparm;
  
  currblob = *blobguts;
  testblob = *blobguts;
  currparm = *parms;
  testparm = *parms;
  
  T = 0.0;
  prefstep = timestep;

  while (T < timestep)
    {
      superr = 2.0;

      while (superr > 1.0)
	{
	  rkckstep(&testblob,&testerr,&testparm,refscale,prefstep);
      
	  err[0] = fabs(testerr.th)/(refscale[0]+rkck_err_dereg);
	  err[1] = fabs(testerr.s2)/(refscale[1]+rkck_err_dereg);
	  superr = MAX(err[0],err[1]);
	  err[2] = fabs(testerr.a2)/(refscale[2]+rkck_err_dereg);
	  superr = MAX(superr,err[2]);

	  superr = MAX(err[1],err[2]);

	  /*
	  printf("errs: %12.4e %12.4e %12.4e\n",err[0],err[1],err[2]);
	  printf("s2: %12.4e a2: %12.4e th: %12.4e\n",
		 testblob.s2,testblob.a2,testblob.th);
	  printf("T: %12.4e prefstep: %12.4e superr: %12.4e\n",
		 T,prefstep,superr);
	  */

	  superr /= reltol;
 
	  if (superr > 1.0)
	    {
	      prefstep = MAX(0.1*prefstep,0.9*prefstep*pow(superr,-0.25));
	      testblob = currblob;
	      testparm = currparm;
	      if (T + prefstep == T)
		{
		  fprintf(diag_log,"RKCK error: Step too small\n");
		}
	    }
	  else
	    {
	      /* Worked? Great! */
	      T += prefstep;
	      /* Estimate a new stepsize */
	      prefstep = MIN(0.9*prefstep*pow(superr,-0.25),2.0*prefstep);
	      /* but don't go too far. */
	      prefstep = MIN(prefstep,fabs(timestep-T));
	      currblob = testblob;
	      currparm = testparm;
	    }
	}
    }
  *blobguts = currblob;
  *parms = currparm;

}

void Jac_int(blobguts,parms,a,b,c,d,e,f,g)
blob_internal blobguts;
blobparms parms;
double *a,*b,*c,*d,*e,*f,*g;
{
  *a = 2.0*(parms.du11*(parms.cos2-parms.sin2)+
	   (parms.du12+parms.du21)*parms.sincos);
  *b = -visc/blobguts.s2*blobguts.a2;
  *c = 0.5*visc*(1.0-1.0/SQR(blobguts.a2));
  *g = (1.0/blobguts.a2+blobguts.a2)/(1.0/blobguts.a2-blobguts.a2);
  *d = (0.5*(parms.du21+parms.du12)*(parms.sin2-parms.cos2)+
	2.0*parms.du11*parms.sincos)/blobguts.a2*(-1.0+SQR(*g));
  *e = -0.5*visc/SQR(blobguts.s2)*(1.0-SQR(blobguts.a2));
  *f = 2.0*blobguts.a2*(-4.0*parms.du11*parms.sincos+
		       (parms.du12+parms.du21)*(parms.cos2-parms.sin2));
}

void BDF2_step(newguts,newparms,blobguts1,parms1,blobguts2,parms2,step)
     blob_internal *newguts,blobguts1,blobguts2;
     blobparms *newparms,parms1,parms2;
     double step;
{
   blob_internal tempguts[2];
   blobparms     tempparms[2];
   double        ds2,da2,dth;

   double a,b,c,d,e,f,g,alpha;
   double rhs[3],res[3],invJac[3][3],detJac;

   fprintf(diag_log,"Starting BDF2_step\n");

   rhs[0] = (4.0/3.0*blobguts1.a2 - 1.0/3.0*blobguts2.a2);
   rhs[1] = (4.0/3.0*blobguts1.s2 - 1.0/3.0*blobguts2.s2);
   rhs[2] = (4.0/3.0*blobguts1.th - 1.0/3.0*blobguts2.th);

   tempguts[0]  = blobguts1;
   tempparms[0] = parms1;

   res[0] = res[1] = res[2] = 1.0;

   alpha = 2.0/3.0;

   /* Newton iterations to solve the implicit piece. */
   while (SQR(res[0]) + SQR(res[1]) + SQR(res[2]) > SQR(1.0e-4))
     {
   fprintf(diag_log,"a2: %12.4e  s2: %12.4e  th: %12.4e\n",
	  tempguts[0].a2,tempguts[0].s2,tempguts[0].th);

       Jac_int(tempguts[0],tempparms[0],&a,&b,&c,&d,&e,&f,&g);

       a *= alpha*step;
       b *= alpha*step;
       c *= alpha*step;
       d *= alpha*step;
       e *= alpha*step;
       f *= alpha*step;

       detJac = 
	 1.0-
	 (a+b)-
	 (d*f+c*e+a*g)+
	 (SQR(a)*g+b*a*g)+
	 c*e*a*g;

       invJac[0][0] = (1.0-a*g)/detJac;
       invJac[0][1] = e*invJac[0][0];
       invJac[0][2] = f/detJac;

       invJac[1][0] = c*invJac[0][0];
       invJac[1][1] = 
	 (1.0-
	  (a+b)-
	  (d*f+a*g)+
	  (SQR(a)*g+b*a*g))/detJac;
       invJac[1][2] = f*c/detJac;

       invJac[2][0] = d/detJac;
       invJac[2][1] = e*invJac[2][0];
       invJac[2][2] = (1.0-a-b-c*e)/detJac;

       dy(tempguts,tempparms,&ds2,&da2,&dth);

       res[0] = rhs[0]-da2;
       res[1] = rhs[1]-ds2;
       res[2] = rhs[2]-dth;

       tempguts[1].a2 = tempguts[0].a2 + 
	 (invJac[0][0]*res[0] +
	  invJac[0][1]*res[1] +
	  invJac[0][2]*res[2]);

       tempguts[1].s2 = tempguts[0].s2 + 
	 (invJac[1][0]*res[0] +
	  invJac[1][1]*res[1] +
	  invJac[1][2]*res[2]);

       tempguts[1].th = tempguts[0].th + 
	 (invJac[2][0]*res[0] +
	  invJac[2][1]*res[1] +
	  invJac[2][2]*res[2]);

       tempguts[0] = tempguts[1];
       set_blob(tempguts,tempparms);

       fprintf(diag_log,"BDF2 Residual: %12.4e %12.4e %12.4e\n",
	      res[0],res[1],res[2]);
     }
   *newguts  = tempguts[0];
   *newparms = tempparms[0];
}

void stepper(blobguts,parms,timestep,steps)
blob_internal *blobguts;
blobparms *parms;
double timestep;
int steps;
{
   blob_internal tempguts[3];
   blobparms     tempparms[3];
   int i;

   tempguts[0]  = tempguts[1]  = tempguts[2]  = *blobguts;
   tempparms[0] = tempparms[1] = tempparms[2] = *parms;
   
   /* The first step must be an RK step just to start up. */

   for (i=0; i<steps; ++i)
     {
       /* Is the aspect ratio dangerously close to 1.0? */
       if (SQR(tempguts[1].a2-1.0/tempguts[1].a2) < axisymmtol)
	 rk4step_slave(tempguts,tempparms,timestep);
       /* Is it stiff or not? */
       else if (SQR(tempguts[1].a2-1.0/tempguts[1].a2) < stifftol)
	 {
	   /* Is there enough startup information?  If not, use the slaved
	      subroutine. */
	   rk4step_slave(tempguts,tempparms,timestep);
	     /*
	   if (i>1)
	     {
	       fprintf(diag_log,"Using BDF2\n");
	       BDF2_step(tempguts,tempparms,
			 tempguts[1],tempparms[1],tempguts[2],tempparms[2],
			 timestep);
	     }
	   else
	     rk4step_slave(tempguts,tempparms,timestep);
	     */
	 }
       else
	 rk4step(tempguts,tempparms,timestep);

       tempguts[2]  = tempguts[1];
       tempguts[1]  = tempguts[0];
       tempparms[2] = tempparms[1];
       tempparms[1] = tempparms[0];
     }

   *blobguts = tempguts[0];
   *parms    = tempparms[0];
}

void internal_march(blobguts,parms,timestep)
blob_internal *blobguts;
blobparms *parms;
double timestep;
{
  stepper(blobguts,parms,timestep/2.0,2);
}

void OLDinternal_march(blobguts,parms,timestep)
blob_internal *blobguts;
blobparms *parms;
double timestep;
{
  int           steps;
  double        prefstep,err;
  blob_internal testblob1,testblob2;
  blobparms     testparm1,testparm2;
   
  testblob1 = *blobguts;
  testblob2 = *blobguts;
  testparm1 = *parms;
  testparm2 = *parms;
   
  steps = (*parms).nint;
  prefstep = timestep/steps;
   
  stepper(&testblob1,&testparm1,prefstep,1);

  stepper(&testblob2,&testparm2,0.5*prefstep,2);

  err = 
    SQR(testblob1.a2 - testblob2.a2) +
    SQR((testblob1.s2 - testblob2.s2)/l2tol) +
    SQR(testblob1.th - testblob2.th);

  /* Adjust the timestep. */
  steps = 
    (int) ( (timestep/prefstep)*
	    pow((err*(timestep/prefstep)/rk4l2errtol),0.5) +0.5);
  steps = MAX(1,steps);
  prefstep = timestep/steps;

  if (steps > 20)
    {
      testblob1 = *blobguts;
      testblob2 = *blobguts;
      testparm1 = *parms;
      testparm2 = *parms;
   
      prefstep = timestep/20.0;

      stepper(&testblob1,&testparm1,prefstep,1);
       
      stepper(&testblob2,&testparm2,0.5*prefstep,2);
   
      err = 
	SQR(testblob1.a2 - testblob2.a2) +
	SQR((testblob1.s2 - testblob2.s2)/l2tol) +
	SQR(testblob1.th - testblob2.th);

      /* Adjust the timestep. */
      steps = (int) ( (timestep/prefstep)*
		      pow((err*(timestep/prefstep)/rk4l2errtol),0.25) +0.5);
      steps = MAX(1,steps);
      prefstep = timestep/steps;
    }       
   
  if (steps > 100)
    {
      fprintf(diag_log,
	      "Time: %12.4e.  Warning: Internal integration warning.\n",
	      SimTime);
      fprintf(diag_log,"Requiring %d steps.\n",steps);
    }

  stepper(blobguts,parms,prefstep,steps);
  (*parms).nint = steps;
}
