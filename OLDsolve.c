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
 * Newark, DE 19715-2553                                                */

#include "global.h"

#define axisymmtol 1.0e-6
#define rk4itertol 14
#define rk4l2errtol 1.0e-14

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
  double dumb,arg;

  dumb = (*parms).du11*((*blobguts).a2+1.0/(*blobguts).a2);
      
  arg = -(dumb +
	  sqrt(SQR(dumb)+
	       ((*parms).du21/(*blobguts).a2+
		(*parms).du12*(*blobguts).a2)*
	       ((*parms).du12/(*blobguts).a2+
		(*parms).du21*(*blobguts).a2)))/
    ((*blobguts).a2*(*parms).du12+(*parms).du21/(*blobguts).a2);


  if (finite(arg))
    {
      (*blobguts).th = atan(arg);
      set_blob(blobguts,parms);
      
      /* Check to see if we chose the stable equilibrium point. */
      /* If not, choose the other one. */
      if ( ((*blobguts).a2-1.0/(*blobguts).a2)*
	   (((*parms).sin2-(*parms).cos2)*(*parms).du11-
	    (*parms).sincos*((*parms).du12+(*parms).du21)) > 0.0)
	(*blobguts).th = 
	  atan(-(dumb -
		 sqrt(SQR(dumb)+
		      ((*parms).du21/(*blobguts).a2+
		       (*parms).du12*(*blobguts).a2)*
		      ((*parms).du12/(*blobguts).a2+
		       (*parms).du21*(*blobguts).a2)))/
	       ((*blobguts).a2*(*parms).du12+
		(*parms).du21/(*blobguts).a2));
      
    }
  else
    /* We should probably be more careful here.  It is possible that
       theta=0 is an unstable equilibrium.  The situation where theta=0
       is an equilibrium at all is a rare circumstance, so we can fix
       this later. */
    (*blobguts).th = 0.0;
  
  set_blob(blobguts,parms);
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

void rk4march(blobguts,parms,prefstep,steps)
blob_internal *blobguts;
blobparms *parms;
double prefstep;
int       steps;
{
   blob_internal tempguts[5];
   blobparms     tempparms[5];
   double     ds2[4],da2[4],dth[4];
   int           i;

   tempguts[0] = *blobguts;
   tempparms[0] = *parms;

   for (i=0; i<steps; ++i)
     {
       /* Check for rapid reorientation. */
       if (fabs(1.0/(*blobguts).a2-(*blobguts).a2) >= axisymmtol)
	 {
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
	
	   tempguts[0] = tempguts[4];
	   tempparms[0] = tempparms[4];
	 }
       else
	 {
	   set_blob(tempguts,tempparms);

	   th_slave(tempguts,tempparms);
	   dy_slave(tempguts,tempparms,ds2,da2);
	
	   /* Forward Euler */
	   tempguts[1].s2 = tempguts[0].s2+0.5*prefstep*ds2[0];
	   tempguts[1].a2 = tempguts[0].a2+0.5*prefstep*da2[0];
	   th_slave(tempguts+1,tempparms+1);
	   tempparms[1] = tempparms[0];
	   set_blob(tempguts+1,tempparms+1);
	
	   dy_slave(tempguts+1,tempparms+1,ds2+1,da2+1);

	   /* Backward Euler */
	   tempguts[2].s2 = tempguts[0].s2+0.5*prefstep*ds2[1];
	   tempguts[2].a2 = tempguts[0].a2+0.5*prefstep*da2[1];
	   th_slave(tempguts+2,tempparms+2);
	   tempparms[2] = tempparms[1];
	   set_blob(tempguts+2,tempparms+2);
	
	   dy_slave(tempguts+2,tempparms+2,ds2+2,da2+2);

	   /* Midpoint rule. */
	   tempguts[3].s2 = tempguts[0].s2+prefstep*ds2[2];
	   tempguts[3].a2 = tempguts[0].a2+prefstep*da2[2];
	   th_slave(tempguts+3,tempparms+3);
	   tempparms[3] = tempparms[2];
	   set_blob(tempguts+3,tempparms+3);
	
	   dy_slave(tempguts+3,tempparms+3,ds2+3,da2+3);

	   /* Simpson's rule corrector. */
	   tempguts[4].s2 = tempguts[0].s2+
	     prefstep*(ds2[0] + 2.0*ds2[1] + 2.0*ds2[2] + ds2[3])/6.0;
	   tempguts[4].a2 = tempguts[0].a2+
	     prefstep*(da2[0] + 2.0*da2[1] + 2.0*da2[2] + da2[3])/6.0;
	   th_slave(tempguts+4,tempparms+4);
	   tempparms[4] = tempparms[3];
	   set_blob(tempguts+4,tempparms+4);
	
	   tempguts[0] = tempguts[4];
	   tempparms[0] = tempparms[4];
	 }
     }
   
   *blobguts = tempguts[4];
   *parms    = tempparms[4];
}
   
void rk4internal(blobguts,parms,timestep)
blob_internal *blobguts;
blobparms *parms;
double timestep;
{
   int           steps,crashtest;
   double        prefstep,err;
   blob_internal testblob1,testblob2;
   blobparms     testparm1,testparm2;
   
   testblob1 = *blobguts;
   testblob2 = *blobguts;
   testparm1 = *parms;
   testparm2 = *parms;
   
   crashtest = 0;
   steps = 4;
   prefstep = timestep/4.0;
   
   rk4march(&testblob1,&testparm1,prefstep,steps);
   
   steps *= 2;
   prefstep /= 2.0;
   rk4march(&testblob2,&testparm2,prefstep,steps);
   
   err = 
     SQR(testblob1.a2 - testblob2.a2) +
     SQR((testblob1.s2 - testblob2.s2)/l2tol) +
     SQR((testblob1.th - testblob2.th)*(testblob2.a2-1.0/testblob2.a2));
   
   while ( (err > rk4l2errtol) &&  (crashtest < rk4itertol) )
     {
	++crashtest;
	
	testblob1 = testblob2;
	testparm1 = testparm2;
	testblob2 = *blobguts;
	testparm2 = *parms;
	
	steps *= 2;
	prefstep /= 2.0;
	rk4march(&testblob2,&testparm2,prefstep,steps);
	
	err = 
	  SQR(testblob1.a2 - testblob2.a2) +
	  SQR((testblob1.s2 - testblob2.s2)/l2tol) +
	  SQR((testblob1.th - testblob2.th)*(testblob2.a2-1.0/testblob2.a2));
     }

   if ( (crashtest == rk4itertol) && (err > rk4l2errtol) )
     {
       fprintf(diag_log,
	       "RK4 failed to converge.  Error: %12.4e %12.4e %12.4e %12.4e\n",
	       SQR(testblob1.a2 - testblob2.a2) +
	       SQR((testblob1.s2 - testblob2.s2)/l2tol) +
	       SQR(testblob1.th - testblob2.th),
	       fabs(testblob1.a2 - testblob2.a2),
	       fabs((testblob1.s2 - testblob2.s2)/l2tol),
	       fabs(testblob1.th - testblob2.th));
       fflush(diag_log);
     }
   
   *blobguts = testblob2;
   *parms    = testparm2;
}
