/* LAMB-OSEEN.C */
/* Copyright (c) 2001 Louis F. Rossi                                    *
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

#include "global_min.h"
#include "particle.h"

#ifdef XANTISYMM
#define CARDINALITY  N/2
#else
#define CARDINALITY  N
#endif

#ifdef MULTIPROC
#include "multiproc.h"
#endif

#define correction1 1.0
#define correction2 1.0
#define correction3 0.0

double fxx(double x,double y,double r2,double s2,
	   double oneminusdecay,double decay)
{
  return((2*x*x*y/r2/r2/s2-0.5*y/r2/s2+0.25*x*x*y/r2/s2/s2)*decay+
	 (-8*x*x*y/r2/r2/r2 + 2*y/r2/r2)*oneminusdecay);
}

double fyy(double x,double y,double r2,double s2,
	   double oneminusdecay,double decay)
{
  return((-1.5*y/r2/s2+2*y*y*y/r2/r2/s2+0.25*y*y*y/r2/s2/s2)*decay +
	 (6*y/r2/r2-8*y*y*y/r2/r2/r2)*oneminusdecay);
}

double fxxx(double x,double y,double r2,double s2,
	    double oneminusdecay,double decay)
{
  return((-12*y*x*x*x/r2/r2/r2/s2+6*x*y/r2/r2/s2-1.5*x*x*x*y/r2/r2/s2/s2+
	  0.75*x*y/r2/s2/s2-0.125*x*x*x*y/r2/s2/s2/s2)*decay+
	 24*(2*x*x*x*y/r2/r2/r2/r2-x*y/r2/r2/r2)*oneminusdecay);
}

double fxxy(double x,double y,double r2,double s2,
	    double oneminusdecay,double decay)
{
  return((-12*x*x*y*y/r2/r2/r2/s2-1.5*x*x*y*y/r2/r2/s2/s2+
	  2/r2/s2-0.5/r2/s2+0.25/s2/s2-0.125*x*x*y*y/r2/s2/s2/s2)*decay+
	 2*(-4/r2/r2+24*x*x*y*y/r2/r2/r2/r2+1/r2/r2)*oneminusdecay);
}

double fyyy(double x,double y,double r2,double s2,
	    double oneminusdecay,double decay)
{
  return((12*y*y/r2/r2/s2-1.5/r2/s2+1.5*y*y/r2/s2/s2-12*y*y*y*y/r2/r2/r2/s2-
	  1.5*y*y*y*y/r2/r2/s2/s2-0.125*y*y*y*y/r2/s2/s2/s2)*decay+
	 6*(1/r2/r2-8*y*y/r2/r2/r2+8*y*y*y*y/r2/r2/r2/r2)*oneminusdecay);
}

double fxxxx(double x,double y,double r2,double s2,
	     double oneminusdecay,double decay)
{
  return((96*x*x*x*x*y/r2/r2/r2/r2/s2-
	  72*y*x*x/r2/r2/r2/s2+
	  12*x*x*x*x*y/r2/r2/r2/s2/s2+
	  6*y/r2/r2/s2-
	  9*x*x*y/r2/r2/s2/s2+
	  x*x*x*x*y/r2/r2/s2/s2/s2+
	  3*y/4/r2/s2/s2-
	  3*x*x*y/4/r2/s2/s2/s2+
	  x*x*x*x*y/16/r2/s2/s2/s2/s2)*decay+
	 24*(-16*x*x*x*x*y/r2/r2/r2/r2/r2+
	     12*x*x*y/r2/r2/r2/r2-
	     y/r2/r2/r2)*oneminusdecay);
}

double fxxxy(double x,double y,double r2,double s2,
	     double oneminusdecay,double decay)
{
  return((96*x*x*x*y*y/r2/r2/r2/r2/s2+
	  12*x*x*x*y*y/r2/r2/r2/s2/s2-
	  12*x*x*x/r2/r2/r2/s2-
	  36*x*y*y/r2/r2/r2/s2+
	  6*x/r2/r2/s2-
	  9*x*y*y/2/r2/r2/s2/s2-
	  3*x*x*x/2/r2/r2/s2/s2+
	  x*x*x*y*y/r2/r2/s2/s2/s2+
	  3*x/4/r2/s2/s2-
	  3*x*y*y/8/r2/s2/s2/s2-
	  x*x*x/8/r2/s2/s2/s2+
	  x*x*x*y*y/16/r2/s2/s2/s2/s2)*decay+
	 24*(2*x*x*x/r2/r2/r2/r2+
	     6*x*y*y/r2/r2/r2/r2-
	     16*x*x*x*y*y/r2/r2/r2/r2/r2-
	     x/r2/r2/r2)*oneminusdecay);
}

double fxxyy(double x,double y,double r2,double s2,
	     double oneminusdecay,double decay)
{
  return((-12*y*y*y/r2/r2/r2/s2-
	  3*y*y*y/2/r2/r2/s2/s2-
	  36*x*x*y/r2/r2/r2/s2+
	  6*y/r2/r2/s2-
	  9*x*x*y/2/r2/r2/s2/s2+
	  3*y/4/r2/s2/s2-
	  3*x*x*y/8/r2/s2/s2/s2+
	  96*x*x*y*y*y/r2/r2/r2/r2/s2+
	  12*x*x*y*y*y/r2/r2/r2/s2/s2+
	  x*x*y*y*y/r2/r2/s2/s2/s2-
	  y*y*y/8/r2/s2/s2/s2+
	  x*x*y*y*y/16/r2/s2/s2/s2/s2)*decay+
	 24*(2*y*y*y/r2/r2/r2/r2+
	     6*x*x*y/r2/r2/r2/r2-
	     y/r2/r2/r2-
	     16*x*x*y*y*y/r2/r2/r2/r2/r2)*oneminusdecay);
}

double fxyyy(double x,double y,double r2,double s2,
	     double oneminusdecay,double decay)
{
  return((-72*x*y*y/r2/r2/r2/s2+
	  6*x/r2/r2/s2-
	  9*x*y*y/r2/r2/s2/s2+
	  96*x*y*y*y*y/r2/r2/r2/r2/s2+
	  12*x*y*y*y*y/r2/r2/r2/s2/s2+
	  x*y*y*y*y/r2/r2/s2/s2/s2+
	  3*x/4/r2/s2/s2-
	  3*x*y*y/4/r2/s2/s2/s2+
	  x*y*y*y*y/16/r2/s2/s2/s2/s2)*decay+
	 24*(-x/r2/r2/r2+
	     12*x*y*y/r2/r2/r2/r2-
	     16*x*y*y*y*y/r2/r2/r2/r2/r2)*oneminusdecay);
}

double fyyyy(double x,double y,double r2,double s2,
	     double oneminusdecay,double decay)
{
  return((-120*y*y*y/r2/r2/r2/s2-
	  15*y*y*y/r2/r2/s2/s2+
	  30*y/r2/r2/s2+
	  15*y/4/r2/s2/s2-
	  5*y*y*y/4/r2/s2/s2/s2+
	  96*y*y*y*y*y/r2/r2/r2/r2/s2+
	  12*y*y*y*y*y/r2/r2/r2/s2/s2+
	  y*y*y*y*y/r2/r2/s2/s2/s2+
	  y*y*y*y*y/16/r2/s2/s2/s2/s2)*decay+
	 24*(20*y*y*y/r2/r2/r2/r2-
	     16*y*y*y*y*y/r2/r2/r2/r2/r2-
	     5*y/r2/r2/r2)*oneminusdecay);
}

/* This is an exact velocity field for a Lamb vortex.  This subroutine 
   is used for diagnostic purposes only. */
/* This one has a circulation of Pi. */

void dpos_vel_linear(int vort)
{
  double x,y,s2,r2,decay,oneminusdecay;
  double Uxx,Uyy,Uxy,Vxx,Vyy,Vxy;
  double Uxxx,Uxyy,Uxxy,Uyyy,Vxxx,Vxyy,Vxxy,Vyyy;
  double Uxxxx,Uxxxy,Uxxyy,Uxyyy,Uyyyy;
  double Vxxxx,Vxxxy,Vxxyy,Vxyyy,Vyyyy;
  double termA,termB,termC;

  x  = mblob[vort].blob0.x;
  y  = mblob[vort].blob0.y;
  s2 = 0.0625+visc*SimTime;

  r2 = SQR(mblob[vort].blob0.x) + SQR(mblob[vort].blob0.y);

  decay = exp(-r2/4.0/s2);

  oneminusdecay = 1.0-decay;

  termA = (tmpparms[vort].cos2*blobguts[vort].a2+
	   tmpparms[vort].sin2/blobguts[vort].a2);

  termB = (tmpparms[vort].sin2*blobguts[vort].a2+
	   tmpparms[vort].cos2/blobguts[vort].a2);

  termC = tmpparms[vort].sincos*(blobguts[vort].a2-1.0/blobguts[vort].a2);

  if (r2 == 0.0)
    {
      mblob[vort].blob0.dx = 0.0;
      mblob[vort].blob0.dy = 0.0;

      tmpparms[vort].du11 = 0.0;
      tmpparms[vort].du12 = -0.125/(0.0625+visc*SimTime);
      tmpparms[vort].du21 = 0.125/(0.0625+visc*SimTime);
    }
  else
    {
      Uxx = 0.5*fxx(x,y,r2,s2,oneminusdecay,decay);
      Uyy = 0.5*fyy(x,y,r2,s2,oneminusdecay,decay);
      Vxx = -0.5*fyy(y,x,r2,s2,oneminusdecay,decay);
      Vyy = -0.5*fxx(y,x,r2,s2,oneminusdecay,decay);

      Uxy = -Vyy;
      Vxy = -Uxx;

      Uxxx = 0.5*fxxx(x,y,r2,s2,oneminusdecay,decay);
      Uxyy = 0.5*fxxx(y,x,r2,s2,oneminusdecay,decay);
      Uxxy = 0.5*fxxy(x,y,r2,s2,oneminusdecay,decay);
      Uyyy = 0.5*fyyy(x,y,r2,s2,oneminusdecay,decay);

      Vxxx = -0.5*fyyy(y,x,r2,s2,oneminusdecay,decay);
      Vxyy = -Uxxy;
      Vxxy = -Uxxx;
      Vyyy = -Uxyy;

      Uxxxx = 0.5*fxxxx(x,y,r2,s2,oneminusdecay,decay);
      Uxxxy = 0.5*fxxxy(x,y,r2,s2,oneminusdecay,decay);
      Uxxyy = 0.5*fxxyy(x,y,r2,s2,oneminusdecay,decay);
      Uxyyy = 0.5*fxyyy(x,y,r2,s2,oneminusdecay,decay);
      Uyyyy = 0.5*fyyyy(x,y,r2,s2,oneminusdecay,decay);

      Vxxxx = -0.5*fyyyy(y,x,r2,s2,oneminusdecay,decay);
      Vxxxy = -Uxxxx;
      Vxxyy = -Uxxxy;
      Vxyyy = -Uxxyy;
      Vyyyy = -Uxyyy;

      mblob[vort].blob0.dx =
	-0.5*mblob[vort].blob0.y/r2*oneminusdecay+

	correction1*blobguts[vort].s2*
	(Uxx*termA + Uyy*termB + 2.0*Uxy*termC)-

	correction3*1.5*SQR(blobguts[vort].s2)*
	(Uxxxx*SQR(termA)+
	 4*Uxxxy*termA*termC+
	 2*Uxxyy*(termA*termB+2*SQR(termC))+
	 4*Uxyyy*termB*termC+
	 Uyyyy*SQR(termB));
      
      mblob[vort].blob0.dy =
	0.5*mblob[vort].blob0.x/r2*oneminusdecay+

	correction1*blobguts[vort].s2*
	(Vxx*termA + Vyy*termB + 2.0*Vxy*termC)-

	correction3*1.5*SQR(blobguts[vort].s2)*
	(Vxxxx*SQR(termA)+
	 4*Vxxxy*termA*termC+
	 2*Vxxyy*(termA*termB+2*SQR(termC))+
	 4*Vxyyy*termB*termC+
	 Vyyyy*SQR(termB));

      tmpparms[vort].du11 = mblob[vort].blob0.x*mblob[vort].blob0.y*
	(oneminusdecay/r2/r2-0.25*decay/r2/(0.0625+visc*SimTime))-

	correction2*blobguts[vort].s2*
	(Uxxx*termA + Uxyy*termB + 2.0*Uxxy*termC);

      tmpparms[vort].du12 = -0.5*oneminusdecay/r2+
	SQR(mblob[vort].blob0.y/r2)*oneminusdecay-
	0.25*SQR(mblob[vort].blob0.y)*decay/r2/(0.0625+visc*SimTime)-

	correction2*blobguts[vort].s2*
	(Uxxy*termA + Uyyy*termB + 2.0*Uxyy*termC);

      tmpparms[vort].du21 = 0.5*oneminusdecay/r2-
	SQR(mblob[vort].blob0.x/r2)*oneminusdecay+
	0.25*SQR(mblob[vort].blob0.x)*decay/r2/(0.0625+visc*SimTime)-

	correction2*blobguts[vort].s2*
	(Vxxx*termA + Vxyy*termB + 2.0*Vxxy*termC);
    }
}

Vector dpos_vel_gen_linear(Vector pos,Blob_internal the_blobguts,
			   Blob_parms parms)
{
  Vector v;

  double s2,r2,decay,oneminusdecay;
  double Uxx,Uyy,Uxy,Vxx,Vyy,Vxy;
  double termA,termB,termC;

  s2 = 0.0625+visc*SimTime;

  r2 = SQR(pos.x) + SQR(pos.y);

  decay = exp(-r2/4.0/s2);

  oneminusdecay = 1.0-decay;

  termA = (parms.cos2*the_blobguts.a2 + parms.sin2/the_blobguts.a2);

  termB = (parms.sin2*the_blobguts.a2 + parms.cos2/the_blobguts.a2);

  termC = parms.sincos*(the_blobguts.a2-1.0/the_blobguts.a2);


  decay = exp(-r2/4.0/(0.0625+visc*SimTime));

  oneminusdecay = (1.0-exp(-r2/4.0/(0.0625+visc*SimTime)));

  Uxx = 0.5*fxx(pos.x,pos.y,r2,s2,oneminusdecay,decay);
  Uyy = 0.5*fyy(pos.x,pos.y,r2,s2,oneminusdecay,decay);
  Vxx = -0.5*fyy(pos.y,pos.x,r2,s2,oneminusdecay,decay);
  Vyy = -0.5*fxx(pos.y,pos.x,r2,s2,oneminusdecay,decay);
  
  Uxy = -Vyy;
  Vxy = -Uxx;

  if (r2 == 0.0)
    {
      v.x = 0.0;
      v.y = 0.0;
    }
  else
    {
      v.x = -0.5*pos.y/r2*oneminusdecay +

	correction1*the_blobguts.s2*
	(Uxx*termA + Uyy*termB + 2.0*Uxy*termC);
     
      v.y =  0.5*pos.x/r2*oneminusdecay+

	correction1*the_blobguts.s2*
	(Vxx*termA + Vyy*termB + 2.0*Vxy*termC);
    }
  return(v);
}
