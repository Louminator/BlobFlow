/* POLYNA.C */
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

#include "global.h"

#ifdef XANTISYMM
#define CARDINALITY  N/2
#else
#define CARDINALITY  N
#endif

#ifdef MULTIPROC
#include "multiproc.h"
#endif

/* This is an exact velocity field for a perturbed linear vortex with a 
   turnover time of 4. */

void dpos_vel_linear(vort)
int vort;
{
  double x,y,r2;
  double U,V,Ux,Uy,Vx,Vy;
  double Uxx,Uyy,Uxy,Vxx,Vyy,Vxy;
  double Uxxx,Uxyy,Uxxy,Uyyy,Vxxx,Vxyy,Vxxy,Vyyy;
  double Uxxxx,Uxxxy,Uxxyy,Uxyyy,Uyyyy;
  double Vxxxx,Vxxxy,Vxxyy,Vxyyy,Vyyyy;
  double termA,termB,termC;
  double correction1=1.0,correction2=1.0,correction3=0.0;
  double xi=0.1;

  x  = mblob[vort].blob0.x;
  y  = mblob[vort].blob0.y;

  r2 = SQR(mblob[vort].blob0.x) + SQR(mblob[vort].blob0.y);

  termA = (tmpparms[vort].cos2*blobguts[vort].a2+
	   tmpparms[vort].sin2/blobguts[vort].a2);

  termB = (tmpparms[vort].sin2*blobguts[vort].a2+
	   tmpparms[vort].cos2/blobguts[vort].a2);

  termC = tmpparms[vort].sincos*(blobguts[vort].a2-1.0/blobguts[vort].a2);

  U = -(M_PI/2 + xi*(4*r2-3))*y;
  V =  (M_PI/2 + xi*(4*r2-3))*x;

  Ux = -8*xi*x*y;
  Uy = -(M_PI/2 + xi*(4*r2-3)) - 8*xi*y*y;

  Vx =  (M_PI/2 + xi*(4*r2-3)) + 8*xi*x*x;
  Vy = -Ux;

  Uxx = -8*xi*y;
  Uxy = -8*xi*x;
  Uyy = -24*xi*y;
  Vxx =  24*xi*x;
  Vxy = -Uxx;
  Vyy = -Uxy;

  Uxxx =  0;
  Uxxy = -8*xi;
  Uxyy =  0;
  Uyyy = -24*xi;

  Vxxx =  24*xi;
  Vxxy = -Uxxx;
  Vxyy = -Uxxy;
  Vyyy = -Uxyy;

  Uxxxx = 0;
  Uxxxy = 0;
  Uxxyy = 0;
  Uxyyy = 0;
  Uyyyy = 0;

  Vxxxx = 0;
  Vxxxy = -Uxxxx;
  Vxxyy = -Uxxxy;
  Vxyyy = -Uxxyy;
  Vyyyy = -Uxyyy;

  mblob[vort].blob0.dx = U +

    correction1*blobguts[vort].s2*
    (Uxx*termA + Uyy*termB + 2.0*Uxy*termC)-

    correction3*1.5*SQR(blobguts[vort].s2)*
    (Uxxxx*SQR(termA)+
     4*Uxxxy*termA*termC+
     2*Uxxyy*(termA*termB+2*SQR(termC))+
     4*Uxyyy*termB*termC+
     Uyyyy*SQR(termB));
      
  mblob[vort].blob0.dy = V +

    correction1*blobguts[vort].s2*
    (Vxx*termA + Vyy*termB + 2.0*Vxy*termC)-

    correction3*1.5*SQR(blobguts[vort].s2)*
    (Vxxxx*SQR(termA)+
     4*Vxxxy*termA*termC+
     2*Vxxyy*(termA*termB+2*SQR(termC))+
     4*Vxyyy*termB*termC+
     Vyyyy*SQR(termB));

  tmpparms[vort].du11 = Ux -
    correction2*blobguts[vort].s2*
    (Uxxx*termA + Uxyy*termB + 2.0*Uxxy*termC);

  tmpparms[vort].du12 = Uy -
    correction2*blobguts[vort].s2*
    (Uxxy*termA + Uyyy*termB + 2.0*Uxyy*termC);

  tmpparms[vort].du21 = Vx -
    correction2*blobguts[vort].s2*
    (Vxxx*termA + Vxyy*termB + 2.0*Vxxy*termC);
}

Vector dpos_vel_gen_linear(pos,the_blobguts,parms)
     Vector    pos;
     Blob_internal the_blobguts;
     Blob_parms parms;
{
  Vector v;

  double x,y,r2;
  double U,V;
  double Uxx,Uyy,Uxy,Vxx,Vyy,Vxy;
  double termA,termB,termC;
  double correction1=1.0;
  double xi=0.1;

  r2 = SQR(pos.x) + SQR(pos.y);

  termA = (parms.cos2*the_blobguts.a2 + parms.sin2/the_blobguts.a2);

  termB = (parms.sin2*the_blobguts.a2 + parms.cos2/the_blobguts.a2);

  termC = parms.sincos*(the_blobguts.a2-1.0/the_blobguts.a2);

  x = pos.x;
  y = pos.y;

  U = -(M_PI/2 + xi*(4*r2-3))*y;
  V =  (M_PI/2 + xi*(4*r2-3))*x;

  Uxx = -8*xi*y;
  Uxy = -8*xi*x;
  Uyy = -24*xi*y;
  Vxx =  24*xi*x;
  Vxy = -Uxx;
  Vyy = -Uxy;

  v.x = U +

    correction1*the_blobguts.s2*
    (Uxx*termA + Uyy*termB + 2.0*Uxy*termC);
     
  v.y = V +

    correction1*the_blobguts.s2*
    (Vxx*termA + Vyy*termB + 2.0*Vxy*termC);

  return(v);
}
