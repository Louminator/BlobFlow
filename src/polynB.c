/* POLYNB.C */
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

#define xi  0.1

#define xi1  3.0*xi
#define xi2 -3.5*xi
#define xi3  1.0*xi

#define correction1 1.0
#define correction2 0.0
#define correction3 0.0

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

  x  = mblob[vort].blob0.x;
  y  = mblob[vort].blob0.y;

  r2 = SQR(x) + SQR(y);

  termA = (tmpparms[vort].cos2*blobguts[vort].a2+
	   tmpparms[vort].sin2/blobguts[vort].a2);

  termB = (tmpparms[vort].sin2*blobguts[vort].a2+
	   tmpparms[vort].cos2/blobguts[vort].a2);

  termC = tmpparms[vort].sincos*(blobguts[vort].a2-1.0/blobguts[vort].a2);

  U = -M_PI/2*y - xi1*2*y - xi2*4*r2*y - xi3*6*r2*r2*y;
  V =  M_PI/2*x + xi1*2*x + xi2*4*r2*x + xi3*6*r2*r2*x;

  Ux =                  -xi2*8*x*y        - xi3*24*r2*x*y;
  Uy = -M_PI/2 - xi1*2 - xi2*4*(r2+2*y*y) - xi3*6*(r2*r2 + 4*r2*y*y);

  Vx =  M_PI/2 + xi1*2 + xi2*4*(r2+2*x*x) + xi3*6*(r2*r2 + 4*r2*x*x);
  Vy = -Ux;

  Uxx = -xi2*8*y  - xi3*24*(r2*y + 2*x*x*y);
  Uxy = -xi2*8*x  - xi3*24*(r2*x + 2*x*y*y);
  Uyy = -xi2*24*y - xi3*24*(3*r2*y + 2*y*y*y);
  Vxx =  xi2*24*x + xi3*24*(3*r2*x + 2*x*x*x);
  Vxy = -Uxx;
  Vyy = -Uxy;

  Uxxx =         - xi3*144*x*y;
  Uxxy = -xi2*8  - xi3*72*r2;
  Uxyy =         - xi3*144*x*y;
  Uyyy = -xi2*24 - xi3*72*(r2 + 4*y*y);

  Vxxx =  xi2*24 + xi3*72*(r2 + 4*x*x);
  Vxxy = -Uxxx;
  Vxyy = -Uxxy;
  Vyyy = -Uxyy;

  Uxxxx = -xi3*144*y;
  Uxxxy = -xi3*144*x;
  Uxxyy = -xi3*144*y;
  Uxyyy = -xi3*144*x;
  Uyyyy = -xi3*720*y;

  Vxxxx =  xi3*720*x;
  Vxxxy = -Uxxxx;
  Vxxyy = -Uxxxy;
  Vxyyy = -Uxxyy;
  Vyyyy = -Uxyyy;

  mblob[vort].blob0.dx = U +

    correction1*blobguts[vort].s2*
    (Uxx*termA + Uyy*termB + 2.0*Uxy*termC)+

    correction3*2.5*SQR(blobguts[vort].s2)*
    (Uxxxx*SQR(termA)+
     4*Uxxxy*termA*termC+
     2*Uxxyy*(termA*termB+2*SQR(termC))+
     4*Uxyyy*termB*termC+
     Uyyyy*SQR(termB));
      
  mblob[vort].blob0.dy = V +

    correction1*blobguts[vort].s2*
    (Vxx*termA + Vyy*termB + 2.0*Vxy*termC)+

    correction3*2.5*SQR(blobguts[vort].s2)*
    (Vxxxx*SQR(termA)+
     4*Vxxxy*termA*termC+
     2*Vxxyy*(termA*termB+2*SQR(termC))+
     4*Vxyyy*termB*termC+
     Vyyyy*SQR(termB));

  tmpparms[vort].du11 = Ux +

    correction2*blobguts[vort].s2*
    (-termA*Uxxx-2*termC*Uxxy-termB*Uxyy);

  tmpparms[vort].du12 = Uy +

    correction2*blobguts[vort].s2*
    (-termA*Uxxy-2*termC*Uxyy-termB*Uyyy);

  tmpparms[vort].du21 = Vx +

    correction2*blobguts[vort].s2*
    (-termA*Vxxx-2*termC*Vxxy-termB*Vxyy);

}

vector dpos_vel_gen_linear(pos,the_blobguts,parms)
     vector    pos;
     blob_internal the_blobguts;
     blobparms parms;
{
  vector v;

  double x,y,r2;
  double U,V;
  double Uxx,Uyy,Uxy,Vxx,Vyy,Vxy;
  double termA,termB,termC;
  double Uxxxx,Uxxxy,Uxxyy,Uxyyy,Uyyyy;
  double Vxxxx,Vxxxy,Vxxyy,Vxyyy,Vyyyy;

  r2 = SQR(pos.x) + SQR(pos.y);

  termA = (parms.cos2*the_blobguts.a2 + parms.sin2/the_blobguts.a2);

  termB = (parms.sin2*the_blobguts.a2 + parms.cos2/the_blobguts.a2);

  termC = parms.sincos*(the_blobguts.a2-1.0/the_blobguts.a2);

  x = pos.x;
  y = pos.y;

  U = -M_PI/2*y - xi1*2*y - xi2*4*r2*y - xi3*6*r2*r2*y;
  V =  M_PI/2*x + xi1*2*x + xi2*4*r2*x + xi3*6*r2*r2*x;

  Uxx = -xi2*8*y  - xi3*24*(r2*y+2*x*x*y);
  Uxy = -xi2*8*x  - xi3*24*(r2*x+2*x*y*y);
  Uyy = -xi2*24*y - xi3*24*(3*r2*y+2*y*y*y);
  Vxx =  xi2*24*x + xi3*24*(3*r2*y+2*y*y*y);
  Vxy = -Uxx;
  Vyy = -Uxy;

  Uxxxx = -xi3*144*y;
  Uxxxy = -xi3*144*x;
  Uxxyy = -xi3*144*y;
  Uxyyy = -xi3*144*x;
  Uyyyy = -xi3*720*y;

  Vxxxx =  xi3*720*x;
  Vxxxy = -Uxxxx;
  Vxxyy = -Uxxxy;
  Vxyyy = -Uxxyy;
  Vyyyy = -Uxyyy;

  v.x = U +

    correction1*the_blobguts.s2*
    (Uxx*termA + Uyy*termB + 2.0*Uxy*termC) -

    correction3*1.5*SQR(the_blobguts.s2)*
    (Uxxxx*SQR(termA)+
     4*Uxxxy*termA*termC+
     2*Uxxyy*(termA*termB+2*SQR(termC))+
     4*Uxyyy*termB*termC+
     Uyyyy*SQR(termB));

  v.y = V +

    correction1*the_blobguts.s2*
    (Vxx*termA + Vyy*termB + 2.0*Vxy*termC)-

    correction3*1.5*SQR(the_blobguts.s2)*
    (Vxxxx*SQR(termA)+
     4*Vxxxy*termA*termC+
     2*Vxxyy*(termA*termB+2*SQR(termC))+
     4*Vxyyy*termB*termC+
     Vyyyy*SQR(termB));

  return(v);
}
