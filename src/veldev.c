/* VELDEV.C */
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

#include <stdio.h>
#include <math.h>
#include "global_min.h"
#include "particle.h"
#include "biot-savart.h"

Tensor induced_veldev(Blob_external *the_blob,Blob_internal *the_blobguts,
		      Blob_parms *parms,double tmpdx,double tmpdy,
		      double r[5],double t[5],double evena[5],double odda[5])
{
   double eps,dx,dy;
   double evenb[5],oddb[5];
   double c0,c1,c2,c3,c4;
   double s2,s4,s6,s8,s10;
   
   double tempa,tempb,tempc;
   
   Tensor result;
   
   /* Change bases to the local axes. */
   
   dx =  (*parms).costh*tmpdx+(*parms).sinth*tmpdy;
   dy = -(*parms).sinth*tmpdx+(*parms).costh*tmpdy;
   
   eps = (sqrt((*the_blobguts).a2)-1.0)/(sqrt((*the_blobguts).a2)+1.0);
   
   s2 = (*the_blobguts).s2;
   s4 = SQR(s2);
   s6 = s2*s4;
   s8 = s4*s4;
   s10 = s4*s6;

   evenb[4] = -224.0+13440*t[1]-32256.0*t[3];
   oddb[4]  = -10752.0*t[0]+28672.0*t[2];

   evenb[3] = eps*(1920.0-38400.0*t[1]+71680.0*t[3])+480.0*t[0]-2240.0*t[2];
   oddb[3]  = eps*(21760.0*t[0]-53760.0*t[2])-320.0+1920.0*t[1];
   
   evenb[2] = SQR(eps)*(-2768.0+36160.0*t[1]-53760.0*t[3])+
     eps*(-1408*t[0]+3840.0*t[2])+
     8.0-160.0*t[1];
   oddb[2]  = SQR(eps)*(-13568.0*t[0]+30720.0*t[2])+
     eps*(512.0-2560.0*t[1])+
     128.0*t[0];
   
   evenb[1] = SQR(eps)*eps*(1248.0-12672.0*t[1]+15360.0*t[3])+
     SQR(eps)*(1020.0*t[0]-1920.0*t[2])+
     eps*(-48.0+192.0*t[1])-12.0*t[0];
   oddb[1]  = SQR(eps)*eps*(2496.0*t[0]-5120.0*t[2])+
     SQR(eps)*(-168.0+768.0*t[1])+
     eps*(-96.0*t[0])+8.0;
   
   evenb[0] = SQR(SQR(eps))*(-160.0+1248.0*t[1]-1280.0*t[3])+
     eps*SQR(eps)*(-168.0*t[0]+256.0*t[2])+
     SQR(eps)*(16.0-48.0*t[1])+
     eps*(8.0*t[0])-1.0;
   oddb[0]  = 0.0;
   
   c4 = SQR(SQR(eps))*(0.5*(evena[4]+odda[4])+SQR(dx)*(evenb[4]+oddb[4])/r[0])/r[4];
   c3 = eps*SQR(eps)*(0.5*(evena[3]+odda[3])+SQR(dx)*(evenb[3]+oddb[3])/r[0])/r[3];
   c2 = SQR(eps)*(0.5*(evena[2]+odda[2])+SQR(dx)*(evenb[2]+oddb[2])/r[0])/r[2];
   c1 = eps*(0.5*(evena[1]+odda[1])+SQR(dx)*(evenb[1]+oddb[1])/r[0])/r[1];
   c0 = (0.5*(evena[0]+odda[0])+SQR(dx)*(evenb[0]+oddb[0])/r[0])/r[0];
   
  if (r[0]/s2<0.009)
    tempc = ((*the_blob).strength/(2.0*s2))*
      (0.5-eps+eps*SQR(eps))*exp(-r[0]/(4.0*s2));
  else
    tempc = -((*the_blob).strength/(2.0*s2))*
      ((r[0]*((c0+s2*(8.0*c1+s2*(96.0*c2+s2*(1536.0*c3+30720.0*c4*s2))))+
	     r[0]*((c1+s2*(12.0*c2+s2*(192.0*c3+3840.0*c4*s2)))+
		  r[0]*((c2+s2*(16.0*c3+320.0*c4*s2))+
		       r[0]*((c3+20.0*c4*s2)+
			    c4*r[0]))))-
	(0.5-eps+eps*SQR(eps)))*exp(-r[0]/(4.0*s2))+
       (-s2)*(4.0*c0+s2*(32.0*c1+s2*(384.0*c2+s2*(6144.0*c3+s2*122880.0*c4))))*
       (1-exp(-r[0]/(4.0*s2))));

   c4 = SQR(SQR(eps))*(0.5*(evena[4]-odda[4])+SQR(dy)*(evenb[4]-oddb[4])/r[0])/r[4];
   c3 = eps*SQR(eps)*(0.5*(evena[3]-odda[3])+SQR(dy)*(evenb[3]-oddb[3])/r[0])/r[3];
   c2 = SQR(eps)*(0.5*(evena[2]-odda[2])+SQR(dy)*(evenb[2]-oddb[2])/r[0])/r[2];
   c1 = eps*(0.5*(evena[1]-odda[1])+SQR(dy)*(evenb[1]-oddb[1])/r[0])/r[1];
   c0 = (0.5*(evena[0]-odda[0])+SQR(dy)*(evenb[0]-oddb[0])/r[0])/r[0];
   
  if (r[0]/s2<0.009)
    tempb = -((*the_blob).strength/(2.0*s2))*
      (0.5+eps-eps*SQR(eps))*exp(-r[0]/(4.0*s2));
  else
    tempb = ((*the_blob).strength/(2.0*s2))*
      ((r[0]*((c0+s2*(8.0*c1+s2*(96.0*c2+s2*(1536.0*c3+30720.0*c4*s2))))+
	     r[0]*((c1+s2*(12.0*c2+s2*(192.0*c3+3840.0*c4*s2)))+
		  r[0]*((c2+s2*(16.0*c3+320.0*c4*s2))+
		       r[0]*((c3+20.0*c4*s2)+
			    c4*r[0]))))-
	(0.5+eps-eps*SQR(eps)))*exp(-r[0]/(4.0*s2))+
       (-s2)*(4.0*c0+s2*(32.0*c1+s2*(384.0*c2+s2*(6144.0*c3+s2*122880.0*c4))))*
       (1-exp(-r[0]/(4.0*s2))));

   evena[4] = -2016.0+24192*t[1]-32256.0*t[3];
   evena[3] = eps*(4480.0-53760.0*t[1]+71680.0*t[3])+
     1120.0*t[0]-2240.0*t[2];
   evena[2] = SQR(eps)*(-3600.0+41280.0*t[1]-53760.0*t[3])+
     eps*(-1920.0*t[0]+3840.0*t[2])+40.0-160.0*t[1];
   evena[1] = evenb[1];
   evena[0] = evenb[0];
   
   c4 = SQR(SQR(eps))*evena[4]/(r[4]*r[0]);
   c3 = eps*SQR(eps)*evena[3]/r[4];
   c2 = SQR(eps)*evena[2]/r[3];
   c1 = eps*evena[1]/r[2];
   c0 = evena[0]/r[1];
   
  if (r[0]/s2<0.009)
    tempa = 0.0;
  else
    tempa = ((*the_blob).strength/(2.0*s2))*dx*dy*
      ((r[0]*((c0+s2*(8.0*c1+s2*(96.0*c2+s2*(1536.0*c3+30720.0*c4*s2))))+
	     r[0]*((c1+s2*(12.0*c2+s2*(192.0*c3+3840.0*c4*s2)))+
		  r[0]*((c2+s2*(16.0*c3+320.0*c4*s2))+
		       r[0]*((c3+20.0*c4*s2)+
			    c4*r[0])))))*exp(-r[0]/(4.0*s2)) +
       (-s2)*(4.0*c0+s2*(32.0*c1+s2*(384.0*c2+s2*(6144.0*c3+s2*122880.0*c4))))*
       (1-exp(-r[0]/(4.0*s2))));

   result.du11 = (tempa*((*parms).cos2-(*parms).sin2)-
		  (tempb+tempc)*(*parms).sincos);
   
   result.du12 = (2.0*tempa*(*parms).sincos+
		  tempb*(*parms).cos2-tempc*(*parms).sin2);
   
   result.du21 = (2.0*tempa*(*parms).sincos-
      tempb*(*parms).sin2+tempc*(*parms).cos2);
   
   return(result);
}
