/* VELOCITY.C */
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
 * or until 15 January 2001                                             *
 * Louis Rossi                                                          *
 * Department of Mathematical Sciences                                  *
 * University of Massachusetts Lowell                                   *
 * One University Ave                                                   *
 * Lowell, MA 01854                                                     *
 * 
 * or after 15 January 2001                                             *
 * Louis Rossi                                                          *
 * Department of Mathematical Sciences                                  *
 * University of Delaware                                               *
 * Newark, DE 19715-2553                                                */

#include <stdio.h>
#include <math.h>
#include "global.h"

vector induced_vel(the_blob,the_blobguts,parms,tmpdx,tmpdy,r,t,even,odd)
blob_external *the_blob;
blob_internal *the_blobguts;
blobparms *parms;
double tmpdx,tmpdy;
double r[5],t[5];
double even[5],odd[5];
{
   double eps,dx,dy;
   double c0,c1,c2,c3,c4;
   double s2,s4,s6,s8,s10;
   
   double tempu,tempv;
   
   vector result;
   
   /* Change bases to the local axes. */
   
   dx =  (*parms).costh*tmpdx+(*parms).sinth*tmpdy;
   
   dy = -(*parms).sinth*tmpdx+(*parms).costh*tmpdy;
   
   eps = (sqrt((*the_blobguts).a2)-1.0)/(sqrt((*the_blobguts).a2)+1.0);
   
   s2 = (*the_blobguts).s2;
   s4 = SQR(s2);
   s6 = s2*s4;
   s8 = s4*s4;
   s10 = s4*s6;
   
   r[0] = SQR(dx*(1-eps)/(1+eps))+SQR(dy*(1+eps)/(1-eps));
   t[0] = (SQR(dx*(1-eps)/(1+eps))-SQR(dy*(1+eps)/(1-eps)))/r[0];
   
   r[1] = SQR(r[0]);
   t[1] = SQR(t[0]);
   
   r[2] = r[1]*r[0];
   t[2] = t[1]*t[0];

   r[3] = SQR(r[1]);
   t[3] = SQR(t[1]);
   
   r[4] = r[2]*r[1];
   t[4] = t[2]*t[1];
   
   
   even[4] = (224.0-2688.0*t[1]+3584.0*t[3]);
   odd[4] =  (896.0*t[0]-1792.0*t[2]);
   
   even[3] = eps*(-640.0+7040.0*t[1]-8960.0*t[3])-160.0*t[0]+320.0*t[2];
   odd[3]  = eps*(-1920.0*t[0]+3840.0*t[2])+40.0-160.0*t[1];
   
   even[2] = SQR(eps)*(656.0-6464.0*t[1]+7680.0*t[3])+
     eps*(352.0*t[0]-640.0*t[2])-8.0+32.0*t[1];
   odd[2]  = SQR(eps)*(1312.0*t[0]-2560.0*t[2])+
     eps*(-64.0+256.0*t[1])-16.0*t[0];
   
   even[1] = SQR(eps)*eps*(-288.0+2400.0*t[1]-2560.0*t[3])+
     SQR(eps)*(-244.0*t[0]+384.0*t[2])+
     eps*(16.0-48.0*t[1])+4.0*t[0];
   odd[1]  = SQR(eps)*eps*(-288.0*t[0]+512.0*t[2])+
     SQR(eps)*(26.0-96.0*t[1])+
     eps*16.0*t[0]-2.0;
   
   even[0] = SQR(eps)*SQR(eps)*(48.0-288.0*t[1]+256.0*t[3])+
     SQR(eps)*eps*(52.0*t[0]-64.0*t[2])+
     SQR(eps)*(-8.0+16.0*t[1])-eps*4.0*t[0]+1.0;
   odd[0] = 0.0;
   

   
   c4 = SQR(SQR(eps))*(even[4]+odd[4])/r[4];
   c3 = SQR(eps)*eps*(even[3]+odd[3])/r[3];
   c2 = SQR(eps)*(even[2]+odd[2])/r[2];
   c1 = eps*(even[1]+odd[1])/r[1];
   c0 = (even[0]+odd[0])/r[0];
   
   if (r[0]/s2<0.009)
     tempv = (dx/2.0)*((*the_blob).strength/(2.0*s2))*
     2.0*(0.5-eps+eps*SQR(eps))*exp(-r[0]/(4.0*s2));
   else
     tempv = -(dx/2.0)*((*the_blob).strength/(2.0*s2))*
     ((r[0]*((c0+s2*(8.0*c1+s2*(96.0*c2+s2*(1536.0*c3+30720.0*c4*s2))))+
	    r[0]*((c1+s2*(12.0*c2+s2*(192.0*c3+3840.0*c4*s2)))+
		 r[0]*((c2+s2*(16.0*c3+320.0*c4*s2))+
		       r[0]*((c3+20.0*c4*s2)+
			    c4*r[0]))))-
	2.0*(0.5-eps+eps*SQR(eps)))*exp(-r[0]/(4.0*s2))+
       (-s2)*(4.0*c0+s2*(32.0*c1+s2*(384.0*c2+s2*(6144.0*c3+s2*122880.0*c4))))*
       (1-exp(-r[0]/(4.0*s2))));

   c4 = SQR(SQR(eps))*(even[4]-odd[4])/r[4];
   c3 = SQR(eps)*eps*(even[3]-odd[3])/r[3];
   c2 = SQR(eps)*(even[2]-odd[2])/r[2];
   c1 = eps*(even[1]-odd[1])/r[1];
   c0 = (even[0]-odd[0])/r[0];
   
  if (r[0]/s2<0.009)
    tempu = -(dy/2.0)*((*the_blob).strength/(2.0*s2))*
      2.0*(0.5+eps-eps*SQR(eps))*exp(-r[0]/(4.0*s2));
  else
    tempu = (dy/2.0)*((*the_blob).strength/(2.0*s2))*
      ((r[0]*((c0+s2*(8.0*c1+s2*(96.0*c2+s2*(1536.0*c3+30720.0*c4*s2))))+
	     r[0]*((c1+s2*(12.0*c2+s2*(192.0*c3+3840.0*c4*s2)))+
		  r[0]*((c2+s2*(16.0*c3+320.0*c4*s2))+
		       r[0]*((c3+20.0*c4*s2)+
			    c4*r[0]))))-
	2.0*(0.5+eps-eps*SQR(eps)))*exp(-r[0]/(4.0*s2))+
       (-s2)*(4.0*c0+s2*(32.0*c1+s2*(384.0*c2+s2*(6144.0*c3+s2*122880.0*c4))))*
       (1-exp(-r[0]/(4.0*s2))));

   /* Rotate back to the global axes. */
   
   result.x = ((*parms).costh*tempu-(*parms).sinth*tempv);
   result.y = ((*parms).sinth*tempu+(*parms).costh*tempv);
   
   return(result);
}
