/* VELDEVDEV.C */
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
 * 
 * Louis Rossi                                                          *
 * Department of Mathematical Sciences                                  *
 * University of Delaware                                               *
 * Newark, DE 19715-2553                                                */

#include <stdio.h>
#include <math.h>
#include "global.h"

Tensor induced_veldev(the_blob,the_blobguts,parms,tmpdx,tmpdy,r,t)
Blob_external *the_blob;
Blob_internal *the_blobguts;
Blob_parms *parms;
double tmpdx,tmpdy;
double r[5],t[5];
{
   double eps,dx,dy;
   double evena[5],odda[5],evenb[5],oddb[5],evenc[5],oddc[5],
     evend[5],oddd[5],evene[5],odde[5],evenf[5],oddf[5],
     eveng[5],oddg[5],evenh[5],oddh[5],eveni[5],oddi[5],evenj[5],oddj[5];
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

   /* Stuff for v_xx or u_yy */

   evena[4] = -18816.0 - 43008.0*t[1] + 645120.0*t[3];
   odda[4]  = 204288.0*t[0] - 774144.0*t[2];

   evenb[4] = -672.0 + 40320.0*t[1] - 96768.0*t[3];
   oddb[4]  = -32256.0*t[0] + 86016.0*t[2];


   evena[3] = 7680.0 + 376320.0*t[1] - 1290240.0*t[3];
   odda[3]  = -422400.0*t[0] + 1290240.0*t[2];

   evenb[3] = 5760.0 - 115200.0*t[1] + 215040.0*t[3];
   oddb[3]  = 65280.0*t[0] - 161280.0*t[2];

   evenc[3] = 1920.0*t[0] + 35840.0*t[2];
   oddc[3]  = 4160.0 - 40320.0*t[1];

   evend[3] = 1440.0*t[0] - 6720.0*t[2];
   oddd[3]  = -960.0 + 5760.0*t[1];


   evena[2] = 24192.0 - 449280.0*t[1] + 860160.0*t[3];
   odda[2]  = 249600.0*t[0] - 645120.0*t[2];

   evenb[2] = -8304.0 + 108480.0*t[1] - 161280.0*t[3];
   oddb[2]  = -40704.0*t[0] + 92160.0*t[2];

   evenc[2] = 11520.0*t[0] - 53760.0*t[2];
   oddc[2]  = -7680.0 + 46080.0*t[1];

   evend[2] = -4224.0*t[0] + 11520.0*t[2];
   oddd[2]  = 1536.0 - 7680.0*t[1];

   evene[2] = 192.0 + 1920.0*t[1];
   odde[2]  = -1920.0*t[0];

   evenf[2] = 24.0 - 480.0*t[1];
   oddf[2]  = 384.0*t[0];


   evena[1] = -13056.0 + 157440.0*t[1] - 215040.0*t[3];
   odda[1]  = -39168.0*t[0] + 92160.0*t[2];

   evenb[1] = 3744.0 - 38016.0*t[1] + 46080.0*t[3];
   oddb[1]  = 7488.0*t[0] - 15360.0*t[2];

   evenc[1] = -10464.0*t[0] + 23040.0*t[2];
   oddc[1]  = 2088.0 - 11520.0*t[1];

   evend[1] = 3060.0*t[0] - 5760.0*t[2];
   oddd[1]  = -504.0 + 2304.0*t[1];

   evene[1] = 384.0 - 1920.0*t[1];
   odde[1]  = 1152.0*t[0];

   evenf[1] = -144.0 + 576.0*t[1];
   oddf[1]  = -288.0*t[0];

   eveng[1] = 96.0*t[0];
   oddg[1]  = -72.0;

   evenh[1] = -36.0*t[0];
   oddh[1]  = 24.0;


   evena[0] = 1344.0 - 13056.0*t[1] + 15360.0*t[3];
   evenb[0] = -480.0 + 3744.0*t[1] -3840.0*t[3];
   evenc[0] = 1392.0*t[0] - 2560.0*t[2];
   evend[0] = -504.0*t[0] + 768.0*t[2];
   evene[0] = -96.0 + 384.0*t[1];
   evenf[0] = 48.0 - 144.0*t[1];
   eveng[0] = -48.0*t[0];
   evenh[0] = 24.0*t[0];
   eveni[0] = 4.0;
   evenj[0] = -3.0;

   c4 = ((evena[4] + odda[4])/r[0]*dx*dx + 
	 (evenb[4] + oddb[4]))/r[4]/r[0]*dx*SQR(eps)*SQR(eps);

   c3 = (((evena[3] + odda[3])*SQR(dx)/r[0] + (evenb[3] + oddb[3]))*eps + 
	 ((evenc[3] + oddc[3])*SQR(dx)/r[0] + (evend[3] + oddd[3]))
	 )*dx*SQR(eps)*eps/r[4];

   c2 = (((evena[2] + odda[2])*SQR(dx)/r[0]+(evenb[2] + oddb[2]))*SQR(eps) +
	 ((evenc[2] + oddc[2])*SQR(dx)/r[0]+(evend[2] + oddd[2]))*eps +
	 ((evene[2] + odde[2])*SQR(dx)/r[0]+(evenf[2] + oddf[2]))
	 )*dx*SQR(eps)/r[3];

   c1 = (((evena[1] + odda[1])*SQR(dx)/r[0]+(evenb[1] + oddb[1]))*SQR(eps)*eps+
	 ((evenc[1] + oddc[1])*SQR(dx)/r[0]+(evend[1] + oddd[1]))*SQR(eps)+
	 ((evene[1] + odde[1])*SQR(dx)/r[0]+(evenf[1] + oddf[1]))*eps+
	 ((eveng[1] + oddg[1])*SQR(dx)/r[0]+(evenh[1] + oddh[1]))
	 )*dx*eps/r[2];

   c0 = ((evena[0]*SQR(dx)/r[0] + evenb[0])*SQR(eps)*SQR(eps)+
	 (evenc[0]*SQR(dx)/r[0] + evend[0])*SQR(eps)*eps+
	 (evene[0]*SQR(dx)/r[0] + evenf[0])*SQR(eps)+
	 (eveng[0]*SQR(dx)/r[0] + evenh[0])*eps+
	 (eveni[0]*SQR(dx)/r[0] + evenj[0])
	 )*dx/r[1];

   /* Fix the last two lines to make it more systematic. */

   /* Stuff for u_xx=-v_yx or v_yy=-u_xy */

   evena[4] = 24192.0 - 387072.0*t[1] + 645120.0*t[3];
   odda[4]  = 96768.0*t[0] - 258048.0*t[2];

   evenb[4] = -2016.0 +24192.0*t[0] - 32256.0*t[1];
   oddb[4]  = 0.0;


   evena[3] = -53760.0 + 806400.0*t[1] - 1290240.0*t[3];
   odda[3]  = -161280.0*t[0] + 430080.0*t[2];

   evenb[3] = 4480.0 - 53760.0*t[1] + 71680.0*t[3];
   oddb[3]  = 0.0;

   evenc[3] = -13440.0*t[0] + 35840.0*t[2];
   oddc[3]  = 2240.0 - 13440.0*t[1];

   evend[3] = 1120.0*t[0] - 2240.0*t[2];
   oddd[3]  = 0.0;

   
   evena[2] = 41600.0 - 572160.0*t[1] + 860160.0*t[3];
   odda[2]  = 83200.0*t[0] - 215040.0*t[2];

   evenb[2] = -3600.0 + 41280.0*t[1] - 53760.0*t[3];
   oddb[2]  = 0.0;

   evenc[2] = 21760.0*t[0] - 53760.0*t[2];
   oddc[2]  = -2560.0 + 15360.0*t[1];

   evend[2] = -1920.0*t[0] + 3840.0*t[1];
   oddd[2]  = 0.0;

   evene[2] = -320.0 - 1920.0*t[1];
   odde[2]  = -640.0*t[0];

   evenf[2] = 40.0 - 160.0*t[1];
   oddf[2]  = 0.0;


   c4 = ((evena[4] + odda[4])*SQR(dx)/r[0] + (evenb[4] + oddb[4])
	 )*dy*SQR(eps)*SQR(eps)/r[4]/r[0];

   c3 = (((evena[3] + odda[3])*SQR(dx)/r[0] + (evenb[3] + oddb[3]))*eps+
	 ((evenc[3] + oddc[3])*SQR(dx)/r[0] + (evend[3] + oddd[3]))
	 )*SQR(eps)*eps*dy/r[4];

   c2 = (((evena[2] + odda[2])*SQR(dx)/r[0] + (evenb[2] + oddb[2]))*SQR(eps)+
	 ((evenc[2] + oddc[2])*SQR(dx)/r[0] + (evend[2] + oddd[2]))*eps+
	 ((evene[2] + odde[2])*SQR(dx)/r[0] + (evenf[2] + oddf[2]))
	 )*SQR(eps)*dy/r[3];

   c1 = (((evena[1] + odda[1])*SQR(dx)/r[0] + 
	  (evenb[1] + oddb[1]))*SQR(eps)*eps +
	 ((evenc[1] + oddc[1])*SQR(dx)/r[0] + (evend[1] + oddd[1]))*SQR(eps) +
	 ((evene[1] + odde[1])*SQR(dx)/r[0] + (evenf[1] + oddf[1]))*eps +
	 ((eveng[1] + oddg[1])*SQR(dx)/r[0] + (evenh[1] + oddh[1]))
	 )*eps*dy/r[2];

   c0 = ((evena[0]*SQR(dx)/r[0] + evenb[0])*SQR(eps)*SQR(eps) +
	 (evenc[0]*SQR(dx)/r[1] + evend[0]/r[0])*SQR(eps)*eps +
	 (evene[0]*SQR(dx)/r[0] + evenf[0])*SQR(eps) +
	 (eveng[0]*SQR(dx)/r[1] + evenh[0]/r[0])*eps +
	 eveni[0]*SQR(dx)/r[0] + evenj[0]
	 )*dy/r[1];

   /* End of work so far */

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
