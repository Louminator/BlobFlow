/* VELOCITY.C */
/* Copyright (c) 2000, 2004 Louis F. Rossi                                    *
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

void induced_v(the_blob,the_blobguts,parms,tmpdx,tmpdy,result)
blob_external *the_blob;
blob_internal *the_blobguts;
blobparms *parms;
double tmpdx,tmpdy;
double result[9];

{
  double r[maxexp],t[maxexp];

  double psi_RT[maxpolyn][maxexp];
  double psi_x_RT[maxpolyn][maxexp],psi_y_RT[maxpolyn][maxexp];
  double psi_xx_RT[maxpolyn][maxexp],psi_yy_RT[maxpolyn][maxexp];
  double psi_xy_RT[maxpolyn][maxexp];
  double psi_xxx_RT[maxpolyn][maxexp],psi_xxy_RT[maxpolyn][maxexp];
  double psi_xyy_RT[maxpolyn][maxexp],psi_yyy_RT[maxpolyn][maxexp];

  double psi_x,psi_y,psi_xx,psi_yy,psi_xy;
  double psi_xxx,psi_xxy,psi_xyy,psi_yyy;

  double psi_c[maxpolyn],psi_x_c[maxpolyn],psi_y_c[maxpolyn];
  double psi_xx_c[maxpolyn],psi_xy_c[maxpolyn],psi_yy_c[maxpolyn];
  double psi_xxx_c[maxpolyn],psi_xxy_c[maxpolyn];
  double psi_xyy_c[maxpolyn],psi_yyy_c[maxpolyn];

  double eps,str,a2,s2,s4,dx,dy;

  /* Change bases to the local axes. */
   
  dx =  (*parms).costh*tmpdx+(*parms).sinth*tmpdy;
  dy = -(*parms).sinth*tmpdx+(*parms).costh*tmpdy;
   
  eps = (sqrt((*the_blobguts).a2)-1.0)/(sqrt((*the_blobguts).a2)+1.0);
   
  s2 = (*the_blobguts).s2;
  s4 = SQR(s2);

  a2  = (*the_blobguts).a2;
  str = (*the_blob).strength;

  build_rt(dx,dy,eps,r,t);

  build_psi(dx,dy,eps,r,t,psi_c,psi_RT);

  psi_x = build_psi_x(dx,dy,str,s2,s4,a2,eps,r,t,
		      psi_c,psi_x_c,psi_RT,psi_x_RT);

  psi_y = build_psi_y(dx,dy,str,s2,s4,a2,eps,r,t,
		      psi_c,psi_y_c,psi_RT,psi_y_RT);
  
  psi_xx = build_psi_xx(dx,dy,str,s2,s4,a2,eps,r,t,
			psi_x_c,psi_xx_c,psi_RT,psi_xx_RT);
  
  psi_yy = build_psi_yy(dx,dy,str,s2,s4,a2,eps,r,t,
			psi_y_c,psi_yy_c,psi_RT,psi_yy_RT);
  
  psi_xy = build_psi_xy(dx,dy,str,s2,s4,a2,eps,r,t,
			psi_y_c,psi_xy_c,psi_RT,psi_xy_RT);
  
  psi_xxx = build_psi_xxx(dx,dy,str,s2,s4,a2,eps,r,t,
			  psi_xx_c,psi_xxx_c,psi_xx_RT,psi_xxx_RT);
  
  psi_xxy = build_psi_xxy(dx,dy,str,s2,s4,a2,eps,r,t,
			  psi_xx_c,psi_xxy_c,psi_xx_RT,psi_xxy_RT);
  
  psi_xyy = build_psi_xyy(dx,dy,str,s2,s4,a2,eps,r,t,
			  psi_yy_c,psi_xyy_c,psi_yy_RT,psi_xyy_RT);
  
  psi_yyy = build_psi_yyy(dx,dy,str,s2,s4,a2,eps,r,t,
			  psi_yy_c,psi_yyy_c,psi_yy_RT,psi_yyy_RT);
  

  /* Rotate everything back into the standard reference frame. */
   
  /* u = -psi_y */
  result[0] = (*parms).costh*(-psi_y)-(*parms).sinth*psi_x;
  /* v = psi_x */
  result[1] = (*parms).sinth*(-psi_y)+(*parms).costh*psi_x;

  /* u_x = -psi_xy */
  result[2] = (-psi_xy*((*parms).cos2-(*parms).sin2)-
	       (-psi_yy+psi_xx)*(*parms).sincos);
   
  /* u_y = -psi_yy */
  result[3] = -((*parms).sin2*psi_xx+
		2.0*(*parms).sincos*psi_xy+
		(*parms).cos2*psi_yy);
   
  /* v_x = psi_xx */
  result[4] = (2.0*(-psi_xy)*(*parms).sincos-
	       (-psi_yy)*(*parms).sin2+psi_xx*(*parms).cos2);

  /* u_xx = -psi_xxy */
  result[5] = -( (*parms).cos2*(*parms).sinth*(psi_xxx-2.0*psi_xyy)+
		 (*parms).cos2*(*parms).costh*psi_xxy+
		 (*parms).sin2*(*parms).costh*(psi_yyy-2.0*psi_xxy)+
		 (*parms).sin2*(*parms).sinth*psi_xyy );

  /* u_xy = -psi_xyy */
  result[6] = -( (*parms).sincos*(*parms).sinth*(psi_xxx-2.0*psi_xyy)+
		 (*parms).cos2*(*parms).costh*psi_xyy-
		 (*parms).sin2*(*parms).sinth*psi_xxy+
		 (*parms).sincos*(*parms).costh*(-psi_yyy+2.0*psi_xxy) );

  /* u_yy = -psi_yyy */
  result[7] = -( (*parms).sin2*(*parms).sinth*psi_xxx+
		 3.0*(*parms).sin2*(*parms).costh*psi_xxy+
		 3.0*(*parms).sincos*(*parms).costh*psi_xyy+
		 (*parms).cos2*(*parms).costh*psi_yyy );

  /* v_xx = psi_xxx */
  result[8] = ( (*parms).cos2*(*parms).costh*psi_xxx -
		3.0*(*parms).cos2*(*parms).sinth*psi_xxy +
		3.0*(*parms).sin2*(*parms).costh*psi_xyy -
		(*parms).sin2*(*parms).sinth*psi_yyy );
}
