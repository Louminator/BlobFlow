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
#include "global_min.h"
#include "particle.h"
#include "biot-savart.h"

void induced_v_asympt(Blob_external *the_blob,Blob_internal *the_blobguts,
		      Blob_parms *parms,double tmpdx,double tmpdy,
		      double result[9])
{
  double r[MAXEXP],t[MAXEXP];

  double psi_RT[MAXPOLYN][MAXEXP];
  double psi_x_RT[MAXPOLYN][MAXEXP],psi_y_RT[MAXPOLYN][MAXEXP];
  double psi_xx_RT[MAXPOLYN][MAXEXP],psi_yy_RT[MAXPOLYN][MAXEXP];
  double psi_xy_RT[MAXPOLYN][MAXEXP];
  double psi_xxx_RT[MAXPOLYN][MAXEXP],psi_xxy_RT[MAXPOLYN][MAXEXP];
  double psi_xyy_RT[MAXPOLYN][MAXEXP],psi_yyy_RT[MAXPOLYN][MAXEXP];

  double psi_x,psi_y,psi_xx,psi_yy,psi_xy;
  double psi_xxx,psi_xxy,psi_xyy,psi_yyy;

  double psi_c[MAXPOLYN],psi_x_c[MAXPOLYN],psi_y_c[MAXPOLYN];
  double psi_xx_c[MAXPOLYN],psi_xy_c[MAXPOLYN],psi_yy_c[MAXPOLYN];
  double psi_xxx_c[MAXPOLYN],psi_xxy_c[MAXPOLYN];
  double psi_xyy_c[MAXPOLYN],psi_yyy_c[MAXPOLYN];

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
  
#ifdef CORRECTVEL4
  psi_xxx = build_psi_xxx(dx,dy,str,s2,s4,a2,eps,r,t,
			  psi_xx_c,psi_xxx_c,psi_xx_RT,psi_xxx_RT);
  
  psi_xxy = build_psi_xxy(dx,dy,str,s2,s4,a2,eps,r,t,
			  psi_xx_c,psi_xxy_c,psi_xx_RT,psi_xxy_RT);
  
  psi_xyy = build_psi_xyy(dx,dy,str,s2,s4,a2,eps,r,t,
			  psi_yy_c,psi_xyy_c,psi_yy_RT,psi_xyy_RT);
  
  psi_yyy = build_psi_yyy(dx,dy,str,s2,s4,a2,eps,r,t,
			  psi_yy_c,psi_yyy_c,psi_yy_RT,psi_yyy_RT);
#endif

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

#ifdef CORRECTVEL4
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
		 3.0*(*parms).sinth*(*parms).cos2*psi_xyy+
		 (*parms).cos2*(*parms).costh*psi_yyy );

  /* v_xx = psi_xxx */
  result[8] = ( (*parms).cos2*(*parms).costh*psi_xxx -
		3.0*(*parms).cos2*(*parms).sinth*psi_xxy +
		3.0*(*parms).sin2*(*parms).costh*psi_xyy -
		(*parms).sin2*(*parms).sinth*psi_yyy );
#endif
}

/* Rodrigo's spectral interp code. */

void induced_v_platte(Blob_external *the_blob,Blob_internal *the_blobguts,
		      Blob_parms *parms,double tmpdx,double tmpdy,
		      double result[9])
{
  double str,a2,s2,a,dx,dy,phi[9];

  s2  = (*the_blobguts).s2;
  a2  = (*the_blobguts).a2;
  a   = sqrt(a2);
  str = (*the_blob).strength*2*M_PI;

  dx = ( (*parms).costh*tmpdx+(*parms).sinth*tmpdy)/sqrt(s2);
  dy = (-(*parms).sinth*tmpdx+(*parms).costh*tmpdy)/sqrt(s2);

  eval_biot(a, dx, dy, phi);

  phi[0] *= str/sqrt(s2);
  phi[1] *= str/sqrt(s2);
  phi[2] *= str/s2;
  phi[3] *= str/s2;
  phi[4] *= str/s2;
#ifdef CORRECTVEL4
  phi[5] *= str/s2/sqrt(s2);
  phi[6] *= str/s2/sqrt(s2);
  phi[7] *= str/s2/sqrt(s2);
  phi[8] *= str/s2/sqrt(s2);
#endif

  /* Rotate everything back into the standard reference frame. */

  /* u = -psi_y */
  result[0] = (*parms).costh*(-phi[1])-(*parms).sinth*phi[0];
  /* v = psi_x */
  result[1] = (*parms).sinth*(-phi[1])+(*parms).costh*phi[0];

  /* u_x = -psi_xy */
  result[2] = (-phi[2]*((*parms).cos2-(*parms).sin2)-
	       (-phi[4]+phi[3])*(*parms).sincos);
   
  /* u_y = -psi_yy */
  result[3] = -((*parms).sin2*phi[3]+
		2.0*(*parms).sincos*phi[2]+
		(*parms).cos2*phi[4]);
   
  /* v_x = psi_xx */
  result[4] = (2.0*(-phi[2])*(*parms).sincos-
	       (-phi[4])*(*parms).sin2+phi[3]*(*parms).cos2);

#ifdef CORRECTVEL4
  /* u_xx = -psi_xxy */
  result[5] = -( (*parms).cos2*(*parms).sinth*(phi[5]-2.0*phi[8])+
		 (*parms).cos2*(*parms).costh*phi[7]+
		 (*parms).sin2*(*parms).costh*(phi[6]-2.0*phi[7])+
		 (*parms).sin2*(*parms).sinth*phi[8] );

  /* u_xy = -psi_xyy */
  result[6] = -( (*parms).sincos*(*parms).sinth*(phi[5]-2.0*phi[8])+
		 (*parms).cos2*(*parms).costh*phi[8]-
		 (*parms).sin2*(*parms).sinth*phi[7]+
		 (*parms).sincos*(*parms).costh*(-phi[6]+2.0*phi[7]) );

  /* u_yy = -psi_yyy */
  result[7] = -( (*parms).sin2*(*parms).sinth*phi[5]+
		 3.0*(*parms).sin2*(*parms).costh*phi[7]+
		 3.0*(*parms).sinth*(*parms).cos2*phi[8]+
		 (*parms).cos2*(*parms).costh*phi[6] );

  /* v_xx = psi_xxx */
  result[8] = ( (*parms).cos2*(*parms).costh*phi[5] -
		3.0*(*parms).cos2*(*parms).sinth*phi[7] +
		3.0*(*parms).sin2*(*parms).costh*phi[8] -
		(*parms).sin2*(*parms).sinth*phi[6] );

#endif
}

void induced_v(Blob_external *the_blob,Blob_internal *the_blobguts,
	       Blob_parms *parms,double tmpdx,double tmpdy,
	       double result[9])
{
  /* induced_v_asympt(the_blob,the_blobguts,parms,tmpdx,tmpdy,result); */
  induced_v_platte(the_blob,the_blobguts,parms,tmpdx,tmpdy,result);
}

