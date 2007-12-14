/* MPCOEFFS.C */
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
 * Newark, DE 19716                                                     */

#include "global.h"

#define R2 256.0

void Calc_Coeffs(int vort, complex z, complex mp[])
{
   double Gamma,a2,s2,dx,dy,tmp1;
   complex ctmp1,ctmp2,bdy_corr;
   int    i,j;
   
   Gamma = mblob[vort].blob0.strength;
   a2 = blobguts[vort].a2;
   s2 = blobguts[vort].s2;
   
   dx = z.re;
   dy = z.im;

   mp[0].re = 0.0;
   mp[0].im = -Gamma*(1.0-exp(-R2));
    
   mp[1].re = -(Gamma)*dy*(1.0-exp(-R2));
   mp[1].im = (Gamma)*dx*(1.0-exp(-R2));

   for (i=2; i<PMAX; ++i)
     {
	/* Calculate the boundary truncation correction. */
	bdy_corr.re = 0.0;
	bdy_corr.im = 0.0;
	for (j=0; j<(i/2)-1; ++j)
	  {
	     tmp1 = C(i-1,j)*C(i-1-j,j+1)*pow(0.5,2.0*j+1.0)*
	       pow(blobguts[vort].a2-1.0/blobguts[vort].a2,j+1.0);
	     ctmp1.re = -dx;
	     ctmp1.im = -dy;
	     ctmp1 = cpowi(ctmp1,i-2*j-2);
	     ctmp2.re = cos((2.0*j+2.0)*blobguts[vort].th);
	     ctmp2.im = sin((2.0*j+2.0)*blobguts[vort].th);
	     ctmp1 = cmult(ctmp1,ctmp2);
	     ctmp1.re *= tmp1;
	     ctmp1.im *= tmp1;
	     
	     bdy_corr.re += ctmp1.re;
	     bdy_corr.im += ctmp1.im;
	  }

	mp[i].re = -dx*mp[i-1].re + dy*mp[i-1].im +
	  2.0*s2*(i-1)*(a2-1.0/a2)*
	  ( (tmpparms[vort].cos2-tmpparms[vort].sin2)*mp[i-2].re -
	   2.0*tmpparms[vort].sincos*mp[i-2].im) -
	  Gamma*exp(-R2)*bdy_corr.im;
	mp[i].im = -dx*mp[i-1].im - dy*mp[i-1].re +
	  2.0*s2*(i-1)*(a2-1.0/a2)*
	  ( (tmpparms[vort].cos2-tmpparms[vort].sin2)*mp[i-2].im +
	   2.0*tmpparms[vort].sincos*mp[i-2].re) +
	  Gamma*exp(-R2)*bdy_corr.re;
     }
}

