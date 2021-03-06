/* DEFORM.C */
/* Copyright (c) 2009 Louis F. Rossi                                    *
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

#include "global_min.h"
#include "particle.h"

double dts2(Blob_internal *the_blobguts)
{
  return((visc/2.0)*( (*the_blobguts).a2 + 1.0/(*the_blobguts).a2 ) );
}

double dta2(Blob_internal *the_blobguts,Blob_parms *parms)
{
  return(2.0*( (*parms).du11*( (*parms).cos2 - (*parms).sin2 ) +
	       ( (*parms).du12+(*parms).du21 )*
	       (*parms).sincos )*(*the_blobguts).a2 +
	 (visc/(2.0*(*the_blobguts).s2))*(1.0 - SQR((*the_blobguts).a2)));
}

double dtth(Blob_internal *the_blobguts,Blob_parms *parms,
	    double prefstep,double axisymmtol)
{
  double retval;

  retval = ((*parms).du21 - (*parms).du12)/2.0 +
    ( (((*parms).du21 + (*parms).du12)/2.0)*
      ((*parms).sin2-(*parms).cos2) +
      2.0*(*parms).du11*(*parms).sincos )*
    (1.0/(*the_blobguts).a2 + (*the_blobguts).a2)/
    (1.0/(*the_blobguts).a2 - (*the_blobguts).a2);

  if (fabs((*the_blobguts).a2-1/(*the_blobguts).a2)/prefstep < axisymmtol)
    return(((*parms).du21 - (*parms).du12)/2.0);
  else
    return(retval);
}
