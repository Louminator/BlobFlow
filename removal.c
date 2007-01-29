/* REMOVAL.C */
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

/* Clip the domain by removing elements that beyond reasonable bounds. */
void remove_vtx(int i)
{
  int k;

  for (k=i; k<N-1; ++k)
    {
      mblob[k]    = mblob[k+1];
      blobguts[k] = blobguts[k+1];
      tmpparms[k] = tmpparms[k+1];
    }
  N--;
}

void clip(double R)
{
  double R2;
  int k;

  R2 = SQR(R);

  for (k=0; k<N; ++k)
    if (SQR(mblob[k].blob0.x)+SQR(mblob[k].blob0.y) > R2)
      remove_vtx(k);
}
