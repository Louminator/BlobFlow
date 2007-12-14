/* REFLECT.C */
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

#include "global.h"

void reflect_X()
{
   int i;
   
   if (2*N > NMax)
     {
	fprintf(comp_log,"The problem size is too big!\n");
	fprintf(comp_log,"Note: If you use XANTISYMMETRY, NMAX must be adjusted to buffer twice\n");
	fprintf(comp_log,"the size of the half-problem.\n");

	stop(-1);
     }
   
   for (i=0; i<N; ++i)
     {
	mblob[i+N] = mblob[i];
	blobguts[i+N] = blobguts[i];

	mblob[i+N].blob0.y *= -1.0;
	mblob[i+N].blob0.strength *= -1.0;
	blobguts[i+N].th *= -1.0;
	set_blob(&(blobguts[i+N]),&(tmpparms[i+N]));
     }
   
   N *= 2;
}
