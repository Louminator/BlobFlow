/* CACHE_RESORT.C */
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

void cache_resort()
{
   int i,j,index,size;
   FineGridLink *trace;
   
   if (N > NMAX/2+1)
     fprintf(comp_log,
	     "Time %12.4e: Cache resorting is not possible.  Problem size is too large.\n",
	     SimTime);
   else
     {
	for (i=0; i<N; ++i)
	  {
	     mblob[i+N] = mblob[i];
	     blobguts[i+N] = blobguts[i];
	  }
	
	size = (int) ldexp(1.0,mplevels);

	index = 0;
	for (i=0; i<size; ++i)
	  for (j=0; j<size; ++j)
	  {
	     trace = FineGridLinks[i+j*size];
	     
	     while (trace != NULL)
	       {
		  mblob[index] = mblob[N+trace->element];
		  tmpparms[index] = tmpparms[N+trace->element];
		  trace->element = index;
		  ++index;
		  trace = trace->next;
	       }
	  }
     }
}
