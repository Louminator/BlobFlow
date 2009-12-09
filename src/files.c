/* FILES.C */
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

#include "global_min.h"
#include "particle.h"
#include "biot-savart.h"

void write_vorts(int frameno)
{
  int         i;
  char        vortex_name[FILENAME_LEN];
  FILE        *vortex_file,*fopen();
   
  sprintf(vortex_name,"%s%04d%s",filename,frameno,".vtx");

  vortex_file = fopen(vortex_name,"w");

  for (i=0; i<N; ++i)
      fprintf(vortex_file,"%10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n",
	      mblob[i].blob0.x,mblob[i].blob0.y,
	      mblob[i].blob0.strength,blobguts[i].s2,
	      blobguts[i].a2,
	      blobguts[i].th);

  fclose(vortex_file);
}

void write_partition(int frameno)
{
   int i,l;
   char grid_name[FILENAME_LEN];
   FILE *gridfile,*fopen();
   
   sprintf(grid_name,"%s%04d%s",filename,frameno,".mp");
   
   gridfile = fopen(grid_name,"w");
   
   fprintf(gridfile,"%12.4e %12.4e\n",minX,minY);
   fprintf(gridfile,"%12.4e %12.4e\n",maxX,minY);
   fprintf(gridfile,"%12.4e %12.4e\n",maxX,maxY);
   fprintf(gridfile,"%12.4e %12.4e\n",minX,maxY);
   fprintf(gridfile,"%12.4e %12.4e\n",minX,minY);
   fprintf(gridfile,"\n");

   for (l=0; l<mplevels; ++l)
     {
	for (i=1; i< (int) ldexp(1.0,l+1); ++i)
	  {
	     fprintf(gridfile,"%12.4e %12.4e\n",
		     minX+distX*i/ldexp(1.0,l+1),minY);
	     fprintf(gridfile,"%12.4e %12.4e\n",
		     minX+distX*i/ldexp(1.0,l+1),maxY);
	     fprintf(gridfile,"%12.4e %12.4e\n",
		     minX+distX*i/ldexp(1.0,l+1),minY);
	  }
	fprintf(gridfile,"\n");
	
	for (i=1; i< (int) ldexp(1.0,l+1); ++i)
	  {
	     fprintf(gridfile,"%12.4e %12.4e\n",
		     minX,minY+distY*i/ldexp(1.0,l+1));
	     fprintf(gridfile,"%12.4e %12.4e\n",
		     maxX,minY+distY*i/ldexp(1.0,l+1));
	     fprintf(gridfile,"%12.4e %12.4e\n",
		     minX,minY+distY*i/ldexp(1.0,l+1));
	  }
	fprintf(gridfile,"\n");
     }
   fclose(gridfile);
}
