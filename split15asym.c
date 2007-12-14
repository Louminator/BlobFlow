/* SPLIT15.C */
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
 * Louis Rossi                                                          *
 * Department of Mathematical Sciences                                  *
 * University of Delaware                                               *
 * Newark, DE 19715-2553                                                */

#include "global.h"

/* Use a lookup table for asymmetric splitting parameters.*/

/* alpha, rA, rB, gammaA, gammaB, aA2. */
/* Central element: gamma <- 1 - gammaA - gammaB */
/* Pair A: x <- x+-rA*s2, a2 <- a2*aA2, gamma <- gammaA/2 */
/* Pair B: x <- x+-rB*s2, gamma <- gammaB/2 */


void split15asymvels(vels,the_blob,the_blobguts,parms,parm_ptr)
     Vector        *vels;
     Blob_external the_blob;
     Blob_internal the_blobguts;
     Blob_parms     parms;
     double        *parm_ptr;
{
     double     rA,rB;
     double     gammaA,gammaB,aA2;
     Vector     pos;

     rA = (*parm_ptr)*sqrt(the_blobguts.s2*the_blobguts.a2);
     rB = *(parm_ptr+1)*sqrt(the_blobguts.s2/the_blobguts.a2);
     gammaA = *(parm_ptr+2);
     gammaB = *(parm_ptr+3);
     aA2 = *(parm_ptr+4);
   
#ifdef XANTISYMM
     /* Double the number of computational elements by reflecting them
      * about the x-axis. */
     reflect_X();
#endif 

     pos.x = the_blob.x;
     pos.y = the_blob.y;

#ifdef LINEAR
     (*vels) = dpos_vel_gen_linear(pos,the_blobguts,parms);
#else
     (*vels) = dpos_vel_gen(pos,the_blobguts,parms);
#endif

     pos.x = the_blob.x + rA*parms.costh;
     pos.y = the_blob.y + rA*parms.sinth;

#ifdef LINEAR
     *(vels+1) = dpos_vel_gen_linear(pos,the_blobguts,parms);
#else
     *(vels+1) = dpos_vel_gen(pos,the_blobguts,parms);
#endif
	
     pos.x = the_blob.x - rA*parms.costh;
     pos.y = the_blob.y - rA*parms.sinth;

#ifdef LINEAR
     *(vels+2) = dpos_vel_gen_linear(pos,the_blobguts,parms);
#else
     *(vels+2) = dpos_vel_gen(pos,the_blobguts,parms);
#endif

     pos.x = the_blob.x + rB*parms.sinth;
     pos.y = the_blob.y - rB*parms.costh;
	
#ifdef LINEAR
     *(vels+3) = dpos_vel_gen_linear(pos,the_blobguts,parms);
#else
     *(vels+3) = dpos_vel_gen(pos,the_blobguts,parms);
#endif

     pos.x = the_blob.x - rB*parms.sinth;
     pos.y = the_blob.y + rB*parms.costh;

#ifdef LINEAR
     *(vels+4) = dpos_vel_gen_linear(pos,the_blobguts,parms);
#else
     *(vels+4) = dpos_vel_gen(pos,the_blobguts,parms);
#endif

#ifdef XANTISYMM 
     /* re-adjust */
     N /= 2;
#endif 
}

void split15asym(the_blob,the_blobguts,parms,
		 new_blob1,new_blobguts1,new_parms1,
		 new_blob2,new_blobguts2,new_parms2,
		 new_blob3,new_blobguts3,new_parms3,
		 new_blob4,new_blobguts4,new_parms4,parm_ptr)
     Blob_external *the_blob,*new_blob1,*new_blob2,*new_blob3,*new_blob4;
     Blob_internal *the_blobguts,*new_blobguts1,*new_blobguts2,
       *new_blobguts3,*new_blobguts4;
     Blob_parms *parms,*new_parms1,*new_parms2,*new_parms3,*new_parms4;
     double        *parm_ptr;
{
   Blob_external old_blob;
   Blob_internal old_blobguts;
   double     rA,rB;
   double     gammaA,gammaB,aA2;
   
   rA = (*parm_ptr)*sqrt((*the_blobguts).s2*(*the_blobguts).a2);
   rB = *(parm_ptr+1)*sqrt((*the_blobguts).s2/(*the_blobguts).a2);
   gammaA = *(parm_ptr+2);
   gammaB = *(parm_ptr+3);
   aA2 = *(parm_ptr+4);
   
   if (N >= NMAX-5)
     {
	fprintf(comp_log,
		"Overload!  The problem size has grown beyond NMAX.\n");
	stop(-1);
     }
   
   old_blob = *the_blob;
   old_blobguts = *the_blobguts;

   *new_blob1 = *the_blob;
   *new_blob2 = *the_blob;
   *new_blob3 = *the_blob;
   *new_blob4 = *the_blob;

   *new_blobguts1 = *the_blobguts;
   *new_blobguts2 = *the_blobguts;
   *new_blobguts3 = *the_blobguts;
   *new_blobguts4 = *the_blobguts;
   
   (*the_blobguts).s2 = old_blobguts.s2*SQR(alpha);
   (*new_blobguts1).s2 = (*the_blobguts).s2;
   (*new_blobguts2).s2 = (*the_blobguts).s2;
   (*new_blobguts3).s2 = (*the_blobguts).s2;
   (*new_blobguts4).s2 = (*the_blobguts).s2;

   (*new_blobguts1).a2 = (*the_blobguts).a2*aA2;
   (*new_blobguts2).a2 = (*the_blobguts).a2*aA2;
   (*new_blobguts3).a2 = (*the_blobguts).a2;
   (*new_blobguts4).a2 = (*the_blobguts).a2;
   
   (*the_blob).strength = (1-gammaA-gammaB)*old_blob.strength;
   (*new_blob1).strength = 0.5*gammaA*old_blob.strength;
   (*new_blob2).strength = 0.5*gammaA*old_blob.strength;
   (*new_blob3).strength = 0.5*gammaB*old_blob.strength;
   (*new_blob4).strength = 0.5*gammaB*old_blob.strength;
   
   (*the_blob).x = old_blob.x;
   (*the_blob).y = old_blob.y;
	
   (*new_blob1).x = old_blob.x + rA*(*parms).costh;
   (*new_blob1).y = old_blob.y + rA*(*parms).sinth;
	
   (*new_blob2).x = old_blob.x - rA*(*parms).costh;
   (*new_blob2).y = old_blob.y - rA*(*parms).sinth;

   (*new_blob3).x = old_blob.x + rB*(*parms).sinth;
   (*new_blob3).y = old_blob.y - rB*(*parms).costh;
	
   (*new_blob4).x = old_blob.x - rB*(*parms).sinth;
   (*new_blob4).y = old_blob.y + rB*(*parms).costh;

   *new_parms1 = *parms;
   *new_parms2 = *parms;
   *new_parms3 = *parms;
   *new_parms4 = *parms;

   set_blob(the_blobguts,parms);
   set_blob(new_blobguts1,new_parms1);
   set_blob(new_blobguts2,new_parms2);
   set_blob(new_blobguts3,new_parms3);
   set_blob(new_blobguts4,new_parms4);

   (*parms).refinecnt = -1;
   (*new_parms1).refinecnt = -1;
   (*new_parms2).refinecnt = -1;
   (*new_parms3).refinecnt = -1;
   (*new_parms4).refinecnt = -1;
}
