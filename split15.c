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

void split15vels(vels,the_blob,the_blobguts,parms)
     vector        *vels;
     blob_external the_blob;
     blob_internal the_blobguts;
     blobparms     parms;
{
     double     r,r2;
     vector     pos;
   
#ifdef XANTISYMM
     /* Double the number of computational elements by reflecting them
      * about the x-axis. */
     reflect_X();
#endif 


     r2 = the_blobguts.s2*the_blobguts.a2*(1-alpha)*
       SQR(4.9034-5.9821*(1-alpha)+0.040012698*(1-alpha)*(1-alpha)+
	   3.650027*(1-alpha)*(1-alpha)*(1-alpha));
     r = sqrt(r2);
	
     pos.x = the_blob.x;
     pos.y = the_blob.y;

     /*(*vels) = dpos_vel_gen(pos);*/
     (*vels) = dpos_vel_gen_linear(pos,the_blobguts,parms);
	
     pos.x = the_blob.x + r*parms.costh;
     pos.y = the_blob.y + r*parms.sinth;

     /*(*vels+1) = dpos_vel_gen(pos);*/
     *(vels+1) = dpos_vel_gen_linear(pos,the_blobguts,parms);
	
     pos.x = the_blob.x - r*parms.costh;
     pos.y = the_blob.y - r*parms.sinth;

     /* *(vels+2) = dpos_vel_gen(pos); */
     *(vels+2) = dpos_vel_gen_linear(pos,the_blobguts,parms);

     r /= the_blobguts.a2;
	
     pos.x = the_blob.x + r*parms.sinth;
     pos.y = the_blob.y - r*parms.costh;
	
     /* *(vels+3) = dpos_vel_gen(pos); */
     *(vels+3) = dpos_vel_gen_linear(pos,the_blobguts,parms);

     pos.x = the_blob.x - r*parms.sinth;
     pos.y = the_blob.y + r*parms.costh;

     /* *(vels+4) = dpos_vel_gen(pos); */
     *(vels+4) = dpos_vel_gen_linear(pos,the_blobguts,parms);

#ifdef XANTISYMM 
     /* re-adjust */
     N /= 2;
#endif 
}

void split15(the_blob,the_blobguts,parms,
	     new_blob1,new_blobguts1,new_parms1,
	     new_blob2,new_blobguts2,new_parms2,
	     new_blob3,new_blobguts3,new_parms3,
	     new_blob4,new_blobguts4,new_parms4)
     blob_external *the_blob,*new_blob1,*new_blob2,*new_blob3,*new_blob4;
     blob_internal *the_blobguts,*new_blobguts1,*new_blobguts2,
       *new_blobguts3,*new_blobguts4;
     blobparms *parms,*new_parms1,*new_parms2,*new_parms3,*new_parms4;
{
   blob_external old_blob;
   blob_internal old_blobguts;
   double     r,r2,a2,amp;
   
   if (N >= NMax-5)
     {
	fprintf(comp_log,
		"Overload!  The problem size has grown beyond NMax.\n");
	stop(-1);
     }
   
   old_blob = *the_blob;
   old_blobguts = *the_blobguts;

   r2 = old_blobguts.s2*old_blobguts.a2*(1-alpha)*
     SQR(4.9034-5.9821*(1-alpha)+0.040012698*(1-alpha)*(1-alpha)+
	 3.650027*(1-alpha)*(1-alpha)*(1-alpha));
   r = sqrt(r2);

   amp = 0.3337 + 0.57049*(1-alpha)+
     0.87065*(1-alpha)*(1-alpha)-1.2181*(1-alpha)*(1-alpha)*(1-alpha);

   a2 = 1 + (1-alpha)*
     (-1.0993-0.30263*(1-alpha)+3.8644*(1-alpha)*(1-alpha)-
      4.0458*(1-alpha)*(1-alpha)*(1-alpha));
	
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

   (*new_blobguts1).a2 = (*the_blobguts).a2*a2;
   (*new_blobguts2).a2 = (*the_blobguts).a2*a2;
   (*new_blobguts3).a2 = (*the_blobguts).a2/a2;
   (*new_blobguts4).a2 = (*the_blobguts).a2/a2;
   
   (*the_blob).strength = (1-amp)*old_blob.strength;
   (*new_blob1).strength = amp*old_blob.strength/4.0;
   (*new_blob2).strength = amp*old_blob.strength/4.0;
   (*new_blob3).strength = amp*old_blob.strength/4.0;
   (*new_blob4).strength = amp*old_blob.strength/4.0;
   
   (*the_blob).x = old_blob.x;
   (*the_blob).y = old_blob.y;
	
   (*new_blob1).x = old_blob.x + r*(*parms).costh;
   (*new_blob1).y = old_blob.y + r*(*parms).sinth;
	
   (*new_blob2).x = old_blob.x - r*(*parms).costh;
   (*new_blob2).y = old_blob.y - r*(*parms).sinth;

   r /= old_blobguts.a2;
	
   (*new_blob3).x = old_blob.x + r*(*parms).sinth;
   (*new_blob3).y = old_blob.y - r*(*parms).costh;
	
   (*new_blob4).x = old_blob.x - r*(*parms).sinth;
   (*new_blob4).y = old_blob.y + r*(*parms).costh;

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
