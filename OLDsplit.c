/* SPLIT.C */
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

void split12(the_blob,the_blobguts,parms,new_blob,new_blobguts,new_parms)
blob_external *the_blob,*new_blob;
blob_internal *the_blobguts,*new_blobguts;
blobparms *parms,*new_parms;
{
   blob_external old_blob;
   blob_internal old_blobguts;
   double     r,r2;
   double     thisalpha;
   
   if (N >= NMax)
     {
	fprintf(comp_log,
		"Overload!  The problem size has grown beyond NMax.\n");
	stop(-1);
     }
   
   old_blob = *the_blob;
   *new_blob = *the_blob;
   old_blobguts = *the_blobguts;
   *new_blobguts = *the_blobguts;
   
   thisalpha = alpha;

   (*the_blobguts).s2 = old_blobguts.s2*thisalpha*thisalpha;
   (*new_blobguts).s2 = (*the_blobguts).s2;
   
   (*the_blob).strength = old_blob.strength/2.0;
   (*new_blob).strength = old_blob.strength/2.0;
   
   if (old_blobguts.a2 > 1.0)
     {
	r2 = old_blobguts.s2*old_blobguts.a2*2.0*(1.0-SQR(SQR(thisalpha)));
	r = sqrt(r2);
	(*the_blobguts).a2 = thisalpha*thisalpha*old_blobguts.a2;
	(*new_blobguts).a2 = (*the_blobguts).a2;
	
	(*the_blob).x = old_blob.x + r*(*parms).costh;
	(*the_blob).y = old_blob.y + r*(*parms).sinth;
	
	(*new_blob).x = old_blob.x - r*(*parms).costh;
	(*new_blob).y = old_blob.y - r*(*parms).sinth;
	
     }
   else
     {
	r2 = 2.0*(1.0-SQR(SQR(thisalpha)))*old_blobguts.s2/old_blobguts.a2;
	r = sqrt(r2);
	(*the_blobguts).a2 = old_blobguts.a2/SQR(thisalpha);
	(*new_blobguts).a2 = (*the_blobguts).a2;
	
	(*the_blob).x = old_blob.x + r*(*parms).sinth;
	(*the_blob).y = old_blob.y - r*(*parms).costh;
	
	(*new_blob).x = old_blob.x - r*(*parms).sinth;
	(*new_blob).y = old_blob.y + r*(*parms).costh;
     }
   
   (*the_blob).order = 1;
   (*new_blob).order = 1;
   set_blob(the_blobguts,parms);
   set_blob(new_blobguts,new_parms);
}

void split14_mom2(the_blob,the_blobguts,parms,
	     new_blob1,new_blobguts1,new_parms1,
	     new_blob2,new_blobguts2,new_parms2,
	     new_blob3,new_blobguts3,new_parms3)
blob_external *the_blob,*new_blob1,*new_blob2,*new_blob3;
blob_internal *the_blobguts,*new_blobguts1,*new_blobguts2,*new_blobguts3;
blobparms *parms,*new_parms1,*new_parms2,*new_parms3;
{
   blob_external old_blob;
   blob_internal old_blobguts;
   double     r,r2;
   
   if (N >= NMax-4)
     {
	fprintf(comp_log,
		"Overload!  The problem size has grown beyond NMax.\n");
	stop(-1);
     }
   
   old_blob = *the_blob;
   *new_blob1 = *the_blob;
   *new_blob2 = *the_blob;
   *new_blob3 = *the_blob;
   old_blobguts = *the_blobguts;
   *new_blobguts1 = *the_blobguts;
   *new_blobguts2 = *the_blobguts;
   *new_blobguts3 = *the_blobguts;
   
   (*the_blobguts).s2 = old_blobguts.s2*SQR(alpha);
   (*new_blobguts1).s2 = (*the_blobguts).s2;
   (*new_blobguts2).s2 = (*the_blobguts).s2;
   (*new_blobguts3).s2 = (*the_blobguts).s2;
   
   (*the_blob).strength = old_blob.strength/4.0;
   (*new_blob1).strength = old_blob.strength/4.0;
   (*new_blob2).strength = old_blob.strength/4.0;
   (*new_blob3).strength = old_blob.strength/4.0;
   
   if (old_blobguts.a2 > 1.0)
     {
	r2 = old_blobguts.s2*old_blobguts.a2*4.0*(1.0-SQR(alpha));
	r = sqrt(r2);
	
	(*the_blob).x = old_blob.x + r*(*parms).costh;
	(*the_blob).y = old_blob.y + r*(*parms).sinth;
	
	(*new_blob1).x = old_blob.x - r*(*parms).costh;
	(*new_blob1).y = old_blob.y - r*(*parms).sinth;

	r /= old_blobguts.a2;
	
	(*new_blob2).x = old_blob.x + r*(*parms).sinth;
	(*new_blob2).y = old_blob.y - r*(*parms).costh;
	
	(*new_blob3).x = old_blob.x - r*(*parms).sinth;
	(*new_blob3).y = old_blob.y + r*(*parms).costh;
     }
   else
     {
	r2 = old_blobguts.s2/old_blobguts.a2*4.0*(1.0-SQR(alpha));
	r = sqrt(r2);
	
	(*the_blob).x = old_blob.x + r*(*parms).costh;
	(*the_blob).y = old_blob.y + r*(*parms).sinth;
	
	(*new_blob1).x = old_blob.x - r*(*parms).costh;
	(*new_blob1).y = old_blob.y - r*(*parms).sinth;

	r *= old_blobguts.a2;
	
	(*new_blob2).x = old_blob.x + r*(*parms).sinth;
	(*new_blob2).y = old_blob.y - r*(*parms).costh;
	
	(*new_blob3).x = old_blob.x - r*(*parms).sinth;
	(*new_blob3).y = old_blob.y + r*(*parms).costh;
     }
   
   (*the_blob).order = 1;
   (*new_blob1).order = 1;
   (*new_blob2).order = 1;
   (*new_blob3).order = 1;
   set_blob(the_blobguts,parms);
   set_blob(new_blobguts1,new_parms1);
   set_blob(new_blobguts2,new_parms2);
   set_blob(new_blobguts3,new_parms3);
}

void splitvels(vels,the_blob,the_blobguts,parms)
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

     r2 = the_blobguts.s2*the_blobguts.a2*(1.0-alpha)*(9.7052*alpha-1.7235);
     r = sqrt(r2);
	
     pos.x = the_blob.x + r*parms.costh;
     pos.y = the_blob.y + r*parms.sinth;

     /*(*vels) = dpos_vel_gen(pos);*/
     (*vels) = dpos_vel_gen_exact_pi(pos,the_blobguts,parms);
	
     pos.x = the_blob.x - r*parms.costh;
     pos.y = the_blob.y - r*parms.sinth;

     /* *(vels+1) = dpos_vel_gen(pos); */
     *(vels+1) = dpos_vel_gen_exact_pi(pos,the_blobguts,parms);

     r /= the_blobguts.a2;
	
     pos.x = the_blob.x + r*parms.sinth;
     pos.y = the_blob.y - r*parms.costh;
	
     /* *(vels+2) = dpos_vel_gen(pos); */
     *(vels+2) = dpos_vel_gen_exact_pi(pos,the_blobguts,parms);

     pos.x = the_blob.x - r*parms.sinth;
     pos.y = the_blob.y + r*parms.costh;

     /* *(vels+3) = dpos_vel_gen(pos); */
     *(vels+3) = dpos_vel_gen_exact_pi(pos,the_blobguts,parms);

#ifdef XANTISYMM 
     /* re-adjust */
     N /= 2;
#endif 
}

void split14(the_blob,the_blobguts,parms,
	     new_blob1,new_blobguts1,new_parms1,
	     new_blob2,new_blobguts2,new_parms2,
	     new_blob3,new_blobguts3,new_parms3)
     blob_external *the_blob,*new_blob1,*new_blob2,*new_blob3;
blob_internal *the_blobguts,*new_blobguts1,*new_blobguts2,*new_blobguts3;
blobparms *parms,*new_parms1,*new_parms2,*new_parms3;
{
   blob_external old_blob;
   blob_internal old_blobguts;
   double     r,r2;
   
   if (N >= NMax-4)
     {
	fprintf(comp_log,
		"Overload!  The problem size has grown beyond NMax.\n");
	stop(-1);
     }
   
   old_blob = *the_blob;
   *new_blob1 = *the_blob;
   *new_blob2 = *the_blob;
   *new_blob3 = *the_blob;
   old_blobguts = *the_blobguts;
   *new_blobguts1 = *the_blobguts;
   *new_blobguts2 = *the_blobguts;
   *new_blobguts3 = *the_blobguts;
   
   (*the_blobguts).s2 = old_blobguts.s2*SQR(alpha);
   (*new_blobguts1).s2 = (*the_blobguts).s2;
   (*new_blobguts2).s2 = (*the_blobguts).s2;
   (*new_blobguts3).s2 = (*the_blobguts).s2;
   
   (*the_blob).strength = old_blob.strength/4.0;
   (*new_blob1).strength = old_blob.strength/4.0;
   (*new_blob2).strength = old_blob.strength/4.0;
   (*new_blob3).strength = old_blob.strength/4.0;
   
   r2 = old_blobguts.s2*old_blobguts.a2*(1.0-alpha)*(9.7052*alpha-1.7235);
   r = sqrt(r2);
	
   (*the_blob).x = old_blob.x + r*(*parms).costh;
   (*the_blob).y = old_blob.y + r*(*parms).sinth;
	
   (*new_blob1).x = old_blob.x - r*(*parms).costh;
   (*new_blob1).y = old_blob.y - r*(*parms).sinth;

   r /= old_blobguts.a2;
	
   (*new_blob2).x = old_blob.x + r*(*parms).sinth;
   (*new_blob2).y = old_blob.y - r*(*parms).costh;
	
   (*new_blob3).x = old_blob.x - r*(*parms).sinth;
   (*new_blob3).y = old_blob.y + r*(*parms).costh;

   set_blob(the_blobguts,parms);
   set_blob(new_blobguts1,new_parms1);
   set_blob(new_blobguts2,new_parms2);
   set_blob(new_blobguts3,new_parms3);

   (*parms).refinecnt = -1;
   (*new_parms1).refinecnt = -1;
   (*new_parms2).refinecnt = -1;
   (*new_parms3).refinecnt = -1;
}

