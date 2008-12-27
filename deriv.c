/* DERIV.C */
/* Copyright (c) 2001 Louis F. Rossi                                    *
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

#include "global.h"

#ifdef XANTISYMM
#define CARDINALITY  N/2
#else
#define CARDINALITY  N
#endif

#ifdef MULTIPROC
#include "multiproc.h"
#endif

void wipe_vort_vel_field()
{
  int i;

  for (i=0; i<N; ++i)
    {
      mblob[i].blob0.dx = 0.0;
      mblob[i].blob0.dy = 0.0;

      tmpparms[i].du11 = 0.0;
      tmpparms[i].du12 = 0.0;
      tmpparms[i].du21 = 0.0;
   
      tmpparms[i].u_xx = 0.0;
      tmpparms[i].u_xy = 0.0;
      tmpparms[i].u_yy = 0.0;
      tmpparms[i].v_xx = 0.0;
    }  
}

void vel_field()
{
   int j;

   for (j=0; j<N; ++j)
     set_blob(&(blobguts[j]),&(tmpparms[j]));
   
#ifdef XANTISYMM
   /* Double the number of computational elements by reflecting them
    * about the x-axis. */
   reflect_X();
#endif 

   wipe_vort_vel_field();

   Create_Hierarchy();

   
#ifdef CACHERESORT
	cache_resort();
#endif

#ifdef MULTIPROC  
   if (total_processes == 2)
     peer();
   else
     {
	if (rank == 0) 
	  master();
	else 
	  slave();
	
	finish();
     }
#else
   /* single processor code */

   for (j=0; j<CARDINALITY; ++j) 
     { 
#ifdef LINEAR
       dpos_vel_linear(j);
#else
#ifdef NOFASTMP
	dpos_vel(j);
#else
	dpos_vel_fast(j);
#endif
#endif
     }
#endif  
   
#ifdef XANTISYMM 
   /* re-adjust */
   N /= 2;
#endif   

#ifdef CORRECTVEL4
   correct_vel_4();
#endif
   Release_Links(mplevels);

   if (B != 0)
     {
       solve_bdy_matrix(walls,Bpiv,BdyMat);
       for (j=0; j<N; ++j)
	 bdy_vel(mblob[j].blob0.x,mblob[j].blob0.y,
		 &(mblob[j].blob0.dx),&(mblob[j].blob0.dy));
     }
}

void dpos_vel(vort)
int vort;
{
   int i;
    
   mblob[vort].blob0.dx = 0.0;
   mblob[vort].blob0.dy = 0.0;

   tmpparms[vort].du11 = 0.0;
   tmpparms[vort].du12 = 0.0;
   tmpparms[vort].du21 = 0.0;
   
   tmpparms[vort].u_xx = 0.0;
   tmpparms[vort].u_xy = 0.0;
   tmpparms[vort].u_yy = 0.0;
   tmpparms[vort].v_xx = 0.0;
   
   for (i=0; i<N; ++i)
     vort_vort_interaction(vort,i);

   potential_flow(vort);
}

void dpos_vel_fast(vort)
{
   /* MP_Sum(vort,mplevels-2); */

   mblob[vort].blob0.dx = 0.0;
   mblob[vort].blob0.dy = 0.0;

   tmpparms[vort].du11 = 0.0;
   tmpparms[vort].du12 = 0.0;
   tmpparms[vort].du21 = 0.0;
   
   tmpparms[vort].u_xx = 0.0;
   tmpparms[vort].u_xy = 0.0;
   tmpparms[vort].u_yy = 0.0;
   tmpparms[vort].v_xx = 0.0;
   
   velsum_cputime_ref = clock();
   MP_Sum(vort,mplevels);

   velsum_cputime += clock()-velsum_cputime_ref;
   
   /* Uses a special finer grid for merging */
   veldirect_cputime_ref = clock();
   MP_Direct3(vort,mplevels);
   veldirect_cputime += clock()-veldirect_cputime_ref;

   potential_flow(vort);
}

Vector dpos_vel_gen(pos,the_blobguts,parms)
     Vector pos;
     Blob_internal the_blobguts;
     Blob_parms parms;
{
   int i;
   double dx,dy;
   double result[9];
   Vector sum;
   
   sum.x = 0.0;
   sum.y = 0.0;

   for (i=0; i<N; ++i)
     {
	dx = pos.x-mblob[i].blob0.x;
	dy = pos.y-mblob[i].blob0.y;

	if ( (dx != 0.0) || (dy != 0.0) )
	  {
	     induced_v(&(mblob[i].blob0),&(blobguts[i]),
		       &(tmpparms[i]),dx,dy,result);

	     sum.x += result[0];
	     sum.y += result[1];
	  }
     }

   return(sum);
}

Vector vel_gen(x,y)
     double x,y;
{
  Vector v,w;
  double dx,dy;
  double result[9];
  int i;

  v.x = 0.0;
  v.y = 0.0;

  for (i=0; i<N; ++i)
    {
      dx = x-mblob[i].blob0.x;
      dy = y-mblob[i].blob0.y;
	       
      if ( (dx != 0.0) || (dy != 0.0) )
	{
	  induced_v(&(mblob[i].blob0),
		    &(blobguts[i]),
		    &(tmpparms[i]),dx,dy,
		    result);

	  v.x += result[0];
	  v.y += result[1];
	}
    }

  w = potential_vel(x,y);

  v.x += w.x;
  v.y += w.y;

  return(v);
}

void correct_vel_4()
{
  int i;
  double termA, termB, termC;

  for (i=0; i<N; ++i)
    {
      termA = (tmpparms[i].cos2*blobguts[i].a2+
	       tmpparms[i].sin2/blobguts[i].a2);
      
      termB = (tmpparms[i].sin2*blobguts[i].a2+
	       tmpparms[i].cos2/blobguts[i].a2);
      
      termC = tmpparms[i].sincos*(blobguts[i].a2-1.0/blobguts[i].a2);
      
      mblob[i].blob0.dx += blobguts[i].s2*
	(tmpparms[i].u_xx*termA + 
	 tmpparms[i].u_yy*termB + 
	 2.0*tmpparms[i].u_xy*termC);
	
      mblob[i].blob0.dy += blobguts[i].s2*
	(tmpparms[i].v_xx*termA + 
	 (-tmpparms[i].u_xy)*termB + 
	 2.0*(-tmpparms[i].u_xx)*termC);
    }
}

double dts2(the_blobguts)
     Blob_internal *the_blobguts;
{
  return((visc/2.0)*( (*the_blobguts).a2 + 1.0/(*the_blobguts).a2 ) );
}

double dta2(the_blobguts,parms)
     Blob_internal *the_blobguts;
     Blob_parms *parms;
{
  double base,tmpa2,thresh,out;

  base = 2.0*( (*parms).du11*( (*parms).cos2 - (*parms).sin2 ) +
	       ( (*parms).du12+(*parms).du21 )*
	       (*parms).sincos )*(*the_blobguts).a2 +
    (visc/(2.0*(*the_blobguts).s2))*(1.0 - SQR((*the_blobguts).a2));

  tmpa2 = (*the_blobguts).a2;
  if (tmpa2<1.0) tmpa2 = 1.0/tmpa2;

  thresh = 0.5*(tanh(-2.0*(tmpa2-3.0))+1.0);

#ifdef A2THRESH
  /* Put in a threshold to prevent the aspect ratio from going berzerk. */
  out = base*thresh;
#else
  out = base;
#endif

  return(out);
  
}

double dtth(the_blobguts,parms)
Blob_internal *the_blobguts;
Blob_parms *parms;
{
  return( ((*parms).du21 - (*parms).du12)/2.0 +
	  ( (((*parms).du21 + (*parms).du12)/2.0)*
	    ((*parms).sin2-(*parms).cos2) +
	    2.0*(*parms).du11*(*parms).sincos )*
	  (1.0/(*the_blobguts).a2 + (*the_blobguts).a2)/
	  (1.0/(*the_blobguts).a2 - (*the_blobguts).a2) );
}
