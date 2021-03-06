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

#include "global_min.h"
#include "particle.h"
#include "boundary.h"
#include "biot-savart.h"
#include "linear_velocity.h"
#include "potential.h"

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

void neg_a2_check()
{
  int j;
  int failure = 0;

  for (j=0; j<N; ++j)
    if (blobguts[j].a2<0.0)
      {
	fprintf(diag_log,"Blob %d has a negative aspect ratio of %lf",
		j,blobguts[j].a2);
	failure = 1;
      }

  if (failure==1)
    {
      printf("BlobFlow is halting because at least one blob has an aspect ratio that is\n");
      printf("negative. This is generally a sign that BlobFlow is not capturing the\n");
      printf("the flow dynamics with the chosen TimeStep, and you will need to decrease it\n");
      printf("to resolve the flow.\n");
      printf("It is time to sleep.\n");
      exit(-1);
    }
}

void vel_field()
{
   int j;

   vel_cputime_ref = clock();
   vel_cputime = 0;
   velsum_cputime = 0;
   veldirect_cputime = 0;
   mp_cputime = 0;

   for (j=0; j<N; ++j)
     set_blob(&(blobguts[j]),&(tmpparms[j]));
   
   if (xantisymm)
     /* Double the number of computational elements by reflecting them
      * about the x-axis. */
     reflect_X();

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

   if (xantisymm)
     for (j=0; j<N/2; ++j) 
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
   else
     for (j=0; j<N; ++j) 
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
   
   if (xantisymm)
   /* re-adjust */
     N /= 2;
   
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

   vel_cputime += clock()-vel_cputime_ref;
   fprintf(cpu_log,"%07d %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e\n",
	   N,
	   ((double)(vel_cputime))/((double)CLOCKS_PER_SEC),
	   ((double)(velsum_cputime))/((double)CLOCKS_PER_SEC),
	   ((double)(veldirect_cputime))/((double)CLOCKS_PER_SEC),
	   ((double)(mp_cputime))/((double)CLOCKS_PER_SEC),
	   ((double)(mp_Init_Fine_Grid))/((double)CLOCKS_PER_SEC),
	   ((double)(mp_Advance_Coeffs))/((double)CLOCKS_PER_SEC)
	   );
   fflush(cpu_log);
}

void dpos_vel(int vort)
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

void dpos_vel_fast(int vort)
{
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

