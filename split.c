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

void chksplit()
{
   int i,j,k,l,m;
 
   for (i=0; i<oldN; ++i)
     {
       if (tmpparms[i].refinecnt == -1)
	 {
	   if ( blobguts[i].s2 + (MaxOrder-1.0)*dts2(&(blobguts[i]))*TimeStep> l2tol )
	     {
	       tmpparms[i].refinecnt = MaxOrder-1;
	       refineindex[refinestack] = i;
	       ++refinestack;
	     }
	 }
       else
	 --tmpparms[i].refinecnt;
     }

   for (j=0; j<refinestack; ++j)
     {
       i = refineindex[j];
       if (tmpparms[i].refinecnt == 0)
	 switch (split_method)
	 {
	 case 1:
	   /* 1 -> 4 splitting.  Second order. */
	   if (blobguts[i].s2 > l2tol)
	     {
	       mblob[N] = mblob[i];
	       mblob[N+1] = mblob[i];
	       mblob[N+2] = mblob[i];
	     
	       split14(&(mblob[i].blob0),&(blobguts[i]),&(tmpparms[i]),
		       &(mblob[N].blob0),&(blobguts[N]),&(tmpparms[N]),
		       &(mblob[N+1].blob0),&(blobguts[N+1]),&(tmpparms[N+1]),
		       &(mblob[N+2].blob0),&(blobguts[N+2]),&(tmpparms[N+2]));
	     
	       /*Watch out.  The new blobs should have the right strength
		 because the merging algorithm uses it to track the motion 
		 of the cluster. */

	       mblob[i].blob1.strength = mblob[i].blob0.strength;
	       mblob[i].blob2.strength = mblob[i].blob0.strength;
	       mblob[i].blob3.strength = mblob[i].blob0.strength;
	     
	       for (k=0; k<3; ++k)
		 {
		   mblob[N+k].blob1.strength = mblob[N+k].blob0.strength;
		   mblob[N+k].blob2.strength = mblob[N+k].blob0.strength;
		   mblob[N+k].blob3.strength = mblob[N+k].blob0.strength;
		 }

	       mblob[i].blob1.dx = refinevels[j][0][0].x;
	       mblob[i].blob1.dy = refinevels[j][0][0].y;
	       mblob[i].blob2.dx = refinevels[j][1][0].x;
	       mblob[i].blob2.dy = refinevels[j][1][0].y;
	       mblob[i].blob3.dx = refinevels[j][2][0].x;
	       mblob[i].blob3.dy = refinevels[j][2][0].y;

	       for (k=0; k<3; ++k)
		 {
		   mblob[N+k].blob1.dx = refinevels[j][0][k+1].x;
		   mblob[N+k].blob1.dy = refinevels[j][0][k+1].y;
		   mblob[N+k].blob2.dx = refinevels[j][1][k+1].x;
		   mblob[N+k].blob2.dy = refinevels[j][1][k+1].y;
		   mblob[N+k].blob3.dx = refinevels[j][2][k+1].x;
		   mblob[N+k].blob3.dy = refinevels[j][2][k+1].y;
		 }
		   
	       for (k=j; k<refinestack-1; ++k)
		 {
		   refineindex[k] = refineindex[k+1];
		   for (l=0; l<3; ++l)
		     for (m=0; m<MaxSplitConf; ++m)
		       refinevels[k][l][m] = refinevels[k+1][l][m];
		 }

	       --refinestack;
	       --j;
	       ++nsplit;
	       N=N+3;
	     }
	   else
	     {
	       /* Hold on.  Not fat enough to split yet. */
	       /* Put off the split until we are fat enough. */
	       tmpparms[i].refinecnt=1;
	       for (l=0; l<2; ++l)
		 for (m=0; m<4; ++m)
		   refinevels[j][2-l][m] = refinevels[j][1-l][m];
	       split14vels(refinevels[j][tmpparms[i].refinecnt-1],
			 mblob[i].blob0,blobguts[i],tmpparms[i]);
	     }
	   break;
	 case 2:
	   /* 1 -> 5 splitting.  Third order. */
	   if (blobguts[i].s2 > l2tol)
	     {
	       mblob[N] = mblob[i];
	       mblob[N+1] = mblob[i];
	       mblob[N+2] = mblob[i];
	       mblob[N+3] = mblob[i];
	     
	       split15(&(mblob[i].blob0),&(blobguts[i]),&(tmpparms[i]),
		       &(mblob[N].blob0),&(blobguts[N]),&(tmpparms[N]),
		       &(mblob[N+1].blob0),&(blobguts[N+1]),&(tmpparms[N+1]),
		       &(mblob[N+2].blob0),&(blobguts[N+2]),&(tmpparms[N+2]),
		       &(mblob[N+3].blob0),&(blobguts[N+3]),&(tmpparms[N+3]));
	     
	       /*Watch out.  The new blobs should have the right strength
		 because the merging algorithm uses it to track the motion 
		 of the cluster. */

	       mblob[i].blob1.strength = mblob[i].blob0.strength;
	       mblob[i].blob2.strength = mblob[i].blob0.strength;
	       mblob[i].blob3.strength = mblob[i].blob0.strength;
	     
	       for (k=0; k<4; ++k)
		 {
		   mblob[N+k].blob1.strength = mblob[N+k].blob0.strength;
		   mblob[N+k].blob2.strength = mblob[N+k].blob0.strength;
		   mblob[N+k].blob3.strength = mblob[N+k].blob0.strength;
		 }

	       mblob[i].blob1.dx = refinevels[j][0][0].x;
	       mblob[i].blob1.dy = refinevels[j][0][0].y;
	       mblob[i].blob2.dx = refinevels[j][1][0].x;
	       mblob[i].blob2.dy = refinevels[j][1][0].y;
	       mblob[i].blob3.dx = refinevels[j][2][0].x;
	       mblob[i].blob3.dy = refinevels[j][2][0].y;

	       for (k=0; k<4; ++k)
		 {
		   mblob[N+k].blob1.dx = refinevels[j][0][k+1].x;
		   mblob[N+k].blob1.dy = refinevels[j][0][k+1].y;
		   mblob[N+k].blob2.dx = refinevels[j][1][k+1].x;
		   mblob[N+k].blob2.dy = refinevels[j][1][k+1].y;
		   mblob[N+k].blob3.dx = refinevels[j][2][k+1].x;
		   mblob[N+k].blob3.dy = refinevels[j][2][k+1].y;
		 }
		   
	       for (k=j; k<refinestack-1; ++k)
		 {
		   refineindex[k] = refineindex[k+1];
		   for (l=0; l<3; ++l)
		     for (m=0; m<MaxSplitConf; ++m)
		       refinevels[k][l][m] = refinevels[k+1][l][m];
		 }

	       --refinestack;
	       --j;
	       ++nsplit;
	       N=N+4;
	     }
	   else
	     {
	       /* Hold on.  Not fat enough to split yet. */
	       /* Put off the split until we are fat enough. */
	       tmpparms[i].refinecnt=1;
	       for (l=0; l<2; ++l)
		 for (m=0; m<5; ++m)
		   refinevels[j][2-l][m] = refinevels[j][1-l][m];
	       split15vels(refinevels[j][tmpparms[i].refinecnt-1],
			 mblob[i].blob0,blobguts[i],tmpparms[i]);
	     }
	   break;
	 case 3:
	   /* 1 -> 5 asymmetric splitting.  Third order. */
	   if (blobguts[i].s2 > l2tol)
	     {
	       mblob[N] = mblob[i];
	       mblob[N+1] = mblob[i];
	       mblob[N+2] = mblob[i];
	       mblob[N+3] = mblob[i];

	       split15asym(&(mblob[i].blob0),&(blobguts[i]),&(tmpparms[i]),
			   &(mblob[N].blob0),&(blobguts[N]),&(tmpparms[N]),
			   &(mblob[N+1].blob0),&(blobguts[N+1]),
			   &(tmpparms[N+1]),
			   &(mblob[N+2].blob0),&(blobguts[N+2]),
			   &(tmpparms[N+2]),
			   &(mblob[N+3].blob0),&(blobguts[N+3]),
			   &(tmpparms[N+3]),split_parm_ptr);
	     
	       /*Watch out.  The new blobs should have the right strength
		 because the merging algorithm uses it to track the motion 
		 of the cluster. */

	       mblob[i].blob1.strength = mblob[i].blob0.strength;
	       mblob[i].blob2.strength = mblob[i].blob0.strength;
	       mblob[i].blob3.strength = mblob[i].blob0.strength;
	     
	       for (k=0; k<4; ++k)
		 {
		   mblob[N+k].blob1.strength = mblob[N+k].blob0.strength;
		   mblob[N+k].blob2.strength = mblob[N+k].blob0.strength;
		   mblob[N+k].blob3.strength = mblob[N+k].blob0.strength;
		 }

	       mblob[i].blob1.dx = refinevels[j][0][0].x;
	       mblob[i].blob1.dy = refinevels[j][0][0].y;
	       mblob[i].blob2.dx = refinevels[j][1][0].x;
	       mblob[i].blob2.dy = refinevels[j][1][0].y;
	       mblob[i].blob3.dx = refinevels[j][2][0].x;
	       mblob[i].blob3.dy = refinevels[j][2][0].y;

	       for (k=0; k<4; ++k)
		 {
		   mblob[N+k].blob1.dx = refinevels[j][0][k+1].x;
		   mblob[N+k].blob1.dy = refinevels[j][0][k+1].y;
		   mblob[N+k].blob2.dx = refinevels[j][1][k+1].x;
		   mblob[N+k].blob2.dy = refinevels[j][1][k+1].y;
		   mblob[N+k].blob3.dx = refinevels[j][2][k+1].x;
		   mblob[N+k].blob3.dy = refinevels[j][2][k+1].y;
		 }
		   
	       for (k=j; k<refinestack-1; ++k)
		 {
		   refineindex[k] = refineindex[k+1];
		   for (l=0; l<3; ++l)
		     for (m=0; m<MaxSplitConf; ++m)
		       refinevels[k][l][m] = refinevels[k+1][l][m];
		 }

	       --refinestack;
	       --j;
	       ++nsplit;
	       N=N+4;
	     }
	   else
	     {
	       /* Hold on.  Not fat enough to split yet. */
	       /* Put off the split until we are fat enough. */
	       tmpparms[i].refinecnt=1;
	       for (l=0; l<2; ++l)
		 for (m=0; m<5; ++m)
		   refinevels[j][2-l][m] = refinevels[j][1-l][m];
	       split15asymvels(refinevels[j][tmpparms[i].refinecnt-1],
			       mblob[i].blob0,blobguts[i],tmpparms[i],
			       split_parm_ptr);
	     }
	   break;
	 default:
	   printf("Error in SplitMethod.\n");
	   exit(-1);
	 }
       else
	 switch (split_method)
	 {
	 case 1:
	   /* 1 -> 4 splitting.  Second order. */
	   split14vels(refinevels[j][tmpparms[i].refinecnt-1],
		       mblob[i].blob0,blobguts[i],tmpparms[i]);
	   break;
	 case 2:
	   /* 1 -> 5 splitting.  Third order. */
	   split15vels(refinevels[j][tmpparms[i].refinecnt-1],
		       mblob[i].blob0,blobguts[i],tmpparms[i]);
	   break;
	 case 3:
	   /* 1 -> 5 asymmetric splitting.  Third order. */
	   split15asymvels(refinevels[j][tmpparms[i].refinecnt-1],
			   mblob[i].blob0,blobguts[i],tmpparms[i],
			   split_parm_ptr);
	   break;
	 default:
	   printf("Error in SplitMethod.\n");
	   exit(-1);
	 }
     }
}
