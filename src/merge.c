/* MERGE.C */
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
 * Newark, DE 19715-2553                                                */

/* MERGECHECK is a compiler flag to dump some diagnostics. */
/* It can be very useful because merging algorithms can be tough to debug. */

/* MERGEDIAG is a compiler flag to write merge information to the
   diagnostic log at runtime. */

#include "global.h"

void rectify()
{
   int i;
   
   /* Make sure all the blobs are oriented so that the prefered axis is the
    * major axis. */
   
  for (i=0; i<N; ++i)
    if (blobguts[i].a2 < 1.0)
     {
	flip_blob(&(blobguts[i]));
	set_blob(&(blobguts[i]),&(tmpparms[i]));
     }
}

double m_xx(int vort)
{
  return(mblob[vort].blob0.strength*(SQR(mblob[vort].blob0.x)+
				     2.0*blobguts[vort].s2*
				     (tmpparms[vort].sin2/
				      blobguts[vort].a2+
				      tmpparms[vort].cos2*
				      blobguts[vort].a2)));
}

double m_yy(int vort)
{
  return(mblob[vort].blob0.strength*(SQR(mblob[vort].blob0.y)+
				     2.0*blobguts[vort].s2*
				     (tmpparms[vort].cos2/
				      blobguts[vort].a2+
				      tmpparms[vort].sin2*
				      blobguts[vort].a2)));
}

double m_xy(int vort)
{
  return(mblob[vort].blob0.strength*(mblob[vort].blob0.x*
				     mblob[vort].blob0.y+
				     2.0*blobguts[vort].s2*
				     tmpparms[vort].sincos*
				     (blobguts[vort].a2-
				      1.0/blobguts[vort].a2)));
}

double m_xxx(int vort)
{
  return(mblob[vort].blob0.strength*(6.0*blobguts[vort].s2*
				     (tmpparms[vort].sin2/blobguts[vort].a2+
				      tmpparms[vort].cos2*blobguts[vort].a2)*
				     mblob[vort].blob0.x+
				     SQR(mblob[vort].blob0.x)*
				     mblob[vort].blob0.x));
}

double m_xxy(int vort)
{
  return(mblob[vort].blob0.strength*
	 (blobguts[vort].s2*(4.0*tmpparms[vort].sincos*mblob[vort].blob0.x*
			     (blobguts[vort].a2-1.0/blobguts[vort].a2)+
			     2.0*mblob[vort].blob0.y*
			     (tmpparms[vort].sin2/blobguts[vort].a2+
			      tmpparms[vort].cos2*blobguts[vort].a2))+
	  SQR(mblob[vort].blob0.x)*mblob[vort].blob0.y));
}

double m_xyy(int vort)
{
  return(mblob[vort].blob0.strength*
	 (blobguts[vort].s2*(4.0*tmpparms[vort].sincos*mblob[vort].blob0.y*
			     (blobguts[vort].a2-1.0/blobguts[vort].a2)+
			     2.0*mblob[vort].blob0.x*
			     (tmpparms[vort].sin2*blobguts[vort].a2+
			      tmpparms[vort].cos2/blobguts[vort].a2))+
	  mblob[vort].blob0.x*SQR(mblob[vort].blob0.y)));
}

double m_yyy(int vort)
{
  return(mblob[vort].blob0.strength*(6.0*blobguts[vort].s2*
				     (tmpparms[vort].sin2*blobguts[vort].a2+
				      tmpparms[vort].cos2/blobguts[vort].a2)*
				     mblob[vort].blob0.y+
				     SQR(mblob[vort].blob0.y)*
				     mblob[vort].blob0.y));
}

double m_xxxx(int vort)
{
  return(mblob[vort].blob0.strength*
	 (12.0*SQR(blobguts[vort].s2)*
	  (SQR(tmpparms[vort].sin2/blobguts[vort].a2)+
	   SQR(tmpparms[vort].cos2*blobguts[vort].a2)+
	   2.0*SQR(tmpparms[vort].sincos))+
	  12.0*blobguts[vort].s2*SQR(mblob[vort].blob0.x)*
	  (tmpparms[vort].sin2/blobguts[vort].a2+
	   tmpparms[vort].cos2*blobguts[vort].a2)+
	  SQR(SQR(mblob[vort].blob0.x))));
}

double m_xxxy(int vort)
{
  return(mblob[vort].blob0.strength*
	 (12.0*SQR(blobguts[vort].s2)*tmpparms[vort].sincos*
	  (-tmpparms[vort].sin2/SQR(blobguts[vort].a2)+
	   tmpparms[vort].cos2*SQR(blobguts[vort].a2)+
	   (tmpparms[vort].sin2-tmpparms[vort].cos2))+
	  6.0*blobguts[vort].s2*
	  (SQR(mblob[vort].blob0.x)*tmpparms[vort].sincos*
	   (blobguts[vort].a2-1.0/blobguts[vort].a2)+
	   mblob[vort].blob0.x*mblob[vort].blob0.y*
	   (tmpparms[vort].sin2/blobguts[vort].a2+
	    tmpparms[vort].cos2*blobguts[vort].a2))+
	  SQR(mblob[vort].blob0.x)*mblob[vort].blob0.x*
	  mblob[vort].blob0.y));
}

double m_xxyy(int vort)
{
  return(mblob[vort].blob0.strength*
	 (4.0*SQR(blobguts[vort].s2)*
	  (SQR(tmpparms[vort].sincos)*
	   (3.0*SQR(blobguts[vort].a2)+3.0/SQR(blobguts[vort].a2)-4.0)+
	   SQR(tmpparms[vort].sin2)+SQR(tmpparms[vort].cos2))+
	  2.0*blobguts[vort].s2*
	  (SQR(mblob[vort].blob0.x)*
	   (tmpparms[vort].sin2*blobguts[vort].a2+
	    tmpparms[vort].cos2/blobguts[vort].a2)+
	   SQR(mblob[vort].blob0.y)*
	   (tmpparms[vort].sin2/blobguts[vort].a2+
	    tmpparms[vort].cos2*blobguts[vort].a2)+
	   4.0*mblob[vort].blob0.x*mblob[vort].blob0.y*
	   tmpparms[vort].sincos*
	   (blobguts[vort].a2-1.0/blobguts[vort].a2))+
	  SQR(mblob[vort].blob0.x*mblob[vort].blob0.y)));
}

double m_xyyy(int vort)
{
  return(mblob[vort].blob0.strength*
	 (12.0*SQR(blobguts[vort].s2)*tmpparms[vort].sincos*
	  (-tmpparms[vort].cos2/SQR(blobguts[vort].a2)+
	   tmpparms[vort].sin2*SQR(blobguts[vort].a2)+
	   (tmpparms[vort].cos2-tmpparms[vort].sin2))+
	  6.0*blobguts[vort].s2*
	  (SQR(mblob[vort].blob0.y)*tmpparms[vort].sincos*
	   (blobguts[vort].a2-1.0/blobguts[vort].a2)+
	   mblob[vort].blob0.x*mblob[vort].blob0.y*
	   (tmpparms[vort].sin2*blobguts[vort].a2+
	    tmpparms[vort].cos2/blobguts[vort].a2))+
	  mblob[vort].blob0.x*SQR(mblob[vort].blob0.y)*
	  mblob[vort].blob0.y));
}

double m_yyyy(int vort)
{
  return(mblob[vort].blob0.strength*
	 (12.0*SQR(blobguts[vort].s2)*
	  (SQR(tmpparms[vort].cos2/blobguts[vort].a2)+
	   SQR(tmpparms[vort].sin2*blobguts[vort].a2)+
	   2.0*SQR(tmpparms[vort].sincos))+
	  12.0*blobguts[vort].s2*SQR(mblob[vort].blob0.y)*
	  (tmpparms[vort].cos2/blobguts[vort].a2+
	   tmpparms[vort].sin2*blobguts[vort].a2)+
	  SQR(SQR(mblob[vort].blob0.y))));
}

double searchlist_distribution(ind,indmax,maxerr,size)
int ind[NMAX],indmax;
double maxerr;
int *size;
{
   int mergelist[NMAX],mergeable,i,j,k,l,m,gotone,listorder;
   double precost,postcost,tmppostcost;
   
   double cum_tot[15],this_tot[15],trial_tot[15],tmpa2,tmpa4,tmps2;
   double err_estimate,curr_err_estimate;
   double tmpC,tmpD,naM2,tmpnaM2,reg_scale;

   /* 0  m_0
      1  m_x
      2  m_y
      3  m_xx
      4  m_yy
      5  m_xy
      6  m_xxx
      7  m_xxy
      8  m_xyy
      9  m_yyy  
     10  m_xxxx
     11  m_xxxy
     12  m_xxyy
     13  m_xyyy
     14  m_yyyy */

   curr_err_estimate = 0.0;
   precost = 0.0;
   postcost = 0.0;

   mergelist[0] = ind[0];
   mergeable = 1;

   if (blobguts[ind[0]].a2>1.0)
     naM2 = blobguts[ind[0]].a2;
   else
     naM2 = 1.0/blobguts[ind[0]].a2;

   cum_tot[0] = mblob[ind[0]].blob0.strength;

   /* First order moments. */

   cum_tot[1] = 
     mblob[ind[0]].blob0.strength*mblob[ind[0]].blob0.x;
   cum_tot[2] = 
     mblob[ind[0]].blob0.strength*mblob[ind[0]].blob0.y;

   /* Second order moments. */

   cum_tot[3] = m_xx(ind[0]);
   cum_tot[4] = m_yy(ind[0]);
   cum_tot[5] = m_xy(ind[0]);

   /* Third order moments. */

   cum_tot[6] = m_xxx(ind[0]);
   cum_tot[7] = m_xxy(ind[0]);
   cum_tot[8] = m_xyy(ind[0]);
   cum_tot[9] = m_yyy(ind[0]);

   /* Fourth order moments. */
   
   cum_tot[10] = m_xxxx(ind[0]);
   cum_tot[11] = m_xxxy(ind[0]);
   cum_tot[12] = m_xxyy(ind[0]);
   cum_tot[13] = m_xyyy(ind[0]);
   cum_tot[14] = m_yyyy(ind[0]);

   gotone = 1;

   precost = l2tol*(1.0-4.0*SQR(alpha))+3.0*blobguts[ind[0]].s2;
   
   while (gotone==1)
     {
	gotone = 0;
	
	for (j=1; j<indmax; ++j)
	  {
	     this_tot[0] = mblob[ind[j]].blob0.strength;

	     /* First order moments */
	     this_tot[1] = 
	       mblob[ind[j]].blob0.strength*mblob[ind[j]].blob0.x;
	     this_tot[2] = 
	       mblob[ind[j]].blob0.strength*mblob[ind[j]].blob0.y;

	     /* Second order moments */

	     this_tot[3] = m_xx(ind[j]);
	     this_tot[4] = m_yy(ind[j]);
	     this_tot[5] = m_xy(ind[j]);

	     /* Third order moments */

	     this_tot[6] = m_xxx(ind[j]);
	     this_tot[7] = m_xxy(ind[j]);
	     this_tot[8] = m_xyy(ind[j]);
	     this_tot[9] = m_yyy(ind[j]);

	     /* Fourth order moments */
   
	     this_tot[10] = m_xxxx(ind[j]);
	     this_tot[11] = m_xxxy(ind[j]);
	     this_tot[12] = m_xxyy(ind[j]);
	     this_tot[13] = m_xyyy(ind[j]);
	     this_tot[14] = m_yyyy(ind[j]);

	     /* Now move all moments to the center or vorticity */

	     trial_tot[0] = cum_tot[0]+this_tot[0];
	     trial_tot[1] = (cum_tot[1]+this_tot[1])/trial_tot[0];
	     trial_tot[2] = (cum_tot[2]+this_tot[2])/trial_tot[0];

	     /* Try doing this in reverse order to simplify the
		calculation by using moments relative to the
		origin. */

	     /* Fourth order moments. */

	     trial_tot[10] = 
	       ((cum_tot[10]+this_tot[10])-
		4.0*trial_tot[1]*(cum_tot[6]+this_tot[6])+
		6.0*SQR(trial_tot[1])*(cum_tot[3]+this_tot[3]))/
	       trial_tot[0]-
	       3.0*SQR(SQR(trial_tot[1]));

	     trial_tot[11] = 
	       ((cum_tot[11]+this_tot[11])-
		3.0*trial_tot[1]*(cum_tot[7]+this_tot[7])+
		3.0*SQR(trial_tot[1])*(cum_tot[5]+this_tot[5])-
		trial_tot[2]*(cum_tot[6]+this_tot[6])+
		3.0*trial_tot[1]*trial_tot[2]*
		(cum_tot[3]+this_tot[3]))/trial_tot[0]-
	       3.0*SQR(trial_tot[1])*trial_tot[1]*trial_tot[2];

	     trial_tot[12] =
	       ((cum_tot[12]+this_tot[12])-
		2.0*trial_tot[2]*(cum_tot[7]+this_tot[7])-
		2.0*trial_tot[1]*(cum_tot[8]+this_tot[8])+
		SQR(trial_tot[2])*(cum_tot[3]+this_tot[3])+
		SQR(trial_tot[1])*(cum_tot[4]+this_tot[4])+
		4.0*trial_tot[1]*trial_tot[2]*(cum_tot[5]+this_tot[5]))/
	       trial_tot[0]-
	       3.0*SQR(trial_tot[1]*trial_tot[2]);
	     
	     trial_tot[13] = 
	       ((cum_tot[13]+this_tot[13])-
		3.0*trial_tot[2]*(cum_tot[8]+this_tot[8])+
		3.0*SQR(trial_tot[2])*(cum_tot[5]+this_tot[5])-
		trial_tot[1]*(cum_tot[9]+this_tot[9])+
		3.0*trial_tot[1]*trial_tot[2]*
		(cum_tot[4]+this_tot[4]))/trial_tot[0]-
	       3.0*SQR(trial_tot[2])*trial_tot[1]*trial_tot[2];
	     
	     trial_tot[14] = 
	       ((cum_tot[14]+this_tot[14])-
		4.0*trial_tot[2]*(cum_tot[9]+this_tot[9])+
		6.0*SQR(trial_tot[2])*(cum_tot[4]+this_tot[4]))/
	       trial_tot[0]-
	       3.0*SQR(SQR(trial_tot[2]));

	     /* Third order moments. */

	     trial_tot[6] = 
	       ((cum_tot[6]+this_tot[6])-
		3.0*trial_tot[1]*(cum_tot[3]+this_tot[3]))/trial_tot[0]+
	       2.0*SQR(trial_tot[1])*trial_tot[1];
	     trial_tot[7] = 
	       ((cum_tot[7]+this_tot[7])-
		2.0*trial_tot[1]*(cum_tot[5]+this_tot[5])-
		trial_tot[2]*(cum_tot[3]+this_tot[3]))/trial_tot[0]+
	       2.0*SQR(trial_tot[1])*trial_tot[2];
	     trial_tot[8] = 
	       ((cum_tot[8]+this_tot[8])-
		2.0*trial_tot[2]*(cum_tot[5]+this_tot[5])-
		trial_tot[1]*(cum_tot[4]+this_tot[4]))/trial_tot[0]+
	       2.0*SQR(trial_tot[2])*trial_tot[1];
	     trial_tot[9] = 
	       ((cum_tot[9]+this_tot[9])-
		3.0*trial_tot[2]*(cum_tot[4]+this_tot[4]))/trial_tot[0]+
	       2.0*SQR(trial_tot[2])*trial_tot[2];

	     /* Second order moments. */

	     trial_tot[3] = 
	       (cum_tot[3]+this_tot[3])/trial_tot[0]-SQR(trial_tot[1]);
	     trial_tot[4] = 
	       (cum_tot[4]+this_tot[4])/trial_tot[0]-SQR(trial_tot[2]);
	     trial_tot[5] = 
	       (cum_tot[5]+this_tot[5])/trial_tot[0]-
	       trial_tot[1]*trial_tot[2];

	     tmpa4 = 
	       (SQR(trial_tot[3])+SQR(trial_tot[4])+2.0*SQR(trial_tot[5])+
		(trial_tot[3]+trial_tot[4])*
		sqrt(SQR(trial_tot[3]-trial_tot[4])+4.0*SQR(trial_tot[5])))/
	       (2.0*(trial_tot[3]*trial_tot[4]-SQR(trial_tot[5])));
	     
	     tmpa2 = sqrt(tmpa4);
	     
	     tmps2 = 
	       (trial_tot[3]+trial_tot[4])*tmpa2/
	       (2.0*(tmpa4+1.0));

	     /* Remember that the moments are scaled by a factor of
		1/2/Pi, so we must correct them here. */

	     reg_scale = -0.5*exp(1.0)*SQR(alpha)/(1.0-SQR(alpha))*log(alpha);

	     err_estimate =
	       fabs(trial_tot[0])*
	       (fabs(trial_tot[10]+
		     2.0*trial_tot[12]+trial_tot[14]-
		     SQR(tmps2)*(12.0*(tmpa4+1.0/tmpa4)+8.0))+
		4.0*sqrt(reg_scale*tmps2)*
		sqrt(SQR(trial_tot[6]+trial_tot[8])+
		     SQR(trial_tot[9]+trial_tot[7])))
	       /pow(reg_scale*tmps2,3.0)/64.0;

	     tmppostcost = l2tol*(1.0-4.0*SQR(alpha))+3.0*tmps2;

	     if (blobguts[ind[j]].a2>naM2)
	       tmpnaM2 = blobguts[ind[j]].a2;
	     else if (1.0/blobguts[ind[j]].a2>naM2)
	       tmpnaM2 = 1.0/blobguts[ind[j]].a2;
	     else
	       tmpnaM2 = naM2;

	     /*
	     if ( (tmpa2 < aM2) &&
		  (tmpa2 > 1.0/aM2) &&
		  (tmps2 < l2tol) &&
		  (tmps2 > alpha*alpha*l2tol) &&
		  (err_estimate <
		   maxerr*mergeable/(mergeable+indmax) ) )
	       */
	     if ( (tmpa2 < tmpnaM2) &&
		  (tmpa2 > 1.0/tmpnaM2) &&
		  (tmps2 < l2tol) &&
		  (tmps2 > alpha*alpha*l2tol) &&
		  (err_estimate < maxerr) )
	       {
		 naM2 = tmpnaM2;

		 precost += l2tol*(1.0-4.0*SQR(alpha))+3.0*blobguts[ind[j]].s2;
		 postcost = tmppostcost;

		 curr_err_estimate=err_estimate;

		 mergelist[mergeable] = ind[j];
		 ++mergeable;
		 
		 for (k=j; k<indmax-1; ++k)
		   ind[k] = ind[k+1];
		 --indmax;
		 --j;
		 
		 for (k=0; k<15; k++)
		   cum_tot[k] += this_tot[k];
		 
		 gotone = 1;
	       }
	  } /*j loop */
     } /* gotone loop */

   /* Assign centroid velocity to post-merger element. */

   if ( (mergeable > 1) && (precost > postcost) ) 
     {
       if (tmpparms[mergelist[0]].refinecnt != -1)
	 {
	   for (i=0; refineindex[i] != mergelist[0]; ++i);
	   
	   for (k=i; k<refinestack-1; ++k)
	     {
	       refineindex[k] = refineindex[k+1];
	       for (l=0; l<3; ++l)
		 for (m=0; m<4; ++m)
		   refinevels[k][l][m] = refinevels[k+1][l][m];
	     }
	   --refinestack;
	 }

       mblob[mergelist[0]].blob1.dx *= mblob[mergelist[0]].blob1.strength;
       mblob[mergelist[0]].blob1.dy *= mblob[mergelist[0]].blob1.strength;

       mblob[mergelist[0]].blob2.dx *= mblob[mergelist[0]].blob2.strength;
       mblob[mergelist[0]].blob2.dy *= mblob[mergelist[0]].blob2.strength;

       mblob[mergelist[0]].blob3.dx *= mblob[mergelist[0]].blob3.strength;
       mblob[mergelist[0]].blob3.dy *= mblob[mergelist[0]].blob3.strength;

       mblob[mergelist[0]].blob4.dx *= mblob[mergelist[0]].blob4.strength;
       mblob[mergelist[0]].blob4.dy *= mblob[mergelist[0]].blob4.strength;

       listorder = mblob[mergelist[0]].blob0.order;

       for (j=1; j<mergeable; ++j)
	 {
	   mblob[mergelist[0]].blob1.dx +=
	     mblob[mergelist[j]].blob1.strength*mblob[mergelist[j]].blob1.dx;
	   mblob[mergelist[0]].blob1.dy +=
	     mblob[mergelist[j]].blob1.strength*mblob[mergelist[j]].blob1.dy;

	   mblob[mergelist[0]].blob2.dx +=
	     mblob[mergelist[j]].blob2.strength*mblob[mergelist[j]].blob2.dx;
	   mblob[mergelist[0]].blob2.dy +=
	     mblob[mergelist[j]].blob2.strength*mblob[mergelist[j]].blob2.dy;

	   mblob[mergelist[0]].blob3.dx +=
	     mblob[mergelist[j]].blob3.strength*mblob[mergelist[j]].blob3.dx;
	   mblob[mergelist[0]].blob3.dy +=
	     mblob[mergelist[j]].blob3.strength*mblob[mergelist[j]].blob3.dy;

	   mblob[mergelist[0]].blob4.dx +=
	     mblob[mergelist[j]].blob4.strength*mblob[mergelist[j]].blob4.dx;
	   mblob[mergelist[0]].blob4.dy +=
	     mblob[mergelist[j]].blob4.strength*mblob[mergelist[j]].blob4.dy;

	   if (mblob[mergelist[j]].blob0.order < listorder)
	     listorder = mblob[mergelist[j]].blob0.order;

	   if (tmpparms[mergelist[j]].refinecnt != -1)
	     {
	       for (i=0; refineindex[i] != mergelist[j]; ++i);

	       for (k=i; k<refinestack-1; ++k)
		 {
		   refineindex[k] = refineindex[k+1];
		   for (l=0; l<3; ++l)
		     for (m=0; m<4; ++m)
		       refinevels[k][l][m] = refinevels[k+1][l][m];
		 }
	       --refinestack;
	     }

	 }

       mblob[mergelist[0]].blob1.dx /= cum_tot[0];
       mblob[mergelist[0]].blob1.dy /= cum_tot[0];
       mblob[mergelist[0]].blob2.dx /= cum_tot[0];
       mblob[mergelist[0]].blob2.dy /= cum_tot[0];
       mblob[mergelist[0]].blob3.dx /= cum_tot[0];
       mblob[mergelist[0]].blob3.dy /= cum_tot[0];
       mblob[mergelist[0]].blob4.dx /= cum_tot[0];
       mblob[mergelist[0]].blob4.dy /= cum_tot[0];

#ifdef MERGECHECK
	 for (j=0; j<mergeable; ++j)
	   printf("%10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n",
		  mblob[mergelist[j]].blob0.x,
		  mblob[mergelist[j]].blob0.y,
		  mblob[mergelist[j]].blob0.strength,
		  blobguts[mergelist[j]].s2,
		  blobguts[mergelist[j]].a2,
		  blobguts[mergelist[j]].th);
	 printf("\n");
#endif

       /* Nullify all of the merged vortices. */

       for (j=1; j<mergeable; ++j)
	 mblob[mergelist[j]].blob0.strength = 0.0;
       
       mblob[mergelist[0]].blob0.strength = cum_tot[0];
       mblob[mergelist[0]].blob0.x        = cum_tot[1]/cum_tot[0];
       mblob[mergelist[0]].blob0.y        = cum_tot[2]/cum_tot[0];
       
       cum_tot[1] /= cum_tot[0];
       cum_tot[2] /= cum_tot[0];
       
       cum_tot[3] = cum_tot[3]/cum_tot[0]-SQR(cum_tot[1]);
       cum_tot[4] = cum_tot[4]/cum_tot[0]-SQR(cum_tot[2]);
       cum_tot[5] = cum_tot[5]/cum_tot[0]-cum_tot[1]*cum_tot[2];
       
       tmpa4 = 
	 (SQR(cum_tot[3])+SQR(cum_tot[4])+2.0*SQR(cum_tot[5])+
	  (cum_tot[3]+cum_tot[4])*
	  sqrt(SQR(cum_tot[3]-cum_tot[4])+4.0*SQR(cum_tot[5])))/
	 (2.0*(cum_tot[3]*cum_tot[4]-SQR(cum_tot[5])));
       
       tmpa2 = sqrt(tmpa4);
       
       tmps2 = 
	 (cum_tot[3]+cum_tot[4])*tmpa2/
	 (2.0*(tmpa4+1.0));
       
       blobguts[mergelist[0]].a2       = tmpa2;
       blobguts[mergelist[0]].s2       = tmps2;
       
       tmpC = cum_tot[3]*(tmpa2-1.0/tmpa2);
       tmpD = cum_tot[4]*(tmpa2-1.0/tmpa2);

       if (fabs(tmpa2*(tmpC+sqrt(SQR(tmpC)-4.0*SQR(cum_tot[5])))-
		1.0/tmpa2*(tmpD+sqrt(SQR(tmpD)-4.0*SQR(cum_tot[5])))) <
	   fabs(tmpa2*(tmpC-sqrt(SQR(tmpC)-4.0*SQR(cum_tot[5])))-
		1.0/tmpa2*(tmpD-sqrt(SQR(tmpD)-4.0*SQR(cum_tot[5])))))
	 blobguts[mergelist[0]].th =
	   atan2(tmpa2*(tmpC+sqrt(SQR(tmpC)-4.0*SQR(cum_tot[5])))/2.0,
		 cum_tot[5]);
       else
	 blobguts[mergelist[0]].th =
	   atan2(tmpa2*(tmpC-sqrt(SQR(tmpC)-4.0*SQR(cum_tot[5])))/2.0,
		 cum_tot[5]);
       
       mblob[mergelist[0]].blob0.order    = listorder;
       tmpparms[mergelist[0]].refinecnt = -1;
       
       set_blob(&(blobguts[mergelist[0]]),&(tmpparms[mergelist[0]]));

#ifdef XANTISYMM
       /* Double the number of computational elements by reflecting them
	* about the x-axis. */
       reflect_X();
#endif 

#ifdef LINEAR
       dpos_vel_linear(mergelist[0]);
#else
       dpos_vel(mergelist[0]);
#endif

#ifdef XANTISYMM 
       /* re-adjust */
       N /= 2;
#endif   


#ifdef MERGECHECK
	 j=0;
	 printf("%10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n",
		mblob[mergelist[j]].blob0.x,
		mblob[mergelist[j]].blob0.y,
		mblob[mergelist[j]].blob0.strength,
		blobguts[mergelist[j]].s2,
		blobguts[mergelist[j]].a2,
		blobguts[mergelist[j]].th);

	 /* exit(3); */
#endif

       ++nmerge;

       *size -= (mergeable-1);

#ifdef MERGEDIAG
       fprintf(diag_log,"<%d> ",mergeable);
       fprintf(diag_log,"Error estimate: %10.4e\n",curr_err_estimate);
#endif
     }
   else
     /* If there is no merger then there is no induced error. */
     curr_err_estimate = 0.0;

   return(curr_err_estimate);
}

double searchlist_uniform(ind,indmax,maxerr,size)
int ind[NMAX],indmax;
double maxerr;
int    *size;
{
   int mergelist[NMAX],mergeable,i,j,k,gotone;
   
   double cum_tot[6],this_tot[6],trial_tot[6],tmpa2,tmpa4,tmps2,
             a2_rescale,s2_rescale,r_rescale,relerr2,relerr3,temp,relerr;
   
   relerr = (aM2/SQR(alpha) - 1.0)*
     pow(aM2/SQR(SQR(alpha)),1.0/(aM2/SQR(alpha)-1.0)) +
     0.86*clusterR/SQR(alpha);

   mergelist[0] = ind[0];
   mergeable = 1;
   
   cum_tot[0] = mblob[ind[0]].blob0.strength;
   cum_tot[1] = 
     mblob[ind[0]].blob0.strength*mblob[ind[0]].blob0.x;
   cum_tot[2] = 
     mblob[ind[0]].blob0.strength*mblob[ind[0]].blob0.y;
   cum_tot[3] = 
     mblob[ind[0]].blob0.strength*(SQR(mblob[ind[0]].blob0.x)+
				   2.0*blobguts[ind[0]].s2*
				   (tmpparms[ind[0]].sin2/
				    blobguts[ind[0]].a2+
				    tmpparms[ind[0]].cos2*
				    blobguts[ind[0]].a2));
   cum_tot[4] = 
     mblob[ind[0]].blob0.strength*(SQR(mblob[ind[0]].blob0.y)+
				   2.0*blobguts[ind[0]].s2*
				   (tmpparms[ind[0]].cos2/
				    blobguts[ind[0]].a2+
				    tmpparms[ind[0]].sin2*
				    blobguts[ind[0]].a2));
   cum_tot[5] = 
     mblob[ind[0]].blob0.strength*(mblob[ind[0]].blob0.x*
				   mblob[ind[0]].blob0.y+
				   2.0*blobguts[ind[0]].s2*
				   tmpparms[ind[0]].sincos*
				   (blobguts[ind[0]].a2-
				    1.0/blobguts[ind[0]].a2));
   
   gotone = 1;
   relerr3 = 0.0;
   
   while (gotone==1)
     {
	gotone = 0;
	
	for (j=1; j<indmax; ++j)
	  {
	     this_tot[0] = mblob[ind[j]].blob0.strength;
	     this_tot[1] = 
	       mblob[ind[j]].blob0.strength*mblob[ind[j]].blob0.x;
	     this_tot[2] = 
	       mblob[ind[j]].blob0.strength*mblob[ind[j]].blob0.y;
	     this_tot[3] = 
	       mblob[ind[j]].blob0.strength*(SQR(mblob[ind[j]].blob0.x)+
					     2.0*blobguts[ind[j]].s2*
					     (tmpparms[ind[j]].sin2/
					      blobguts[ind[j]].a2+
					      tmpparms[ind[j]].cos2*
					      blobguts[ind[j]].a2));
	     this_tot[4] = 
	       mblob[ind[j]].blob0.strength*(SQR(mblob[ind[j]].blob0.y)+
					     2.0*blobguts[ind[j]].s2*
					     (tmpparms[ind[j]].cos2/
					      blobguts[ind[j]].a2+
					      tmpparms[ind[j]].sin2*
					      blobguts[ind[j]].a2));
	     this_tot[5] = 
	       mblob[ind[j]].blob0.strength*(mblob[ind[j]].blob0.x*
					     mblob[ind[j]].blob0.y+
					     2.0*blobguts[ind[j]].s2*
					     tmpparms[ind[j]].sincos*
					     (blobguts[ind[j]].a2-
					      1.0/blobguts[ind[j]].a2));
	     
	     trial_tot[0] = cum_tot[0]+this_tot[0];
	     trial_tot[1] = (cum_tot[1]+this_tot[1])/trial_tot[0];
	     trial_tot[2] = (cum_tot[2]+this_tot[2])/trial_tot[0];
	     trial_tot[3] = 
	       (cum_tot[3]+this_tot[3])/trial_tot[0]-SQR(trial_tot[1]);
	     trial_tot[4] = 
	       (cum_tot[4]+this_tot[4])/trial_tot[0]-SQR(trial_tot[2]);
	     trial_tot[5] = 
	       (cum_tot[5]+this_tot[5])/trial_tot[0]-trial_tot[1]*trial_tot[2];
	     
	     tmpa4 = 
	       (SQR(trial_tot[3])+SQR(trial_tot[4])+2.0*SQR(trial_tot[5])+
		(trial_tot[3]+trial_tot[4])*
		sqrt(SQR(trial_tot[3]-trial_tot[4])+4.0*SQR(trial_tot[5])))/
	       (2.0*(trial_tot[3]*trial_tot[4]-SQR(trial_tot[5])));
	     
	     tmpa2 = sqrt(tmpa4);
	     
	     tmps2 = 
	       (trial_tot[3]+trial_tot[4])*tmpa2/
	       (2.0*(tmpa4+1.0));
	     
	     if ( (tmpa2 < aM2) &&
		 (tmpa2 > 1.0/aM2) &&
		 (tmps2 < l2tol) )
	       {
		  r_rescale = 0.0;
		  a2_rescale = 1.0;
		  s2_rescale = 1.0;
		  
		  for (i=0; i<mergeable; ++i)
		    {
		       temp = 
			 sqrt((SQR(trial_tot[1]-mblob[ind[i]].blob0.x) +
			       SQR(trial_tot[2]-mblob[ind[i]].blob0.y))/
			      tmps2)/2.0;
		       if (temp > r_rescale) r_rescale = temp;
		       
		       temp = blobguts[ind[i]].s2/tmps2;
		       if (temp<1.0) temp=1.0/temp;
		       if (temp > s2_rescale) s2_rescale=temp;
		       
		       temp = blobguts[ind[i]].a2/tmpa2;
		       if (temp<1.0) temp=1.0/temp;
		       if (temp > a2_rescale) a2_rescale=temp;
		    }
		  
		  relerr2 = (s2_rescale*a2_rescale - 1.0)*
		    pow(a2_rescale*SQR(s2_rescale),
			1.0/(1.0/(s2_rescale*a2_rescale)-1.0)) +
		    0.86*r_rescale;

		  if (relerr2*fabs(trial_tot[0])/2.0/l2tol < 
		      maxerr)
		    {
		       relerr3 = relerr2;
		       
		       mergelist[mergeable] = ind[j];
		       ++mergeable;
		       
		       for (k=j; k<indmax-1; ++k)
			 ind[k] = ind[k+1];
		       --indmax;
		       --j;
		       
		       cum_tot[0] += this_tot[0];
		       cum_tot[1] += this_tot[1];
		       cum_tot[2] += this_tot[2];
		       cum_tot[3] += this_tot[3];
		       cum_tot[4] += this_tot[4];
		       cum_tot[5] += this_tot[5];
		  
		       gotone = 1;
		    }
	       }
	  } /*j loop */
     } /* gotone loop */
   
   if (mergeable > 1) 
     {
	/* Nullify all of the merged vortices. */
	for (j=1; j<mergeable; ++j)
	  mblob[mergelist[j]].blob0.strength = 0.0;
	
	mblob[mergelist[0]].blob0.strength = cum_tot[0];
	mblob[mergelist[0]].blob0.x        = cum_tot[1]/cum_tot[0];
	mblob[mergelist[0]].blob0.y        = cum_tot[2]/cum_tot[0];
	
	cum_tot[1] /= cum_tot[0];
	cum_tot[2] /= cum_tot[0];
	
	cum_tot[3] = cum_tot[3]/cum_tot[0]-SQR(cum_tot[1]);
	cum_tot[4] = cum_tot[4]/cum_tot[0]-SQR(cum_tot[2]);
	cum_tot[5] = cum_tot[5]/cum_tot[0]-cum_tot[1]*cum_tot[2];
	     
	tmpa4 = 
	  (SQR(cum_tot[3])+SQR(cum_tot[4])+2.0*SQR(cum_tot[5])+
	   (cum_tot[3]+cum_tot[4])*
	   sqrt(SQR(cum_tot[3]-cum_tot[4])+4.0*SQR(cum_tot[5])))/
	  (2.0*(cum_tot[3]*cum_tot[4]-SQR(cum_tot[5])));
	
	tmpa2 = sqrt(tmpa4);
	
	tmps2 = 
	  (cum_tot[3]+cum_tot[4])*tmpa2/
	  (2.0*(tmpa4+1.0));
	
	blobguts[mergelist[0]].a2       = tmpa2;
	blobguts[mergelist[0]].s2       = tmps2;
	
	if (cum_tot[3]-2.0*tmps2/tmpa2 == 0.0)
	  blobguts[mergelist[0]].th = M_PI_2;
	else
	  blobguts[mergelist[0]].th = 
	  atan(cum_tot[5]/(cum_tot[3]-2.0*tmps2/tmpa2));
	
	mblob[mergelist[0]].blob0.order    = 1;
	set_blob(&(blobguts[mergelist[0]]),&(tmpparms[mergelist[0]]));
#ifdef XANTISYMM
       /* Double the number of computational elements by reflecting them
	* about the x-axis. */
       reflect_X();
#endif 
       dpos_vel(mergelist[0]);
#ifdef XANTISYMM 
       /* re-adjust */
       N /= 2;
#endif   
	++nmerge;
	
#ifdef MERGEDIAG
	fprintf(diag_log,"<%d> ",mergeable);
#endif
     }
   else
     cum_tot[0] = 0.0;

   /* If there is no merger then there is no induced error. */
   
   return(relerr3*fabs(cum_tot[0])/2.0/l2tol);
}

double clusterlist(ind,indmax,maxerr,searchlist)
int ind[NMAX],indmax;
double maxerr;
double (*searchlist)(int *, int, double,int *);
{
   int clusterlist[NMAX],clustersize,i,j,size;
   
   double scale1,scale2,dx,dy,dx_rescale,dy_rescale,budget;

   /* Some sort variables. */
   double clusterdist[NMAX],currR2,tmpR2;
   int k,l,currind;
   
   budget = 0.0;

   size = indmax;

   /* Pass 1. - Be cautious and spread the error around. */
   
   for (i=0; i<indmax-1; ++i)
     if (mblob[ind[i]].blob0.strength != 0.0)
     {	
	clusterlist[0] = ind[i];
	clustersize = 1;
	
	scale1 = 2.0*sqrt(blobguts[ind[i]].s2*blobguts[ind[i]].a2);
	scale2 = 2.0*sqrt(blobguts[ind[i]].s2/blobguts[ind[i]].a2);
	
	for (j=0; j<indmax; ++j)
	  if ((i != j) && 
	      (mblob[ind[j]].blob0.strength != 0.0) )
	  {
	     dx =  tmpparms[ind[i]].costh*
	       (mblob[ind[j]].blob0.x-mblob[ind[i]].blob0.x)+
	       tmpparms[ind[i]].sinth*
	       (mblob[ind[j]].blob0.y-mblob[ind[i]].blob0.y);
	     
	     dy = -tmpparms[ind[i]].sinth*
	       (mblob[ind[j]].blob0.x-mblob[ind[i]].blob0.x)+
	       tmpparms[ind[i]].costh*
	       (mblob[ind[j]].blob0.y-mblob[ind[i]].blob0.y);
	     
	     dx_rescale = dx/scale1;
	     dy_rescale = dy/scale2;
	     
	     tmpR2 = SQR(dx_rescale)+SQR(dy_rescale);
	     if ( tmpR2 < SQR(clusterR))
	       {
		  clusterlist[clustersize] = ind[j];
		  clusterdist[clustersize] = tmpR2;
		  clustersize++;
	       }
	  }
	if (clustersize > 1)
	  {
	    /* Insertion sort */
	    clusterdist[0] = 0.0;

	    for (k=1; k<clustersize; ++k)
	      {
		currR2 = clusterdist[k];
		currind = clusterlist[k];
		for (l=k-1; l>=0; --l)
		  {
		    if (clusterdist[l] <= currR2) 
		      break;
		    clusterdist[l+1] = clusterdist[l];
		    clusterlist[l+1] = clusterlist[l];
		  }
		clusterdist[l+1] = currR2;
		clusterlist[l+1] = currind;
	      }
	    /* Try to merge. */
	    /*
	    budget += (*searchlist)(clusterlist,clustersize,
				    (maxerr-budget)*clustersize/size,&size);
	    */
	    budget += (*searchlist)(clusterlist,clustersize,
				    (maxerr-budget)*clustersize/indmax,&size);
	}
     }

   /* Pass 2. - Be greedy. */
   
   for (i=0; i<indmax-1; ++i)
     if (mblob[ind[i]].blob0.strength != 0.0)
     {	
	clusterlist[0] = ind[i];
	clustersize = 1;
	
	scale1 = 2.0*sqrt(blobguts[ind[i]].s2*blobguts[ind[i]].a2);
	scale2 = 2.0*sqrt(blobguts[ind[i]].s2/blobguts[ind[i]].a2);
	
	for (j=0; j<indmax; ++j)
	  if ((i != j) && 
	      (mblob[ind[j]].blob0.strength != 0.0) )
	  {
	     dx =  tmpparms[ind[i]].costh*
	       (mblob[ind[j]].blob0.x-mblob[ind[i]].blob0.x)+
	       tmpparms[ind[i]].sinth*
	       (mblob[ind[j]].blob0.y-mblob[ind[i]].blob0.y);
	     
	     dy = -tmpparms[ind[i]].sinth*
	       (mblob[ind[j]].blob0.x-mblob[ind[i]].blob0.x)+
	       tmpparms[ind[i]].costh*
	       (mblob[ind[j]].blob0.y-mblob[ind[i]].blob0.y);
	     
	     dx_rescale = dx/scale1;
	     dy_rescale = dy/scale2;
	     
	     tmpR2 = SQR(dx_rescale)+SQR(dy_rescale);
	     if (tmpR2 < SQR(clusterR))
	       {
		  clusterlist[clustersize] = ind[j];
		  clusterdist[clustersize] = tmpR2;
		  clustersize++;
	       }
	  }
	if (clustersize > 1)
	    /* Insertion sort */
	    clusterdist[0] = 0.0;

	    for (k=1; k<clustersize; ++k)
	      {
		currR2 = clusterdist[k];
		currind = clusterlist[k];
		for (l=k-1; l>=0; --l)
		  {
		    if (clusterdist[l] <= currR2) 
		      break;
		    clusterdist[l+1] = clusterdist[l];
		    clusterlist[l+1] = clusterlist[l];
		  }
		clusterdist[l+1] = currR2;
		clusterlist[l+1] = currind;
	      }
	    /* Try to merge. */
	  budget += (*searchlist)(clusterlist,clustersize,
				  (maxerr-budget),&size);
     }

   return(budget);
}

double posneglists(double avail_budget,int k, int l, int size)
{
   int indp[NMAX],indpmax,indn[NMAX],indnmax;
   
   FineGridLink *trace;
   
   double cum_error;
   
   indpmax = 0;
   indnmax = 0;
   trace = FineGridLinks[k+l*size];
   
   while (trace != NULL)
     {
	if (mblob[trace->element].blob0.strength>0.0)
	  {
		  indp[indpmax] = trace->element;
		  ++indpmax;
	  }
	else
	  if (mblob[trace->element].blob0.strength<0.0)
	  {
	     indn[indnmax] = trace->element;
	     ++indnmax;
	  }
	trace = trace->next;
     }
	
   cum_error = 0.0;
   
   if (indpmax>1) 
     switch (merge_estimator)
       {
       case 1: 
	 cum_error = 
	   clusterlist(indp,indpmax,
		       avail_budget*indpmax/(indpmax+indnmax),
		       searchlist_distribution);
	 break;
       case 2: 
	   clusterlist(indp,indpmax,
		       avail_budget*indpmax/(indpmax+indnmax),
		       searchlist_uniform);
	 break;
       default:
	 printf("Error in MergeErrorEstimator.\n");
	 exit(-1);
       }
   if (indnmax>1) 
     switch (merge_estimator)
       {
       case 1: 
	 cum_error = 
	   clusterlist(indn,indnmax,
		       avail_budget*indnmax/(indpmax+indnmax),
		       searchlist_distribution);
	 break;
       case 2: 
	   clusterlist(indn,indnmax,
		       avail_budget*indnmax/(indpmax+indnmax),
		       searchlist_uniform);
	 break;
       default:
	 printf("Error in MergeErrorEstimator.\n");
	 exit(-1);
       }
   return(cum_error);
}

double calcaM2()
{
   int i;
   double temp;

   temp = 1.0;
   
   for (i=0; i<N; ++i)
     if (blobguts[i].a2 > temp) temp=blobguts[i].a2;
   
   return(temp);
}

void merge()
{
   int k,l,size;
   double cum_error,budget;
   
   rectify();
   
   aM2 = calcaM2();
   
   size = (int) (ldexp(1.0,mplevels));
   
   if (merge_budget != 0.0)
     budget = merge_budget;
   else
     budget = merge_c*exp(SimTime*merge_growth_rate);
   
#ifdef MERGEDIAG
   fprintf(diag_log,"Merging (time=%12.4e): ",SimTime);
#endif

   /*
   cum_error = posneglists(budget,8,8,size);

   printf("Cum_error %12.4e %d\n",cum_error,size);
   */

   for (k=0; k<size; ++k)
     for (l=0; l<size; ++l)
       cum_error = posneglists(budget,k,l,size);
   
#ifdef MERGEDIAG
   fprintf(diag_log,"\n");
   fflush(diag_log);
#endif
}

void resort()
{
    int  i,j,removed,tmpindex[NMAX];

    removed = 0;

    for (j=0; j<refinestack; ++j)
      tmpindex[j] = refineindex[j];

    for (i=0; i<N; ++i)
	{
	    if (removed !=0)
	     {
		mblob[i-removed] = mblob[i];
		blobguts[i-removed] = blobguts[i];
		tmpparms[i-removed] = tmpparms[i];
	     }
	    if (mblob[i].blob0.strength == 0.0)
	      {
		for (j=0; j<refinestack; ++j)
		  if (refineindex[j]>i) 
		    --tmpindex[j];
		++removed;
	      }
	}
    for (j=0; j<refinestack; ++j)
      refineindex[j] = tmpindex[j];
 
    N -= removed;
}
