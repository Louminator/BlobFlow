/* MPSUM.C */
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
 * Louis Rossi                                                          *
 * Department of Mathematical Sciences                                  *
 * University of Delaware                                               *
 * Newark, DE 19716                                                     */

#undef DEBUG_MP


#include "global.h"

void vort_vort_interaction(i,j)
     int i,j;
{
   double dx,dy;
   tensor a;
   double result[9];

   dx = mblob[i].blob0.x-mblob[j].blob0.x;
   dy = mblob[i].blob0.y-mblob[j].blob0.y;
	       
   if ( (dx != 0.0) || (dy != 0.0) )
     {
       induced_v(&(mblob[j].blob0),
		 &(blobguts[j]),
		 &(tmpparms[j]),dx,dy,
		 result);

       mblob[i].blob0.dx += result[0];
       mblob[i].blob0.dy += result[1];

       tmpparms[i].du11 += result[2];
       tmpparms[i].du12 += result[3];
       tmpparms[i].du21 += result[4];

       tmpparms[i].u_xx += result[5];
       tmpparms[i].u_xy += result[6];
       tmpparms[i].u_yy += result[7];
       tmpparms[i].v_xx += result[8];
     }
   else
     {
       a.du11 = 0.0;
       a.du12 = -(mblob[j].blob0.strength/(2.0*blobguts[j].s2))/
	 (1.0 + 1.0/blobguts[j].a2);
       a.du21 = (mblob[j].blob0.strength/(2.0*blobguts[j].s2))/
	 (1.0 + blobguts[j].a2);
		    
       tmpparms[i].du11 += (a.du11*
			    (tmpparms[j].cos2-tmpparms[j].sin2)-
			    (a.du12+a.du21)*tmpparms[j].sincos);
		    
       tmpparms[i].du12 += (2.0*a.du11*tmpparms[j].sincos+
			    a.du12*tmpparms[j].cos2-
			    a.du21*tmpparms[j].sin2);
		    
       tmpparms[i].du21 += (2.0*a.du11*tmpparms[j].sincos-
			    a.du12*tmpparms[j].sin2+
			    a.du21*tmpparms[j].cos2);
     }
}

double maxC()
{
   int i;
   double max,tmp;
   
   max = 2.0*blobguts[0].s2*fabs(blobguts[0].a2-1.0/blobguts[0].a2);
   
   for (i=1; i<N; ++i)
     {
	tmp = 2.0*blobguts[i].s2*fabs(blobguts[i].a2-1.0/blobguts[i].a2);
	if (tmp>max) max = tmp;
     }
   return(max);
}

void MP_Sum(int vort, int levels)
{
   int i,j,p,l,size,mini,maxi,minj,maxj,Pallowed;
   complex *Coeff_Array,z[PMax+1],tmpz;
   double R2,mC;
   
   mC = maxC();
   
   for (l=1; l<levels-2; ++l)
     {
	size = ldexp(1.0,l+1);
	Coeff_Array = Level_Ptr[l];

	maxi = 2*(gridx[l-1][vort]+2);
	maxj = 2*(gridy[l-1][vort]+2);
	mini = 2*(gridx[l-1][vort]-1);
	minj = 2*(gridy[l-1][vort]-1);

	if (maxi>size) maxi=size;
	if (maxj>size) maxj=size;
	if (mini<0) mini = 0;
	if (minj<0) minj = 0;
	     
	for (i=mini; i<maxi; ++i)
	  for (j=minj; j<maxj; ++j)
	  {
	     /* Ignore nearest neighbors. */
	     
	     if ((i < gridx[l][vort]-1) || (i > gridx[l][vort]+1) ||
		 (j < gridy[l][vort]-1) || (j > gridy[l][vort]+1) )
	       {
#ifdef DEBUG_MP
		  printf("%d summing over level %d (%d, %d)\n",vort,l,i,j);
#endif
		  
		  if (mC == 0.0)
		    Pallowed = PMax;
		  else
		    {
		       R2 = SQR((i-gridx[l][vort])*distX/size)+
			 SQR((j-gridy[l][vort])*distY/size);
		       Pallowed = ((int) (R2/mC))-1;
		       /*printf("Pallowed: %d maxC: %12.4e\n",Pallowed,mC);*/
		       if (Pallowed>PMax) Pallowed=PMax;
		    }
		  
		  tmpz.re = (mblob[vort].blob0.x-(minX+(distX/size)*(i+0.5)));
		  tmpz.im = (mblob[vort].blob0.y-(minY+(distY/size)*(j+0.5)));
		  
		  z[0].re =  tmpz.re/(SQR(tmpz.re)+SQR(tmpz.im));
		  z[0].im = -tmpz.im/(SQR(tmpz.re)+SQR(tmpz.im));
		  
		  for (p=1; p<=Pallowed; ++p)
		    {
		       if (p % 2 == 1)
			 z[p] = cmult(z[p/2],z[p/2]);
		       else
			 z[p] = cmult(z[p/2],z[p/2-1]);
		    }
		  
		  for (p=0; p<Pallowed; ++p)
		    {
		       tmpz = cmult(z[p],*(Coeff_Array+(i+j*size)*PMax+p));
			
		       mblob[vort].blob0.dx += tmpz.re;
		       mblob[vort].blob0.dy -= tmpz.im;
		       
		       tmpz = cmult(z[p+1],*(Coeff_Array+(i+j*size)*PMax+p));
		       /* or */
		       /* tmpz = cmult(tmpz,z[0]); */
		       
		       tmpparms[vort].du11 += -(p+1)*tmpz.re;
		       tmpparms[vort].du12 -= -(p+1)*tmpz.im;
		       tmpparms[vort].du21 -= -(p+1)*tmpz.im;

		       tmpz = cmult(tmpz,z[0]);
		       
		       tmpparms[vort].u_xx += (p+1)*(p+2)*tmpz.re;
		       tmpparms[vort].u_xy -= (p+1)*(p+2)*tmpz.im;
		       tmpparms[vort].u_yy -= (p+1)*(p+2)*tmpz.re;
		       tmpparms[vort].v_xx -= (p+1)*(p+2)*tmpz.im;

		    } /* for p */
		  
	       } /* if not nearest neighbor */
	  } /* for i,j */
     } /* for l */

   if (levels > 2)
     {
	l = levels-2;

	size = ldexp(1.0,l+1);
	Coeff_Array = Level_Ptr[l];

	maxi = 2*(gridx[l-1][vort]+2);
	maxj = 2*(gridy[l-1][vort]+2);
	mini = 2*(gridx[l-1][vort]-1);
	minj = 2*(gridy[l-1][vort]-1);

	if (maxi>size) maxi=size;
	if (maxj>size) maxj=size;
	if (mini<0) mini = 0;
	if (minj<0) minj = 0;
	     
	for (i=mini; i<maxi; ++i)
	  for (j=minj; j<maxj; ++j)
	  {
	     /* Ignore nearest neighbors. */
	     
	     if ((i < gridx[l][vort]-2) || (i > gridx[l][vort]+2) ||
		 (j < gridy[l][vort]-2) || (j > gridy[l][vort]+2) )
	       {
#ifdef DEBUG_MP
		  printf("%d summing over level %d (%d, %d)\n",vort,l,i,j);
#endif
		  
		  if (mC == 0.0)
		    Pallowed = PMax;
		  else
		    {
		       R2 = SQR((i-gridx[l][vort])*distX/size)+
			 SQR((j-gridy[l][vort])*distY/size);
		       Pallowed = ((int) (R2/mC))-1;
		       /*printf("Pallowed: %d maxC: %12.4e\n",Pallowed,mC);*/
		       if (Pallowed>PMax) Pallowed=PMax;
		    }
		  
		  tmpz.re = (mblob[vort].blob0.x-(minX+(distX/size)*(i+0.5)));
		  tmpz.im = (mblob[vort].blob0.y-(minY+(distY/size)*(j+0.5)));
		  
		  z[0].re =  tmpz.re/(SQR(tmpz.re)+SQR(tmpz.im));
		  z[0].im = -tmpz.im/(SQR(tmpz.re)+SQR(tmpz.im));
		  
		  for (p=1; p<=Pallowed; ++p)
		    {
		       if (p % 2 == 1)
			 z[p] = cmult(z[p/2],z[p/2]);
		       else
			 z[p] = cmult(z[p/2],z[p/2-1]);
		    }
		  
		  for (p=0; p<Pallowed; ++p)
		    {
		       tmpz = cmult(z[p],*(Coeff_Array+(i+j*size)*PMax+p));
		       
		       mblob[vort].blob0.dx += tmpz.re;
		       mblob[vort].blob0.dy -= tmpz.im;
		       
		       tmpz = cmult(z[p+1],*(Coeff_Array+(i+j*size)*PMax+p));
		       /* or */
		       /* tmpz = cmult(tmpz,z[0]); */
		       
		       tmpparms[vort].du11 += -(p+1)*tmpz.re;
		       tmpparms[vort].du12 -= -(p+1)*tmpz.im;
		       tmpparms[vort].du21 -= -(p+1)*tmpz.im;

		       tmpz = cmult(tmpz,z[0]);
		       
		       tmpparms[vort].u_xx += (p+1)*(p+2)*tmpz.re;
		       tmpparms[vort].u_xy -= (p+1)*(p+2)*tmpz.im;
		       tmpparms[vort].u_yy -= (p+1)*(p+2)*tmpz.re;
		       tmpparms[vort].v_xx -= (p+1)*(p+2)*tmpz.im;
		    } /* for p */
	       } /* if not nearest neighbor */
	  } /* for i,j */
     } /* if levels>2 */

   if (levels > 1)
     {
	l = levels-1;

	size = ldexp(1.0,l+1);
	Coeff_Array = Level_Ptr[l];

	maxi = 2*(gridx[l-1][vort]+3);
	maxj = 2*(gridy[l-1][vort]+3);
	mini = 2*(gridx[l-1][vort]-2);
	minj = 2*(gridy[l-1][vort]-2);

	if (maxi>size) maxi=size;
	if (maxj>size) maxj=size;
	if (mini<0) mini = 0;
	if (minj<0) minj = 0;
	     
	for (i=mini; i<maxi; ++i)
	  for (j=minj; j<maxj; ++j)
	  {
	     /* Ignore nearest neighbors. */
	     
	     if ((i < gridx[l][vort]-4) || (i > gridx[l][vort]+4) ||
		 (j < gridy[l][vort]-4) || (j > gridy[l][vort]+4) )
	       {
#ifdef DEBUG_MP
		  printf("%d summing over level %d (%d, %d)\n",vort,l,i,j);
#endif
		  
		  if (mC == 0.0)
		    Pallowed = PMax;
		  else
		    {
		       R2 = SQR((i-gridx[l][vort])*distX/size)+
			 SQR((j-gridy[l][vort])*distY/size);
		       Pallowed = ((int) (R2/mC))-1;
		       /*printf("Pallowed: %d maxC: %12.4e\n",Pallowed,mC);*/
		       if (Pallowed>PMax) Pallowed=PMax;
		    }
		  
		  tmpz.re = (mblob[vort].blob0.x-(minX+(distX/size)*(i+0.5)));
		  tmpz.im = (mblob[vort].blob0.y-(minY+(distY/size)*(j+0.5)));
		  
		  z[0].re =  tmpz.re/(SQR(tmpz.re)+SQR(tmpz.im));
		  z[0].im = -tmpz.im/(SQR(tmpz.re)+SQR(tmpz.im));
		  
		  for (p=1; p<=Pallowed; ++p)
		    {
		       if (p % 2 == 1)
			 z[p] = cmult(z[p/2],z[p/2]);
		       else
			 z[p] = cmult(z[p/2],z[p/2-1]);
		    }
		  
		  for (p=0; p<Pallowed; ++p)
		    {
		       tmpz = cmult(z[p],*(Coeff_Array+(i+j*size)*PMax+p));
		       
		       mblob[vort].blob0.dx += tmpz.re;
		       mblob[vort].blob0.dy -= tmpz.im;
		       
		       tmpz = cmult(z[p+1],*(Coeff_Array+(i+j*size)*PMax+p));
		       
		       tmpparms[vort].du11 += -(p+1)*tmpz.re;
		       tmpparms[vort].du12 -= -(p+1)*tmpz.im;
		       tmpparms[vort].du21 -= -(p+1)*tmpz.im;

		       tmpz = cmult(tmpz,z[0]);
		       
		       tmpparms[vort].u_xx += (p+1)*(p+2)*tmpz.re;
		       tmpparms[vort].u_xy -= (p+1)*(p+2)*tmpz.im;
		       tmpparms[vort].u_yy -= (p+1)*(p+2)*tmpz.re;
		       tmpparms[vort].v_xx -= (p+1)*(p+2)*tmpz.im;
		    } /* for p */
	       } /* if not nearest neighbor */
	  } /* for i,j */
     } /* if levels>1 */
}

void MP_Direct(int vort, int levels)
{
   int i,j,size,mini,maxi,minj,maxj;
   FineGridLink *trace;
   
   /* Now use direct summation at the finest level.  (Boo hoo) */
   
   size = (int) ldexp(1.0,levels);
   maxi = gridx[levels-1][vort]+2;
   maxj = gridy[levels-1][vort]+2;
   mini = gridx[levels-1][vort]-1;
   minj = gridy[levels-1][vort]-1;
   if (maxi>size) maxi=size;
   if (maxj>size) maxj=size;
   if (mini<0) mini = 0;
   if (minj<0) minj = 0;
  

   for (i=mini; i<maxi; ++i)
     for (j=minj; j<maxj; ++j)
     {
	trace = FineGridLinks[i+j*size];

	while (trace != NULL)
	  {
	    vort_vort_interaction(vort,trace->element);
	    trace = trace->next;
	  }
     }
}

/* Uses a special finer grid for merging. */

void MP_Direct2(int vort, int levels)
{
   int i,j,size,mini,maxi,minj,maxj;
   FineGridLink *trace;
   
   /* Now use direct summation at the finest level.  (Boo hoo) */
   
   size = (int) ldexp(1.0,levels);
   
   maxi = 4*(gridx[levels-1][vort]/4 + 2);
   maxj = 4*(gridy[levels-1][vort]/4 + 2);
   mini = 4*(gridx[levels-1][vort]/4 - 1);
   minj = 4*(gridy[levels-1][vort]/4 - 1);
   
   if (maxi>size) maxi=size;
   if (maxj>size) maxj=size;
   if (mini<0) mini = 0;
   if (minj<0) minj = 0;
  
   for (i=mini; i<maxi; ++i)
     for (j=minj; j<maxj; ++j)
       {
	  trace = FineGridLinks[i+j*size];

	  while (trace != NULL)
	    {	    
	      vort_vort_interaction(vort,trace->element);
	      trace = trace->next;
	    }
       }
}

void MP_Direct3(int vort, int levels)
{
   int i,j,size,mini,maxi,minj,maxj;
   FineGridLink *trace;
   
   /* Now use direct summation at the finest level.  (Boo hoo) */
   
   size = (int) ldexp(1.0,levels);
   
   maxi = gridx[levels-1][vort] + 5;
   maxj = gridy[levels-1][vort] + 5;
   mini = gridx[levels-1][vort] - 4;
   minj = gridy[levels-1][vort] - 4;
   
   if (maxi>size) maxi=size;
   if (maxj>size) maxj=size;
   if (mini<0) mini = 0;
   if (minj<0) minj = 0;
  
   for (i=mini; i<maxi; ++i)
     for (j=minj; j<maxj; ++j)
       {
#ifdef DEBUG_MP
	  printf("%d DIRECT summing over %d (%d, %d)\n",vort,levels-1,i,j);
#endif

	  trace = FineGridLinks[i+j*size];

	  while (trace != NULL)
	    {
	      ++numk2;

	      vort_vort_interaction(vort,trace->element);
	      trace = trace->next;
	    }
       }
}
