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

#undef DEBUG_MP

#include "global.h"
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
		       /*
		       if (p==0)
			 printf("Err: %12.4e ",
				sqrt(SQR(z[p].re)+SQR(z[p].im)));
		       if (p==PMax-1)
			 printf("%12.4e\n",
				sqrt(SQR(z[p].re)+SQR(z[p].im)));
		       */
			
		       mblob[vort].blob0.dx += tmpz.re;
		       mblob[vort].blob0.dy -= tmpz.im;
		       
		       tmpz = cmult(z[p+1],*(Coeff_Array+(i+j*size)*PMax+p));
		       
		       tmpparms[vort].du11 += -(p+1)*tmpz.re;
		       tmpparms[vort].du12 -= -(p+1)*tmpz.im;
		       tmpparms[vort].du21 -= -(p+1)*tmpz.im;
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
		       
		       /*
		       if (p==0)
			 printf("Err: %12.4e ",
				sqrt(SQR(z[p].re)+SQR(z[p].im)));
		       if (p==PMax-1)
			 printf("%12.4e\n",
				sqrt(SQR(z[p].re)+SQR(z[p].im)));
		       */
		       
		       mblob[vort].blob0.dx += tmpz.re;
		       mblob[vort].blob0.dy -= tmpz.im;
		       
		       tmpz = cmult(z[p+1],*(Coeff_Array+(i+j*size)*PMax+p));
		       
		       tmpparms[vort].du11 += -(p+1)*tmpz.re;
		       tmpparms[vort].du12 -= -(p+1)*tmpz.im;
		       tmpparms[vort].du21 -= -(p+1)*tmpz.im;
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
		       
		       /*
		       if (p==0)
			 printf("Err: %12.4e ",
				sqrt(SQR(tmpz.re)+SQR(tmpz.im)));
		       if (p==PMax-1)
			 printf("%12.4e\n",
				sqrt(SQR(tmpz.re)+SQR(tmpz.im)));
			*/
		       
		       mblob[vort].blob0.dx += tmpz.re;
		       mblob[vort].blob0.dy -= tmpz.im;
		       
		       tmpz = cmult(z[p+1],*(Coeff_Array+(i+j*size)*PMax+p));
		       
		       tmpparms[vort].du11 += -(p+1)*tmpz.re;
		       tmpparms[vort].du12 -= -(p+1)*tmpz.im;
		       tmpparms[vort].du21 -= -(p+1)*tmpz.im;
		    } /* for p */
	       } /* if not nearest neighbor */
	  } /* for i,j */
     } /* if levels>1 */
}

void MP_Direct(int vort, int levels)
{
   int i,j,size,mini,maxi,minj,maxj;
   FineGridLink *trace;
   double dx,dy,eps;
   double r[5],t[5],even[5],odd[5];
   vector v;   
   tensor a;   
   
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
	     dx = mblob[vort].blob0.x-mblob[trace->element].blob0.x;
	     dy = mblob[vort].blob0.y-mblob[trace->element].blob0.y;

	     if ( (dx != 0.0) || (dy != 0.0) )
	       {

		   v = induced_vel(&(mblob[trace->element].blob0),
				  &(blobguts[trace->element]),
				  &(tmpparms[trace->element]),dx,dy,
				  r,t,even,odd);
		  mblob[vort].blob0.dx += v.x;
		  mblob[vort].blob0.dy += v.y;
		  
		  a = induced_veldev(&(mblob[trace->element].blob0),
				     &(blobguts[trace->element]),
				     &(tmpparms[trace->element]),dx,dy,
				     r,t,even,odd);
	     
		  tmpparms[vort].du11 += a.du11;
		  tmpparms[vort].du12 += a.du12;
		  tmpparms[vort].du21 += a.du21;
	       }
	     else
	       {
		  eps = (1.0-sqrt(blobguts[trace->element].a2))/
		    (1.0+sqrt(blobguts[trace->element].a2));
		  a.du11 = 0.0;
		  a.du12 = -(mblob[trace->element].blob0.strength/
			     (2.0*blobguts[trace->element].s2))*
		    (0.5+eps-eps*SQR(eps));
		  a.du21 = (mblob[trace->element].blob0.strength/
			    (2.0*blobguts[trace->element].s2))*
		    (0.5-eps+eps*SQR(eps));
		  
		  tmpparms[vort].du11 += (a.du11*
					  (tmpparms[trace->element].cos2-
					   tmpparms[trace->element].sin2)-
					  (a.du12+a.du21)*
					  tmpparms[trace->element].sincos);
		  
		  tmpparms[vort].du12 += (2.0*a.du11*
					  tmpparms[trace->element].sincos+
					  a.du12*
					  tmpparms[trace->element].cos2-
					  a.du21*
					  tmpparms[trace->element].sin2);
		  
		  tmpparms[vort].du21 += (2.0*a.du11*
					  tmpparms[trace->element].sincos-
					  a.du12*
					  tmpparms[trace->element].sin2+
					  a.du21*
					  tmpparms[trace->element].cos2);
		  
	       }
	     trace = trace->next;
	  }
     }
}

/* Uses a special finer grid for merging. */

void MP_Direct2(int vort, int levels)
{
   int i,j,size,mini,maxi,minj,maxj;
   FineGridLink *trace;
   double dx,dy,eps;
   double r[5],t[5],even[5],odd[5];
   vector v;
   tensor a;
   
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
	       dx = mblob[vort].blob0.x-mblob[trace->element].blob0.x;
	       dy = mblob[vort].blob0.y-mblob[trace->element].blob0.y;

	       if ( (dx != 0.0) || (dy != 0.0) )
		 {
		    v = induced_vel(&(mblob[trace->element].blob0),
				    &(blobguts[trace->element]),
				    &(tmpparms[trace->element]),dx,dy,
				    r,t,even,odd);
		    mblob[vort].blob0.dx += v.x;
		    mblob[vort].blob0.dy += v.y;
		    
		    a = induced_veldev(&(mblob[trace->element].blob0),
				       &(blobguts[trace->element]),
				       &(tmpparms[trace->element]),dx,dy,
				       r,t,even,odd);
		    
		    tmpparms[vort].du11 += a.du11;
		    tmpparms[vort].du12 += a.du12;
		    tmpparms[vort].du21 += a.du21;
		 }
	       else
		 {
		    eps = (1.0-sqrt(blobguts[trace->element].a2))/
		      (1.0+sqrt(blobguts[trace->element].a2));
		    a.du11 = 0.0;
		    a.du12 = -(mblob[trace->element].blob0.strength/
			       (2.0*blobguts[trace->element].s2))*
		      (0.5+eps-eps*SQR(eps));
		    a.du21 = (mblob[trace->element].blob0.strength/
			      (2.0*blobguts[trace->element].s2))*
		      (0.5-eps+eps*SQR(eps));
		    
		    tmpparms[vort].du11 += (a.du11*
					    (tmpparms[trace->element].cos2-
					     tmpparms[trace->element].sin2)-
					    (a.du12+a.du21)*
					    tmpparms[trace->element].sincos);
		    
		    tmpparms[vort].du12 += (2.0*a.du11*
					    tmpparms[trace->element].sincos+
					    a.du12*
					    tmpparms[trace->element].cos2-
					    a.du21*
					    tmpparms[trace->element].sin2);
		    
		    tmpparms[vort].du21 += (2.0*a.du11*
					    tmpparms[trace->element].sincos-
					    a.du12*
					    tmpparms[trace->element].sin2+
					    a.du21*
					    tmpparms[trace->element].cos2);
		 }
	       trace = trace->next;
	    }
       }
}

void MP_Direct3(int vort, int levels)
{
   int i,j,size,mini,maxi,minj,maxj;
   FineGridLink *trace;
   double dx,dy,eps;
   double r[5],t[5],even[5],odd[5];
   vector v;
   tensor a;
   
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

	       dx = mblob[vort].blob0.x-mblob[trace->element].blob0.x;
	       dy = mblob[vort].blob0.y-mblob[trace->element].blob0.y;
	       
	       if ( (dx != 0.0) || (dy != 0.0) )
		 {
		    v = induced_vel(&(mblob[trace->element].blob0),
				    &(blobguts[trace->element]),
				    &(tmpparms[trace->element]),dx,dy,
				    r,t,even,odd);
		    mblob[vort].blob0.dx += v.x;
		    mblob[vort].blob0.dy += v.y;
		    
		    a = induced_veldev(&(mblob[trace->element].blob0),
				       &(blobguts[trace->element]),
				       &(tmpparms[trace->element]),dx,dy,
				       r,t,even,odd);
		    
		    tmpparms[vort].du11 += a.du11;
		    tmpparms[vort].du12 += a.du12;
		    tmpparms[vort].du21 += a.du21;
		 }
	       else
		 {
		    eps = (1.0-sqrt(blobguts[trace->element].a2))/
		      (1.0+sqrt(blobguts[trace->element].a2));
		    a.du11 = 0.0;
		    a.du12 = -(mblob[trace->element].blob0.strength/
			       (2.0*blobguts[trace->element].s2))*
		      (0.5+eps-eps*SQR(eps));
		    a.du21 = (mblob[trace->element].blob0.strength/
			      (2.0*blobguts[trace->element].s2))*
		      (0.5-eps+eps*SQR(eps));
		    
		    tmpparms[vort].du11 += (a.du11*
					    (tmpparms[trace->element].cos2-
					     tmpparms[trace->element].sin2)-
					    (a.du12+a.du21)*
					    tmpparms[trace->element].sincos);
		    
		    tmpparms[vort].du12 += (2.0*a.du11*
					    tmpparms[trace->element].sincos+
					    a.du12*
					    tmpparms[trace->element].cos2-
					    a.du21*
					    tmpparms[trace->element].sin2);
		    
		    tmpparms[vort].du21 += (2.0*a.du11*
					    tmpparms[trace->element].sincos-
					    a.du12*
					    tmpparms[trace->element].sin2+
					    a.du21*
					    tmpparms[trace->element].cos2);
		 }
	       trace = trace->next;
	    }
       }
}
