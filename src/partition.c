/* PARTITION.C */
/* Copyright (c) 2000,2004 Louis F. Rossi                                    *
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

#include "global_min.h"
#include "particle.h"
#include "biot-savart.h"
#ifdef MULTIPROC
#include "multiproc.h"
#endif

int C(int n, int k)
{
   int i,sum;
   
   sum = 1;
   if (n != 0)
     {
	for (i=1; i<=n; ++i)
	  sum *= i;
	for (i=1; i<=k; ++i)
	  sum /= i;
	for (i=1; i<=n-k; ++i)
	  sum /= i;
     }
   return(sum);
}

Complex cpowi(Complex c,int n)
{
   int     i;
   Complex result;
   
   result = c;
   
   for (i=1; i<n; ++i)
     result = cmult(result,c);
 
   return(result);
}

Complex cmult(Complex c1,Complex c2)
{
   Complex result;
   
   result.re = c1.re*c2.re-c1.im*c2.im;
   result.im = c1.re*c2.im+c1.im*c2.re;
   
   return(result);
}

int Set_Level()
{
   int i;
   int level1,level2,level;
   double maxa2,maxs2,mingrid;
   double cX,cY,MGminX,MGmaxX,MGminY,MGmaxY;
   
   minX = mblob[0].blob0.x;
   maxX = mblob[0].blob0.x;
   minY = mblob[0].blob0.y;
   maxY = mblob[0].blob0.y;
   for (i=1; i<N; ++i)
     {
	if (mblob[i].blob0.x<minX) minX = mblob[i].blob0.x; 
	if (mblob[i].blob0.x>maxX) maxX = mblob[i].blob0.x; 
	if (mblob[i].blob0.y<minY) minY = mblob[i].blob0.y; 
	if (mblob[i].blob0.y>maxY) maxY = mblob[i].blob0.y; 
     }
   
   distX = maxX-minX;
   distY = maxY-minY;

   cX = (maxX+minX)/2.0;
   cY = (maxY+minY)/2.0;

   if (distY > distX)
     {
	minX -= (distY-distX)/2.0;
	maxX += (distY-distX)/2.0;
	distX = distY;
     }
   else
     {
	minY -= (distX-distY)/2.0;
	maxY += (distX-distY)/2.0;
	distY = distX;
     }
   
   maxa2 = 1.0;
   maxs2 = blobguts[0].s2;
   
   for (i=1; i<N; ++i)
     {
	if (1.0/blobguts[i].a2 > maxa2)
	  maxa2 = 1.0/blobguts[i].a2;
	else
	  if (blobguts[i].a2 > maxa2)
	    maxa2 = blobguts[i].a2;
	if (blobguts[i].s2 > maxs2)
	    maxs2 = blobguts[i].s2;
     }
   
   /* Minimum 5 elements per cell on average. */

   level1 = (int) ceil(log(N/5.0)/log(4.0));
   
   /* Using a finer mesh for merging */

   /*
   level2 = (int) ceil(log(distX/2.0/sqrt(maxs2*maxa2))/log(2.0));

   if (level1>level2)
     level = level2;
   else
     level = level1;
   
   if (level > LMAX) level = LMAX;
   if (level < 1) level = 1;
   */

   /* Choose two levels deeper than the minimum fine mesh size because 
      the mp-sum algorithm will piece together the sub-sized regions.
      This yields some significant savings.*/
   level2 = (int) (2+MAX(ceil(log(distX/8.0/sqrt(maxs2*maxa2))/log(2)),
			 ceil(log(distY/8.0/sqrt(maxs2*maxa2))/log(2))));

   if (level1>level2)
     level = level2;
   else
     level = level1;
   
   if (level > LMAX) level = LMAX;
   if (level < 1) level = 1;

   /* Calculate the smallest prudent grid cell. */
   mingrid = 8.0*sqrt(maxs2*maxa2);

   /* Since we are going two levels deeper, we need to reduce the size of 
      mingrid.*/
   mingrid /= 4.0;

   MGminX = cX - ldexp(1.0,level-1)*mingrid;
   MGmaxX = cX + ldexp(1.0,level-1)*mingrid;
   MGminY = cY - ldexp(1.0,level-1)*mingrid;
   MGmaxY = cY + ldexp(1.0,level-1)*mingrid;

   if (MGminX<minX) minX=MGminX;
   if (MGmaxX>maxX) maxX=MGmaxX;

   if (MGminY<minY) minY=MGminY;
   if (MGmaxY>maxY) maxY=MGmaxY;

   distX = distY = maxX-minX;

   return(level);
}

void partition(int levels)
{
   int i,l,size;
   FineGridLink *addone;
   
   /* Create a membership list for each vortex. */
   for (l=0; l<levels; ++l)
     {
	for (i=0; i<N; ++i)
	  {
	     gridx[l][i] = (int) 
	       ( (mblob[i].blob0.x-minX)*ldexp(1.0,l+1)/distX );
	     
	     gridy[l][i] = (int) 
	       ( (mblob[i].blob0.y-minY)*ldexp(1.0,l+1)/distY );
	     
	     /* Just in case the position is at the extreme. */
	     if (gridx[l][i] == (int) ldexp(1.0,l+1)) --gridx[l][i];
	     if (gridy[l][i] == (int) ldexp(1.0,l+1)) --gridy[l][i];
	  }
     }
   
   size = (int) ldexp(1.0,levels);
   
   FineGridLinks = malloc(sizeof(*addone)*SQR(size));
   
   for (i=0; i<SQR(size); ++i)
     FineGridLinks[i] = NULL;

   for (i=0; i<N; ++i)
     {
	addone = malloc(sizeof(*addone));
	addone->element = i;
	if ((FineGridLinks[gridx[levels-1][i]+size*gridy[levels-1][i]]) 
	    == NULL)
	  {
	    FineGridLinks[gridx[levels-1][i]+size*gridy[levels-1][i]] = 
	      addone;
	    addone->next = NULL;
	  }
	else
	  {   /* insert */
	    addone->next = 
	      FineGridLinks[gridx[levels-1][i]+size*gridy[levels-1][i]];
	    FineGridLinks[gridx[levels-1][i]+size*gridy[levels-1][i]]  = 
	      addone;
	  }
     }
}

/* This subroutine is parallelized crudely assuming a symmetric system.  
   The particles are divided evenly among available CPUs. Since the
   multipole coefficients are needed for all particles even with
   xantisymm. */

void Init_Fine_Grid(int levels)
{
   int i,j,p,size;
   Complex *Coeff_Array,tmpz,dz,mp[PMAX];

#ifdef MULTIPROC
   int start,end,rank,total_processes,buffsize;
   Complex *Coeff_buff;
#endif
   
   /* These need to be explicitly determined here for some reason.*/
   /* Check into this later.  Without explicitly finding rank and
      total_processes, they get wolloped for some reason */

#ifdef MULTIPROC
   MPI_Comm_size(MPI_COMM_WORLD, &total_processes);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

   Coeff_Array = Level_Ptr[levels-1];
   
   /* Wipe the finest level mesh */
   tmpz.re = 0.0; tmpz.im = 0.0;
   size = (int) ldexp(1.0,levels);
   
   for (i=0; i < size; ++i)
     for (j=0; j < size; ++j)
       for (p=0; p<PMAX; ++p)
	 *(Coeff_Array+(i+j*size)*PMAX+p) = tmpz;

#ifdef MULTIPROC
   MPI_Comm_size(MPI_COMM_WORLD, &total_processes);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   buffsize = PMAX*((int) ldexp(1.0,2*(levels)));

   Coeff_buff = malloc(sizeof(Complex)*buffsize);

   start = rank*((int) N/total_processes);
   end   = (rank+1)*((int) N/total_processes);
   
   if (rank==total_processes-1)
     end = N;

   for (i=start; i<end; ++i)
     {
	dz.re=(minX+(gridx[levels-1][i]+0.5)*(distX/size))-mblob[i].blob0.x;
	dz.im=(minY+(gridy[levels-1][i]+0.5)*(distY/size))-mblob[i].blob0.y;

	Calc_Coeffs(i,dz,mp);
	
	for (p=0; p<PMAX; ++p)
	  {
	     (*(Coeff_Array+
	       (gridx[levels-1][i]+gridy[levels-1][i]*size)*PMAX+p)).re +=
	       mp[p].re;
	     (*(Coeff_Array+
	       (gridx[levels-1][i]+gridy[levels-1][i]*size)*PMAX+p)).im +=
	       mp[p].im;
	  }
     }

   MPI_Allreduce(Coeff_Array,Coeff_buff,2*buffsize,MPI_DOUBLE,
	      MPI_SUM,MPI_COMM_WORLD);

   /*  This is the equivalent of an MPI_All_Reduce and should be
       removed if the line above works properly. */
   /*
   MPI_Reduce(Coeff_Array,Coeff_buff,2*buffsize,MPI_DOUBLE,
	      MPI_SUM,0,MPI_COMM_WORLD);

   MPI_Bcast (Coeff_buff, 2*buffsize, MPI_DOUBLE, 0, MPI_COMM_WORLD );
   */

   for (j=0; j<PMAX*size*size; ++j)
     *(Coeff_Array+j) = *(Coeff_buff+j);

   free(Coeff_buff);
#else
   for (i=0; i<N; ++i)
     {
	dz.re=(minX+(gridx[levels-1][i]+0.5)*(distX/size))-mblob[i].blob0.x;
	dz.im=(minY+(gridy[levels-1][i]+0.5)*(distY/size))-mblob[i].blob0.y;

	Calc_Coeffs(i,dz,mp);
	
	for (p=0; p<PMAX; ++p)
	  {
	     (*(Coeff_Array+
	       (gridx[levels-1][i]+gridy[levels-1][i]*size)*PMAX+p)).re +=
	       mp[p].re;
	     (*(Coeff_Array+
	       (gridx[levels-1][i]+gridy[levels-1][i]*size)*PMAX+p)).im +=
	       mp[p].im;
	  }
     }
#endif
}

/* Use coefficients at the finest level and pass them up to the coarsest 
 level */

void Advance_Coeffs(int levels)
{
   int i,j,p,p1,l,size;
   Complex *Coeff_Array,*Coeff_Array_Finer,tmpz,dz[PMAX+1];
   
#ifdef MULTIPROC
   int start,end,rank,total_processes,buffsize;
   Complex *Coeff_buff;
#endif
   
   /* These need to be explicitly determined here for some reason.*/
   /* Check into this later.  Without explicitly finding rank and
      total_processes, they get wolloped for some reason */

#ifdef MULTIPROC
   MPI_Comm_size(MPI_COMM_WORLD, &total_processes);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

   for (l=levels-2; l>=0; --l)
     {
#ifdef MULTIPROC
       buffsize = PMAX*((int) ldexp(1.0,2*(l+1)));
       Coeff_buff = malloc(sizeof(Complex)*buffsize);
#endif

	Coeff_Array = Level_Ptr[l];
	Coeff_Array_Finer = Level_Ptr[l+1];
	size = ldexp(1.0,l+1);
	
	/* Wipe the current layer */
	tmpz.re = 0.0; tmpz.im = 0.0;
	for (i=0; i<size; ++i)
	  for (j=0; j<size; ++j)
	    for (p=0; p<PMAX; ++p)
	      *(Coeff_Array+(i+j*size)*PMAX+p) = tmpz;

	dz[0].re = 1.0;
	dz[0].im = 0.0;
	
	
	/* Lower left children */
	
	dz[1].re = -distX/size/4.0;
	dz[1].im = -distY/size/4.0;
	
	for (p=2; p<=PMAX; ++p)
	  {
	     if (p % 2 == 0)
	       dz[p] = cmult(dz[p/2],dz[p/2]);
	     else
	       dz[p] = cmult(dz[p/2],dz[(p/2)+1]);
	  }

	for (i=0; i<size; ++i)
	  for (j=0; j<size; ++j)
#ifdef MULTIPROC
	    /* Split up the work. */
	    if ( ( (i*size+j) % total_processes) == rank)
	      {
#endif
	    for (p=0; p<PMAX; ++p)
	      for (p1=0; p1<=p; ++p1)
	  {
	     tmpz = cmult(dz[p-p1],
			  *(Coeff_Array_Finer+(2*i+2*j*2*size)*PMAX+p1));
	     (*(Coeff_Array+(i+j*size)*PMAX+p)).re +=
	       tmpz.re*C(p,p1);
	     (*(Coeff_Array+(i+j*size)*PMAX+p)).im +=
	       tmpz.im*C(p,p1);
	  }

#ifdef MULTIPROC
	      }
#endif
     	/* Lower right children */
	
	dz[1].re =  distX/size/4.0;
	dz[1].im = -distY/size/4.0;
	
	for (p=2; p<=PMAX; ++p)
	  {
	     if (p % 2 == 0)
	       dz[p] = cmult(dz[p/2],dz[p/2]);
	     else
	       dz[p] = cmult(dz[p/2],dz[(p/2)+1]);
	  }

	for (i=0; i<size; ++i)
	  for (j=0; j<size; ++j)
#ifdef MULTIPROC
	    /* Split up the work. */
	    if ( ( (i*size+j) % total_processes) == rank)
	      {
#endif
	    for (p=0; p<PMAX; ++p)
	      for (p1=0; p1<=p; ++p1)
	  {
	     tmpz = cmult(dz[p-p1],
			  *(Coeff_Array_Finer+(2*i+1+2*j*2*size)*PMAX+p1));
	     (*(Coeff_Array+(i+j*size)*PMAX+p)).re +=
	       tmpz.re*C(p,p1);
	     (*(Coeff_Array+(i+j*size)*PMAX+p)).im +=
	       tmpz.im*C(p,p1);
	  }
#ifdef MULTIPROC
	      }
#endif

     	/* Upper left children */
	
	dz[1].re = -distX/size/4.0;
	dz[1].im =  distY/size/4.0;
	
	for (p=2; p<=PMAX; ++p)
	  {
	     if (p % 2 == 0)
	       dz[p] = cmult(dz[p/2],dz[p/2]);
	     else
	       dz[p] = cmult(dz[p/2],dz[(p/2)+1]);
	  }

	for (i=0; i<size; ++i)
	  for (j=0; j<size; ++j)
#ifdef MULTIPROC
	    /* Split up the work. */
	    if ( ( (i*size+j) % total_processes) == rank)
	      {
#endif
	    for (p=0; p<PMAX; ++p)
	      for (p1=0; p1<=p; ++p1)
	  {
	     tmpz = cmult(dz[p-p1],
			  *(Coeff_Array_Finer+(2*i+(2*j+1)*2*size)*PMAX+p1));
	     (*(Coeff_Array+(i+j*size)*PMAX+p)).re +=
	       tmpz.re*C(p,p1);
	     (*(Coeff_Array+(i+j*size)*PMAX+p)).im +=
	       tmpz.im*C(p,p1);
	  }
#ifdef MULTIPROC
	      }
#endif

     	/* Upper right children */
	
	dz[1].re =  distX/size/4.0;
	dz[1].im =  distY/size/4.0;
	
	for (p=2; p<=PMAX; ++p)
	  {
	     if (p % 2 == 0)
	       dz[p] = cmult(dz[p/2],dz[p/2]);
	     else
	       dz[p] = cmult(dz[p/2],dz[(p/2)+1]);
	  }

	for (i=0; i<size; ++i)
	  for (j=0; j<size; ++j)
#ifdef MULTIPROC
	    /* Split up the work. */
	    if ( ( (i*size+j) % total_processes) == rank)
	      {
#endif
	    for (p=0; p<PMAX; ++p)
	      for (p1=0; p1<=p; ++p1)
	  {
	     tmpz = cmult(dz[p-p1],
			  *(Coeff_Array_Finer+(2*i+1+(2*j+1)*2*size)*PMAX+p1));
	     (*(Coeff_Array+(i+j*size)*PMAX+p)).re +=
	       tmpz.re*C(p,p1);
	     (*(Coeff_Array+(i+j*size)*PMAX+p)).im +=
	       tmpz.im*C(p,p1);
	  }
#ifdef MULTIPROC
	      }
#endif

#ifdef MULTIPROC
	MPI_Allreduce(Coeff_Array,Coeff_buff,2*buffsize,MPI_DOUBLE,
		      MPI_SUM,MPI_COMM_WORLD);
	
	for (j=0; j<PMAX*size*size; ++j)
	  *(Coeff_Array+j) = *(Coeff_buff+j);

	free(Coeff_buff);
#endif
     }
}

void Release_Links(int levels)
{
   int i,size;
   FineGridLink *takeone,*taketwo;
   
   size = ldexp(1.0,levels);
   
   for (i=0; i<SQR(size); ++i)
     {
	takeone = *(FineGridLinks+i);
	
	while (takeone != NULL)
	  {
	     taketwo = takeone;
	     takeone = takeone->next;
	     free(taketwo);
	  }
     }
   free(FineGridLinks);
}

void Create_Hierarchy()
{
   mp_cputime_ref = clock();
   mplevels = Set_Level();
   partition(mplevels);
   Init_Fine_Grid(mplevels);
   Advance_Coeffs(mplevels);
   mp_cputime += clock()-mp_cputime_ref;
}
