/* BOUNDARY.C */
/* Copyright (c) 2005 Louis F. Rossi                                    *
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
 * Newark, DE 19716-2553                                                */

#include "global.h"

#ifdef MULTIPROC
#include "multiproc.h"
#endif

double normtol = 1.0e-10;

double f23(a,b,c,x)
     double a,b,c,x;
{
  return((x+b/2/a)*log(a*x*x+b*x+c) + 
         sqrt(4*a*c-b*b)/a*atan((2*a*x+b)/sqrt(4*a*c-b*b)) - 2*x);
}

double f26(A,B,C,x)
     double A,B,C,x;
{
  return((x+A/B)*atan((A+B*x)/C)-C/2/B*log(1+SQR((B*x+A)/C)));
}

double f31(a,b,c,A,B,C,D,x,edgetol)
     double a,b,c,A,B,C,D,x,edgetol;
{
  double q,y;

  q = 4*a*c-SQR(b);
  
  if (fabs(C+D*x) <edgetol)
    y = -(B*C-A*D)/D*(C-D*b/2/a)*2/sqrt(q)*atan((2*a*x+b)/sqrt(q)) - 
      (B*C-A*D)/2/a*log(a*x*x+b*x+c);
  else
    y = (x+C/D)*atan((A+B*x)/(C+D*x)) - 
      (B*C-A*D)/D*(C-D*b/2/a)*2/sqrt(q)*atan((2*a*x+b)/sqrt(q)) - 
      (B*C-A*D)/2/a*log(a*x*x+b*x+c);

  return(y);
}

double intu(a,ba,ca,bb,cb,l1,edgetol)
     double a,ba,ca,bb,cb,l1,edgetol;
{
  double y;
  /* First half of u integral. */
  if ((fabs(a*SQR(l1/2)+ba*l1/2+ca) < SQR(edgetol)) || 
      (fabs(a*SQR(l1/2)-ba*l1/2+ca) < SQR(edgetol)))
    y = 2*l1*(log(sqrt(a)*l1)-1);
  else
    y = f23(a,ba,ca,l1/2) - f23(a,ba,ca,-l1/2);
  
  /*  Second half of u integral. */
  if ((fabs(a*SQR(l1/2)+bb*l1/2+cb) < SQR(edgetol)) || 
      (fabs(a*SQR(l1/2)-bb*l1/2+cb) < SQR(edgetol)))
    y = y - 2*l1*(log(sqrt(a)*l1)-1);
  else
    y = y - f23(a,bb,cb,l1/2) + f23(a,bb,cb,-l1/2);

  return(y);
}



double intv(a,ba,ca,bb,cb,Aa,Ab,B,C,D,l1,edgetol)
     double a,ba,ca,bb,cb,Aa,Ab,B,C,D,l1,edgetol;
{
  double y;

  if ( (fabs(D) < normtol) && (fabs(C) < edgetol) )
    y = 0;
  else
    {
      if (fabs(D) < normtol)
	y = f26(Aa,B,C,l1/2) - f26(Aa,B,C,-l1/2);
      else
	if (SQR(Aa*D-B*C) < SQR(edgetol))
	  y = l1*atan(Aa/C);
	else
	  y = f31(a,ba,ca,Aa,B,C,D,l1/2,edgetol) - 
	    f31(a,ba,ca,Aa,B,C,D,-l1/2,edgetol);

      if (fabs(D) < normtol)
	y = y - f26(Ab,B,C,l1/2) + f26(Ab,B,C,-l1/2);
      else
	if (SQR(Ab*D-B*C) < SQR(edgetol))
	  y = y - l1*atan(Ab/C);
	else
	    y = y - 
	      f31(a,bb,cb,Ab,B,C,D,l1/2,edgetol) + 
	      f31(a,bb,cb,Ab,B,C,D,-l1/2,edgetol);
    }

  return(y);
}


/* Calculate the effect of wall 0 on wall 1. */
double wallflux(x1,y1,n1x,n1y,l1,x0,y0,n0x,n0y,l0)
     double x1,y1,n1x,n1y,l1,x0,y0,n0x,n0y,l0;
{
  double edgetol,uflux,vflux,flux;
  double dx,dy,A,B,C,D,Aa,Ab,a,ba,ca,bb,cb;

  edgetol = 1.0e-4*MIN(l0,l1);

  dx = x1-x0;
  dy = y1-y0;

  A = dx*n0y  - dy*n0x;
  B = n0x*n1x + n0y*n1y;
  C = dx*n0x  + dy*n0y;
  D = n0x*n1y - n0y*n1x;
  
  Aa = A+l0/2;
  Ab = A-l0/2;
  
  a = B*B + D*D;
  ba = 2*(C*D+Aa*B);
  ca = Aa*Aa + C*C;
  bb = 2*(C*D+Ab*B);
  cb = Ab*Ab + C*C;

  uflux = intu(a,ba,ca,bb,cb,l1,edgetol);

  vflux = intv(a,ba,ca,bb,cb,Aa,Ab,B,C,D,l1,edgetol);

  flux = -D*uflux/4/M_PI+B*vflux/2/M_PI;
  
  return(flux);
}

/* Build matrix expressing the flux into the walls induced by other
   wall elements. */
void buildflux(walls,A,n)
     Panel  walls[BMAX];
     double A[BMAX][BMAX];
     int n;
{
  int k,j;

  for (k=0; k<n; ++k)
    for (j=0; j<n; ++j)
      if (k==j)
	A[j][k] = -walls[k].l/2.0;
      else
	A[j][k] = 
	  -wallflux(walls[k].x,walls[k].y,walls[k].nx,walls[k].ny,walls[k].l,
		    walls[j].x,walls[j].y,walls[j].nx,walls[j].ny,walls[j].l);
}

/* Compute the influence of element j on position (x,y). */
void vort_pos_interaction(j,x,y,u,v)
     int    j;
     double x,y,*u,*v;
{
   double dx,dy;
   double result[9];

   dx = x - mblob[j].blob0.x;
   dy = y - mblob[j].blob0.y;
	       
   /* if ( (dx != 0.0) || (dy != 0.0) ) */
   if ((SQR(dx)/blobguts[j].s2 >1.0e-2) || (SQR(dy)/blobguts[j].s2 > 1.0e-2))
     {
       induced_v(&(mblob[j].blob0),
		 &(blobguts[j]),
		 &(tmpparms[j]),dx,dy,
		 result);

       *u += result[0];
       *v += result[1];
     }
}

void vort_bdy_vel(double x, double y, double* u, double* v)
{
  int i;
  *u = 0.0;
  *v = 0.0;
  for (i=0; i<N; ++i)
    vort_pos_interaction(i,x,y,u,v);
}

/* A little subroutine to calculate the phase of the vortex. */

/*
void calc_phase(mxx,myy,mzz)
{
  int j;
  double *mxx,*mxy,*myy;

  *mxx = *myy = *mxy = 0.0;

  for (j=0; j<N; ++j)
    {
      *mxx += mblob0[j].strength*
	(2.0*blobguts[j].s2*(tmpparms[j].cos2*blobguts[j].a2+
			     tmpparms[j].sin2/blobguts[j].a2) +
	 SQR(mblob0[j].x));
      *myy += mblob0[j].strength*
	(2.0*blobguts[j].s2*(tmpparms[j].sin2*blobguts[j].a2+
			     tmpparms[j].cos2/blobguts[j].a2) +
	 SQR(mblob0[j].y));
      *mxy += mblob0[j].strength*
	(2.0*blobguts[j].s2*tmpparms[j].sincos*
	 (blobguts[j].a2-1.0/blobguts[j].a2) +
	 mblob0[j].x*mblob0[j].y);
    }
}
*/

/* Calculate the flux out of the walls from the fluid flow. */
void vort_panel_flux(wall,b,n,k)
     Panel  wall[BMAX];
     double b[BMAX];
     int    n,k;
{
  double u,v,tmpu,tmpv,x,y;
  int    j;
  /* Gauss interpolation points. 11th order.*/
  /*
  double g_x[] = {0.23861918,0.66120939,0.93246951};
  double g_w[] = {0.46791393,0.36076157,0.17132449};
  int    gpts = 3;
  */
  /* Gauss interpolation points. 7th order.*/
  double g_x[] = {0.33998104,0.86113631};
  double g_w[] = {0.65214515,0.34785485};
  int    gpts=2;

  tmpu = 0.0;
  tmpv = 0.0;

  for (j=0; j<gpts; ++j)
    {
      x = wall[k].x + wall[k].ny*g_x[j]*wall[k].l/2.0;
      y = wall[k].y - wall[k].nx*g_x[j]*wall[k].l/2.0;
	  
      vort_bdy_vel(x,y,&u,&v);
      tmpu += g_w[j]*u;
      tmpv += g_w[j]*v;

      x = wall[k].x - wall[k].ny*g_x[j]*wall[k].l/2.0;
      y = wall[k].y + wall[k].nx*g_x[j]*wall[k].l/2.0;
	  
      vort_bdy_vel(x,y,&u,&v);
      tmpu += g_w[j]*u;
      tmpv += g_w[j]*v;
    }
  tmpu *= wall[k].l/2.0;
  tmpv *= wall[k].l/2.0;

  b[k] = tmpu*wall[k].nx+tmpv*wall[k].ny;
      
  /* Code to include forcing for Travis' experiments. */
  /*
  int m=3;
  double phase;
  
  b[k] = += 
  */

}

#ifdef MULTIPROC  
void buildrhs_master(wall,b,n)
     Panel  wall[BMAX];
     double b[BMAX];
     int    n;
{
  int flag,k,proc=1,recvd;
  
  recvd = 0;
  while ( recvd < B )
    {
      MPI_Iprobe(proc,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,&mpistatus);
      if (flag==1)
	{
	  MPI_Recv(&k,1,MPI_INT,proc,
		   MPI_ANY_TAG,MPI_COMM_WORLD,&mpistatus);
	  MPI_Recv(&(b[k]),1,MPI_DOUBLE,proc,
		   MPI_ANY_TAG,MPI_COMM_WORLD,&mpistatus);
	  ++recvd;
	}

      ++proc;

      /* Consider adding some code to make proc 0 do some work too */
      if (proc==total_processes)
	{
	  proc = 1;
	}
    }
  MPI_Bcast (b, B, MPI_DOUBLE, 0, MPI_COMM_WORLD );
}

void buildrhs_slave(wall,b,n)
     Panel  wall[BMAX];
     double b[BMAX];
     int    n;
{
  int iters,leftover,worksize,k,index;

  iters    = B/(total_processes-1);
  leftover = B % (total_processes-1);

  for (k=0; k<iters; ++k)
    {
      index = k*(total_processes-1) + (rank-1);
      vort_panel_flux(wall,b,n,index);
      MPI_Send(&index,1,MPI_INT,0,0,MPI_COMM_WORLD);
      MPI_Send(&(b[index]),1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
    }
  if (rank-1 < leftover)
    {
      index = iters*(total_processes-1) + (rank-1);
      vort_panel_flux(wall,b,n,index);
      MPI_Send(&index,1,MPI_INT,0,0,MPI_COMM_WORLD);
      MPI_Send(&(b[index]),1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
    }
  MPI_Bcast (b, B, MPI_DOUBLE, 0, MPI_COMM_WORLD );
}
#endif

/* Calculate the flux out of the walls from the fluid flow. */
void buildrhs(wall,b,n)
     Panel  wall[BMAX];
     double b[BMAX];
     int    n;
{
  int    k;

#ifdef MULTIPROC  
  if (rank == 0) 
    buildrhs_master(wall,b,n);
  else 
    buildrhs_slave(wall,b,n);
  
  finish();
#else
  for (k=0; k<n; ++k)
    vort_panel_flux(wall,b,n,k);
#endif
}

void factor_bdy_matrix(walls,ipiv,A)
     Panel  walls[BMAX];
     double A[BMAX][BMAX];
     int    ipiv[BMAX];
{
#ifndef NOBOUNDARY
  int info,lda=BMAX,n=B;
#endif

  buildflux(walls,A,B);

#ifndef NOBOUNDARY
  dgetrf_(&n,&n,A,&lda,ipiv,&info);
#endif
}

void solve_bdy_matrix(walls,ipiv,A)
     Panel  walls[BMAX];
     double A[BMAX][BMAX];
     int    ipiv[BMAX];
{
  double b[BMAX];
  int k;
#ifndef NOBOUNDARY
  int    info,lda=BMAX,ldb=BMAX,n=B,nrhs=1;
  char   trans[]={'N'};
#endif

  buildrhs(walls,b,B);
  
  /* dgesv_(&n,&nrhs,A,&lda,ipiv,b,&ldb,&info); */
#ifndef NOBOUNDARY
  dgetrs_(trans,&B,&nrhs,A,&lda,ipiv,b,&ldb,&info);
#endif
  
  for (k=0; k<B; ++k)
    walls[k].m = b[k];
}     

void bdy_vel(x,y,u,v)
     double x, y, *u, *v;
{
  int k;
  double dx,dy,xt,yt,ut,vt;

  for (k=0; k<B; ++k)
    {
      dx = x-walls[k].x;
      dy = y-walls[k].y;
      xt = dx*walls[k].ny - dy*walls[k].nx;
      yt = dx*walls[k].nx + dy*walls[k].ny;
      ut = walls[k].m/4/M_PI*log((SQR(xt+walls[k].l/2) + SQR(yt))/
				 (SQR(xt-walls[k].l/2) + SQR(yt)));
      if (yt != 0.0)
	vt = walls[k].m/2/M_PI*
	  (atan((xt+walls[k].l/2)/yt)-atan((xt-walls[k].l/2)/yt));
      else
	vt = walls[k].m/4*(SGN(xt+walls[k].l/2)-SGN(xt-walls[k].l/2));
      *u +=  ut*walls[k].ny + vt*walls[k].nx;
      *v += -ut*walls[k].nx + vt*walls[k].ny; 
    }
}
