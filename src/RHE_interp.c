/* RHE_interp.c */
/* Copyright (c) 2009 Louis F. Rossi                                    *
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

#include "global.h"

/* Interpolate blobs onto the finest multipole grid with a gridn X
   gridn regular grid. */

void interp_fine_grid (int FineGridi, int FineGridj, int levels,
		       int gridn, double *field)
{
  int i,j,k,l,size,mini,maxi,minj,maxj;
  FineGridLink *trace;
  
  double h,gridsize;
  double posx,posy,circ,c,s,cos2,sin2,var,a2;
  double x,y,dx,dy;
  double temp;

  /* This assumes that the finest grid levels are >= 2 sigma a. */
  /* gridn^2 is the number of particles created per fine grid cell. */
  /* s2 is the core size of the new particles.  */

  size = (int) ldexp(1.0,levels);
   
  gridsize = (maxX-minX)/size;
  h = gridsize/gridn;

  maxi = FineGridi + 5;
  maxj = FineGridj + 5;
  mini = FineGridi - 4;
  minj = FineGridj - 4;
  
  if (maxi>size) maxi=size;
  if (maxj>size) maxj=size;
  if (mini<0)    mini = 0;
  if (minj<0)    minj = 0;

  for (k=0; k<gridn; ++k)
    for (l=0; l<gridn; ++l)
      {
	x = minX+FineGridi*gridsize+(k+0.5)*h;
	y = minY+FineGridj*gridsize+(l+0.5)*h;
	
	*(field+k+gridn*l) = 0.0;
	for (i=mini; i<maxi; ++i)
	  for (j=minj; j<maxj; ++j)
	    {
	      trace = FineGridLinks[i+j*size];
	      
	      while (trace != NULL)
		{
		  posx = mblob[trace->element].blob0.x;
		  posy = mblob[trace->element].blob0.y;
		  circ = mblob[trace->element].blob0.strength;
		  c    = tmpparms[trace->element].costh;
		  s    = tmpparms[trace->element].sinth;
		  cos2 = tmpparms[trace->element].cos2;
		  sin2 = tmpparms[trace->element].sin2;
		  a2   = blobguts[trace->element].a2;
		  var  = blobguts[trace->element].s2;
		  
		  dx = x - posx;
		  dy = y - posy;
		  
		  /* Evaluation of omega here. */
		  
		  temp = -(SQR(dx)*(cos2/a2+sin2*a2)+
			   SQR(dy)*(sin2/a2+cos2*a2)+
			   2.0*dx*dy*c*s*(1.0/a2-a2))/
		    (4.0*var);
		  /* A little code to avoid underflow handling. */
		  if (-temp<16)
		    *(field+k+gridn*l) += 
		      circ*exp(temp)/(2.0*var);
		  
		  trace = trace->next;
		}
	    }
	*(field+k+gridn*l) *= h*h;
      }
}

double smart_read(int i, int j, int gridn, 
		  double *left,double *center,double *right,
		  double *down, double *up)
{
  double *ptr;
  int flag=0;

  ptr = center;

  if (i<0)
    {
      ptr = left;
      i += gridn;
      flag = 1;
    }
  if (i>=gridn)
    {
      ptr = right;
      i -= gridn;
      flag = 1;
    }

  if (j<0)
    {
      ptr = down;
      j += gridn;
      flag = 1;
    }
  if (j>=gridn)
    {
      ptr = up;
      j -= gridn;
      flag = 1;
    }

  if (ptr != NULL)
    return(*(ptr+i+gridn*j));
  else
    return(0.0);
}

void fill_vec(double v[13], int k, int l, int gridn, 
	      double *left, double *center, double *right,
	      double *down, double *up)
{
  v[0]  = smart_read(k-3,l,  gridn,left,center,right,down,up);
  v[1]  = smart_read(k-2,l,  gridn,left,center,right,down,up);
  v[2]  = smart_read(k-1,l,  gridn,left,center,right,down,up);
  v[3]  = smart_read(k,  l,  gridn,left,center,right,down,up);
  v[4]  = smart_read(k+1,l,  gridn,left,center,right,down,up);
  v[5]  = smart_read(k+2,l,  gridn,left,center,right,down,up);
  v[6]  = smart_read(k+3,l,  gridn,left,center,right,down,up);
  v[7]  = smart_read(k,  l-3,gridn,left,center,right,down,up);
  v[8]  = smart_read(k,  l-2,gridn,left,center,right,down,up);
  v[9]  = smart_read(k,  l-1,gridn,left,center,right,down,up);
  v[10] = smart_read(k,  l+1,gridn,left,center,right,down,up);
  v[11] = smart_read(k,  l+2,gridn,left,center,right,down,up);
  v[12] = smart_read(k,  l+3,gridn,left,center,right,down,up);
}

double RHE_RK4(double **circs, int gridn, double s2)

{
  int size,i,j,k,l,levels;
  double gridsize,h;
  double **g1,**g2,**g3,**k1,**k2,**k3;
  double *field1,*field2,*field3,*field4;
  double *circ,v[13];
  double alpha;

  /* CAUTION: right now the code only works at the finest level. */
  levels = mplevels;  
  size = (int) ldexp(1.0,levels);

  g1    = malloc(sizeof(field1)*SQR(size));
  if (!g1)
    fprintf(diag_log,"g2 malloc failed.\n");
  g2    = malloc(sizeof(field1)*SQR(size));
  if (!g2)
    fprintf(diag_log,"g2 malloc failed.\n");
  g3    = malloc(sizeof(field1)*SQR(size));
  if (!g3)
    fprintf(diag_log,"g3 malloc failed.\n");
  k1    = malloc(sizeof(field1)*SQR(size));
  if (!k1)
    fprintf(diag_log,"k1 malloc failed.\n");
  k2    = malloc(sizeof(field1)*SQR(size));
  if (!k2)
    fprintf(diag_log,"k2 malloc failed.\n");
  k3    = malloc(sizeof(field1)*SQR(size));
  if (!k3)
    fprintf(diag_log,"k3 malloc failed.\n");

  gridsize = (maxX-minX)/size;

  h = gridsize/gridn;

  alpha = s2/h/h;

  /* Initialize the circulations. */

  for (i=0; i<SQR(size); ++i)
    {
      circs[i] = NULL;
      g1[i]    = NULL;
      g2[i]    = NULL;
      g3[i]    = NULL;
      k1[i]    = NULL;
      k2[i]    = NULL;
      k3[i]    = NULL;
    }

  /* Walk through the fine grid and grab the circs. */

  for (j=0; j<size; ++j)
    for (i=0; i<size; ++i)
      {
	/* You only need to create a new field if you want to keep the old one. */
	field1 = malloc(sizeof(h)*SQR(gridn));
	if (!field1)
	  fprintf(diag_log,"field1 malloc failed.\n");
	interp_fine_grid(i, j, levels, gridn, field1);
	circs[j*size+i] = field1;
      }

  /* v_1-6   = v_i-3,j...v_i+3,j */
  /* v_7-9   = v_i,j-3...v_i,j-1 */
  /* v_10-12 = v_i,j+1...v_i,j+3 */

  /* Half step FCE */

  for (j=0; j<size; ++j)
    for (i=0; i<size; ++i)
      {
	circ = circs[j*size+i];
	field1 = malloc(sizeof(h)*SQR(gridn));
	field2 = malloc(sizeof(h)*SQR(gridn));
	if (!field1)
	  fprintf(diag_log,"field1 malloc failed.\n");
	if (!field2)
	  fprintf(diag_log,"field2 malloc failed.\n");
	/* printf("%d %d\n",i,j); */
	for (l=0; l<gridn; ++l)
	  for (k=0; k<gridn; ++k)
	    {
	      if ((i==0) && (j == 0)) /* Bottom left*/
		fill_vec(v,k,l,gridn,
			 NULL,
			 circs[j*size+i],
			 circs[j*size+i+1],
			 NULL,
			 circs[(j+1)*size+i]);
	      else if ((i == size-1) && (j == 0)) /* Bottom right*/
		fill_vec(v,k,l,gridn,
			 circs[j*size+i-1],
			 circs[j*size+i],
			 NULL,
			 NULL,
			 circs[(j+1)*size+i]);
	      else if ((i==0) && (j == size-1)) /* Top left*/
		fill_vec(v,k,l,gridn,
			 NULL,
			 circs[j*size+i],
			 circs[j*size+i+1],
			 circs[(j-1)*size+i],
			 NULL);
	      else if ((i == size-1) && (j == size-1)) /* Top right*/
		fill_vec(v,k,l,gridn,
			 circs[j*size+i-1],
			 circs[j*size+i],
			 NULL,
			 circs[(j-1)*size+i],
			 NULL);
	      else if (j == 0) /* Bottom */
		fill_vec(v,k,l,gridn,
			 circs[j*size+i-1],
			 circs[j*size+i],
			 circs[j*size+i+1],
			 NULL,
			 circs[(j+1)*size+i]);
	      else if (j == size-1) /* Top */
		fill_vec(v,k,l,gridn,
			 circs[j*size+i-1],
			 circs[j*size+i],
			 circs[j*size+i+1],
			 circs[(j-1)*size+i],
			 NULL);
	      else if (i == 0) /* Left */
		fill_vec(v,k,l,gridn,
			 NULL,
			 circs[j*size+i],
			 circs[j*size+i+1],
			 circs[(j-1)*size+i],
			 circs[(j+1)*size+i]);
	      else if (i == size-1) /* Right */
		fill_vec(v,k,l,gridn,
			 circs[j*size+i-1],
			 circs[j*size+i],
			 NULL,
			 circs[(j-1)*size+i],
			 circs[(j+1)*size+i]);
	      else
		fill_vec(v,k,l,gridn,
			 circs[j*size+i-1],
			 circs[j*size+i],
			 circs[j*size+i+1],
			 circs[(j-1)*size+i],
			 circs[(j+1)*size+i]);

	      *(field1+k+gridn*l) = 
		-alpha/180*
		(2*v[0] - 27*v[1] + 
		 270*v[2] - 980*v[3] + 270*v[4] -27*v[5] + 
		 2*v[6] + 2*v[7] - 27*v[8] + 270*v[9] + 
		 270*v[10] - 27*v[11] + 2*v[12]);

	      *(field2+k+gridn*l) = *(circ+k+gridn*l) + 
		0.5*(*(field1+k+gridn*l));

	    }

	k1[j*size+i] = field1;
	g1[j*size+i] = field2;
      }

  /* Half step BCE */

  for (j=0; j<size; ++j)
    for (i=0; i<size; ++i)
      {
	circ = circs[j*size+i];
	field1 = malloc(sizeof(h)*SQR(gridn));
	field2 = malloc(sizeof(h)*SQR(gridn));
	if (!field1)
	  fprintf(diag_log,"field1 malloc failed.\n");
	if (!field2)
	  fprintf(diag_log,"field2 malloc failed.\n");
	for (l=0; l<gridn; ++l)
	  for (k=0; k<gridn; ++k)
	    {
	      if ((i==0) && (j == 0)) /* Bottom left*/
		fill_vec(v,k,l,gridn,
			 NULL,
			 g1[j*size+i],
			 g1[j*size+i+1],
			 NULL,
			 g1[(j+1)*size+i]);
	      else if ((i == size-1) && (j == 0)) /* Bottom right*/
		fill_vec(v,k,l,gridn,
			 g1[j*size+i-1],
			 g1[j*size+i],
			 NULL,
			 NULL,
			 g1[(j+1)*size+i]);
	      else if ((i==0) && (j == size-1)) /* Top left*/
		fill_vec(v,k,l,gridn,
			 NULL,
			 g1[j*size+i],
			 g1[j*size+i+1],
			 g1[(j-1)*size+i],
			 NULL);
	      else if ((i == size-1) && (j == size-1)) /* Top right*/
		fill_vec(v,k,l,gridn,
			 g1[j*size+i-1],
			 g1[j*size+i],
			 NULL,
			 g1[(j-1)*size+i],
			 NULL);
	      else if (j == 0) /* Bottom */
		fill_vec(v,k,l,gridn,
			 g1[j*size+i-1],
			 g1[j*size+i],
			 g1[j*size+i+1],
			 NULL,
			 g1[(j+1)*size+i]);
	      else if (j == size-1) /* Top */
		fill_vec(v,k,l,gridn,
			 g1[j*size+i-1],
			 g1[j*size+i],
			 g1[j*size+i+1],
			 g1[(j-1)*size+i],
			 NULL);
	      else if (i == 0) /* Left */
		fill_vec(v,k,l,gridn,
			 NULL,
			 g1[j*size+i],
			 g1[j*size+i+1],
			 g1[(j-1)*size+i],
			 g1[(j+1)*size+i]);
	      else if (i == size-1) /* Right */
		fill_vec(v,k,l,gridn,
			 g1[j*size+i-1],
			 g1[j*size+i],
			 NULL,
			 g1[(j-1)*size+i],
			 g1[(j+1)*size+i]);
	      else
		fill_vec(v,k,l,gridn,
			 g1[j*size+i-1],
			 g1[j*size+i],
			 g1[j*size+i+1],
			 g1[(j-1)*size+i],
			 g1[(j+1)*size+i]);

	      *(field1+k+gridn*l) = 
		-alpha/180*
		(2*v[0] - 27*v[1] + 
		 270*v[2] - 980*v[3] + 270*v[4] -27*v[5] + 
		 2*v[6] + 2*v[7] - 27*v[8] + 270*v[9] + 
		 270*v[10] - 27*v[11] + 2*v[12]);

	      *(field2+k+gridn*l) = *(circ+k+gridn*l) + 
		0.5*(*(field1+k+gridn*l));

	    }
	k2[j*size+i] = field1;
	g2[j*size+i] = field2;
      }

  /* Midpoint */

  for (j=0; j<size; ++j)
    for (i=0; i<size; ++i)
      {
	circ = circs[j*size+i];
	field1 = malloc(sizeof(h)*SQR(gridn));
	field2 = malloc(sizeof(h)*SQR(gridn));
	if (!field1)
	  fprintf(diag_log,"field1 malloc failed.\n");
	if (!field2)
	  fprintf(diag_log,"field2 malloc failed.\n");
	for (l=0; l<gridn; ++l)
	  for (k=0; k<gridn; ++k)
	    {
	      if ((i==0) && (j == 0)) /* Bottom left*/
		fill_vec(v,k,l,gridn,
			 NULL,
			 g2[j*size+i],
			 g2[j*size+i+1],
			 NULL,
			 g2[(j+1)*size+i]);
	      else if ((i == size-1) && (j == 0)) /* Bottom right*/
		fill_vec(v,k,l,gridn,
			 g2[j*size+i-1],
			 g2[j*size+i],
			 NULL,
			 NULL,
			 g2[(j+1)*size+i]);
	      else if ((i==0) && (j == size-1)) /* Top left*/
		fill_vec(v,k,l,gridn,
			 NULL,
			 g2[j*size+i],
			 g2[j*size+i+1],
			 g2[(j-1)*size+i],
			 NULL);
	      else if ((i == size-1) && (j == size-1)) /* Top right*/
		fill_vec(v,k,l,gridn,
			 g2[j*size+i-1],
			 g2[j*size+i],
			 NULL,
			 g2[(j-1)*size+i],
			 NULL);
	      else if (j == 0) /* Bottom */
		fill_vec(v,k,l,gridn,
			 g2[j*size+i-1],
			 g2[j*size+i],
			 g2[j*size+i+1],
			 NULL,
			 g2[(j+1)*size+i]);
	      else if (j == size-1) /* Top */
		fill_vec(v,k,l,gridn,
			 g2[j*size+i-1],
			 g2[j*size+i],
			 g2[j*size+i+1],
			 g2[(j-1)*size+i],
			 NULL);
	      else if (i == 0) /* Left */
		fill_vec(v,k,l,gridn,
			 NULL,
			 g2[j*size+i],
			 g2[j*size+i+1],
			 g2[(j-1)*size+i],
			 g2[(j+1)*size+i]);
	      else if (i == size-1) /* Right */
		fill_vec(v,k,l,gridn,
			 g2[j*size+i-1],
			 g2[j*size+i],
			 NULL,
			 g2[(j-1)*size+i],
			 g2[(j+1)*size+i]);
	      else
		fill_vec(v,k,l,gridn,
			 g2[j*size+i-1],
			 g2[j*size+i],
			 g2[j*size+i+1],
			 g2[(j-1)*size+i],
			 g2[(j+1)*size+i]);


	      *(field1+k+gridn*l) = 
		-alpha/180*
		(2*v[0] - 27*v[1] + 
		 270*v[2] - 980*v[3] + 270*v[4] -27*v[5] + 
		 2*v[6] + 2*v[7] - 27*v[8] + 270*v[9] + 
		 270*v[10] - 27*v[11] + 2*v[12]);

	      *(field2+k+gridn*l) = *(circ+k+gridn*l) + 
		(*(field1+k+gridn*l));
	    }
	k3[j*size+i] = field1;
	g3[j*size+i] = field2;
      }

  /* Simpsons corrector */

  for (j=0; j<size; ++j)
    for (i=0; i<size; ++i)
      {
	circ = circs[j*size+i];
	field4 = malloc(sizeof(h)*SQR(gridn));
	if (!field4)
	  fprintf(diag_log,"field4 malloc failed.\n");
	field1 = k1[j*size+i];
	field2 = k2[j*size+i];
	field3 = k3[j*size+i];

	for (l=0; l<gridn; ++l)
	  for (k=0; k<gridn; ++k)
	    {
	      if ((i==0) && (j == 0)) /* Bottom left*/
		fill_vec(v,k,l,gridn,
			 NULL,
			 g3[j*size+i],
			 g3[j*size+i+1],
			 NULL,
			 g3[(j+1)*size+i]);
	      else if ((i == size-1) && (j == 0)) /* Bottom right*/
		fill_vec(v,k,l,gridn,
			 g3[j*size+i-1],
			 g3[j*size+i],
			 NULL,
			 NULL,
			 g3[(j+1)*size+i]);
	      else if ((i==0) && (j == size-1)) /* Top left*/
		fill_vec(v,k,l,gridn,
			 NULL,
			 g3[j*size+i],
			 g3[j*size+i+1],
			 g3[(j-1)*size+i],
			 NULL);
	      else if ((i == size-1) && (j == size-1)) /* Top right*/
		fill_vec(v,k,l,gridn,
			 g3[j*size+i-1],
			 g3[j*size+i],
			 NULL,
			 g3[(j-1)*size+i],
			 NULL);
	      else if (j == 0) /* Bottom */
		fill_vec(v,k,l,gridn,
			 g3[j*size+i-1],
			 g3[j*size+i],
			 g3[j*size+i+1],
			 NULL,
			 g3[(j+1)*size+i]);
	      else if (j == size-1) /* Top */
		fill_vec(v,k,l,gridn,
			 g3[j*size+i-1],
			 g3[j*size+i],
			 g3[j*size+i+1],
			 g3[(j-1)*size+i],
			 NULL);
	      else if (i == 0) /* Left */
		fill_vec(v,k,l,gridn,
			 NULL,
			 g3[j*size+i],
			 g3[j*size+i+1],
			 g3[(j-1)*size+i],
			 g3[(j+1)*size+i]);
	      else if (i == size-1) /* Right */
		fill_vec(v,k,l,gridn,
			 g3[j*size+i-1],
			 g3[j*size+i],
			 NULL,
			 g3[(j-1)*size+i],
			 g3[(j+1)*size+i]);
	      else
		fill_vec(v,k,l,gridn,
			 g3[j*size+i-1],
			 g3[j*size+i],
			 g3[j*size+i+1],
			 g3[(j-1)*size+i],
			 g3[(j+1)*size+i]);

	      *(field4+k+gridn*l) = 
		-alpha/180*
		(2*v[0] - 27*v[1] + 
		 270*v[2] - 980*v[3] + 270*v[4] -27*v[5] + 
		 2*v[6] + 2*v[7] - 27*v[8] + 270*v[9] + 
		 270*v[10] - 27*v[11] + 2*v[12]);

	      /* Update the circulations */
	      *(circ+k+gridn*l) +=
		*(field1+k+gridn*l)/6 +
		*(field2+k+gridn*l)/3 +
		*(field3+k+gridn*l)/3 +
		*(field4+k+gridn*l)/6;
	    }
      }

  /* Don't forget to free up all these pointers. */

  for (i=0; i<SQR(size); ++i)
    {
      free(g1[i]);
      free(g2[i]);
      free(g3[i]);
      free(k1[i]);
      free(k2[i]);
      free(k3[i]);
    }

  return(h);
}

double RHE_RK4_grd_init(double X0, double X1, double Y0, double Y1,
			int gridn, double *circs)

{
  int    i,j,k,l;
  double h,s2;
  double *g1,*g2,*g3,*k1,*k2,*k3;
  double *field1,*field2,*field3,*field4;
  double v[13];
  double alpha;

  /* CAUTION: right now the code only works at the finest level. */
  /*  levels = mplevels;  
      size = (int) ldexp(1.0,levels); */

  h = (X1-X0)/gridn;

  if (fabs((h-(Y1-Y0)/gridn)/h) > 1.0e-6)
    {
      fprintf(diag_log,"The initialization grid must be square.\n");
      exit(-1);
    }

  s2 = h*h;

  alpha = s2/h/h;

  /* Initialize the circulations. */

  /* v_1-6   = v_i-3,j...v_i+3,j */
  /* v_7-9   = v_i,j-3...v_i,j-1 */
  /* v_10-12 = v_i,j+1...v_i,j+3 */

  /* Half step FCE */

  field1 = malloc(sizeof(h)*SQR(gridn));
  field2 = malloc(sizeof(h)*SQR(gridn));
  if (!field1)
    fprintf(diag_log,"field1 malloc failed.\n");
  if (!field2)
    fprintf(diag_log,"field2 malloc failed.\n");
  for (l=0; l<gridn; ++l)
    for (k=0; k<gridn; ++k)
      {
	fill_vec(v,k,l,gridn,NULL,circs,NULL,NULL,NULL);

	*(field1+k+gridn*l) = 
	  -alpha/180*
	  (2*v[0] - 27*v[1] + 
	   270*v[2] - 980*v[3] + 270*v[4] -27*v[5] + 
	   2*v[6] + 2*v[7] - 27*v[8] + 270*v[9] + 
	   270*v[10] - 27*v[11] + 2*v[12]);
	
	*(field2+k+gridn*l) = *(circs+k+gridn*l) + 
	  0.5*(*(field1+k+gridn*l));
      }
  
  k1 = field1;
  g1 = field2;

  /* Half step BCE */

  field1 = malloc(sizeof(h)*SQR(gridn));
  field2 = malloc(sizeof(h)*SQR(gridn));
  if (!field1)
    fprintf(diag_log,"field1 malloc failed.\n");
  if (!field2)
    fprintf(diag_log,"field2 malloc failed.\n");
  for (l=0; l<gridn; ++l)
    for (k=0; k<gridn; ++k)
      {
	fill_vec(v,k,l,gridn,NULL,circs,NULL,NULL,NULL);

	*(field1+k+gridn*l) = 
	  -alpha/180*
	  (2*v[0] - 27*v[1] + 
	   270*v[2] - 980*v[3] + 270*v[4] -27*v[5] + 
	   2*v[6] + 2*v[7] - 27*v[8] + 270*v[9] + 
	   270*v[10] - 27*v[11] + 2*v[12]);
	
	*(field2+k+gridn*l) = *(circs+k+gridn*l) + 
	  0.5*(*(field1+k+gridn*l));
	
      }
  k2 = field1;
  g2 = field2;

  /* Midpoint */

  field1 = malloc(sizeof(h)*SQR(gridn));
  field2 = malloc(sizeof(h)*SQR(gridn));
  if (!field1)
    fprintf(diag_log,"field1 malloc failed.\n");
  if (!field2)
    fprintf(diag_log,"field2 malloc failed.\n");
  for (l=0; l<gridn; ++l)
    for (k=0; k<gridn; ++k)
      {
	fill_vec(v,k,l,gridn,NULL,circs,NULL,NULL,NULL);
	
	*(field1+k+gridn*l) = 
	  -alpha/180*
	  (2*v[0] - 27*v[1] + 
	   270*v[2] - 980*v[3] + 270*v[4] -27*v[5] + 
	   2*v[6] + 2*v[7] - 27*v[8] + 270*v[9] + 
	   270*v[10] - 27*v[11] + 2*v[12]);
	
	*(field2+k+gridn*l) = *(circs+k+gridn*l) + 
	  (*(field1+k+gridn*l));
      }
  k3 = field1;
  g3 = field2;

  /* Simpsons corrector */

  field4 = malloc(sizeof(h)*SQR(gridn));
  if (!field4)
    fprintf(diag_log,"field4 malloc failed.\n");
  field1 = k1;
  field2 = k2;
  field3 = k3;

  for (l=0; l<gridn; ++l)
    for (k=0; k<gridn; ++k)
      {
	fill_vec(v,k,l,gridn,NULL,circs,NULL,NULL,NULL);

	*(field4+k+gridn*l) = 
	  -alpha/180*
	  (2*v[0] - 27*v[1] + 
	   270*v[2] - 980*v[3] + 270*v[4] -27*v[5] + 
	   2*v[6] + 2*v[7] - 27*v[8] + 270*v[9] + 
	   270*v[10] - 27*v[11] + 2*v[12]);
	
	/* Update the circulations */
	*(circs+k+gridn*l) +=
	  *(field1+k+gridn*l)/6 +
	  *(field2+k+gridn*l)/3 +
	  *(field3+k+gridn*l)/3 +
	  *(field4+k+gridn*l)/6;
      }

  /* Don't forget to free up all these pointers. */

  free(g1);
  free(g2);
  free(g3);
  free(k1);
  free(k2);
  free(k3);
  free(field4);

  return(h);
}

void RHE_interp(double s2, double pop_control)
{
  double **circs,*circ;
  double x,y,h,gridsize;
  int size,levels,gridn;
  int i,j,k,l,count;

  mplevels = Set_Level();
  partition(mplevels);

  /* CAUTION: right now the code only works at the finest level. */
  levels = mplevels;  

  size = (int) ldexp(1.0,levels);

  circs = malloc(sizeof(circ)*SQR(size));
  if (!circs)
    fprintf(diag_log,"circ malloc failed.\n");

  gridsize = (maxX-minX)/size;

  gridn = (int) ((gridsize/sqrt(s2))+0.5);

  h = RHE_RK4(circs,gridn,s2);

  count = 0;

  for (i=0; i<size; ++i)
    for (j=0; j<size; ++j)
      {
	circ = circs[j*size+i];

	for (k=0; k<gridn; ++k)
	  for (l=0; l<gridn; ++l)
	    {
	      if (fabs(circ[l*gridn+k])/h/h>pop_control)
		{
		  x = minX+i*gridsize+(k+0.5)*h;
		  y = minY+j*gridsize+(l+0.5)*h;
		  
		  mblob[count].blob0.x = x;
		  mblob[count].blob0.y = y;
		  mblob[count].blob0.strength = circ[l*gridn+k]/2/M_PI;
		  blobguts[count].s2 = s2;
		  blobguts[count].a2 = 1.0;
		  blobguts[count].th = 0.0;
		  
		  set_blob(&(blobguts[count]),&(tmpparms[count]));

		  ++count;

		  if (count==NMAX)
		    {
		      fprintf(diag_log,"Out of memory in RHE remesh.\n");
		      fprintf(diag_log,"Time to sleep.\n");
		    }
		}
	    }
      }

  for (i=0; i<SQR(size); ++i)
    if (circs[i])
      free(circs[i]);

  Release_Links(mplevels);

  N = count;
}

