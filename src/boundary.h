/* BOUNDARY.H */
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
 *                                                                      *
 * Louis Rossi                                                          *
 * Department of Mathematical Sciences                                  *
 * University of Delaware                                               *
 * Newark, DE 19716                                                     */

/* If you are not able to install the BLAS/LAPACK libraries, you can
   set this flag which will remove all references to the library calls. */

#undef NOBOUNDARY

#define BMAX 250

typedef struct
{
  double m,x,y,nx,ny,l;
}
Panel;

/* Boundaries variables */
extern int    B,Bpiv[BMAX];
extern Panel  walls[BMAX];
extern double BdyMat[BMAX][BMAX];

extern double BoundaryStep;
extern int    BoundaryFrame;

/* Boundary headers */

extern void BoundaryConstrain();

extern void solve_bdy_matrix(Panel[],int[],double[][BMAX]);
extern void factor_bdy_matrix(Panel[],int[],double[][BMAX]);
extern void bdy_vel(double, double, double*, double*);
extern void clip(double);

/* Lapack headers */
#ifndef NOBOUNDARY
extern int dgesv_(int*,int*,double[][BMAX],int*,int[],double[],int*,int*);
extern int dgetrf_( int*, int*, double[][BMAX], int*, int[], int*);
void dgetrs_( char*, int*, int*, double[][BMAX], int*, int[], double[], int*, int*);
#endif
