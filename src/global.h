/* GLOBAL.H */
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
 * or after 15 January 2001                                             *
 * Louis Rossi                                                          *
 * Department of Mathematical Sciences                                  *
 * University of Delaware                                               *
 * Newark, DE 19716                                                     */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <ctype.h>

#define SQR(X) ((X)*(X))
#define a2diff(X) (1.0/(X).a2-(X).a2)
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define SGN(X) ((X) < (0.0) ? (-1) : (1))

#define CORRECTVEL4
#undef  NOFASTMP
#define NOBOUNDARY

/* Set LINEAR to impose your own velocity field. */
#undef  LINEAR
#undef  cashkarp

#define NMAX 65000       /* Maximum number of computational elements. */
#define BMAX 250
#define FILENAME_LEN 500 /* Maximum file name length */

typedef struct
{
    double x,y;
}
Vector;

typedef struct
{
   double du11,du12,du21;  
}
Tensor;

typedef struct
{
  double m,x,y,nx,ny,l;
}
Panel;

typedef struct
{
   double re,im;
}
Complex;

extern FILE          *comp_log,*diag_log,*mpi_log,*cpu_log;

extern int           xantisymm;

/* Boundaries */
extern int    B,Bpiv[BMAX];
extern Panel  walls[BMAX];
extern double BdyMat[BMAX][BMAX];

/* Time integration parameters */
extern double TimeStep,PrefStep,FrameStep,EndTime,SimTime;
extern int    Frame,MaxOrder;

/*Fusion parameters */
extern double merge_budget,merge_p,merge_a2tol,merge_thtol;
extern double merge_mom3wt,merge_mom4wt,clusterR,aM2;
extern double MergeStep,merge_c,merge_growth_rate;
extern int    merge_estimator,nsplit,nmerge,totsplit,totmerge,MergeFrame;

extern double BoundaryStep;
extern int    BoundaryFrame;

/*Splitting parameters */
extern int    split_method;
extern double *split_parm_ptr;

/*Fast multipole variables*/

/*Data export variables */
extern char filename[];

extern clock_t tot_cputime_ref,tot_cputime,
  vel_cputime_ref,vel_cputime,velsum_cputime_ref,velsum_cputime,
  veldirect_cputime_ref,veldirect_cputime,
  mp_cputime,mp_cputime_ref;

/* Function headers */

extern void init(int, char**);
extern int  main(int, char**);
extern void rectify();

extern void flip_blob();
extern void set_blob();

extern Vector induced_vel();
extern Tensor induced_veldev();

extern void rk4(double);

/* Streamfunction approximation routines. */

extern void resort();
extern void write_vorts(int);
extern void write_pmvorts(int);
extern void write_partition(int);
extern void potential_flow(int);
extern Vector potential_vel(double,double);
extern void reflect_X();
extern void cache_resort();

extern void RHE_interp(double,double);

extern void BoundaryConstrain();

extern void stop(int);

extern void chk_nan(int);

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
