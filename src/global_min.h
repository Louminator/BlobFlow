/* GLOBAL_MIN.H */
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
 * or                                                                   *
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

typedef struct
{
   double re,im;
}
Complex;

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

/* Parameters */
extern int           xantisymm;

/* Time integration parameters */
extern double TimeStep,PrefStep,FrameStep,EndTime,SimTime;
extern int    Frame,MaxOrder;

/*Data export variables */

#define FILENAME_LEN 500 /* Maximum file name length */

extern FILE          *comp_log,*diag_log,*mpi_log,*cpu_log;
extern char filename[];


/* Splitting parameters DEPRECATED */
extern int    split_method;
extern double *split_parm_ptr;


/* Fusion parameters DEPRECATED */
extern double merge_budget,merge_p,merge_a2tol,merge_thtol;
extern double merge_mom3wt,merge_mom4wt,clusterR,aM2;
extern double MergeStep,merge_c,merge_growth_rate;
extern int    merge_estimator,nsplit,nmerge,totsplit,totmerge,MergeFrame;

