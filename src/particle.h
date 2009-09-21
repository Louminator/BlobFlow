/* PARTICLE.H */
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
 * or                                                                   *
 * Louis Rossi                                                          *
 * Department of Mathematical Sciences                                  *
 * University of Delaware                                               *
 * Newark, DE 19716                                                     */

/* This header describes particles. */

#define CORRECTVEL4

#define NMAX 50000       /* Maximum number of computational elements. */

typedef struct
{
    double x,y;
    double strength;
    double dx,dy;
}
Blob_external;

typedef struct
{
    double a2,s2,th;
}
Blob_internal;

typedef struct
{
  double sin2,cos2,sincos,costh,sinth,du11,du12,du21;
  double u_xx,u_xy,u_yy,v_xx;
}
Blob_parms;

typedef struct
{
    Blob_external blob0,blob1;
}
Metablob;

extern int           N;
extern Metablob      mblob[NMAX];
extern Blob_internal blobguts[NMAX];
extern Blob_parms    tmpparms[NMAX];
extern double        visc;                       /* Physical constants */
extern double        alpha,l2tol,dtth_delta;     /* Numerical parameters */

extern void vel_field();

extern double dta2(Blob_internal*,Blob_parms*);
extern double dts2(Blob_internal*);
extern double dtth(Blob_internal*,Blob_parms*,double,double);
