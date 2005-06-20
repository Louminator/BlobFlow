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

/* A little routine for a 'cylinder' */
void BoundaryConstrain()
{
  vector v;
  double normu,theta,circ;

  double beta;

  beta = 1.0;

  v = vel_gen(0.0,0.0);
  normu = sqrt(SQR(v.x)+SQR(v.y));
  theta = atan2(v.y,v.x);
  circ = normu*beta*alpha*l2tol/(1.0-exp(-SQR(beta)));

  mblob[N].blob0.strength = -circ;
  mblob[N].blob0.x = -2.0*beta*alpha*l2tol*sin(theta);
  mblob[N].blob0.y =  2.0*beta*alpha*l2tol*cos(theta);
  blobguts[N].s2 = alpha*l2tol;
  blobguts[N].a2 = 1.0;
  blobguts[N].th = 0.0;

  mblob[N].blob0.order = 1;
  tmpparms[N].refinecnt = -1;
  tmpparms[N].nint = 1;
  set_blob(&(blobguts[N]),&(tmpparms[N]));

  mblob[N].blob1 = mblob[N].blob0;
  mblob[N].blob2 = mblob[N].blob0;
  mblob[N].blob3 = mblob[N].blob0;
  mblob[N].blob4 = mblob[N].blob0;

  ++N;

  mblob[N].blob0.strength = circ;
  mblob[N].blob0.x =  2.0*beta*alpha*l2tol*sin(theta);
  mblob[N].blob0.y = -2.0*beta*alpha*l2tol*cos(theta);
  blobguts[N].s2 = alpha*l2tol;
  blobguts[N].a2 = 1.0;
  blobguts[N].th = 0.0;
  mblob[N].blob0.order = 1;
  tmpparms[N].refinecnt = -1;
  tmpparms[N].nint = 1;
  set_blob(&(blobguts[N]),&(tmpparms[N]));

  mblob[N].blob1 = mblob[N].blob0;
  mblob[N].blob2 = mblob[N].blob0;
  mblob[N].blob3 = mblob[N].blob0;
  mblob[N].blob4 = mblob[N].blob0;

  ++N;
}
