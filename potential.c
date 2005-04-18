/* POTENTIAL.C */
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
 * Louis Rossi                                                          *
 * Department of Mathematical Sciences                                  *
 * University of Delaware                                               *
 * Newark, DE 19716                                                     */

#undef DEBUG_MP


#include "global.h"

/* Enter potential flow field information here. */
vector potential_vel(x,y)
     double x,y;
{
  vector v;

  v.x = 1.0;
  v.y = 0.0;

  return(v);
}

void potential_flow(vort)
     int vort;

{
  vector v;

  v = potential_vel(mblob[vort].blob0.x,mblob[vort].blob0.y);

  mblob[vort].blob0.dx += v.x;
  mblob[vort].blob0.dy += v.y;

  tmpparms[vort].du11 += 0.0;
  tmpparms[vort].du12 += 0.0;
  tmpparms[vort].du21 += 0.0;

  tmpparms[vort].u_xx += 0.0;
  tmpparms[vort].u_xy += 0.0;
  tmpparms[vort].u_yy += 0.0;
  tmpparms[vort].v_xx += 0.0;
}

