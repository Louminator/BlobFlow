/* MULTIPROC.H */ 
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

#ifndef MPI_HEADER
#define MPI_HEADER
#include <mpi.h>

#define DONE 0
#define WORK 1
#define MORE 2
#define RESERVED_TAGS 5
#define PARTICLE_DATA_PACKET_SIZE 5

int total_processes, rank;
MPI_Status mpistatus;
MPI_Request mpireq;

/* CAUTION: Maximum number of processes is hardwired at 30. */

void master ( void );
void master2 ( int levels );
void slave   ( void );
void peer    ( void );
void finish  ( void );
void oldfinish  ( void );
void Setup_mpi   ( int, char ** );
void Cleanup_mpi ( void );
#endif
















