/* EVOLVE-MS.C */
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
 * Louis Rossi                                                          *
 * Department of Mathematical Sciences                                  *
 * University of Delaware                                               *
 * Newark, DE 19715-2553                                                */

#include "global_min.h"
#include "particle.h"
#include "biot-savart.h"
#include "boundary.h"
#include "field_interp.h"

#ifdef MULTIPROC
#include "multiproc.h"
#endif

/* Global variables */

FILE *comp_log,*diag_log,*mpi_log,*cpu_log;

int       N;
Metablob  mblob[NMAX];
Blob_internal blobguts[NMAX];
double    visc;
int       xantisymm;
double    alpha,l2tol,dtth_delta;
double    TimeStep,PrefStep,FrameStep,EndTime,SimTime;
int       Frame,MaxOrder;
double    merge_budget,merge_p,merge_a2tol,merge_thtol,MergeStep;
double    merge_mom3wt,merge_mom4wt,clusterR,aM2;
double    merge_c,merge_growth_rate;

double    InterpStep,InterpPopulationControl,InterpVar;
int       Interps;

int       split_method;
double    *split_parm_ptr;

int       merge_estimator,nsplit,nmerge,totmerge,totsplit,MergeFrame;

double    BoundaryStep;
int       BoundaryFrame;

char      filename[FILENAME_LEN];
int       write_vtx,write_partitions,write_vel;

Blob_parms tmpparms[NMAX];

/*Fast multipole variables*/

double minX,maxX,minY,maxY,distX,distY;

Complex *Level_Ptr[LMAX];

int gridx[LMAX][NMAX], gridy[LMAX][NMAX], mplevels;

int numk2;

double aM2;

FineGridLink **FineGridLinks,*trace;

/* Boundaries */

int    B,Bpiv[BMAX];
Panel  walls[BMAX];
double BdyMat[BMAX][BMAX];

double phase,mxx,mxy,myy,old_mxx,old_mxy,old_myy;

clock_t tot_cputime_ref,tot_cputime,
  vel_cputime_ref,vel_cputime,velsum_cputime_ref,velsum_cputime,
  veldirect_cputime_ref,veldirect_cputime,
  mp_cputime,mp_cputime_ref,
  mp_Init_Fine_Grid,mp_Init_Fine_Grid_ref,
  mp_Advance_Coeffs,mp_Advance_Coeffs_ref;


#ifdef MULTIPROC
int total_processes, rank;
MPI_Status mpistatus;
MPI_Request mpireq;

double wicked_big_vectah[NMAX * PARTICLE_DATA_PACKET_SIZE];
#endif


/* Subroutine to stop everything if there is a problem. */
void stop(int retval)
{
  fflush(comp_log);
  fflush(diag_log);
  fflush(mpi_log);
   
#ifdef MULTPROC
  Cleanup_mpi();
#endif
   
  exit(retval);
}

void run()
{
  int j,k;
  double Rblob,temp,lambdaM;
 
  fprintf(comp_log,"Time     N     Split Merge MPlev lambdaM\n");
   
  while (SimTime < EndTime)
    {

      tot_cputime_ref = clock();

      numk2 = 0;

      rk4(TimeStep);

      SimTime += TimeStep;

      for (j=0; j<N; ++j)
	if (blobguts[j].a2 < 0.0)
	  {
	    fprintf(diag_log,
		    "Time: %12.4e. WARNING: Negative a2 in element %d\n",
		    SimTime,j);
	    for (k=0; k<N; ++k)
	      fprintf(diag_log,
		      "str=%10.3e x=%10.3e y=%10.3e s2=%10.3e a2=%10.3e\n",
		      mblob[k].blob0.strength,
		      mblob[k].blob0.x,mblob[k].blob0.y,
		      blobguts[k].s2,blobguts[k].a2);
	    exit(0); 
	  }

      /* Constrain the bdy conditions */
      /*
      if ( (BoundaryStep > 0.0) &&
	   ((SimTime+0.499*TimeStep) >= BoundaryStep*BoundaryFrame) )
	{
	  ++BoundaryFrame;
	  BoundaryConstrain();
	}
      */
      if ( (InterpStep > 0.0) &&
	   ((SimTime+0.499*TimeStep) >= InterpStep*Interps) )
	{
	  RHE_interp2(InterpVar,InterpPopulationControl);
	  ++Interps;
	}

      if ((SimTime+0.499*TimeStep) >= FrameStep*Frame)
	{
	  if (write_vel) write_vels(Frame);
#ifdef MULTIPROC
	  /* Have only the 'root' node do the writes */
	  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	  if (rank == 0) 
	    {
	      if (write_vtx) write_vorts(Frame);
	      if (write_partitions) write_partition(Frame);
	    }
#else 
	  /* single processor */
	  if (write_vel) write_vels(Frame);
	  if (write_partitions) write_partition(Frame);	    
#endif

	  ++Frame; 

	  Rblob = 0.0;

	  lambdaM = 0.0;

	  for (j=0; j<N; ++j)
	    {
	      temp = 2.0*
		(tmpparms[j].du11*(tmpparms[j].cos2-tmpparms[j].sin2)+
		 (tmpparms[j].du12+tmpparms[j].du21)*tmpparms[j].sincos);

	      if (lambdaM < temp) 
		lambdaM = temp;

	      temp = 2.0*
		(tmpparms[j].du11*(tmpparms[j].cos2-tmpparms[j].sin2)-
		 (tmpparms[j].du12+tmpparms[j].du21)*tmpparms[j].sincos);

	      if (lambdaM < temp) 
		lambdaM = temp;
	    }

	  fprintf(comp_log,"%8.2e %-5d %-5d %-5d %-5d %8.2e\n",
		  SimTime,N,totsplit+nsplit,totmerge+nmerge,
		  mplevels,lambdaM);
	  fflush(comp_log);
	  fflush(diag_log);
	  fflush(mpi_log);
	}

    }  /*   end while loop  */

  fprintf(comp_log,"Simulation complete.\n");
  fprintf(comp_log,"Total splitting events: %d\n",totsplit);
  fprintf(comp_log,"Total merging events: %d\n",totmerge);
  fclose(cpu_log);
}



#ifdef MULTIPROC
int main(int argc, char* argv[])
{
  double starttime, endtime;
   
  Setup_mpi(argc, argv);
  init(argc, argv);
   
  starttime = MPI_Wtime();   
   
  fprintf(diag_log,"Running.\n");
  run();
   
  fclose(comp_log);
  fclose(diag_log);
  endtime = MPI_Wtime();
  fclose(mpi_log);
  Cleanup_mpi();
   
  return(1); 
}

#else

int main(int argc, char* argv[])
{
   
  init(argc, argv);
  fprintf(diag_log,"Running.\n");
  run();
   
  fclose(comp_log);
  fclose(diag_log);

  return(1); 
}
#endif
