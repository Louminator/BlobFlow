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

#include "global.h"

#ifdef MULTIPROC
#include "multiproc.h"
#endif

#undef MERGEDIAG

/* Global variables */

FILE *comp_log,*diag_log,*mpi_log,*cpu_log;

int N,oldN;
metablob mblob[NMax];
blob_internal blobguts[NMax];
double    visc;
double    alpha,l2tol,dtth_delta;
double    TimeStep,PrefStep,FrameStep,EndTime,SimTime;
int       Frame,MaxOrder;
double    merge_budget,merge_p,merge_a2tol,merge_thtol,MergeStep;
double    merge_mom3wt,merge_mom4wt,clusterR,aM2;
double    merge_c,merge_growth_rate;
int       split_method;
double    *split_parm_ptr;

int       merge_estimator,nsplit,nmerge,totmerge,totsplit,MergeFrame;

double    BoundaryStep;
int       BoundaryFrame;

char      filename[Title];

blobparms tmpparms[NMax];

/* Dynamic memory allocation might be better here. Nah! */
vector       refinevels[NMax][3][MaxSplitConf];
int          refineindex[NMax],refinestack;

/*Fast multipole variables*/

double minX,maxX,minY,maxY,distX,distY;

complex *Level_Ptr[LMax];

int gridx[LMax][NMax], gridy[LMax][NMax], mplevels;

int numk2;

double aM2;

FineGridLink **FineGridLinks,*trace;

/* Boundaries */

int    B,Bpiv[BMax];
panel  walls[BMax];
double BdyMat[BMax][BMax];

double phase,mxx,mxy,myy,old_mxx,old_mxy,old_myy;

clock_t tot_cputime_ref,tot_cputime,
  vel_cputime_ref,vel_cputime,velsum_cputime_ref,velsum_cputime,
  veldirect_cputime_ref,veldirect_cputime,
  mp_cputime,mp_cputime_ref;

#ifdef MULTIPROC
int total_processes, rank;
MPI_Status mpistatus;
MPI_Request mpireq;

double wicked_big_vectah[NMax * PARTICLE_DATA_PACKET_SIZE];
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

      /* clip(2.5); */
      tot_cputime_ref = clock();
      vel_cputime = 0;
      velsum_cputime = 0;
      veldirect_cputime = 0;
      mp_cputime = 0;

      numk2 = 0;
	
      nsplit = 0;
      nmerge = 0;
	
      /* Track second moments for Travis' forcing */
      /*
      old_mxx = mxx;
      old_mxy = mxy;
      old_myy = myy;
      calc_moments(&mxx,&myy,&mxy);
      */
      /* */

      oldN = N;

      /* All computational elements march forward with split step. */

      /* Take a half step externally. */
      for (j=0; j<N; ++j)
	{
	  push(&mblob[j]);
	     
	  if (mblob[j].blob0.order == 1)
	    {
	      push(&mblob[j]);
	      ab2half(&mblob[j],&tmpparms[j],TimeStep);
	    }
	     
	  if (mblob[j].blob0.order == 2)
	    ab2half(&mblob[j],&tmpparms[j],TimeStep);
	     
	  if (mblob[j].blob0.order == 3)
	    ab3half(&mblob[j],&tmpparms[j],TimeStep);
	     
	  if (mblob[j].blob0.order == 4)
	    ab4half(&mblob[j],&tmpparms[j],TimeStep);
	     
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
	}

      SimTime += TimeStep/2;

      vel_cputime_ref = clock();
      vel_field();
      vel_cputime += clock()-vel_cputime_ref;

      /* Take a full step internally */
	
      for (j=0; j<N; ++j)
	{
#ifdef cashkarp
	  rkckmarch(&(blobguts[j]),&(tmpparms[j]),TimeStep,1.0e-5);
#else
	  internal_march(&(blobguts[j]),&(tmpparms[j]),TimeStep);
#endif
	}

      /* Take a full step externally. */
      for (j=0; j<N; ++j)
	{
	  if (mblob[j].blob0.order == 1)
	    {
	      ab2(&mblob[j],&tmpparms[j],TimeStep);
	    }
	     
	  if (mblob[j].blob0.order == 2)
	    ab2(&mblob[j],&tmpparms[j],TimeStep);
	     
	  if (mblob[j].blob0.order == 3)
	    ab3(&mblob[j],&tmpparms[j],TimeStep);
	     
	  if (mblob[j].blob0.order == 4)
	    ab4(&mblob[j],&tmpparms[j],TimeStep);
	     
	  if (mblob[j].blob0.order < MaxOrder)
	    ++mblob[j].blob0.order;
	     
	  if (blobguts[j].a2 < 0.0)
	    {
	      fprintf(diag_log,
		      "WARNING: Negative a2 in element %d\n",j);
	      for (k=0; k<N; ++k)
		fprintf(diag_log,
			"str=%10.3e x=%10.3e y=%10.3e s2=%10.3e a2=%10.3e\n",
			mblob[k].blob0.strength,
			mblob[k].blob0.x,mblob[k].blob0.y,
			blobguts[k].s2,blobguts[k].a2);
	      exit(0); 
	    }	     
	}

      /* Set commonly used parameters here. Beyond this point, if any  *
       * subroutine messes with computational elements, it should also *
       * reset the parameters.                                         */
      for (j=0; j<N; ++j)
	set_blob(&(blobguts[j]),&(tmpparms[j]));
	
      /*SimTime += TimeStep;*/
      SimTime += TimeStep/2;
	
      /* Split horizontally challenged elements. */
      chksplit();

      /* Constrain the bdy conditions */
      /*
      if ( (BoundaryStep > 0.0) &&
	   ((SimTime+0.499*TimeStep) >= BoundaryStep*BoundaryFrame) )
	{
	  ++BoundaryFrame;
	  BoundaryConstrain();
	}
      */

      if ((SimTime+0.499*TimeStep) >= MergeStep*MergeFrame)
	{
	  mplevels = Set_Level();
	     
	  partition(mplevels);
	     
#ifdef MERGEDIAG
#ifdef MULTIPROC
	  /* Have only the 'root' node do the writes */
	  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	  if (rank == 0) 
	    write_pmvorts(Frame);
#else
	  write_pmvorts(Frame);
#endif
#endif
	  merge();
	  resort();
	     
	  Release_Links(mplevels);
	     
	  ++MergeFrame;
	}
	
      if ((SimTime+0.499*TimeStep) >= FrameStep*Frame)
	{
#ifdef MULTIPROC
	  /* Have only the 'root' node do the writes */
	  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	  if (rank == 0) {
	    write_vorts(Frame);
	    write_partition(Frame);
	  }
#else 
	  /* single processor */
	  write_vorts(Frame);
	  write_partition(Frame);	    
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

      vel_cputime_ref = clock();
      vel_field();
      vel_cputime += clock()-vel_cputime_ref;

      tot_cputime = clock()-tot_cputime_ref;

      fprintf(cpu_log,"%07d %12.4e %12.4e %12.4e %12.4e %12.4e\n",
	      N,
	      ((double)(tot_cputime))/((double)CLOCKS_PER_SEC),
	      ((double)(vel_cputime))/((double)CLOCKS_PER_SEC),
	      ((double)(velsum_cputime))/((double)CLOCKS_PER_SEC),
	      ((double)(veldirect_cputime))/((double)CLOCKS_PER_SEC),
	      ((double)(mp_cputime))/((double)CLOCKS_PER_SEC)
	      );
      fflush(cpu_log);

      totsplit += nsplit;
      totmerge += nmerge;
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
