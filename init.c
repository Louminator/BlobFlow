/* INIT-MS.C */
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

/* Initialization subroutines for eflow. */

#include "global.h"

#ifdef MULTIPROC
#include "multiproc.h"
#endif

#ifdef XANTISYMM
#define CARDINALITY  N/2
#else
#define CARDINALITY  N
#endif

#define MERGEDIAG

/* Script constants. */
static char *inputs[] =
{
  "FrameStep:",
  "EndTime:",
  "Viscosity:",
  "VtxInit:",
  "TimeStep:",

  "MaxOrder:",
  "DtThDelta",
  "l2Tol:",
  "alpha:",
  "SplitMethod:",

  "MergeErrorEstimator:",
  "MergeBudget:",
  "MergeStep:",
  "ClusterRadius:",
  "Mergea2Tol:",

  "MergeThTol:",
  "MergeMom3Wt:",
  "MergeMom4Wt:",
  "MergeC:",
  "MergeGrowthRate:"
};

#include "split15asym.h"

/* Use a lookup table for 15 asymmetric splitting parameters.*/

enum ScriptItems 
{FRAME_STEP,END_TIME,VISCOSITY,VTX_INIT,TIME_STEP,
 MAX_ORDER,DTTH_DELTA,L2_TOL,ALPHA,SPLIT_METHOD,MERGE_ERROR_ESTIMATOR,
 MERGE_BUDGET,MERGE_STEP,CLUSTER_RADIUS,MERGE_A2_TOL,MERGE_TH_TOL,
 MERGE_MOM3_WT,MERGE_MOM4_WT,MERGE_C,MERGE_GROWTH_RATE
};

void check_sim(char *vtxfilename)
{
  if ( (FrameStep <= 0.0) ||
       (EndTime <= 0.0) ||
       (visc < 0.0) ||
       (strcmp(vtxfilename,"") == 0) )
    {
      printf("Simulation file error.\n");

      if (FrameStep <= 0.0)
	printf("FrameStep not set properly.\n");

      if (EndTime <= 0.0)
	printf("EndTime not set properly.\n");

      if (visc < 0.0)
	printf("Viscosity not set properly.\n");
      
      if (strcmp(vtxfilename,"") == 0)
	printf("VtxInit not set properly.\n");

      exit(-1);
    }
}

void read_sim()
{
  char      sim_name[Title],vtxfilename[Title],temp[Title];
  FILE      *sim_file;
  char *word,*p1;
  int  i,inputs_size;
   
  sprintf(sim_name,"%s%s",filename,".sim");
  fprintf(diag_log,"Reading simfile...\n");
  sim_file = fopen(sim_name,"r");
   
  /* ASSUME a 32-bit word. */
  inputs_size = sizeof(inputs)/4;
  word = malloc(100*sizeof(char));
   
  fscanf(sim_file,"%s",word);

  while (!feof(sim_file))
    {
	
      /*Keep this just in case.  On gcc, the memcmp automatically skips
       * over whitespace.  I am not sure if everyone does this.*/
      while ((isspace(*word)) && (word != '\0'))
	++word;
	
      if (*word == '#')
	{
	  while (*word != '\n')
	    *word = fgetc(sim_file);
	}
      else	
	if (*word != '\0')
	  {
	    for (i=0; i<inputs_size; ++i)
	      if (!memcmp(word,inputs[i],strlen(inputs[i])))
		{
		  switch (i)
		    {
		    case FRAME_STEP:
		      fscanf(sim_file,"%lf",&FrameStep);
		      i=inputs_size+1;
		      break;
		    case END_TIME:  
		      fscanf(sim_file,"%lf",&EndTime);
		      i=inputs_size+1;
		      break;
		    case VISCOSITY:
		      fscanf(sim_file,"%lf",&visc);
		      i=inputs_size+1;
		      break;
		    case VTX_INIT:
		      fscanf(sim_file,"%s",vtxfilename);
		      i=inputs_size+1;
		      break;
		    default:
		      printf("Uh oh.  There is an error in the scripting code.\n");
		    }
		}
	    if (i==inputs_size)
	      printf("Unrecognized input variable: %s\n",word);
	  }
      fscanf(sim_file,"%s",word);
    }
  fclose(sim_file);

  check_sim(vtxfilename);
   
  p1 = getenv("ECCSVM_HOME");
   
  if (p1 == NULL) 
    p1 = getenv("PWD");

  sprintf(temp,"%s/%s",p1,vtxfilename);
   
  sim_file = fopen(temp,"r");

  fscanf(sim_file,"%lf%lf%lf%lf%lf%lf",
	 &mblob[0].blob0.x,
	 &mblob[0].blob0.y,
	 &mblob[0].blob0.strength,&blobguts[0].s2,
	 &blobguts[0].a2,&blobguts[0].th);
   
  mblob[0].blob0.order = 1;
  tmpparms[0].refinecnt = -1;
  tmpparms[0].nint = 1;
  set_blob(&(blobguts[0]),&(tmpparms[0]));

  N = 1;
   
  while (!feof(sim_file))
    {
      fscanf(sim_file,"%lf%lf%lf%lf%lf%lf",
	     &mblob[N].blob0.x,
	     &mblob[N].blob0.y,
	     &mblob[N].blob0.strength,&blobguts[N].s2,
	     &blobguts[N].a2,&blobguts[N].th);
	
      mblob[N].blob0.order = 1;
      tmpparms[N].refinecnt = -1;
      tmpparms[N].nint = 1;
      set_blob(&(blobguts[N]),&(tmpparms[N]));
      ++N;
    }
  --N;
   
  fclose(sim_file);
   
  fprintf(comp_log,"Simulation parameters.\n");
  fprintf(comp_log,"%s %12.4e\n",inputs[FRAME_STEP],FrameStep);
  fprintf(comp_log,"%s %12.4e\n",inputs[END_TIME],EndTime);
  fprintf(comp_log,"%s %12.4e\n",inputs[VISCOSITY],visc);
  fprintf(comp_log,"%s %s\n",inputs[VTX_INIT],vtxfilename);
  fprintf(comp_log,"%d vortices read from %s.\n",N,temp);
  fflush(comp_log);
}

void check_ctl()
{
  int i;

  if ( (PrefStep <= 0.0) ||
       (MergeStep<FrameStep) ||
       (MaxOrder <= 0) ||
       (dtth_delta <= 0.0) ||
       (l2tol <= 0.0) ||
       (alpha <= 0.0) || (alpha >= 1.0) ||
       (split_method == 0) ||
       (merge_estimator == 0) ||
       (merge_budget < 0.0) ||
       (clusterR < 0.0) ||
       (MergeStep <= 0.0) ||
       ( (merge_estimator == 1) &&
	 ( (merge_mom3wt < 0.0) || (merge_mom4wt < 0.0) ) ) ||
       ( (merge_estimator == 2) &&
	 ( (merge_a2tol < 0.0) || (merge_thtol < 0.0) ) ) )
    {
      printf("Control file error.\n");

      if (PrefStep <= 0.0)
	printf("PrefStep not set properly.\n");
      if (MergeStep < FrameStep)
	printf("Warning!  Mergestep should not be less than FrameStep!\n");
      if (MaxOrder <= 0)
	printf("MaxOrder not set properly.\n");
      if (dtth_delta <= 0.0)
	printf("DtThDelta not set properly.\n");
      if ((alpha <= 0.0) || (alpha >= 1.0))
	printf("Alpha not set properly.\n");
      if (split_method == 0)
	printf("Split method not set properly.\n");
      if (merge_estimator == 0)
	printf("MergeErrorEstimator not properly set.\n");
      if ( (merge_estimator == 1) &&
	   ( (merge_mom3wt < 0.0) || (merge_mom4wt < 0.0) ) )
	printf("Parameters for MergeErrorEstimator 1 not properly set.\n");
      if ( (merge_estimator == 2) &&
	   ( (merge_a2tol < 0.0) || (merge_thtol < 0.0) ) )
	printf("Parameters for MergeErrorEstimator 2 not properly set.\n");
      if (clusterR < 0.0)
	printf("ClusterR not set properly.\n");
      if (MergeStep <= 0.0)
	printf(" not set properly.\n");

      exit(-1);
    }

  if (split_method == 3)
    {
      fprintf(diag_log,
	      "Using lookup table for asymmetric 1 -> 5 splitting.\n");
      i = (int) (100*(alpha - 0.6));
      alpha = 0.6 + i/100.0;
      fprintf(comp_log,
	      "Resetting alpha to the nearest hundreth: %3.1e\n",
	      alpha);
      split_parm_ptr = *split15asym_parms+i*6+1;
    }
}

void read_ctl()
{
  FILE *control_file;
  char *word,control_name[Title];
  int  i,inputs_size;

  sprintf(control_name,"%s%s",filename,".ctl");
   
  /* ASSUME a 32-bit word. */
  inputs_size = sizeof(inputs)/4;
   
  word = malloc(100*sizeof(char));

  control_file = fopen(control_name,"r");
   
  fscanf(control_file,"%s",word);

  while (!feof(control_file))
    {
	
      /*Keep this just in case.  On gcc, the memcmp automatically skips
       * over whitespace.  I am not sure if everyone does this.*/
      while ((isspace(*word)) && (word != '\0'))
	++word;
	
      if (*word == '#')
	{
	  while (*word != '\n')
	    *word = fgetc(control_file);
	}
      else	
	if (*word != '\0')
	  {
	    for (i=0; i<inputs_size; ++i)
	      if (!memcmp(word,inputs[i],strlen(inputs[i])))
		{
		  switch (i)
		    {
		    case TIME_STEP:
		      fscanf(control_file,"%lf",&PrefStep);
		      i=inputs_size+1;
		      break;
		    case MAX_ORDER:
		      fscanf(control_file,"%d",&MaxOrder);
		      i=inputs_size+1;
		      break;
		    case DTTH_DELTA:
		      fscanf(control_file,"%lf",&dtth_delta);
		      i=inputs_size+1;
		      break;
		    case L2_TOL:
		      fscanf(control_file,"%lf",&l2tol);
		      i=inputs_size+1;
		      break;
		    case ALPHA:
		      fscanf(control_file,"%lf",&alpha);
		      i=inputs_size+1;
		      break;
		    case SPLIT_METHOD:
		      fscanf(control_file,"%d",&split_method);
		      i=inputs_size+1;
		      break;
		    case MERGE_ERROR_ESTIMATOR:
		      fscanf(control_file,"%d",&merge_estimator);
		      i=inputs_size+1;
		      break;
		    case MERGE_BUDGET:
		      fscanf(control_file,"%lf",&merge_budget);
		      i=inputs_size+1;
		      break;
		    case CLUSTER_RADIUS:
		      fscanf(control_file,"%lf",&clusterR);
		      i=inputs_size+1;
		      break;
		    case MERGE_STEP:
		      fscanf(control_file,"%lf",&MergeStep);
		      i=inputs_size+1;
		      break;
		    case MERGE_A2_TOL:
		      fscanf(control_file,"%lf",&merge_a2tol);
		      i=inputs_size+1;
		      break;
		    case MERGE_TH_TOL:
		      fscanf(control_file,"%lf",&merge_thtol);
		      i=inputs_size+1;
		      break;
		    case MERGE_MOM3_WT:
		      fscanf(control_file,"%lf",&merge_mom3wt);
		      i=inputs_size+1;
		      break;
		    case MERGE_MOM4_WT:
		      fscanf(control_file,"%lf",&merge_mom4wt);
		      i=inputs_size+1;
		      break;
		    case MERGE_C:
		      fscanf(control_file,"%lf",&merge_c);
		      i=inputs_size+1;
		      break;
		    case MERGE_GROWTH_RATE:
		      fscanf(control_file,"%lf",&merge_growth_rate);
		      i=inputs_size+1;
		      break;
		    default:
		      printf("Uh oh.  There is an error in the scripting code.\n");
		    }
		}
	    if (i==inputs_size)
	      printf("Unrecognized input variable: %s\n",word);
	  }
      fscanf(control_file,"%s",word);
    }
  fclose(control_file);

  check_ctl();

  fprintf(comp_log,"Approximation control parameters.\n");
  fprintf(comp_log,"%s %12.4e\n",inputs[TIME_STEP],PrefStep);
  fprintf(comp_log,"%s %d\n",inputs[MAX_ORDER],MaxOrder);
  fprintf(comp_log,"%s %12.4e\n",inputs[MERGE_STEP],MergeStep);
  fprintf(comp_log,"%s %12.4e\n",inputs[L2_TOL],l2tol);
  fprintf(comp_log,"%s %12.4e\n",inputs[ALPHA],alpha);
  fprintf(comp_log,"%s %12.4e\n",inputs[MERGE_BUDGET],merge_budget);
  fprintf(comp_log,"%s %12.4e\n",inputs[CLUSTER_RADIUS],clusterR);
  /*  fprintf(comp_log,"%s %12.4e\n",inputs[MERGE_A2_TOL],merge_a2tol);
      fprintf(comp_log,"%s %12.4e\n",inputs[MERGE_TH_TOL],merge_thtol); */
  fflush(comp_log);
}

void startup_new()
{
  int i,j;
  blob_external startup_blobs[3][NMax];

  /* Storage for RK4 */
  /* blob1 - temporary use */
  /* blob2 - FCE */
  /* blob3 - BCE */
  /* blob4 - midpoint */

  for (j=0; j<N; ++j)
    startup_blobs[0][j] = mblob[j].blob0;

  for (i=1; i<3; ++i)
    {
      /* Take a half step of RK4 */

      /* Half step of FCE predictor. */
      for (j=0; j<N; ++j)
	mblob[j].blob1 = mblob[j].blob0;

      /* Note: blob1 stores the initial velocity field at t=0. */

      for (j=0; j<N; ++j)
	{
	  mblob[j].blob0.x = mblob[j].blob1.x + 
	    0.25*TimeStep*mblob[j].blob1.dx;
	  mblob[j].blob0.y = mblob[j].blob1.y + 
	    0.25*TimeStep*mblob[j].blob1.dy;
	}

      SimTime += TimeStep/4;

      vel_field();

      for (j=0; j<N; ++j)
	mblob[j].blob2 = mblob[j].blob0;

      /* Half step of BCE predictor. */
      for (j=0; j<N; ++j)
	{
	  mblob[j].blob0.x = mblob[j].blob1.x + 
	    0.25*TimeStep*mblob[j].blob2.dx;
	  mblob[j].blob0.y = mblob[j].blob1.y + 
	    0.25*TimeStep*mblob[j].blob2.dy;
	}

      vel_field();

      for (j=0; j<N; ++j)
	mblob[j].blob3 = mblob[j].blob0;

      /* Full step of midpoint rule predictor. */

      for (j=0; j<N; ++j)
	{
	  mblob[j].blob0.x = mblob[j].blob1.x + 
	    0.5*TimeStep*mblob[j].blob3.dx;
	  mblob[j].blob0.y = mblob[j].blob1.y + 
	    0.5*TimeStep*mblob[j].blob3.dy;
	}

      SimTime += TimeStep/4;

      vel_field();

      for (j=0; j<N; ++j)
	mblob[j].blob4 = mblob[j].blob0;

      /* Simpson's rule corrector. */

      for (j=0; j<N; ++j)
	{
	  mblob[j].blob0.x = mblob[j].blob1.x + 
	    0.5*TimeStep*(mblob[j].blob1.dx+2.0*mblob[j].blob2.dx+
			  2.0*mblob[j].blob3.dx+mblob[j].blob4.dx)/6.0;
	  mblob[j].blob0.y = mblob[j].blob1.y + 
	    0.5*TimeStep*(mblob[j].blob1.dy+2.0*mblob[j].blob2.dy+
			  2.0*mblob[j].blob3.dy+mblob[j].blob4.dy)/6.0;
	}

      vel_field();

      /* Take a full step internally */
      for (j=0; j<N; ++j)
	{
#ifdef cashkarp
	  rkckmarch(&(blobguts[j]),&(tmpparms[j]),TimeStep,1.0e-6);
#else
	  internal_march(&(blobguts[j]),&(tmpparms[j]),TimeStep);
#endif
	}

      /* Take a half step of RK4 */

      /* Half step of FCE predictor. */
      for (j=0; j<N; ++j)
	mblob[j].blob1 = mblob[j].blob0;

      /* Note: blob1 stores the initial velocity field at t=0. */

      for (j=0; j<N; ++j)
	{
	  mblob[j].blob0.x = mblob[j].blob1.x + 
	    0.25*TimeStep*mblob[j].blob1.dx;
	  mblob[j].blob0.y = mblob[j].blob1.y + 
	    0.25*TimeStep*mblob[j].blob1.dy;
	}

      SimTime += TimeStep/4;
      vel_field();

      for (j=0; j<N; ++j)
	mblob[j].blob2 = mblob[j].blob0;

      /* Half step of BCE predictor. */
      for (j=0; j<N; ++j)
	{
	  mblob[j].blob0.x = mblob[j].blob1.x + 
	    0.25*TimeStep*mblob[j].blob2.dx;
	  mblob[j].blob0.y = mblob[j].blob1.y + 
	    0.25*TimeStep*mblob[j].blob2.dy;
	}

      vel_field();

      for (j=0; j<N; ++j)
	mblob[j].blob3 = mblob[j].blob0;

      /* Full step of midpoint rule predictor. */

      for (j=0; j<N; ++j)
	{
	  mblob[j].blob0.x = mblob[j].blob1.x + 
	    0.5*TimeStep*mblob[j].blob3.dx;
	  mblob[j].blob0.y = mblob[j].blob1.y + 
	    0.5*TimeStep*mblob[j].blob3.dy;
	}

      SimTime += TimeStep/4;
      vel_field();

      for (j=0; j<N; ++j)
	mblob[j].blob4 = mblob[j].blob0;

      /* Simpson's rule corrector. */

      for (j=0; j<N; ++j)
	{
	  mblob[j].blob0.x = mblob[j].blob1.x + 
	    0.5*TimeStep*(mblob[j].blob1.dx+2.0*mblob[j].blob2.dx+
			  2.0*mblob[j].blob3.dx+mblob[j].blob4.dx)/6.0;
	  mblob[j].blob0.y = mblob[j].blob1.y + 
	    0.5*TimeStep*(mblob[j].blob1.dy+2.0*mblob[j].blob2.dy+
			  2.0*mblob[j].blob3.dy+mblob[j].blob4.dy)/6.0;
	}

      vel_field();

      for (j=0; j<N; ++j)
	startup_blobs[i][j] = mblob[j].blob0;

      /* Dump the startup positions. */

      if (FrameStep == TimeStep)
	{
	  write_vorts(i);
	  if (i == 3)
	    Frame = 4;
	}
      else
	write_vorts(9996+i);
    }

  /* Store startup elements back into the global variables. */

  for (j=0; j<N; ++j)
    {
      mblob[j].blob0 = startup_blobs[2][j];
      mblob[j].blob0.order = MaxOrder;
      mblob[j].blob1 = startup_blobs[2][j];
      mblob[j].blob1.order = MaxOrder;
    }

  for (j=0; j<N; ++j)
     mblob[j].blob2 = startup_blobs[1][j];

  for (j=0; j<N; ++j)
     mblob[j].blob3 = startup_blobs[0][j];

  /* Take a half step externally. */
  for (j=0; j<N; ++j)
    ab3half(&mblob[j],&tmpparms[j],TimeStep); 
  
  SimTime += TimeStep/2;
  vel_field();

  /* Take a full step internally */
  
  for (j=0; j<N; ++j)
    {
#ifdef cashkarp
	  rkckmarch(&(blobguts[j]),&(tmpparms[j]),TimeStep,1.0e-6);
#else
	  internal_march(&(blobguts[j]),&(tmpparms[j]),TimeStep);
#endif
    }
  
  /* Take a full step externally. */
  for (j=0; j<N; ++j)
    ab3(&mblob[j],&tmpparms[j],TimeStep);
  
  SimTime += TimeStep/2;
  vel_field();
  
  if (FrameStep == TimeStep)
    {
      write_vorts(3);
      Frame = 4;
    }
  else
    write_vorts(9999);

}

void startup()
{
  int i,j;
  double saveab4dx[NMax],saveab4dy[NMax];

  /* Take eight steps at first order. */
  for (i=0; i<8; ++i)
    {
      /* All computational elements march forward with split step. */
	
      /* Take a half step externally. */
      for (j=0; j<N; ++j)
	{
	  push(&mblob[j]);
	     
	  push(&mblob[j]);
	  ab2half(&mblob[j],&tmpparms[j],TimeStep/8.0);
	  /* ab2(&mblob[j],&tmpparms[j],TimeStep/8.0); */
	}
	
      vel_field();
	
      /* Take a full step internally */
	
      for (j=0; j<N; ++j)
	internal_march(&(blobguts[j]),&(tmpparms[j]),TimeStep/8.0);

      /* Take a full step externally. */
      for (j=0; j<N; ++j)
	ab2(&mblob[j],&tmpparms[j],TimeStep/8.0);

      SimTime += TimeStep/8.0;

      vel_field();
    }

  if (FrameStep == TimeStep)
    write_vorts(1);
  else
    write_vorts(9997);

  /* Take four steps at second order. */

  /* Advance the old velocity data because we have just doubled the
     step size. */
  for (j=0; j<N; ++j)
    {
      mblob[j].blob1.dx = mblob[j].blob2.dx;
      mblob[j].blob1.dy = mblob[j].blob2.dy;
    }

  for (i=0; i<4; ++i)
    {
      /* All computational elements march forward with split step. */
	
      /* Take a half step externally. */
      for (j=0; j<N; ++j)
	{
	  push(&mblob[j]);
	  ab2half(&mblob[j],&tmpparms[j],TimeStep/4.0);
	  /* ab2(&mblob[j],&tmpparms[j],TimeStep/4.0); */
	}
	
      vel_field();
	
      /* Take a full step internally */
	
      for (j=0; j<N; ++j)
	{
	  internal_march(&(blobguts[j]),&(tmpparms[j]),TimeStep/4.0);
	}

      /* Take a full step externally. */
      for (j=0; j<N; ++j)
	ab2(&mblob[j],&tmpparms[j],TimeStep/4.0);

      SimTime += TimeStep/4.0;

      vel_field();
    }

  if (FrameStep == TimeStep)
    write_vorts(2);
  else
    write_vorts(9998);

  /* Take two steps at third order. */

  /* Advance the old velocity data because we have just doubled the
     step size. */
  for (j=0; j<N; ++j)
    {
      mblob[j].blob1.dx = mblob[j].blob2.dx;
      mblob[j].blob1.dy = mblob[j].blob2.dy;
      mblob[j].blob2.dx = mblob[j].blob4.dx;
      mblob[j].blob2.dy = mblob[j].blob4.dy;

      saveab4dx[j] = mblob[j].blob4.dx;
      saveab4dy[j] = mblob[j].blob4.dy;
    }

  for (i=0; i<2; ++i)
    {
      /* All computational elements march forward with split step. */
	
      /* Take a half step externally. */
      for (j=0; j<N; ++j)
	{
	  push(&mblob[j]);
	  ab3half(&mblob[j],&tmpparms[j],TimeStep/2.0); 
	  /* ab3(&mblob[j],&tmpparms[j],TimeStep/2.0); */
	}
	
      vel_field();
	
      /* Take a full step internally */
	
      for (j=0; j<N; ++j)
	{
	  internal_march(&(blobguts[j]),&(tmpparms[j]),TimeStep/2.0);
	}
	
      /* Take a full step externally. */
      for (j=0; j<N; ++j)
	ab3(&mblob[j],&tmpparms[j],TimeStep/2.0);

      SimTime += TimeStep/2.0;

      numk2 = 0;
      vel_field();
    }

  if (FrameStep == TimeStep)
    {
      write_vorts(3);
      Frame = 4;
    }
  else
    write_vorts(9999);

  /* Advance the old velocity data because we have just doubled the
     step size. */
  for (j=0; j<N; ++j)
    {
      mblob[j].blob1.dx = mblob[j].blob2.dx;
      mblob[j].blob1.dy = mblob[j].blob2.dy;
      mblob[j].blob2.dx = mblob[j].blob4.dx;
      mblob[j].blob2.dy = mblob[j].blob4.dy;
      mblob[j].blob3.dx = saveab4dx[j];
      mblob[j].blob3.dy = saveab4dy[j];

      mblob[j].blob0.order = MaxOrder;
    }

  /* Now we can hit the trail at fourth order. */
}

void init(int argc, char *argv[]) 
{
  char      sim_name[Title],control_name[Title],temp[Title],*p1;
  int       i;

  /* Defaults */
  dtth_delta = 1.0e-3;
  MaxOrder   = 3;
   
  /*Allocate memory for multipole coefficients.*/
  for (i=0; i<LMax; ++i)
    Level_Ptr[i] = 
      malloc(sizeof(complex)*PMax*((int) ldexp(1.0,2*(i+1))));
 
  p1 = getenv("ECCSVM_HOME");
  if (p1 == NULL)
    p1 = getenv("PWD");

  else strcpy(temp,p1);
   
  p1 = getenv("ECCSVM_BASENAME");
  if (p1 == NULL) 
    sprintf(filename,"%s/%s",temp,"dipoles");
  else 
    sprintf(filename,"%s/%s",temp,p1);

#ifdef MULTIPROC
  /* set up log files with process number embedded in names */
  sprintf(temp,"%s_comp.%d.log",filename, rank);   
  comp_log = fopen(temp,"w");   
  sprintf(temp,"%s_diag.%d.log",filename, rank);
  diag_log = fopen(temp,"w");   
  sprintf(temp,"%s_mpi.%d.log",filename, rank);
  mpi_log = fopen(temp,"w");   
  sprintf(temp,"%s_cpu.%d.log",filename, rank);
  cpu_log = fopen(temp,"w");   
#else 
  /*  log file names will NOT have process rank # embedded  */
  sprintf(temp,"%s_comp.log",filename);
  comp_log = fopen(temp,"w");
  sprintf(temp,"%s_diag.log",filename);
  diag_log = fopen(temp,"w");
  sprintf(temp,"%s_cpu.log",filename);
  cpu_log = fopen(temp,"w");
#endif

  fprintf(diag_log,"Initializing.\n");
   
  fprintf(comp_log,"ECCSVM base: %s\n",filename);
   
  sprintf(sim_name,"%s%s",filename,".sim");
  sprintf(control_name,"%s%s",filename,".ctl");

  /* Read in the simulation parameters. */

  fprintf(diag_log,"Reading simfile...\n");
  read_sim();
   
  fprintf(diag_log,"Reading ctlfile...\n");
  read_ctl();

  fprintf(diag_log,"Finished reading files.\n");
   
  fprintf(diag_log,"Initializing elements...\n");

  /* Go ahead and merge redundant elements from the initial conditions. */
  mplevels = Set_Level();

  fprintf(diag_log,"Initial mplevel: %d\n",mplevels);

  partition(mplevels);
   
#ifdef MERGEDIAG
  write_pmvorts(0);
#endif

  merge();
  resort();
   
  Release_Links(mplevels);

#ifdef XANTISYMM
  fprintf(comp_log,"X anti-symmetry imposed.\n");
#endif 

#ifdef CORRECTVEL4
  fprintf(comp_log,"Using fourth order curvature corrections.\n");
#endif 

   
#ifdef LINEAR
  fprintf(comp_log,"Operating as a particle tracking code with a known velocity field\n");
#endif 
   
#ifdef MULTIPROC
  fprintf(mpi_log,"Detected %d processes.\n",total_processes);
  if (total_processes == 2)
    fprintf(mpi_log,"Running peer algorithm.\n");
  else
    {
      fprintf(mpi_log,"Running master/slave algorithm.\n");
      if (rank == 0)
	fprintf(mpi_log,"The master is alive and kicking.\n");
      else
	fprintf(mpi_log,"Slave #%d is alive and kicking.\n",rank);
    }
  
#endif

  fflush(mpi_log);
  fflush(comp_log);
  fflush(diag_log);

#ifndef MERGECHECK
  vel_field();
#endif

  for (i=0; i<N; ++i) 
    { 
      mblob[i].blob1 = mblob[i].blob0;
      mblob[i].blob2 = mblob[i].blob0;
      mblob[i].blob3 = mblob[i].blob0;
      mblob[i].blob4 = mblob[i].blob0;
    }

  /*chkvel2();*/
   
#ifdef MULTIPROC
  /* Have only the 'root' node do the writes */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {   
    write_vorts(0);
    write_partition(0);
  }
#else 
  /* there is only one processor anyway */
  write_vorts(0);
  write_partition(0);
#endif

#ifdef MERGECHECK
  printf("Merge events: %d\n",nmerge);
  exit(1);
#endif

  Frame = 1;
  MergeFrame=1;
  TimeStep = PrefStep;
  SimTime = 0.0;

  startup_new();

  nsplit   = 0;
  refinestack = 0;
  nmerge   = 0;
  totsplit = 0;
  totmerge = 0;
}
