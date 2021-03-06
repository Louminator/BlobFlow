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
 * Louis Rossi                                                          *
 * Department of Mathematical Sciences                                  *
 * University of Delaware                                               *
 * Newark, DE 19716-2553                                                */

/* Initialization subroutines for eflow. */

#include "global_min.h"
#include "particle.h"
#include "biot-savart.h"
#include "boundary.h"
#include "field_interp.h"
#include "global_matrices.h"

#include <sys/time.h>
#include <sys/resource.h>

#ifdef MULTIPROC
#include "multiproc.h"
#endif

/* This is a band-aid for the PGI compiler which seems to hit a stack
   limitation at least on our opteron system.  Using -Mbounds
   illustrates the problem.

   Efforts to use dynamic memory allocation also fail though it is not
   at all clear why.  */

/* Script constants. */
static char *inputs[] =
{
  "FrameStep:",
  "EndTime:",
  "Viscosity:",

  "VtxInit:",
  "GrdInit:",
  "GrdX0:",
  "GrdX1:",
  "GrdY0:",
  "GrdY1:",
  "GrdNumPts:",
  "TimeStep:",
  "BdyInit:",
  "XAntiSymmetry",

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
  "MergeGrowthRate:",

  "InterpStep:",
  "InterpPopulationControl:",
  "InterpVar:",
  "BoundaryStep:",

  "WriteVtx:",
  "WritePartitions:",
  "WriteVel:"
};

enum ScriptItems
  {FRAME_STEP,END_TIME,VISCOSITY,VTX_INIT,
   GRD_INIT,GRD_X0,GRD_X1,GRD_Y0,GRD_Y1,GRD_NUM_PTS,
   TIME_STEP,BDY_INIT,XANTISYMM,
   MAX_ORDER,DTTH_DELTA,L2_TOL,ALPHA,SPLIT_METHOD,MERGE_ERROR_ESTIMATOR,
   MERGE_BUDGET,MERGE_STEP,CLUSTER_RADIUS,MERGE_A2_TOL,MERGE_TH_TOL,
   MERGE_MOM3_WT,MERGE_MOM4_WT,MERGE_C,MERGE_GROWTH_RATE,
   INTERP_STEP,INTERP_POPULATION_CONTROL,INTERP_VAR,
   BOUNDARY_STEP,WRITE_VTX,WRITE_PARTITIONS,WRITE_VEL
  };

/* A routine to remove stack size limits. */

void unlimit_stack(void)
{

  struct rlimit rlim = { RLIM_INFINITY, RLIM_INFINITY };

  if ( setrlimit(RLIMIT_STACK, &rlim) == -1 )
    {
      perror("setrlimit error");
      exit(1);
    }
}

void check_sim(char *vtxfilename,char *grdfilename)
{
  if ( (FrameStep <= 0.0) ||
       (EndTime <= 0.0) ||
       (visc < 0.0) ||
       ( (strcmp(vtxfilename,"") == 0) &&
	 (strcmp(grdfilename,"") == 0) ) )
    {
      printf("Simulation file error.\n");

      if (FrameStep <= 0.0)
	printf("FrameStep not set properly.\n");

      if (EndTime <= 0.0)
	printf("EndTime not set properly.\n");

      if (visc < 0.0)
	printf("Viscosity not set properly.\n");

      if ( (strcmp(vtxfilename,"") == 0) &&
	   (strcmp(grdfilename,"") == 0) )
	printf("No initialization file set.\n");

      exit(-1);
    }
}

void bailout_simfile(int i)
{
  printf("Fatal error: Cannot read sim file %s value.\n",inputs[i]);
  exit(-1);
};

void read_sim(char inputdir[])
{
  char      sim_name[FILENAME_LEN],vtxfilename[FILENAME_LEN]="",
    grdfilename[FILENAME_LEN]="",
    bdyfilename[FILENAME_LEN]="",temp[FILENAME_LEN],
    xantisymmyn[5];
  FILE      *sim_file,*bdy_file;
  char      *word;
  int       i,k,l,count,scan_test,inputs_size;
  double    grdX0=0.0,grdX1=0.0,grdY0=0.0,grdY1=0.0,h,*circs,x,y;
  int       gridn=0;

  sprintf(sim_name,"%s%s",filename,".sim");
  sim_file = fopen(sim_name,"r");

  if(sim_file == NULL)
    {
      printf("Fatal error: Sim file cannot be found.\n");
      exit(-1);
    };
   
  /* ASSUME a 32-bit word. */
  inputs_size = sizeof(inputs)/4;
  word = malloc(100*sizeof(char));

  if (fscanf(sim_file,"%s",word) == 0)
    {
      printf("Fatal error: Cannot read sim file keyword.\n");
      exit(-1);
    };

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
		      if (fscanf(sim_file,"%lf",&FrameStep) != 1)
			bailout_simfile(i);
		      i=inputs_size+1;
		      break;
		    case END_TIME:
		      if (fscanf(sim_file,"%lf",&EndTime) != 1)
			bailout_simfile(i);
		      i=inputs_size+1;
		      break;
		    case VISCOSITY:
		      if (fscanf(sim_file,"%lf",&visc) != 1)
			bailout_simfile(i);
		      i=inputs_size+1;
		      break;
		    case VTX_INIT:
		      if (fscanf(sim_file,"%s",vtxfilename) != 1)
			bailout_simfile(i);
		      i=inputs_size+1;
		      break;
		    case GRD_INIT:
		      if (fscanf(sim_file,"%s",grdfilename) != 1)
			bailout_simfile(i);
		      i=inputs_size+1;
		      break;
		    case GRD_X0:
		      if (fscanf(sim_file,"%lf",&grdX0) != 1)
			bailout_simfile(i);
		      i=inputs_size+1;
		      break;
		    case GRD_X1:
		      if (fscanf(sim_file,"%lf",&grdX1) != 1)
			bailout_simfile(i);
		      i=inputs_size+1;
		      break;
		    case GRD_Y0:
		      if (fscanf(sim_file,"%lf",&grdY0) != 1)
			bailout_simfile(i);
		      i=inputs_size+1;
		      break;
		    case GRD_Y1:
		      if (fscanf(sim_file,"%lf",&grdY1) != 1)
			bailout_simfile(i);
		      i=inputs_size+1;
		      break;
		    case GRD_NUM_PTS:
		      if (fscanf(sim_file,"%d",&gridn) != 1)
			bailout_simfile(i);
		      i=inputs_size+1;
		      break;
		    case BDY_INIT:
		      if (fscanf(sim_file,"%s",bdyfilename) != 1)
			bailout_simfile(i);
		      i=inputs_size+1;
		      break;
		    case XANTISYMM:
		      if (fscanf(sim_file,"%s",xantisymmyn) != 1)
			bailout_simfile(i);
		      if (strcmp(xantisymmyn,"y") ||
			  strcmp(xantisymmyn,"yes") ||
			  strcmp(xantisymmyn,"Y") ||
			  strcmp(xantisymmyn,"Yes"))
			xantisymm = 1;
		      i=inputs_size+1;
		      break;
		    default:
		      printf("Uh oh.  There is an error in the scripting code.\n");
		    }
		}
	    if (i==inputs_size)
	      printf("Unrecognized input variable: %s\n",word);
	  }
      if (fscanf(sim_file,"%s",word) == 0)
	{
	  printf("Fatal error: Cannot read sim file keyword.\n");
	  exit(-1);
	};
    }
  fclose(sim_file);

  check_sim(vtxfilename,grdfilename);

  if (strcmp(vtxfilename,"") != 0)
    {
      sprintf(temp,"%s/%s",inputdir,vtxfilename);

      sim_file = fopen(temp,"r");

      if (fscanf(sim_file,"%lf%lf%lf%lf%lf%lf",
		 &mblob[0].blob0.x,
		 &mblob[0].blob0.y,
		 &mblob[0].blob0.strength,&blobguts[0].s2,
		 &blobguts[0].a2,&blobguts[0].th) != 6)
	{
	  printf("Fatal error reading initial data.\n");
	  exit(-1);
	}

      set_blob(&(blobguts[0]),&(tmpparms[0]));

      N = 1;

      while (!feof(sim_file))
	{
	  scan_test = fscanf(sim_file,"%lf%lf%lf%lf%lf%lf",
			     &mblob[N].blob0.x,
			     &mblob[N].blob0.y,
			     &mblob[N].blob0.strength,&blobguts[N].s2,
			     &blobguts[N].a2,&blobguts[N].th);
	  if ( (scan_test != 6) && (scan_test > 0) )
	    {
	      printf("Fatal error reading initial data. ");
	      printf("Read %d numbers.\n",scan_test);
	      exit(-1);
	    }
	  set_blob(&(blobguts[N]),&(tmpparms[N]));
	  ++N;
	}
      --N;

      fclose(sim_file);
    }

  if (strcmp(grdfilename,"") != 0)
    {
      h = (grdX1-grdX0)/gridn;

      if (fabs((h-(grdY1-grdY0)/gridn)/h) > 1.0e-6)
	{
	  fprintf(diag_log,"The initialization grid must be square.\n");
	  exit(-1);
	}

      circs = malloc(sizeof(double)*SQR(gridn));

      sprintf(temp,"%s/%s",inputdir,grdfilename);

      sim_file = fopen(temp,"r");

      for (i=0; i<gridn*gridn; ++i)
	{
	  if (fscanf(sim_file,"%lf",(circs+i)) != 1)
	    {
	      printf("Error: Expected %d values in %s, but only read %d.",
		     SQR(gridn),grdfilename,i);
	      exit(-1);
	    }
	  *(circs+i) *= h*h;
	}

      fclose(sim_file);

      h = RHE_RK4_grd_init(grdX0,grdX1,grdY0,grdY1,gridn,circs);
      count = 0;

      for (k=0; k<gridn; ++k)
	for (l=0; l<gridn; ++l)
	  {
	    if (fabs(circs[l*gridn+k])/h/h>InterpPopulationControl)
	      {
		x = grdX0+(k+0.5)*h;
		y = grdY0+(l+0.5)*h;

		mblob[count].blob0.x = x;
		mblob[count].blob0.y = y;
		mblob[count].blob0.strength = circs[l*gridn+k]/2/M_PI;
		blobguts[count].s2 = h*h;
		blobguts[count].a2 = 1.0;
		blobguts[count].th = 0.0;

		set_blob(&(blobguts[count]),&(tmpparms[count]));

		++count;

		if (count==NMAX)
		  {
		    fprintf(diag_log,"Out of memory in RHE remesh.\n");
		    fprintf(diag_log,"Time to sleep.\n");
		  }
	      }
	    N = count;
	  }
      free(circs);
    }

  if (!(strcmp(bdyfilename,"") == 0))
    {
      sprintf(temp,"%s/%s",inputdir,bdyfilename);

      bdy_file = fopen(temp,"r");

      if (fscanf(bdy_file,"%lf%lf%lf%lf%lf",
	     &walls[0].x,
	     &walls[0].y,
	     &walls[0].nx,
	     &walls[0].ny,
	     &walls[0].l)!= 5)
	{
	  printf("Fatal error reading boundary data.\n");
	  exit(-1);
	};

      B = 1;

      while (!feof(bdy_file))
	{
	  scan_test = fscanf(bdy_file,"%lf%lf%lf%lf%lf",
			     &walls[B].x,
			     &walls[B].y,
			     &walls[B].nx,
			     &walls[B].ny,
			     &walls[B].l);
	  if ( (scan_test !=5 ) && (scan_test > 0) )
	    {
	      printf("Fatal error reading boundary data.\n");
	      exit(-1);
	    };
	  ++B;
	}
      --B;

      fclose(bdy_file);
    }

  fprintf(comp_log,"Simulation parameters.\n");
  fprintf(comp_log,"%s %12.4e\n",inputs[FRAME_STEP],FrameStep);
  fprintf(comp_log,"%s %12.4e\n",inputs[END_TIME],EndTime);
  fprintf(comp_log,"%s %12.4e\n",inputs[VISCOSITY],visc);
  if (strcmp(vtxfilename,"") == 0)
    fprintf(comp_log,"%s %s\n",inputs[VTX_INIT],vtxfilename);
  if (strcmp(grdfilename,"") == 0)
    fprintf(comp_log,"%s %s\n",inputs[GRD_INIT],grdfilename);
  fprintf(comp_log,"%d vortices read from %s.\n",N,temp);
  if (B != 0)
    fprintf(comp_log,"%d wall elements read from %s.\n",B,temp);
  fflush(comp_log);
}

void check_ctl()
{
  if ( (PrefStep <= 0.0) ||
       /* These are no longer relevant. */
       /*
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
	 ( (merge_a2tol < 0.0) || (merge_thtol < 0.0) ) ) ||
       */
       ( (InterpStep > 0.0) && (InterpVar <= 0.0 ) ) )
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
      if ( (InterpStep > 0.0) && (InterpVar <= 0.0) )
	printf("InterpVar not set properly.\n");
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
}

void bailout_ctlfile(int i)
{
  printf("Fatal error: Cannot read ctl file %s value.\n",inputs[i]);
  exit(-1);
};

void read_ctl()
{
  FILE *control_file;
  char *word,control_name[FILENAME_LEN];
  int  i,inputs_size;

  BoundaryStep = -1.0;
  InterpStep   = -1.0;
  Interps      = 0;

  write_vtx = 1;
  write_partitions = 0;
  write_vel=0;

  sprintf(control_name,"%s%s",filename,".ctl");

  /* ASSUME a 32-bit word. */
  inputs_size = sizeof(inputs)/4;

  word = malloc(100*sizeof(char));

  control_file = fopen(control_name,"r");

  if(control_file == NULL)
    {
      printf("Fatal error: ctl file cannot be found.\n");
      exit(-1);
    };

  if (fscanf(control_file,"%s",word) == 0)
    {
      printf("Fatal error: Cannot read ctl file keyword.\n");
      exit(-1);
    };

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
		      if (fscanf(control_file,"%lf",&PrefStep) != 1)
			bailout_ctlfile(i);
		      i=inputs_size+1;
		      break;
		    case MAX_ORDER:
		      if (fscanf(control_file,"%d",&MaxOrder) != 1)
			bailout_ctlfile(i);
		      fprintf(diag_log,
			      "Deprecation: This version does not use MAX_ORDER.\n");
		      fprintf(diag_log,
			      "Deprecation: If you are using a set of parameter files for the older BlobFlow\n");
		      fprintf(diag_log,
			      "Adams family integrators, you should consider increasing the timestep by a\n");
		      fprintf(diag_log,
			      "factor of 4 or more.  A factor of 8 should achieve similar precision.\n");
		      i=inputs_size+1;
		      break;
		    case DTTH_DELTA:
		      if (fscanf(control_file,"%lf",&dtth_delta) != 1)
			bailout_ctlfile(i);
		      fprintf(diag_log,
			      "Deprecated: %s not used in this version.\n",
			      inputs[i]);
		      i=inputs_size+1;
		      break;
		    case L2_TOL:
		      if (fscanf(control_file,"%lf",&l2tol) != 1)
			bailout_ctlfile(i);
		      i=inputs_size+1;
		      break;
		    case ALPHA:
		      if (fscanf(control_file,"%lf",&alpha) != 1)
			bailout_ctlfile(i);
		      fprintf(diag_log,
			      "Deprecated: %s not used in this version.\n",
			      inputs[i]);
		      i=inputs_size+1;
		      break;
		    case SPLIT_METHOD:
		      if (fscanf(control_file,"%d",&split_method) != 1)
			bailout_ctlfile(i);
		      i=inputs_size+1;
		      break;
		    case MERGE_ERROR_ESTIMATOR:
		      if (fscanf(control_file,"%d",&merge_estimator) != 1)
			bailout_ctlfile(i);
		      fprintf(diag_log,
			      "Deprecated: %s not used in this version.\n",
			      inputs[i]);
		      i=inputs_size+1;
		      break;
		    case MERGE_BUDGET:
		      if (fscanf(control_file,"%lf",&merge_budget) != 1)
			bailout_ctlfile(i);
		      fprintf(diag_log,
			      "Deprecated: %s not used in this version.\n",
			      inputs[i]);
		      i=inputs_size+1;
		      break;
		    case CLUSTER_RADIUS:
		      if (fscanf(control_file,"%lf",&clusterR) != 1)
			bailout_ctlfile(i);
		      fprintf(diag_log,
			      "Deprecated: %s not used in this version.\n",
			      inputs[i]);
		      i=inputs_size+1;
		      break;
		    case MERGE_STEP:
		      if (fscanf(control_file,"%lf",&MergeStep) != 1)
			bailout_ctlfile(i);
		      fprintf(diag_log,
			      "Deprecated: %s not used in this version.\n",
			      inputs[i]);
		      i=inputs_size+1;
		      break;
		    case MERGE_A2_TOL:
		      if (fscanf(control_file,"%lf",&merge_a2tol) != 1)
			bailout_ctlfile(i);
		      fprintf(diag_log,
			      "Deprecated: %s not used in this version.\n",
			      inputs[i]);
		      i=inputs_size+1;
		      break;
		    case MERGE_TH_TOL:
		      if (fscanf(control_file,"%lf",&merge_thtol) != 1)
			bailout_ctlfile(i);
		      fprintf(diag_log,
			      "Deprecated: %s not used in this version.\n",
			      inputs[i]);
		      i=inputs_size+1;
		      break;
		    case MERGE_MOM3_WT:
		      if (fscanf(control_file,"%lf",&merge_mom3wt) != 1)
			bailout_ctlfile(i);
		      fprintf(diag_log,
			      "Deprecated: %s not used in this version.\n",
			      inputs[i]);
		      i=inputs_size+1;
		      break;
		    case MERGE_MOM4_WT:
		      if (fscanf(control_file,"%lf",&merge_mom4wt) != 1)
			bailout_ctlfile(i);
		      fprintf(diag_log,
			      "Deprecated: %s not used in this version.\n",
			      inputs[i]);
		      i=inputs_size+1;
		      break;
		    case MERGE_C:
		      if (fscanf(control_file,"%lf",&merge_c) != 1)
			bailout_ctlfile(i);
		      fprintf(diag_log,
			      "Deprecated: %s not used in this version.\n",
			      inputs[i]);
		      i=inputs_size+1;
		      break;
		    case MERGE_GROWTH_RATE:
		      if (fscanf(control_file,"%lf",&merge_growth_rate) != 1)
			bailout_ctlfile(i);
		      fprintf(diag_log,
			      "Deprecated: %s not used in this version.\n",
			      inputs[i]);
		      i=inputs_size+1;
		      break;
		    case INTERP_STEP:
		      if (fscanf(control_file,"%lf",&InterpStep) != 1)
			bailout_ctlfile(i);
		      i=inputs_size+1;
		      break;
		    case INTERP_POPULATION_CONTROL:
		      if (fscanf(control_file,"%lf",
				 &InterpPopulationControl) != 1)
			bailout_ctlfile(i);
		      i=inputs_size+1;
		      break;
		    case INTERP_VAR:
		      if (fscanf(control_file,"%lf",&InterpVar) != 1)
			bailout_ctlfile(i);
		      i=inputs_size+1;
		      break;
		    case BOUNDARY_STEP:
		      if (fscanf(control_file,"%lf",&BoundaryStep) != 1)
			bailout_ctlfile(i);
		      i=inputs_size+1;
		      break;
		    case WRITE_VTX:
		      if (fscanf(control_file,"%d",&write_vtx) != 1)
			bailout_ctlfile(i);
		      i=inputs_size+1;
		      break;
		    case WRITE_PARTITIONS:
		      if (fscanf(control_file,"%d",&write_partitions) != 1)
			bailout_ctlfile(i);
		      i=inputs_size+1;
		      break;
		    case WRITE_VEL:
		      if (fscanf(control_file,"%d",&write_vel) != 1)
			bailout_ctlfile(i);
		      i=inputs_size+1;
		      break;
		    default:
		      printf("Uh oh.  There is an error in the scripting code.\n");
		    }
		}
	    if (i==inputs_size)
	      printf("Unrecognized input variable: %s\n",word);
	  }
      if (fscanf(control_file,"%s",word) == 0)
	{
	  printf("Fatal error: Cannot read ctl file keyword.\n");
	  exit(-1);
	};

    }
  fclose(control_file);

  check_ctl();

  fprintf(comp_log,"Approximation control parameters.\n");
  fprintf(comp_log,"%s %12.4e\n",inputs[TIME_STEP],PrefStep);
  fflush(comp_log);
}

void init(int argc, char *argv[])
{
  char      sim_name[FILENAME_LEN],control_name[FILENAME_LEN],
    temp[FILENAME_LEN],inputdir[FILENAME_LEN],config[FILENAME_LEN];
  int       i,inputdirread=0,configread=0,domdirread=0;

  unlimit_stack();

  xantisymm = 0;

#ifdef MULTIPROC
  if (rank == 0)
    {
#endif
      if (argc > 1)
	{
	  for (i=1; i<argc; ++i)
	    {
	      if (strcmp(argv[i],"-inputdir") == 0)
		{
		  if (i < argc-1)
		    {
		      ++i;
		      strcpy(inputdir,argv[i]);
		      if(inputdir[(strlen(inputdir)-1)] != '/')
			{
			  inputdir[strlen(inputdir)+1] = 
			    inputdir[strlen(inputdir)];
			  inputdir[(strlen(inputdir))] = '/';
			}
		      inputdirread = 1;
		    }
		}

	      if (strcmp(argv[i],"-config") == 0)
		{
		  if (i < argc-1)
		    {
		      ++i;
		      strcpy(config,argv[i]);
		      configread = 1;
		    }
		}

	      if (strcmp(argv[i],"-domdir") == 0)
		{
		  if (i < argc-1)
		    {
		      ++i;
		      strcpy(datarootname,argv[i]);
		      if(datarootname[(strlen(datarootname)-1)] != '/')
			{
			  datarootname[strlen(datarootname)+1] = 
			    datarootname[strlen(datarootname)];
			  datarootname[(strlen(datarootname))] = '/';
			}
		      domdirread = 1;
		    }
		}
	    }

	  if (inputdirread*configread*domdirread == 0)
	    {
	      if (inputdirread == 0)
		printf("No input directory provided.\n");
	      if (configread == 0)
		printf("No configuration name provided.\n");
	      if (domdirread == 0)
		printf("No dom data directory provided.\n");
	      exit(-1);
	    }
	}
      else
	{
	  printf("Usage: eflow -inputdir <dir> -config <config name> -domdir <dir>\n");
	  exit(0);
	}
#ifdef MULTIPROC
    }

  MPI_Bcast(inputdir,FILENAME_LEN,MPI_CHAR,0,MPI_COMM_WORLD);
  MPI_Bcast(config,FILENAME_LEN,MPI_CHAR,0,MPI_COMM_WORLD);
  MPI_Bcast(datarootname,FILENAME_LEN,MPI_CHAR,0,MPI_COMM_WORLD);
#endif

  /* Defaults */
  dtth_delta = 1.0e-3;
  MaxOrder   = 3;

  /*Allocate memory for multipole coefficients.*/
  for (i=0; i<LMAX; ++i)
    Level_Ptr[i] =
      malloc(sizeof(Complex)*PMAX*((int) ldexp(1.0,2*(i+1))));

  sprintf(filename,"%s%s",inputdir,config);

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

  fprintf(diag_log,"Reading spectral coefficients\n");

  read_data();

  fprintf(comp_log,"ECCSVM base: %s\n",filename);

  sprintf(sim_name,"%s%s",filename,".sim");
  sprintf(control_name,"%s%s",filename,".ctl");

  /* Read in the simulation parameters. */

  fprintf(diag_log,"Reading ctlfile...\n");
  read_ctl();

  /* Read in the simulation parameters. */

  fprintf(diag_log,"Reading simfile...\n");
  read_sim(inputdir);

  fprintf(diag_log,"Finished reading files.\n");

  fprintf(diag_log,"Initializing elements...\n");

  /* Go ahead and merge redundant elements from the initial conditions. */
  mplevels = Set_Level();

  fprintf(diag_log,"Initial mplevel: %d\n",mplevels);

  /*
  partition(mplevels);

  RHE_interp(3.125e-4,1.0e-6);

  write_vorts(9990);

  exit(-1);

  Release_Links(mplevels);
  */

  if (B != 0)
    {
      fprintf(comp_log,"Building boundary matrix.\n");
      factor_bdy_matrix(walls,Bpiv,BdyMat);
    }

  if (xantisymm)
    fprintf(comp_log,"X anti-symmetry imposed.\n");

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

  if (write_vel) write_vels(0);
#ifdef MULTIPROC
  /* Have only the 'root' node do the writes */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) 
    {
      if (write_vtx) write_vorts(0);
      if (write_partitions) write_partition(0);
    }
#else
  /* there is only one processor anyway */
  if (write_vtx) write_vorts(0);
  if (write_partitions) write_partition(0);
#endif

  Frame = 1;
  MergeFrame=1;
  BoundaryFrame=1;
  TimeStep = PrefStep;
  SimTime = 0.0;

  /*
  if (BoundaryStep>0.0)
    {
      BoundaryConstrain();
    }
  else
    fprintf(comp_log,"No boundary conditions.\n");
  */

  nsplit   = 0;
  nmerge   = 0;
  totsplit = 0;
  totmerge = 0;
  Interps  = 1;

  fprintf(comp_log,"Initialization complete.\n");

  fflush(mpi_log);
  fflush(comp_log);
  fflush(diag_log);
}
