/* MASTER-SLAVE.C */
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
 * Newark, DE 19716                                                    */

#ifdef MULTIPROC

#include "global_min.h"
#include "particle.h"
#include "biot-savart.h"
#include "multiproc.h"
#define WorkSize 25

/*  the master-slave approach yields dynamic load balancing, so if certain 
    nodes on the multicomputer are busier, they are given less work  */

void blob_to_buffer(Blob_external *blob,Blob_parms *parm,double *buffer)
{
  *buffer     = (*blob).dx;
  *(buffer+1) = (*blob).dy;
  *(buffer+2) = (*parm).du11;
  *(buffer+3) = (*parm).du12;
  *(buffer+4) = (*parm).du21;
  *(buffer+5) = (*parm).u_xx;
  *(buffer+6) = (*parm).u_xy;
  *(buffer+7) = (*parm).u_yy;
  *(buffer+8) = (*parm).v_xx;
}

void buffer_to_blob(double *buffer,Blob_external *blob,Blob_parms *parm)
{
  (*blob).dx   = *buffer;
  (*blob).dy   = *(buffer+1);
  (*parm).du11 = *(buffer+2);
  (*parm).du12 = *(buffer+3);
  (*parm).du21 = *(buffer+4);
  (*parm).u_xx = *(buffer+5);
  (*parm).u_xy = *(buffer+6);
  (*parm).u_yy = *(buffer+7);
  (*parm).v_xx = *(buffer+8);
}

/* It would be nice to dynamically allocate buffers based on
   total_processes.  Also, it would be nice to do something
   intelligent if the number of cpus is huge.  */

void master ( void ) 
{
   int proc,vort,job,flag1;
   int workbuf[WorkSize],i;
   double rbuf[MAX_CPUS][PARTICLE_DATA_PACKET_SIZE*WorkSize];
   MPI_Request mpireqs[MAX_CPUS];
   
   int position,msgsize,membersize,packsize;
   char *buffer[MAX_CPUS];
   int joblimit;
   
   /* Calculate packed buffer size */
   
   MPI_Pack_size(WorkSize,MPI_INT,MPI_COMM_WORLD,&membersize);
   packsize=membersize;
   MPI_Pack_size(PARTICLE_DATA_PACKET_SIZE*WorkSize,
		 MPI_DOUBLE,MPI_COMM_WORLD,&membersize);
   packsize += membersize;
   
   for (i=0; i<total_processes; ++i)
     buffer[i] = malloc(packsize);
   
  /*  INITIALIZE EACH SLAVE WITH FIRST CALCULATIONS */

  /* I assume that N>2*total_processes*WorkSize.  */
  job=0;
  for (proc=1; proc<total_processes; proc++)
    {
       for (i=0; i<WorkSize; ++i)
	 workbuf[i] = job+i;
       MPI_Send (workbuf,WorkSize,MPI_INT,proc,WORK,MPI_COMM_WORLD);
       job += WorkSize;
    }

  for (proc=1; proc<total_processes; proc++)
    {
       for (i=0; i<WorkSize; ++i)
	 workbuf[i] = job+i;
       MPI_Send (workbuf,WorkSize,MPI_INT,proc,WORK,MPI_COMM_WORLD);
       job += WorkSize;
    }

  /* Set up a pipeline of persistent send/recv's. */

   for (proc=1; proc<total_processes; proc++)
     MPI_Recv_init(buffer[proc],packsize,MPI_PACKED,proc,MPI_ANY_TAG,
		   MPI_COMM_WORLD,&(mpireqs[proc-1]));

   /* Start listening */
   MPI_Startall(total_processes-1,mpireqs);

  /* RECEIVE AND DISPATCH */

   proc = 0;

   if (xantisymm)
     joblimit = N/2;
   else
     joblimit = N;

   while (job < joblimit) 
     {
	/* MPI_Waitany(total_processes-1,mpireqs[0],&proc,&mpistatus); */
	
	flag1 = 0;
	while (!(flag1)) 
	  {
	     ++proc;
	     if (proc==total_processes)
	       proc = 1;
	     
	     MPI_Test(&(mpireqs[proc-1]),&flag1,&mpistatus);
	  }
	
	position = 0;
	MPI_Get_count(&mpistatus,MPI_PACKED,&msgsize);
	MPI_Unpack(buffer[proc],msgsize,&position,workbuf,WorkSize,
		   MPI_INT,MPI_COMM_WORLD);
	MPI_Unpack(buffer[proc],msgsize,&position,rbuf[proc],
		   WorkSize*PARTICLE_DATA_PACKET_SIZE,
		   MPI_DOUBLE,MPI_COMM_WORLD);
	
	for (i=0; i<WorkSize; ++i)
	  {
	     vort = workbuf[i];
	     
	     if (vort != -1)
	       buffer_to_blob(&(rbuf[proc][PARTICLE_DATA_PACKET_SIZE*i]),
			      &(mblob[vort].blob0),&(tmpparms[vort]));
	     if (job < joblimit)
	       workbuf[i] = job;
	     else
	       workbuf[i] = -1;
	     ++job;
	  }
	
	MPI_Start(&(mpireqs[proc-1]));
	MPI_Send(workbuf,WorkSize,MPI_INT,proc,WORK,
		 MPI_COMM_WORLD);
     }
   
   /* RECEIVE OUTSTANDING RESULTS */

   for (proc=1; proc<total_processes; ++proc)
     {
	MPI_Wait(&(mpireqs[proc-1]),&mpistatus);
	
	position = 0;
	MPI_Get_count(&mpistatus,MPI_PACKED,&msgsize);
	MPI_Unpack(buffer[proc],msgsize,&position,workbuf,WorkSize,
		   MPI_INT,MPI_COMM_WORLD);
	MPI_Unpack(buffer[proc],msgsize,&position,rbuf[proc],
		   WorkSize*PARTICLE_DATA_PACKET_SIZE,
		   MPI_DOUBLE,MPI_COMM_WORLD);
	
	for (i=0; i<WorkSize; ++i)
	  {
	    vort = workbuf[i];
	     
	     if (vort != -1)
	       buffer_to_blob(&(rbuf[proc][PARTICLE_DATA_PACKET_SIZE*i]),
			      &(mblob[vort].blob0),&(tmpparms[vort]));
	  }
	
	MPI_Start(&(mpireqs[proc-1]));
     }
   
   for (proc=1; proc<total_processes; ++proc)
     {
	MPI_Wait(&(mpireqs[proc-1]),&mpistatus);
	
	position = 0;
	MPI_Get_count(&mpistatus,MPI_PACKED,&msgsize);
	MPI_Unpack(buffer[proc],msgsize,&position,workbuf,WorkSize,
		   MPI_INT,MPI_COMM_WORLD);
	MPI_Unpack(buffer[proc],msgsize,&position,rbuf[proc],
		   WorkSize*PARTICLE_DATA_PACKET_SIZE,
		   MPI_DOUBLE,MPI_COMM_WORLD);
	
	for (i=0; i<WorkSize; ++i)
	  {
	     vort = workbuf[i];
	     
	     if (vort != -1)
	       buffer_to_blob(&(rbuf[proc][PARTICLE_DATA_PACKET_SIZE*i]),
			      &(mblob[vort].blob0),&(tmpparms[vort]));
	  }
     }
   
   for (proc=1; proc<total_processes; ++proc)
     {
	MPI_Request_free(&(mpireqs[proc-1]));
	free(buffer[proc]);
     }
   
   /* KILL THE SLAVES */  
   for (proc=1; proc<total_processes; proc++) 
     {
	MPI_Send ( 0, 0, MPI_INT, proc, DONE, MPI_COMM_WORLD ); 
     }
}

void slave ( void ) 
{
  int         vort,workbuf1[WorkSize],i,first=1;  
  double      sbuf[PARTICLE_DATA_PACKET_SIZE*WorkSize];
  MPI_Request master_req,slave_resp;
  int         position,membersize,packsize;
  char        *buffer;
   
  /* Calculate packed buffer size */
   
  MPI_Pack_size(WorkSize,MPI_INT,MPI_COMM_WORLD,&membersize);
  packsize=membersize;
  MPI_Pack_size(PARTICLE_DATA_PACKET_SIZE*WorkSize,
		MPI_DOUBLE,MPI_COMM_WORLD,&membersize);
  packsize += membersize;
   
  buffer = malloc(packsize);
   
  MPI_Recv_init (workbuf1, WorkSize, MPI_INT, 0, MPI_ANY_TAG, 
		 MPI_COMM_WORLD, &master_req );
  MPI_Send_init (buffer, packsize, MPI_PACKED, 0, 
		 rank+RESERVED_TAGS, MPI_COMM_WORLD,&slave_resp);
   
  while (1)
    {
      MPI_Start(&master_req);
      MPI_Wait(&master_req,&mpistatus);
      if  ( mpistatus.MPI_TAG == DONE )
	{
	  /* Wait for the last send to complete. */
	  MPI_Wait(&slave_resp,&mpistatus);
	  break;
	}


      for (i=0; i<WorkSize; ++i)
	{
	  vort = workbuf1[i];
	    
	  /* If vort == -1, 
	     then you do not really need to do anything. 
	     The master will ignore any data sent.  Just send
	     whatever is in sbuf.*/
	  if (vort != -1)
	    {
#ifdef NOFASTMP
	      dpos_vel(vort);
#else
	      dpos_vel_fast(vort);
#endif
	    }

	  if (vort != -1)
	    blob_to_buffer(&(mblob[vort].blob0),&(tmpparms[vort]),
			   &(sbuf[i*PARTICLE_DATA_PACKET_SIZE]));
	}

	    
      /* Better wait here just in case I compute faster than I thought. */
      /* buffer must be protected. */

      if (!first)
	{
	  MPI_Wait(&slave_resp,&mpistatus);
	  first = 0;
	}

      position=0;
      MPI_Pack(workbuf1,WorkSize,MPI_INT,buffer,packsize,&position,
	       MPI_COMM_WORLD);
      MPI_Pack(sbuf,PARTICLE_DATA_PACKET_SIZE*WorkSize,
	       MPI_DOUBLE,buffer,packsize,&position,
	       MPI_COMM_WORLD);
       
      MPI_Start(&slave_resp);
    }
   
  MPI_Request_free(&master_req);
  MPI_Request_free(&slave_resp);
  free(buffer);
}

void finish ( void ) {

  int j;
  int joblimit;

  if (xantisymm)
    joblimit = N/2;
  else
    joblimit = N;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  if (rank == 0) {
    for (j=0; j<joblimit; j++) 
      blob_to_buffer(&(mblob[j].blob0),&(tmpparms[j]),
		     &(wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*j]));
  }

  MPI_Bcast ( wicked_big_vectah, N*PARTICLE_DATA_PACKET_SIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD );

  /* If you do not trust Bcast, try this. */
  /*
    if (rank == 0)
    {
    for (j=1; j<total_processes; ++j)
    MPI_Send(wicked_big_vectah,N*PARTICLE_DATA_PACKET_SIZE,MPI_DOUBLE,j,
    RESERVED_TAGS+j,MPI_COMM_WORLD);
    }
    else
    MPI_Recv(wicked_big_vectah, N*PARTICLE_DATA_PACKET_SIZE, MPI_DOUBLE, 0, RESERVED_TAGS+rank,
    MPI_COMM_WORLD,&mpistatus);
  */
  
  if (rank != 0)  
    {
      for (j=0; j<joblimit; j++)
	buffer_to_blob(&(wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*j]),
		       &(mblob[j].blob0),&(tmpparms[j]));
    }
  
  /*  explicit synchronization before next timestep  */

  /* MPI_Barrier (MPI_COMM_WORLD); */

}

void Setup_mpi (int argc, char *argv[]) {

  /* initialize mpi */
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &total_processes);   
  if (total_processes < 2) 
    {       
      fprintf(comp_log,"Expecting at least 2 processes for execution.\n");
      stop(-100);
    }
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
}


void Cleanup_mpi (void) 
{
  MPI_Finalize();
}	

#endif
