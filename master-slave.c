/* MS.C */
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

#include "global.h"
#include "multiproc.h"
#define WorkSize 100

#ifdef XANTISYMM
#define CARDINALITY  N/2
#else
#define CARDINALITY  N
#endif

double wicked_big_vectah[NMax * PARTICLE_DATA_PACKET_SIZE];

/*  the master-slave approach yields dynamic load balancing, so if certain 
    nodes on the multicomputer are busier, they are given less work  */

void master ( void ) 
{
   int proc,vort,job,flag1;
   int workbuf[30][WorkSize],i;
   double rbuf[30][PARTICLE_DATA_PACKET_SIZE*WorkSize];
   MPI_Request mpireqs[30];
   
   int position,msgsize,membersize,packsize;
   char *buffer[30];
   
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
	 workbuf[proc][i] = job+i;
       MPI_Send (workbuf[proc],WorkSize,MPI_INT,proc,WORK,MPI_COMM_WORLD);
       job += WorkSize;
    }

  for (proc=1; proc<total_processes; proc++)
    {
       for (i=0; i<WorkSize; ++i)
	 workbuf[proc][i] = job+i;
       MPI_Send (workbuf[proc],WorkSize,MPI_INT,proc,WORK,MPI_COMM_WORLD);
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
   while (job < CARDINALITY) 
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
	MPI_Unpack(buffer[proc],msgsize,&position,workbuf[proc],WorkSize,
		   MPI_INT,MPI_COMM_WORLD);
	MPI_Unpack(buffer[proc],msgsize,&position,rbuf[proc],
		   WorkSize*PARTICLE_DATA_PACKET_SIZE,
		   MPI_DOUBLE,MPI_COMM_WORLD);
	
	for (i=0; i<WorkSize; ++i)
	  {
	     vort = workbuf[proc][i];
	     
	     if (vort != -1)
	       {
		  mblob[vort].blob0.dx = 
		    rbuf[proc][PARTICLE_DATA_PACKET_SIZE*i+0];
		  mblob[vort].blob0.dy = 
		    rbuf[proc][PARTICLE_DATA_PACKET_SIZE*i+1];
		  tmpparms[vort].du11  = 
		    rbuf[proc][PARTICLE_DATA_PACKET_SIZE*i+2];
		  tmpparms[vort].du12  = 
		    rbuf[proc][PARTICLE_DATA_PACKET_SIZE*i+3];
		  tmpparms[vort].du21  = 
		    rbuf[proc][PARTICLE_DATA_PACKET_SIZE*i+4];
		  tmpparms[vort].u_xx  = 
		    rbuf[proc][PARTICLE_DATA_PACKET_SIZE*i+5];
		  tmpparms[vort].u_xy  =
		    rbuf[proc][PARTICLE_DATA_PACKET_SIZE*i+6];
		  tmpparms[vort].u_yy  = 
		    rbuf[proc][PARTICLE_DATA_PACKET_SIZE*i+7];
		  tmpparms[vort].v_xx  = 
		    rbuf[proc][PARTICLE_DATA_PACKET_SIZE*i+8];
	       }
	     
	     if (job < CARDINALITY)
	       workbuf[proc][i] = job;
	     else
	       workbuf[proc][i] = -1;
	     ++job;
	  }
	
	MPI_Start(&(mpireqs[proc-1]));
	MPI_Send(workbuf[proc],WorkSize,MPI_INT,proc,WORK,
		 MPI_COMM_WORLD);
     }
   
   /* RECEIVE OUTSTANDING RESULTS */

   for (proc=1; proc<total_processes; ++proc)
     {
	MPI_Wait(&(mpireqs[proc-1]),&mpistatus);
	
	position = 0;
	MPI_Get_count(&mpistatus,MPI_PACKED,&msgsize);
	MPI_Unpack(buffer[proc],msgsize,&position,workbuf[proc],WorkSize,
		   MPI_INT,MPI_COMM_WORLD);
	MPI_Unpack(buffer[proc],msgsize,&position,rbuf[proc],
		   WorkSize*PARTICLE_DATA_PACKET_SIZE,
		   MPI_DOUBLE,MPI_COMM_WORLD);
	
	for (i=0; i<WorkSize; ++i)
	  {
	     vort = workbuf[proc][i];
	     
	     if (vort != -1)
	       {
		  mblob[vort].blob0.dx = 
		    rbuf[proc][PARTICLE_DATA_PACKET_SIZE*i+0];
		  mblob[vort].blob0.dy = 
		    rbuf[proc][PARTICLE_DATA_PACKET_SIZE*i+1];
		  tmpparms[vort].du11  = 
		    rbuf[proc][PARTICLE_DATA_PACKET_SIZE*i+2];
		  tmpparms[vort].du12  = 
		    rbuf[proc][PARTICLE_DATA_PACKET_SIZE*i+3];
		  tmpparms[vort].du21  = 
		    rbuf[proc][PARTICLE_DATA_PACKET_SIZE*i+4];
		  tmpparms[vort].u_xx  = 
		    rbuf[proc][PARTICLE_DATA_PACKET_SIZE*i+5];
		  tmpparms[vort].u_xy  =
		    rbuf[proc][PARTICLE_DATA_PACKET_SIZE*i+6];
		  tmpparms[vort].u_yy  = 
		    rbuf[proc][PARTICLE_DATA_PACKET_SIZE*i+7];
		  tmpparms[vort].v_xx  = 
		    rbuf[proc][PARTICLE_DATA_PACKET_SIZE*i+8];
	       }
	  }
	
	MPI_Start(&(mpireqs[proc-1]));
     }
   
   for (proc=1; proc<total_processes; ++proc)
     {
	MPI_Wait(&(mpireqs[proc-1]),&mpistatus);
	
	position = 0;
	MPI_Get_count(&mpistatus,MPI_PACKED,&msgsize);
	MPI_Unpack(buffer[proc],msgsize,&position,workbuf[proc],WorkSize,
		   MPI_INT,MPI_COMM_WORLD);
	MPI_Unpack(buffer[proc],msgsize,&position,rbuf[proc],
		   WorkSize*PARTICLE_DATA_PACKET_SIZE,
		   MPI_DOUBLE,MPI_COMM_WORLD);
	
	for (i=0; i<WorkSize; ++i)
	  {
	     vort = workbuf[proc][i];
	     
	     if (vort != -1)
	      {
		 mblob[vort].blob0.dx = 
		   rbuf[proc][PARTICLE_DATA_PACKET_SIZE*i+0];
		 mblob[vort].blob0.dy = 
		   rbuf[proc][PARTICLE_DATA_PACKET_SIZE*i+1];
		 tmpparms[vort].du11  = 
		   rbuf[proc][PARTICLE_DATA_PACKET_SIZE*i+2];
		 tmpparms[vort].du12  = 
		   rbuf[proc][PARTICLE_DATA_PACKET_SIZE*i+3];
		 tmpparms[vort].du21  = 
		   rbuf[proc][PARTICLE_DATA_PACKET_SIZE*i+4];
		 tmpparms[vort].u_xx  = 
		   rbuf[proc][PARTICLE_DATA_PACKET_SIZE*i+5];
		 tmpparms[vort].u_xy  =
		   rbuf[proc][PARTICLE_DATA_PACKET_SIZE*i+6];
		 tmpparms[vort].u_yy  = 
		   rbuf[proc][PARTICLE_DATA_PACKET_SIZE*i+7];
		 tmpparms[vort].v_xx  = 
		   rbuf[proc][PARTICLE_DATA_PACKET_SIZE*i+8];
	      }
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
   int vort,start,workbuf1[WorkSize],workbuf2[WorkSize],i;  
   double sbuf[PARTICLE_DATA_PACKET_SIZE*WorkSize];
   MPI_Request master_req,slave_resp;
   int position,membersize,packsize;
   char *buffer;
   
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
   
  start = 1;

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
	    
#ifdef NOFASTMP
	    dpos_vel(vort);
#else
	    dpos_vel_fast(vort);
#endif
	    
	    /* Better wait here just in case I compute faster than I thought. */
	    /* sbuf must be protected. */
	    if (!start)
	      {
		 MPI_Wait(&slave_resp,&mpistatus);
	      }
	    else
	      start = 0;

	    sbuf[i*PARTICLE_DATA_PACKET_SIZE]   = mblob[vort].blob0.dx; 
	    sbuf[i*PARTICLE_DATA_PACKET_SIZE+1] = mblob[vort].blob0.dy;
	    sbuf[i*PARTICLE_DATA_PACKET_SIZE+2] = tmpparms[vort].du11 ;
	    sbuf[i*PARTICLE_DATA_PACKET_SIZE+3] = tmpparms[vort].du12 ;
	    sbuf[i*PARTICLE_DATA_PACKET_SIZE+4] = tmpparms[vort].du21 ;
	    sbuf[i*PARTICLE_DATA_PACKET_SIZE+5] = tmpparms[vort].u_xx ;
	    sbuf[i*PARTICLE_DATA_PACKET_SIZE+6] = tmpparms[vort].u_xy ;
	    sbuf[i*PARTICLE_DATA_PACKET_SIZE+7] = tmpparms[vort].u_yy ;
	    sbuf[i*PARTICLE_DATA_PACKET_SIZE+8] = tmpparms[vort].v_xx ;
	    workbuf2[i] = workbuf1[i];
	 }

       position=0;
       MPI_Pack(workbuf2,WorkSize,MPI_INT,buffer,packsize,&position,
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

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  if (rank == 0) {
    for (j=0; j<CARDINALITY; j++) 
      {
	wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*j]   = 
	  mblob[j].blob0.dx; 
	wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*j+1] = 
	  mblob[j].blob0.dy;
	wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*j+2] = 
	  tmpparms[j].du11 ; 
	wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*j+3] = 
	  tmpparms[j].du12 ;   
	wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*j+4] = 
	  tmpparms[j].du21 ;
	wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*j+5] = 
	  tmpparms[j].u_xx ;
	wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*j+6] = 
	  tmpparms[j].u_xy ;
	wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*j+7] = 
	  tmpparms[j].u_yy ;
	wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*j+8] = 
	  tmpparms[j].v_xx ;
      }
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
      for (j=0; j<CARDINALITY; j++)
      {	    
	mblob[j].blob0.dx = wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*j];
	mblob[j].blob0.dy = wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*j+1];
	tmpparms[j].du11  = wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*j+2];
	tmpparms[j].du12  = wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*j+3];
	tmpparms[j].du21  = wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*j+4];
	tmpparms[j].u_xx  = wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*j+5];
	tmpparms[j].u_xy  = wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*j+6];
	tmpparms[j].u_yy  = wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*j+7];
	tmpparms[j].v_xx  = wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*j+8];
      }
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
