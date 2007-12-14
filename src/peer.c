/* PEER.C */

#ifdef MULTIPROC

#include "global.h"
#include "multiproc.h"
#define SmallWork 5

#ifdef XANTISYMM
#define CARDINALITY  N/2
#else
#define CARDINALITY  N
#endif

void peer(void)
{
   int workqueue[NMax],worksize,donelist[NMax];
   int i,peer_proc,dummy,more_flag,msg,count,outcount,done_flag;
   
   MPI_Request mpireq;
   
   for (i=0; i<CARDINALITY; ++i)
     donelist[i] = 0;
   
   if (rank == 0)
     {
	for (i=0; i<CARDINALITY/2; ++i)
	  workqueue[i] = i;
	worksize = CARDINALITY/2;
	peer_proc=1;
     }
   
   if (rank == 1)
     {
	for (i=CARDINALITY/2; i<CARDINALITY; ++i)
	  workqueue[i-CARDINALITY/2] = i;
	worksize = CARDINALITY-CARDINALITY/2;
	peer_proc=0;
     }
   
   done_flag = 0;
   
   while (worksize>0)
     {
	--worksize;
	dpos_vel_fast(workqueue[worksize]);
	donelist[workqueue[worksize]] = 1;

	MPI_Iprobe(peer_proc,MORE,MPI_COMM_WORLD,&more_flag,&mpistatus);
	
	if (more_flag)
	  {
	     MPI_Recv(&dummy,1,MPI_INT,peer_proc,MORE,
		      MPI_COMM_WORLD,&mpistatus);
	     if (worksize>SmallWork)
	       {
		  msg = worksize/2; /* Here's some stuff */
		  MPI_Send(&msg,1,MPI_INT,peer_proc,WORK,MPI_COMM_WORLD);
		  MPI_Send(&(workqueue[worksize-msg]),msg,MPI_INT,peer_proc,
			   WORK,MPI_COMM_WORLD);
		  worksize -= msg;
	       }
	     else
	       {
		  MPI_Send(&dummy,1,MPI_INT,peer_proc,MORE,MPI_COMM_WORLD);
		  done_flag = 1;
	       }
	  }
	
	if ( (worksize == 0) && (!done_flag) )
	  {
	     MPI_Send(&dummy,1,MPI_INT,peer_proc,MORE,MPI_COMM_WORLD);
	     MPI_Recv(&msg,1,MPI_INT,peer_proc,MPI_ANY_TAG,
		      MPI_COMM_WORLD,&mpistatus);
	     if (mpistatus.MPI_TAG == WORK)
	       {
		  MPI_Recv(workqueue,msg,MPI_INT,peer_proc,WORK,
			   MPI_COMM_WORLD,&mpistatus);
		  worksize = msg;
	       }
	  }
     }
   
   count = 0;
   
   for (i=0; i<CARDINALITY; ++i)
     {
	if (donelist[i])
	  {
	    blob_to_buffer(&(mblob[i].blob0),&(tmpparms[i]),
			   &(wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*count]));
	     ++count;
	  }
     }
   outcount=count;
   MPI_Send(&outcount,1,MPI_INT,peer_proc,DONE,MPI_COMM_WORLD);
   MPI_Recv(&count,1,MPI_INT,peer_proc,DONE,MPI_COMM_WORLD,&mpistatus);
   MPI_Isend(wicked_big_vectah,outcount*PARTICLE_DATA_PACKET_SIZE,
	     MPI_DOUBLE,peer_proc,DONE,MPI_COMM_WORLD,&mpireq);
   MPI_Recv(&(wicked_big_vectah[outcount*PARTICLE_DATA_PACKET_SIZE]),
	    count*PARTICLE_DATA_PACKET_SIZE,MPI_DOUBLE,peer_proc,
	    DONE,MPI_COMM_WORLD,&mpistatus);
   
   dummy = 0;
   for (i=outcount; i<outcount+count; ++i)
     {
	while (donelist[dummy])
	  ++dummy;
	
	buffer_to_blob(&(wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*i]),
		       &(mblob[dummy].blob0),
		       &(tmpparms[dummy]));
	++dummy;
     }
}


#endif
