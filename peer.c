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

double wicked_big_vectah[NMax * PARTICLE_DATA_PACKET_SIZE];

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
	     wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*count]   = 
	       mblob[i].blob0.dx; 
	     wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*count+1] = 
	       mblob[i].blob0.dy;
	     wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*count+2] = 
	       tmpparms[i].du11 ; 
	     wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*count+3] = 
	       tmpparms[i].du12 ;   
	     wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*count+4] = 
	       tmpparms[i].du21 ;
	     wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*count+5] = 
	       tmpparms[i].u_xx ;
	     wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*count+6] = 
	       tmpparms[i].u_xy ;
	     wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*count+7] = 
	       tmpparms[i].u_yy ;
	     wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*count+8] = 
	       tmpparms[i].v_xx ;
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
	
	mblob[dummy].blob0.dx = 
	  wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*i];
	mblob[dummy].blob0.dy = 
	  wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*i+1];
	tmpparms[dummy].du11  = 
	  wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*i+2];
	tmpparms[dummy].du12  = 
	  wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*i+3];
	tmpparms[dummy].du21  = 
	  wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*i+4];
	tmpparms[dummy].u_xx  = 
	  wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*i+5];
	tmpparms[dummy].u_xy  = 
	  wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*i+6];
	tmpparms[dummy].u_yy  = 
	  wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*i+7];
	tmpparms[dummy].v_xx  = 
	  wicked_big_vectah[PARTICLE_DATA_PACKET_SIZE*i+8];

	++dummy;
     }
}


#endif
