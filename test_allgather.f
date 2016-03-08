C-----------------------------------------------------------------------------
C MPI tutorial example code: Collective Communications
C FILE: mpi_scatter.f
C AUTHOR: Blaise Barney
C LAST REVISED:
C-----------------------------------------------------------------------------
      program scatter
      include 'mpif.h'

      integer size
      parameter(size = 4)
      integer numtasks, rank, sendcount, recvcount, source, ierr
      integer count
      real a
      real sendbuf, recvbuf(size)


      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)

 
c      a = 1
      a = rank
      write(*,*) 'my rank, a...',rank,a
      

c      write(*,*) 'Barrier...',rank,a
c      call MPI_Barrier(MPI_COMM_WORLD,ierr)
c      write(*,*) 'All processors synchronized...',rank,a

      source = 0
      sendcount = 1
      recvcount = 1
      count = 1
      if (numtasks .eq. size) then
         call MPI_ALLGATHER(a,1,MPI_REAL,recvbuf,1,MPI_REAL,
     x        MPI_COMM_WORLD,ierr)
c     call MPI_ALLREDUCE(a,recvbuf,count,
c     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
         write(*,*) 'my rank, recvbuf...',rank,recvbuf
      else
         print *, 'Must specify',SIZE,' processors.  Terminating.' 
      endif

c      if (numtasks .eq. SIZE) then
c         source = 1
c         sendcount = SIZE
c         recvcount = SIZE
c         call MPI_SCATTER(sendbuf, sendcount, MPI_REAL, recvbuf, 
c     &        recvcount, MPI_REAL, source, MPI_COMM_WORLD, ierr)
c         print *, 'rank= ',rank,' Results: ',recvbuf 
c      else
c         print *, 'Must specify',SIZE,' processors.  Terminating.' 
c      endif

      call MPI_FINALIZE(ierr)

      end
