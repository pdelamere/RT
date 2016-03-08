C *****************************************************************************
C MPI Example - Array Assignment - Fortran Version
C FILE: mpi_array.f
C OTHER FILES: make.mpi_array.f
C DESCRIPTION:
C   In this simple example, the master task initiates numtasks-1 number of
C   worker tasks.  It then distributes an equal portion of an array to each
C   worker task.  Each worker task receives its portion of the array, and
C   performs a simple value assignment to each of its elements. The value
C   assigned to each element is simply that element's index in the array+1.
C   Each worker task then sends its portion of the array back to the master
C   task.  As the master receives back each portion of the array, selected
C   elements are displayed.
C AUTHOR: Blaise Barney. Converted to MPI: George L. Gusciora (1/23/95)
C LAST REVISED: 12/14/95 Blaise Barney
C **************************************************************************

      program array 
      include 'mpif.h'

      integer   ARRAYSIZE, MASTER
      parameter (ARRAYSIZE = 600000)
      parameter (MASTER = 0)

      integer  numtasks, numworkers, taskid, dest, index, i, ierr,
     &         arraymsg, indexmsg, source, chunksize
      integer status(MPI_STATUS_SIZE)
      real*4   data(ARRAYSIZE), result(ARRAYSIZE)

C ************************ initializations ***********************************
C Find out how many tasks are in this partition and what my task id is.  Then
C define the number of worker tasks and the array partition size as chunksize.
C Note:  For this example, the MP_PROCS environment variable should be set
C to an odd number...to insure even distribution of the array to numtasks-1
C worker tasks.
C *****************************************************************************

      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, taskid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numtasks, ierr )

      write(*,*)'taskid =',taskid
      numworkers = numtasks-1
      chunksize = (ARRAYSIZE / numworkers)
      arraymsg = 1
      indexmsg = 2

C *************************** master task *************************************
      if (taskid .eq. MASTER) then
      print *, '*********** Starting MPI Example 1 ************'

C     Initialize the array
      do 20 i=1, ARRAYSIZE 
        data(i) =  0.0
 20   continue

C     Send each worker task its portion of the array
      index = 1
      do 30 dest=1, numworkers
        write(*,*) 'Sending to worker task', dest
        call MPI_SEND( index, 1, MPI_INTEGER, dest, indexmsg, 
     &                 MPI_COMM_WORLD, ierr )
c        call MPI_SEND( data(index), chunksize, MPI_REAL, dest, arraymsg,
c     &                 MPI_COMM_WORLD, ierr )
        index = index + chunksize
 30   continue

C     Now wait to receive back the results from each worker task and print 
C     a few sample values 
      do 40 i=1, numworkers
        source = i
        call MPI_RECV( index, 1, MPI_INTEGER, source, indexmsg,
     &                 MPI_COMM_WORLD, status, ierr )
        call MPI_RECV( data(index), chunksize, MPI_REAL, source,
     &                 arraymsg, MPI_COMM_WORLD, status, ierr )
        print *, '---------------------------------------------------'
        print *, 'MASTER: Sample results from worker task ', source
        print *, '   result[', index, ']=', data(index)
        print *, '   result[', index+100, ']=', data(index+100)
        print *, '   result[', index+1000, ']=', data(index+1000)
        print *, ' '
 40   continue

      print *, 'MASTER: All Done!' 
      endif

C *************************** worker task ************************************
      if (taskid .gt. MASTER) then
C     Receive my portion of array from the master task */

        call MPI_RECV( index, 1, MPI_INTEGER, MASTER, indexmsg,
     &                 MPI_COMM_WORLD, status, ierr )
c        call MPI_RECV( data(index), chunksize, MPI_REAL, MASTER,
c     &                 arraymsg, MPI_COMM_WORLD, status, ierr )

C     Do a simple value assignment to each of my array elements
        do 50 i=index, index + chunksize
          data(i) = i + 1
 50     continue

C     Send my results back to the master
        call MPI_SEND( index, 1, MPI_INTEGER, MASTER, indexmsg,
     &                 MPI_COMM_WORLD, ierr )
        call MPI_SEND( data(index), chunksize, MPI_REAL, MASTER,
     &                 arraymsg, MPI_COMM_WORLD, ierr)

      endif
      call MPI_FINALIZE(ierr)
      end
