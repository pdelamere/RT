c----------------------------------------------------------------------
c converts 3d scalar variables from cray to IDL f77_unformatted
c compatable arrays
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      subroutine convert_np()
c----------------------------------------------------------------------

      integer*4 m
      real np(10)

      open(10,file='test.dat',
     x        form='unformatted',status='unknown')
     
      open(20,file='test1.dat',status='unknown',
     x        form='unformatted')

         write(*,*) 'reading m...'
         read(10) m
         write(20) m
         read(10) np
         write(20) np
 30      continue

      return
      end
c----------------------------------------------------------------------




c----------------------------------------------------------------------
      program main
c----------------------------------------------------------------------

      call convert_np()

      stop
      end
c----------------------------------------------------------------------






