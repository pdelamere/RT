      PROGRAM MAIND


      integer m
      real np(10)

c----------------------------------------------------------------------
c Initialize diagnostic output files
c----------------------------------------------------------------------
      open(110,file='test.dat',status='unknown',
     x         form='unformatted')


      m = 5
      do i = 1,10
         np(i) = i
      enddo

      write(110) m
      write(110) np

      close(110)

      stop 
      end



















