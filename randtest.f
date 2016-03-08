      program randtest
      
      real rnum

      do 10 i=1,10
         rnum = ranf()
         write(*,*) rnum
 10      continue


      stop
      end
