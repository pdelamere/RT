      PROGRAM random
     
      integer*2 n
      character filenum
      character(30) filename
      character(2) a(16) 


      a = (/'1 ','2 ','3 ','4 ','5 ','6 ','7 ','8 ',
     x     '9 ','10','11','12','13','14','15','16'/)


      filename = 'c.xp_'

      do i = 1,16
         filename = 'c.xp_'
         write(*,*) trim(a(i))
         filename = trim(filename)//trim(a(i))//".dat"
         write(*,*) 'filename...',filename
      enddo


      stop 
      end




