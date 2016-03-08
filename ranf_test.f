      PROGRAM random
     
      real ranf,time
      integer t1,t2,cnt_rt

      call system_clock(t1,cnt_rt)
      call random_seed()
      do i = 0,1000
      call random_number(ranf)
      write(*,*) 'random number...',ranf
      enddo
      call system_clock(t2,cnt_rt)

      time = (real(t2) - real(t1))/real(cnt_rt)

      write(*,*) 'elapsed time...',(t2-t1)/cnt_rt,cnt_rt,time


      stop 
      end



















