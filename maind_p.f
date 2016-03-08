      PROGRAM MAIND
     
c----------------------------------------------------------------------
c maind.f
c Main driver for an ionospheric chemical release hybrid code
c simulation.  This version adds the electron pressure term to the
c electron momentum equation.  Diffusion is based on density
c gradients at a fixed temperature.  Temperature is specified
c in the simulation parameters.
c----------------------------------------------------------------------

      include 'incurv.h'

c----------------------------------------------------------------------
c Listing of all declared variables
c
c Note that the format for specifying variables at different time
c levels is based upon the number of 1/2 time steps that the varible
c is behind the current time level.  For example, 
c uf2 is 1 time step behind uf, ufp2 is 1 time step ahead of uf,
c and b12 is 1 full time step behind b1 (not 12 1/2 steps behind b 
c since b does not exist...right). b1p2 is an exception...this is a
c temporary holder for b1 at m+1 in the predictor/corrector update
c of the magnetic field.
c----------------------------------------------------------------------
      integer time, t1, t2    !keep track of run time
      external time

      real b0(nz),            !ambient magnetic field
     x     b1(nx,ny,nz,3),    !1st order magnetic field
     x     b12(nx,ny,nz,3),   !b1 at previous time step
     x     b1p2(nx,ny,nz,3),  !temporary b1 at time level m+1
     x     bt(nx,ny,nz,3),    !total magnetic field..mc covarient
     x     btmf(nx,ny,nz,3),  !main cell contravarient bt field
     x     btc(nx,ny,nz,3),   !btmf at cell center for particle move
     x     nf(nx,ny,nz),      !ambient fixed fluid density
     x     nf1(nx,ny,nz),     !nf at n-1/2
     x     nn(nx,ny,nz),      !neutral cloud density
     x     nnd(nx,ny,nz),     !neutral cloud density decrement
     x     np(nx,ny,nz),      !particle ion den at time level n, n+1/2
     x     vp(Ni_max,3),      !particle velocity at t level n+1/2
     x     vp1(Ni_max,3),     !particle velocity at t level n
     x     vplus(Ni_max,3),   !v+ used in velocity update
     x     vminus(Ni_max,3),  !v- used in velocity update
     x     up(nx,ny,nz,3),    !particle flow at time level n, n+1/2
     x     xp(Ni_max,3),      !coordinates of ion particles
     x     uf(nx,ny,nz,3),    !fluid velocity
     x     uf2(nx,ny,nz,3),   !fluid velcity at time level n-1
     x     ufp1(nx,ny,nz,3),  !fluid velocity at time level n+1/2
     x     ufp2(nx,ny,nz,3),  !fluid velocity at time level n+1
     x     ui(nx,ny,nz,3),    !total ion flow velocity
     x     aj(nx,ny,nz,3),    !curlB/(alpha*n) 
     x     nu(nx,ny,nz),      !collision frequency
     x     Ep(Ni_max,3),      !Ion particle electric field
     x     Ef(nx,ny,nz,3),    !fluid electric field
     x     E(nx,ny,nz,3),     !electric field from electron mom eqn
     x     uplus(nx,ny,nz,3), !u plus used in velocity update
     x     uminus(nx,ny,nz,3) !u minus used in velocity update

      real Evp,       !total particle kinetic energy
     x     Euf,       !total fluid kinetic energy
     x     EB1,       !total magnetic field energy
     x     EB1x,      !total b1x energy
     x     EB1y,      !total b1y energy
     x     EB1z,      !total b1z energy
     x     EE,        !total electric field energy
     x     EeP        !total electron pressure energy

      real pup(3),      !total particle momentum
     x     puf(3),      !total fluid momentum
     x     peb(3),      !total momentum carried by E and B fields
     x     input_p(3)   !input momentum

      real chex_rate
      real satnp
      real gradP(nx,ny,nz,3)
      real etemp(nx,ny,nz)
      real ugradu(nx,ny,nz,3)
      real minnf,maxnf
      real divu(nx,ny,nz)
c----------------------------------------------------------------------

      t1 = time()
      seed = float(t1)

c----------------------------------------------------------------------
c Initialize all variables
c----------------------------------------------------------------------
      Ni_tot = 0
      mstart = 0
      ndiag = 0
      prev_Etot = 1.0
      bndry_Eflux = 0.0
      chex_rate = 0.0
      minnf = 10.0e20
      maxnf = 10.0e20

      if (.not.(restart)) then
         do 66 i=1,nx
            do 66 j=1,ny
               do 66 k=1,nz
                  nf(i,j,k) = 10.0e20
                  nf1(i,j,k) = 10.0e20
 66               continue
         endif

c      call grd4()
c      call setup(b0,bt,nf)        !initialize b0, nf
  
      call pgrd_reg()
c      call pgrd_irr()
      call psetup(nf,b0,bt,nu)
c      call grd6()
c      call grd7()
c      call grd6_setup(nf,b0,bt,nu)
c      call cov_to_contra(bt,btmf) 
c      call fgrd_reg()
c      call fsetup(nf,b0,bt)
c      call f_init(b0,b1,b12,uf,uf2,nf) 
c----------------------------------------------------------------------


c----------------------------------------------------------------------
c check for restart flag
c----------------------------------------------------------------------
      write(*,*) 'restart status....',restart
      if (restart) then 
         write(*,*) 'opening restart.vars......'
         open(210,file='restart.vars',status='unknown',
     x            form='unformatted')
         write(*,*) 'reading restart.vars......'
         read(210) b1,b12,b1p2,bt,btmf,nn,np,nf,vp,vp1,vplus,vminus,
     x            up,xp,uf,uf2,ufp2,aj,Ep,Ef,E,uplus,uminus,Evp,Euf,
     x            EB1,EB1x,EB1y,EB1z,EE,EeP,input_E,Ni_tot,
     x            ijkp,mstart,input_p,input_EeP,prev_Etot
         write(*,*) 'restarting hybrid.....'
         endif
c----------------------------------------------------------------------


c----------------------------------------------------------------------
c Initialize diagnostic output files
c----------------------------------------------------------------------
      call assign('assign -F system -N ultrix f:' //'c.npall.dat')
      open(110,file='c.npall.dat',status='unknown',access='sequential',
     x         form='unformatted')

      call assign('assign -F system -N ultrix f:' //'c.nfall.dat')
      open(115,file='c.nfall.dat',status='unknown',access='sequential',
     x         form='unformatted')

      call assign('assign -F system -N ultrix f:' //'c.ufall.dat')
      open(120,file='c.ufall.dat',status='unknown',access='sequential',
     x         form='unformatted')

      call assign('assign -F system -N ultrix f:' //'c.b1all.dat')
      open(130,file='c.b1all.dat',status='unknown',access='sequential',
     x         form='unformatted')

      call assign('assign -F system -N ultrix f:' //'c.ajall.dat')
      open(140,file='c.ajall.dat',status='unknown',access='sequential',
     x         form='unformatted')

      call assign('assign -F system -N ultrix f:' //'c.Eall.dat')
      open(150,file='c.Eall.dat',status='unknown',access='sequential',
     x         form='unformatted')

      call assign('assign -F system -N ultrix f:' //'c.energy.dat')
      open(160,file='c.energy.dat',status='unknown',
     x         access='sequential',form='unformatted')

      call assign('assign -F system -N ultrix f:' //'c.chex.dat')
      open(170,file='c.chex.dat',status='unknown',
     x         access='sequential',form='unformatted')

      call assign('assign -F system -N ultrix f:' //'c.satnp.dat')
      open(175,file='c.satnp.dat',status='unknown',
     x         access='sequential',form='unformatted')

      call assign('assign -F system -N ultrix f:' //'c.up.dat')
      open(180,file='c.up.dat',status='unknown',
     x         access='sequential',form='unformatted')

      call assign('assign -F system -N ultrix f:' //'c.momentum.dat')
      open(190,file='c.momentum.dat',status='unknown',
     x         access='sequential',form='unformatted')

      call assign('assign -F system -N ultrix f:' //'c.p_conserve.dat') 
      open(192,file='c.p_conserve.dat',status='unknown',                
     x         access='sequential',form='unformatted')                  

c      call assign('assign -F system -N ultrix f:' //'c.ugradu.dat')
c      open(300,file='c.ugradu.dat',status='unknown',
c     x         access='sequential',form='unformatted')

      call assign('assign -F system -N ultrix f:' //'c.Ef.dat')
      open(305,file='c.Ef.dat',status='unknown',
     x         access='sequential',form='unformatted')

c      call assign('assign -F system -N ultrix f:' //'c.ui.dat')
c      open(310,file='c.ui.dat',status='unknown',
c     x         access='sequential',form='unformatted')

c      call assign('assign -F system -N ultrix f:' //'c.uf2.dat')
c      open(320,file='c.uf2.dat',status='unknown',
c     x         access='sequential',form='unformatted')

c      call assign('assign -F system -N ultrix f:' //'c.ufp2.dat')
c      open(330,file='c.ufp2.dat',status='unknown',
c     x         access='sequential',form='unformatted')

      call assign('assign -F system -N ultrix f:' //'c.divu.dat')
      open(340,file='c.divu.dat',status='unknown',
     x         access='sequential',form='unformatted')

c----------------------------------------------------------------------


c======================================================================
c  MAIN LOOP!
c======================================================================

      do 1 m = mstart+1, nt

         write(*,*) ' '
         write(*,*) 'time...', m, m*dt

         !Calculate neutral density
c         call Neut_Den(nn,nnd,m)
c         call charge_exchange(np,xp,vp,vp1,m,chex_rate,
c     x                        input_p)
         
         !Ionize cloud and calculate ion density
         call Ionize2(np,vp,vp1,xp,m,input_p,up)
         write(*,*) 'Ni_tot...',Ni_tot !/beta
         write(*,*) ' '
         write(*,*) 'nu....',nu(1,1,1)

c         call get_interp_weights(xp)
c         call update_np(np)

c         call get_gradP(gradP,np,nf,m,etemp)

c         call check_np(np,nf)
c         call fix_b1z(b1)

         !energy diagnostics
c         call p_diag(up,np,vp,vn,Ep)
         call Energy_diag(vp,uf,nf,b1,E,Evp,Euf,EB1,EB1x,EB1y,EB1z,EE,
     x                    EeP,etemp,nu)

         call curlB(b1,nf,np,aj)
         call cov_to_contra(bt,btmf) 
         call face_to_center(btmf,btc)       !interp bt to cell center
         call extrapol_up(up,vp,vp1,np)    
         call get_Ep(Ep,aj,np,nf,up,uf,btc,nu,gradP)
         call get_vplus_vminus(Ep,btc,vp,vplus,vminus)
         call improve_up(vp1,vplus,vminus,up,np)

         call get_Ep(Ep,aj,np,nf,up,uf,btc,nu,gradP)
         call get_vplus_vminus(Ep,btc,vp,vplus,vminus)
         call get_vp_final(Ep,vp,vp1,vplus)

         call move_ion_half(xp,vp,vp1)  !1/2 step ion move to n+1/2
         call get_interp_weights(xp)
         call update_np(np)             !np at n+1/2
         call update_up(vp,np,up)       !up at n+1/2
c**********************************************************************
c SUBCYCLING LOOP!
c**********************************************************************

      do 2 n = 1, ntsub

c         write(*,*) 'subcycle step...',n

c         !convert main cell covarient bt to main cell contravarient
cc         call f_diag(uf,nf,b1,E)
c         call cov_to_contra(bt,btmf) 
c         call curlB(b1,nf,np,aj)     

cc         !update collision frequency
cc         call Col_Freq(nu,aj)

c         !update fluid velocity, uf 

cc only need predict_uf when calculating ugradu

c         call predict_uf(Ef,b0,b1,b12,uf,uf2,ufp2,nu,np,nf,uplus, 
c     x                   uminus,ugradu,up)
c         call correct_uf(Ef,btmf,uf,uf2,ufp2,nu,np,nf,uplus,uminus, 
c     x                   ugradu,aj,up,ufp1)
c         call trans_nf_1(nf,nf1,uf,divu)  
c         call trans_nf_2(nf,nf1,ufp1)

c         !update magnetic field, b1
c         call predict_B(b1,b12,b1p2,bt,btmf,E,aj,up,uf,uf2,np,nf,nu) 

c         call correct_B(b0,b1,b1p2,E,aj,up,uf,np,nf,nu)
c         call f_update_tlev(uf,uf2,b1,b12,b1p2,bt,b0)

         call Momentum_diag(up,uf,np,nf,E,b1,pup,puf,peb,input_p)
c         call check_momentum(uf,nf,b1,ugradu)
c         write(192) m
c         write(192) n
c         write(192) surf_tot, graduu_tot, ugradu_tot

c         call get_bndry_Eflux(uf,b1,E,nf)

 2     continue
c**********************************************************************

         call move_ion_half(xp,vp,vp1)  !final ion move to n+1
c         call get_interp_weights(xp)
c         call update_np(np)
c         call Neut_Den(nn,m)

c         call get_chex_rate(np,nn,up,m,chex_rate)
         write(170) chex_rate
         write(*,*) ' '
         write(*,*) 'chex_rate...',chex_rate

         call get_np_at_sat(np,m,satnp)
         write(175) satnp
         write(*,*) ' '
         write(*,*) 'np at sat...',satnp

         call get_ui(uf,nf,up,np,ui)

         write(*,*) 'Momentum conservation...'
         write(*,*) '  Particles.............',pup(1),pup(2),pup(3)
         write(*,*) '  Fluid.................',puf(1),puf(2),puf(3)
         write(*,*) '  ExB...................',peb(1),peb(2),peb(3)
         write(*,*) '  Normalized............',
     x                        (pup(1)+puf(1)+peb(1))/input_p(1),
     x                        (pup(2)+puf(2)+peb(2))/input_p(2),
     x                        (pup(3)+puf(3)+peb(3))/input_p(3)


c----------------------------------------------------------------------
c diagnostic output
c----------------------------------------------------------------------
         write(160) m
         write(160) input_E,input_EeP,Evp,Euf,EB1,EB1x,EB1y,EB1z,EE,EeP
         write(190) m
         write(190) pup, puf, peb, input_p


         ndiag = ndiag + 1
         if (ndiag .eq. nout) then
c         if (m .ge. 275) then
            write(110) m
            write(110) np
c            write(110) nn
            write(115) m
            write(115) nf
            write(120) m
            write(120) uf
            write(130) m
            write(130) b1
c            write(130) btmf
            write(140) m
            write(140) aj
            write(150) m
            write(150) E
c            write(150) Ef
            write(180) m
            write(180) up
c            write(300) m
c            write(300) ugradu
            write(305) m
            write(305) Ef
c            write(310) m
c            write(310) ui
c            write(320) m
c            write(320) uf2
c            write(330) m
c            write(330) ufp2
            write(340) m
            write(340) divu
            ndiag = 0
            endif

c----------------------------------------------------------------------


c----------------------------------------------------------------------
c Write restart file
c----------------------------------------------------------------------
       if (m .eq. mrestart) then
          open(220,file='restart.vars.new',status='unknown',
     x             form='unformatted')
          write(220) b1,b12,b1p2,bt,btmf,nn,np,nf,vp,vp1,vplus,vminus,
     x            up,xp,uf,uf2,ufp2,aj,Ep,Ef,E,uplus,uminus,Evp,Euf,
     x            EB1,EB1x,EB1y,EB1z,EE,EeP,input_E,Ni_tot,
     x            ijkp,mrestart,input_p,input_EeP,prev_Etot
          endif
c----------------------------------------------------------------------


 1     continue
c======================================================================


       close(110)
       close(115)
       close(120)
       close(130)
       close(140)
       close(150)
       close(160)
       close(170)
       close(175)
       close(180)
       close(190)
       close(192)
       close(210)
       close(220)
c       close(300)
       close(305)
c       close(310)
c       close(320)
c       close(330)
       close(340)

c       call poutput(nn,np,Ep,up)
c       call foutput(uf,b1,bt,aj)
c       call Output(nn,np,vn,Ep,up)

       t2 = time()
       write(*,*) 
       write(*,*) 'Elapsed time....',t2-t1,' sec'
       write(*,*)

       stop 
       end



















