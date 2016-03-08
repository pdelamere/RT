c----------------------------------------------------------------------
      SUBROUTINE f_update_tlev(uf,uf2,b1,b12,b1p2,bt,b0,bdp)
c loops run 1 to n since values are only being copied
c----------------------------------------------------------------------
      PARAMETER (nx = 79, ny = 19, nz = 37)
 
c grid parameters
      PARAMETER (dx = 10000.0, dy = 10000.0)   !units in km
      PARAMETER (delz = 10000.0)          !dz at release coordinates
                                       !this is not the final dz
                                       !d_lamba determines dz
c time stepping parameters
      PARAMETER (tstart = 0.2)   !starting time
      PARAMETER (dt = 1.)     !main time step
      PARAMETER (nt = 5.0)        !number of time steps
c      PARAMETER (dtsub = 0.01) !subcycle time step
      PARAMETER (dtsub_init = 0.1) !subcycle time step
      PARAMETER (ntsub = 10.0)   !number of subcycle time steps
      PARAMETER (nout = 50)      !number of steps to diagnosic output
 
c logical variable for restart
C...Translated by Crescent Bay Software VAST-F    7.0K8 15:17:00   4/ 1/04    
      logical restart
      PARAMETER (restart = .false.)
      PARAMETER (mrestart = 500)      ! use -1 for no save
      PARAMETER (mbegin = 0)      !mstart
 
c inital ambient density at the release point.
      real nf_init,b0_init,nn_coef
      PARAMETER (nf_init = 4e15)   !was 0.01
      PARAMETER (b0_init = 2e-9)    !was 0.2
      PARAMETER (nn_coef = 1e5)
 
c electron ion collision frequency
      real nu_init, eta_init
      PARAMETER (nu_init = 0.001)
      PARAMETER (eta_init = 0.0)
 
c temporary directory for output
      character*9 temp_dir
      PARAMETER (temp_dir = '/tmp/pad/')
 
c total number of neutrals released
      real No
      PARAMETER (No=1.0e28)
 
c neutral cloud expansion characteristics
      PARAMETER(vo = 20.0)
      PARAMETER(vth = 19.3)
      PARAMETER(vsat = 0.0)
      PARAMETER(vsw = 340.0)
 
c total kinetic energy of released neutrals
      PARAMETER(Etot = 2.696e7)       !units of kg m^2 / s^2
 
c density scaling parameter, alpha, and ion particle array dims
 
      real alpha
 
c      PARAMETER (beta =  8.5227271e-18)
      PARAMETER (alpha = 1.9263418e-20) !mH
 
c max number of ion particles to be produced.  This parameter
c is used to dimension the particle arrays.
      integer*4 Ni_max
      PARAMETER (Ni_max = 10000)
 
c earth's magnetic dipole field geometry
c earth radius and equatorial B alt (defining L shell)
      PARAMETER (r_earth=6.38e3)                !earth radius
      PARAMETER (B_max_alt=2.0e3+r_earth)       !ro defining L shell
 
c misc constants
      real mu0,epsilon,pi,rtod,q,mO,mBa,O_to_Ba,km_to_m,kboltz,melec
      real mproton,tempf0
      PARAMETER (pi = 3.14159)
      PARAMETER (rtod = 180.0/pi)  !radians to degreesc
      PARAMETER (mu0 = pi*4.0e-7)  !magnetic permeability of free space
      PARAMETER (epsilon = 8.85e-12) !dielectric constant
      PARAMETER (q = 1.6e-19)        !electron charge (Coulombs)
      PARAMETER (melec = 9.1e-31)  !mass of electron (kg)
      PARAMETER (mproton = 1.67e-27)
      PARAMETER (mO = 1.67e-27)    !mass of H (kg)
c      PARAMETER (mO = 2.3e-25)    !mass of Ba (kg)
      PARAMETER (mBa = 17.0*mO)    !mass of Ba (kg)
      PARAMETER (O_to_Ba = mO/mBa) !convert E and B for particle move
      PARAMETER (km_to_m = 1e3)    !km to meters conversion
      PARAMETER (kboltz = 1.38e-29)   !kg km^2 / s^2 / K
      PARAMETER (tempf0 = 50*11600.)     !K
 
      real Qo, vrad, N_o, Rp, tau_photo
      integer S_radius
      PARAMETER (Qo = 1.0e28)  !neutral source rate
      PARAMETER (vrad = 1.0)   !escape velocity
c      PARAMETER (N_o = 5e34)   !Steady state neutral particle constant
      PARAMETER (N_o = 1e33)   !Steady state neutral particle constant
      PARAMETER (Rp = 20000.0)  !Pluto radius
      PARAMETER (tau_photo = 1.2e6)
c      PARAMETER (dNi = 2500)
      PARAMETER (S_radius = 3) !units of dx
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
c raw grid coordinate data
 
      real qx(nx), qy(ny), qz(nz)
                            !cartisian coordinates of grid points
      real lambda(nz)       !latitude of given grid point
      integer ri, rj, rk    !release coordinates
      real dz_grid(nz)      !dz for grid, variable in z direction
      real dz_cell(nz)      !dz for cell, variable in z direction
 
      common /coords/ qx,qy,qz,lambda,ri,rj,rk,dz_grid,dz_cell
 
 
c Total number of ions produced at a given time during the simulation
      integer*4 Ni_tot, Nn_tot
 
c Location (indices) of particles with the grid
      integer ijkp(Ni_max,3)
      logical ionized(Ni_max)
      common /ion_info/ Ni_tot, Nn_tot, ijkp, ionized
 
      real seed
      common /rndseed/ seed
 
c Indices for boundary cells
      integer ip,im,jp,jm,kp,km
 
c Total input energy (ions)
      real input_E,input_EeP,prev_Etot, bndry_Eflux
      real input_chex, input_bill
      common /energy/ input_E, input_EeP, prev_Etot, bndry_Eflux,
     x                input_chex, input_bill
 
c Momentum conservation
      real surf_tot(3),graduu_tot(3),ugradu_tot(3)
      common /momentum/ surf_tot,graduu_tot,ugradu_tot
 
c Dummy variables for doing vector operations and particle flux
c calculation....etc.
 
      real a(nx,ny,nz,3),c(nx,ny,nz,3),ct(nx,ny,nz,3)
      common /dummy/ a,c,ct
 
c Variable time step for fluid velocity update used in both
c particle update and fluid update
 
      real dtf,ntf,dtsub
      common /timestep/ ntf,dtsub
 
c Weight variables for trilinear interpolation
 
      real wght(Ni_max,8), wquad(Ni_max,3)
      common /weights/ wght, wquad
 
c variable for anomlous resistivity
 
      real etar(nx,ny,nz,3)
      real nuei
      common /resist/ etar,nuei
 
c variable for particle scale
 
      real beta, dNi
      common /scale/ beta, dNi
 
 
c mass array for multi-ion species
      real mrat,m_arr
      common /mass/ mrat(Ni_max),m_arr(Ni_max)
 
 
 
      real uf(nx,ny,nz,3),
     x     uf1(nx,ny,nz,3),
     x     b1(nx,ny,nz,3),
     x     b12(nx,ny,nz,3),
     x     b1p2(nx,ny,nz,3),
     x     bt(nx,ny,nz,3),
     x     b0(ny),
     x     bdp(nx,ny,nz,3)
 
      EXTERNAL VP_F_UPDATE_TLEV1
      INTEGER J3, J4, J5, J6, J7, J8
      CALL VP_PCALL(VP_F_UPDATE_TLEV1,2,6,B1P2,BDP,BT,B0,B1,B12)
 
      return
      end
      SUBROUTINE VP_F_UPDATE_TLEV1(VP_MYTHREAD, VP_NMTHREADS, B1P2, BDP
     1   , BT, B0, B1, B12)
      REAL B1P2(79,19,37,3), BDP(79,19,37,3), BT(79,19,37,3), B0(19), B1
     1   (79,19,37,3), B12(79,19,37,3)
      INTEGER J, K, M, I, J3, J4, J5, J6, J7, J8, VP_MYTHREAD, 
     1   VP_NMTHREADS, J1, J2
      J3 = (19 + VP_NMTHREADS - 1)/VP_NMTHREADS
      J4 = VP_MYTHREAD*J3 + 1
      J5 = MIN(19 - J4 + 1,J3)
      DO J = J4, J5 + J4 - 1
         DO K = 1, 37
            DO I = 1, 79
               BT(I,J,K,1) = B1P2(I,J,K,1) + BDP(I,J,K,1)
               BT(I,J,K,2) = B1P2(I,J,K,2) + BDP(I,J,K,2) + B0(J)
               BT(I,J,K,3) = B1P2(I,J,K,3) + BDP(I,J,K,3)
            END DO
         END DO
      END DO
      J6 = (19 + VP_NMTHREADS - 1)/VP_NMTHREADS
      J7 = VP_MYTHREAD*J6 + 1
      J8 = MIN(19 - J7 + 1,J6)
      DO J = J7, J8 + J7 - 1
         DO K = 1, 37
            DO I = 1, 79
c                  uf2(i,j,k,m) = uf(i,j,k,m)
               B12(I,J,K,1) = B1(I,J,K,1)
               B12(I,J,K,2) = B1(I,J,K,2)
               B12(I,J,K,3) = B1(I,J,K,3)
               B1(I,J,K,1) = B1P2(I,J,K,1)
               B1(I,J,K,2) = B1P2(I,J,K,2)
               B1(I,J,K,3) = B1P2(I,J,K,3)
            END DO
         END DO
      END DO
      RETURN 
      END 
c----------------------------------------------------------------------
 
c----------------------------------------------------------------------
      SUBROUTINE crossf(aa,bbmf,cc)
c The cross product is formed at the main cell center.  a is assumed
c be main cell contravarient (cell face) and b is assumed to be
c main cell covarient (cell edge).  The result is main cell
c contravarient (cell face).
 
c Can only center vectors on main cell for loops 3 to n...and can
c only extrapolate back for loops 2 to n-1.  Must handle other cases
c separately.
 
c The magnetic field does not require extrapolation to cell centers
c on boundaries since dB/dx = 0 is the boundary condition.  That is
c just copy the interior values to the boundary.
c----------------------------------------------------------------------
      PARAMETER (nx = 79, ny = 19, nz = 37)
 
c grid parameters
      PARAMETER (dx = 10000.0, dy = 10000.0)   !units in km
      PARAMETER (delz = 10000.0)          !dz at release coordinates
                                       !this is not the final dz
                                       !d_lamba determines dz
c time stepping parameters
      PARAMETER (tstart = 0.2)   !starting time
      PARAMETER (dt = 1.)     !main time step
      PARAMETER (nt = 5.0)        !number of time steps
c      PARAMETER (dtsub = 0.01) !subcycle time step
      PARAMETER (dtsub_init = 0.1) !subcycle time step
      PARAMETER (ntsub = 10.0)   !number of subcycle time steps
      PARAMETER (nout = 50)      !number of steps to diagnosic output
 
c logical variable for restart
C...Translated by Crescent Bay Software VAST-F    7.0K8 15:17:00   4/ 1/04    
      logical restart
      PARAMETER (restart = .false.)
      PARAMETER (mrestart = 500)      ! use -1 for no save
      PARAMETER (mbegin = 0)      !mstart
 
c inital ambient density at the release point.
      real nf_init,b0_init,nn_coef
      PARAMETER (nf_init = 4e15)   !was 0.01
      PARAMETER (b0_init = 2e-9)    !was 0.2
      PARAMETER (nn_coef = 1e5)
 
c electron ion collision frequency
      real nu_init, eta_init
      PARAMETER (nu_init = 0.001)
      PARAMETER (eta_init = 0.0)
 
c temporary directory for output
      character*9 temp_dir
      PARAMETER (temp_dir = '/tmp/pad/')
 
c total number of neutrals released
      real No
      PARAMETER (No=1.0e28)
 
c neutral cloud expansion characteristics
      PARAMETER(vo = 20.0)
      PARAMETER(vth = 19.3)
      PARAMETER(vsat = 0.0)
      PARAMETER(vsw = 340.0)
 
c total kinetic energy of released neutrals
      PARAMETER(Etot = 2.696e7)       !units of kg m^2 / s^2
 
c density scaling parameter, alpha, and ion particle array dims
 
      real alpha
 
c      PARAMETER (beta =  8.5227271e-18)
      PARAMETER (alpha = 1.9263418e-20) !mH
 
c max number of ion particles to be produced.  This parameter
c is used to dimension the particle arrays.
      integer*4 Ni_max
      PARAMETER (Ni_max = 10000)
 
c earth's magnetic dipole field geometry
c earth radius and equatorial B alt (defining L shell)
      PARAMETER (r_earth=6.38e3)                !earth radius
      PARAMETER (B_max_alt=2.0e3+r_earth)       !ro defining L shell
 
c misc constants
      real mu0,epsilon,pi,rtod,q,mO,mBa,O_to_Ba,km_to_m,kboltz,melec
      real mproton,tempf0
      PARAMETER (pi = 3.14159)
      PARAMETER (rtod = 180.0/pi)  !radians to degreesc
      PARAMETER (mu0 = pi*4.0e-7)  !magnetic permeability of free space
      PARAMETER (epsilon = 8.85e-12) !dielectric constant
      PARAMETER (q = 1.6e-19)        !electron charge (Coulombs)
      PARAMETER (melec = 9.1e-31)  !mass of electron (kg)
      PARAMETER (mproton = 1.67e-27)
      PARAMETER (mO = 1.67e-27)    !mass of H (kg)
c      PARAMETER (mO = 2.3e-25)    !mass of Ba (kg)
      PARAMETER (mBa = 17.0*mO)    !mass of Ba (kg)
      PARAMETER (O_to_Ba = mO/mBa) !convert E and B for particle move
      PARAMETER (km_to_m = 1e3)    !km to meters conversion
      PARAMETER (kboltz = 1.38e-29)   !kg km^2 / s^2 / K
      PARAMETER (tempf0 = 50*11600.)     !K
 
      real Qo, vrad, N_o, Rp, tau_photo
      integer S_radius
      PARAMETER (Qo = 1.0e28)  !neutral source rate
      PARAMETER (vrad = 1.0)   !escape velocity
c      PARAMETER (N_o = 5e34)   !Steady state neutral particle constant
      PARAMETER (N_o = 1e33)   !Steady state neutral particle constant
      PARAMETER (Rp = 20000.0)  !Pluto radius
      PARAMETER (tau_photo = 1.2e6)
c      PARAMETER (dNi = 2500)
      PARAMETER (S_radius = 3) !units of dx
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
c raw grid coordinate data
 
      real qx(nx), qy(ny), qz(nz)
                            !cartisian coordinates of grid points
      real lambda(nz)       !latitude of given grid point
      integer ri, rj, rk    !release coordinates
      real dz_grid(nz)      !dz for grid, variable in z direction
      real dz_cell(nz)      !dz for cell, variable in z direction
 
      common /coords/ qx,qy,qz,lambda,ri,rj,rk,dz_grid,dz_cell
 
 
c Total number of ions produced at a given time during the simulation
      integer*4 Ni_tot, Nn_tot
 
c Location (indices) of particles with the grid
      integer ijkp(Ni_max,3)
      logical ionized(Ni_max)
      common /ion_info/ Ni_tot, Nn_tot, ijkp, ionized
 
      real seed
      common /rndseed/ seed
 
c Indices for boundary cells
      integer ip,im,jp,jm,kp,km
 
c Total input energy (ions)
      real input_E,input_EeP,prev_Etot, bndry_Eflux
      real input_chex, input_bill
      common /energy/ input_E, input_EeP, prev_Etot, bndry_Eflux,
     x                input_chex, input_bill
 
c Momentum conservation
      real surf_tot(3),graduu_tot(3),ugradu_tot(3)
      common /momentum/ surf_tot,graduu_tot,ugradu_tot
 
c Dummy variables for doing vector operations and particle flux
c calculation....etc.
 
      real a(nx,ny,nz,3),c(nx,ny,nz,3),ct(nx,ny,nz,3)
      common /dummy/ a,c,ct
 
c Variable time step for fluid velocity update used in both
c particle update and fluid update
 
      real dtf,ntf,dtsub
      common /timestep/ ntf,dtsub
 
c Weight variables for trilinear interpolation
 
      real wght(Ni_max,8), wquad(Ni_max,3)
      common /weights/ wght, wquad
 
c variable for anomlous resistivity
 
      real etar(nx,ny,nz,3)
      real nuei
      common /resist/ etar,nuei
 
c variable for particle scale
 
      real beta, dNi
      common /scale/ beta, dNi
 
 
c mass array for multi-ion species
      real mrat,m_arr
      common /mass/ mrat(Ni_max),m_arr(Ni_max)
 
 
 
      real aa(nx,ny,nz,3)        !main cell contravarient vector
      real bbmf(nx,ny,nz,3)      !main cell contravarient vector
      real cc(nx,ny,nz,3)        !cross product result, main cell
                                 !contravarient (cell face)
 
      real ax,ay,az,bx,by,bz    !dummy vars
      real temp                 !used to vectorize loop
      real zfrc(nz)             !0.5*dz_grid(k)/dz_cell(k)
c      real ct(nx,ny,nz,3)       !temp main cell center cross product
      real aac(3),bbc(3)
 
 
c extrapolate(/interpolate) to main cell center and do cross product
 
 
      external vp_crossf1
      external vp_crossf2
      integer j3, j4, j5, j6, j7, j8
      real r1, r2, r3, r4, r5, r6
      call periodic(aa)
      call periodic(bbmf)
 
      do 5 k=1,nz
         zfrc(k) = 0.5*dz_grid(k)/dz_cell(k)
 5       continue
 
      call vp_pcall (vp_crossf1, 2, 3, aa, bbmf, zfrc)
 
       call periodic(ct)
 
c extrapolate back to main cell contravarient positions.
c ...just average across cells since cell edges are centered
c about the grid points.
 
      call vp_pcall (vp_crossf2, 2, 1, cc)
 
      call periodic(cc)
 
      return
      end
      subroutine vp_crossf1(vp_mythread, vp_nmthreads, aa, bbmf, zfrc)
      common /dummy/a, c, ct
      REAL aa(79,19,37,3), bbmf(79,19,37,3), zfrc(37), ct(79,19,37,3), a
     1   (79,19,37,3), c(79,19,37,3)
      REAL ax, bx, ay, by, az, bz
      INTEGER k,km,j,jm,i,im,j3,j4,j5,vp_mythread,vp_nmthreads,j1,j2
 
      j3 = (35 + vp_nmthreads - 1)/vp_nmthreads
      j4 = vp_mythread*j3 + 1
      j5 = min(35 - j4 + 1,j3)
      do k = j4, j5 + j4 - 1
         do j = 1, 17
            do i = 1, 77
 
               ax = 0.5*(aa(1+i,j+1,k+1,1)+aa(i,j+1,k+1,1))
               bx = 0.5*(bbmf(1+i,j+1,k+1,1)+bbmf(i,j+1,k+1,1))
 
               ay = 0.5*(aa(1+i,j+1,k+1,2)+aa(1+i,j,k+1,2))
               by = 0.5*(bbmf(1+i,j+1,k+1,2)+bbmf(1+i,j,k+1,2))
 
               az = zfrc(k+1)*(aa(1+i,j+1,k+1,3)-aa(1+i,j+1,k,3)) + aa(1
     1            +i,j+1,k,3)
               bz = zfrc(k+1)*(bbmf(1+i,j+1,k+1,3)-bbmf(1+i,j+1,k,3)) + 
     1            bbmf(1+i,j+1,k,3)
 
               ct(1+i,j+1,k+1,1) = ay*bz - az*by
               ct(1+i,j+1,k+1,2) = az*bx - ax*bz
               ct(1+i,j+1,k+1,3) = ax*by - ay*bx
            end do
 
         end do
      end do
      return 
      end 
      subroutine vp_crossf2(vp_mythread, vp_nmthreads, cc)
      common /dummy/a, c, ct
      REAL ct(79,19,37,3), cc(79,19,37,3), a(79,19,37,3), c(79,19,37,3)
      REAL r1, r2, r3, r4, r5, r6
      INTEGER k,kp,j,jp,i,ip,j6,j7,j8,vp_mythread,vp_nmthreads,j1,j2
 
c extrapolate back to main cell contravarient positions.
c ...just average across cells since cell edges are centered
c about the grid points.
 
      j6 = (35 + vp_nmthreads - 1)/vp_nmthreads
      j7 = vp_mythread*j6 + 1
      j8 = min(35 - j7 + 1,j6)
      do k = j7, j8 + j7 - 1
         j = 1
         do i = 1, 77
 
            cc(1+i,2,k+1,1) = 0.5*(ct(1+i,2,k+1,1)+ct(2+i,2,k+1,1))
            cc(1+i,2,k+1,2) = 0.5*(ct(1+i,2,k+1,2)+ct(1+i,3,k+1,2))
            cc(1+i,2,k+1,3) = 0.5*(ct(1+i,2,k+1,3)+ct(1+i,2,k+2,3))
         end do
 
         do j = 1, 8
            do i = 1, 77
               r1 = 0.5*(ct(1+i,j*2+1,k+1,1)+ct(2+i,j*2+1,k+1,1))
               r2 = 0.5*(ct(1+i,(j+1)*2,k+1,1)+ct(2+i,(j+1)*2,k+1,1))
 
               cc(1+i,j*2+1,k+1,1) = r1
               cc(1+i,(j+1)*2,k+1,1) = r2
               r3 = 0.5*(ct(1+i,j*2+1,k+1,2)+ct(1+i,(j+1)*2,k+1,2))
               r4 = 0.5*(ct(1+i,(j+1)*2,k+1,2)+ct(1+i,j*2+3,k+1,2))
               cc(1+i,j*2+1,k+1,2) = r3
               cc(1+i,(j+1)*2,k+1,2) = r4
               r5 = 0.5*(ct(1+i,j*2+1,k+1,3)+ct(1+i,j*2+1,k+2,3))
               r6 = 0.5*(ct(1+i,(j+1)*2,k+1,3)+ct(1+i,(j+1)*2,k+2,3))
               cc(1+i,j*2+1,k+1,3) = r5
               cc(1+i,(j+1)*2,k+1,3) = r6
            end do
 
         end do
      end do
      return 
      end 
c----------------------------------------------------------------------
 
 
c----------------------------------------------------------------------
      SUBROUTINE cov_to_contra(bt,btmf)
c Converts total magnetic field from main cell covarient positions
c to main cell contravarient positions.  This is then used in the
c fluid velocity update routines.  This routine assumes that cell
c edges and cell centers are "topologically centered".  So the grid
c points do not reside at the cell centers...rather they are offset
c a little so that the cell edges are equidistant from the k and k-1
c grid points.  In extrapolating the coventient vectors to the
c contravarient vector positions, this assymetry is accounted for
c using a linear interpolation of the k and k-1 values to the grid
c point location.
c----------------------------------------------------------------------
      PARAMETER (nx = 79, ny = 19, nz = 37)
 
c grid parameters
      PARAMETER (dx = 10000.0, dy = 10000.0)   !units in km
      PARAMETER (delz = 10000.0)          !dz at release coordinates
                                       !this is not the final dz
                                       !d_lamba determines dz
c time stepping parameters
      PARAMETER (tstart = 0.2)   !starting time
      PARAMETER (dt = 1.)     !main time step
      PARAMETER (nt = 5.0)        !number of time steps
c      PARAMETER (dtsub = 0.01) !subcycle time step
      PARAMETER (dtsub_init = 0.1) !subcycle time step
      PARAMETER (ntsub = 10.0)   !number of subcycle time steps
      PARAMETER (nout = 50)      !number of steps to diagnosic output
 
c logical variable for restart
C...Translated by Crescent Bay Software VAST-F    7.0K8 15:17:00   4/ 1/04    
      logical restart
      PARAMETER (restart = .false.)
      PARAMETER (mrestart = 500)      ! use -1 for no save
      PARAMETER (mbegin = 0)      !mstart
 
c inital ambient density at the release point.
      real nf_init,b0_init,nn_coef
      PARAMETER (nf_init = 4e15)   !was 0.01
      PARAMETER (b0_init = 2e-9)    !was 0.2
      PARAMETER (nn_coef = 1e5)
 
c electron ion collision frequency
      real nu_init, eta_init
      PARAMETER (nu_init = 0.001)
      PARAMETER (eta_init = 0.0)
 
c temporary directory for output
      character*9 temp_dir
      PARAMETER (temp_dir = '/tmp/pad/')
 
c total number of neutrals released
      real No
      PARAMETER (No=1.0e28)
 
c neutral cloud expansion characteristics
      PARAMETER(vo = 20.0)
      PARAMETER(vth = 19.3)
      PARAMETER(vsat = 0.0)
      PARAMETER(vsw = 340.0)
 
c total kinetic energy of released neutrals
      PARAMETER(Etot = 2.696e7)       !units of kg m^2 / s^2
 
c density scaling parameter, alpha, and ion particle array dims
 
      real alpha
 
c      PARAMETER (beta =  8.5227271e-18)
      PARAMETER (alpha = 1.9263418e-20) !mH
 
c max number of ion particles to be produced.  This parameter
c is used to dimension the particle arrays.
      integer*4 Ni_max
      PARAMETER (Ni_max = 10000)
 
c earth's magnetic dipole field geometry
c earth radius and equatorial B alt (defining L shell)
      PARAMETER (r_earth=6.38e3)                !earth radius
      PARAMETER (B_max_alt=2.0e3+r_earth)       !ro defining L shell
 
c misc constants
      real mu0,epsilon,pi,rtod,q,mO,mBa,O_to_Ba,km_to_m,kboltz,melec
      real mproton,tempf0
      PARAMETER (pi = 3.14159)
      PARAMETER (rtod = 180.0/pi)  !radians to degreesc
      PARAMETER (mu0 = pi*4.0e-7)  !magnetic permeability of free space
      PARAMETER (epsilon = 8.85e-12) !dielectric constant
      PARAMETER (q = 1.6e-19)        !electron charge (Coulombs)
      PARAMETER (melec = 9.1e-31)  !mass of electron (kg)
      PARAMETER (mproton = 1.67e-27)
      PARAMETER (mO = 1.67e-27)    !mass of H (kg)
c      PARAMETER (mO = 2.3e-25)    !mass of Ba (kg)
      PARAMETER (mBa = 17.0*mO)    !mass of Ba (kg)
      PARAMETER (O_to_Ba = mO/mBa) !convert E and B for particle move
      PARAMETER (km_to_m = 1e3)    !km to meters conversion
      PARAMETER (kboltz = 1.38e-29)   !kg km^2 / s^2 / K
      PARAMETER (tempf0 = 50*11600.)     !K
 
      real Qo, vrad, N_o, Rp, tau_photo
      integer S_radius
      PARAMETER (Qo = 1.0e28)  !neutral source rate
      PARAMETER (vrad = 1.0)   !escape velocity
c      PARAMETER (N_o = 5e34)   !Steady state neutral particle constant
      PARAMETER (N_o = 1e33)   !Steady state neutral particle constant
      PARAMETER (Rp = 20000.0)  !Pluto radius
      PARAMETER (tau_photo = 1.2e6)
c      PARAMETER (dNi = 2500)
      PARAMETER (S_radius = 3) !units of dx
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
c raw grid coordinate data
 
      real qx(nx), qy(ny), qz(nz)
                            !cartisian coordinates of grid points
      real lambda(nz)       !latitude of given grid point
      integer ri, rj, rk    !release coordinates
      real dz_grid(nz)      !dz for grid, variable in z direction
      real dz_cell(nz)      !dz for cell, variable in z direction
 
      common /coords/ qx,qy,qz,lambda,ri,rj,rk,dz_grid,dz_cell
 
 
c Total number of ions produced at a given time during the simulation
      integer*4 Ni_tot, Nn_tot
 
c Location (indices) of particles with the grid
      integer ijkp(Ni_max,3)
      logical ionized(Ni_max)
      common /ion_info/ Ni_tot, Nn_tot, ijkp, ionized
 
      real seed
      common /rndseed/ seed
 
c Indices for boundary cells
      integer ip,im,jp,jm,kp,km
 
c Total input energy (ions)
      real input_E,input_EeP,prev_Etot, bndry_Eflux
      real input_chex, input_bill
      common /energy/ input_E, input_EeP, prev_Etot, bndry_Eflux,
     x                input_chex, input_bill
 
c Momentum conservation
      real surf_tot(3),graduu_tot(3),ugradu_tot(3)
      common /momentum/ surf_tot,graduu_tot,ugradu_tot
 
c Dummy variables for doing vector operations and particle flux
c calculation....etc.
 
      real a(nx,ny,nz,3),c(nx,ny,nz,3),ct(nx,ny,nz,3)
      common /dummy/ a,c,ct
 
c Variable time step for fluid velocity update used in both
c particle update and fluid update
 
      real dtf,ntf,dtsub
      common /timestep/ ntf,dtsub
 
c Weight variables for trilinear interpolation
 
      real wght(Ni_max,8), wquad(Ni_max,3)
      common /weights/ wght, wquad
 
c variable for anomlous resistivity
 
      real etar(nx,ny,nz,3)
      real nuei
      common /resist/ etar,nuei
 
c variable for particle scale
 
      real beta, dNi
      common /scale/ beta, dNi
 
 
c mass array for multi-ion species
      real mrat,m_arr
      common /mass/ mrat(Ni_max),m_arr(Ni_max)
 
 
 
      real bt(nx,ny,nz,3),   !main cell covarient
     x     btmf(nx,ny,nz,3)  !main cell contravarient
 
      real bx1, bx2, by1, by2, bz1, bz2  !main cell center fields
      real zrat           !ratio for doing linear interpolation
                          !to grid point position.
      real zplus, zminus  !position of main cell edges up and down
      real b_j, b_jm, b_i, b_im !intermediate step in average process
 
      external vp_cov_to_contra1
      integer j3, j4, j5
      real r1, r2, r3, r4
      call vp_pcall (vp_cov_to_contra1, 2, 2, bt, btmf)
 
c      call boundaries(btmf)
      call periodic(btmf)
 
      return
      end
      subroutine vp_cov_to_contra1(vp_mythread, vp_nmthreads, bt, btmf)
      common /coords/qx, qy, qz, lambda, ri, rj, rk, dz_grid, dz_cell
      INTEGER ri, rj, rk
      REAL qz(37), bt(79,19,37,3), btmf(79,19,37,3), qx(79), qy(19), 
     1   lambda(37), dz_grid(37), dz_cell(37)
      REAL zplus, zminus, zrat, b_j, b_jm, bx1, bx2, b_i, b_im, by1, by2
     1   , bz1, bz2, r1, r2, r3, r4
      INTEGER k, kp, km, j, jp, jm, i, ip, im, j3, j4, j5, vp_mythread, 
     1   vp_nmthreads, j1, j2
 
      j3 = (35 + vp_nmthreads - 1)/vp_nmthreads
      j4 = vp_mythread*j3 + 1
      j5 = min(35 - j4 + 1,j3)
      do k = j4, j5 + j4 - 1
         do j = 1, 17
            r1 = 1./2.0
            r2 = 1./2.0
            r3 = 1./2.0
            r4 = 1./2.0
 
c The x component of B resides at the k and k-1 edges, so this
c requires the non-uniform grid interpolation
 
            zplus = (qz(k+2)+qz(k+1))/2.0
            zminus = (qz(k+1)+qz(k))/2.0
            zrat = (qz(k+1)-zminus)/(zplus - zminus)
            do i = 1, 77
               b_jm=bt(1+i,j,k,1)+zrat*(bt(1+i,j,k+1,1)-bt(1+i,j,k,1))
               bx1 = (bt(1+i,j+1,k,1)+zrat*(bt(1+i,j+1,k+1,1)-bt(1+i,j+1
     1            ,k,1))+b_jm)*r1
               b_jm=bt(2+i,j,k,1)+zrat*(bt(2+i,j,k+1,1)-bt(2+i,j,k,1))
               bx2 = (bt(2+i,j+1,k,1)+zrat*(bt(2+i,j+1,k+1,1)-bt(2+i,j+1
     1            ,k,1))+b_jm)*r2
               b_im=bt(i,j+1,k,2)+zrat*(bt(i,j+1,k+1,2)-bt(i,j+1,k,2))
               by1 = (bt(1+i,j+1,k,2)+zrat*(bt(1+i,j+1,k+1,2)-bt(1+i,j+1
     1            ,k,2))+b_im)*r3
               b_im=bt(i,j+2,k,2)+zrat*(bt(i,j+2,k+1,2)-bt(i,j+2,k,2))
               by2 = (bt(1+i,j+2,k,2)+zrat*(bt(1+i,j+2,k+1,2)-bt(1+i,j+2
     1            ,k,2))+b_im)*r4
 
 
               bz1 = 0.25*(bt(1+i,j+1,k+1,3)+bt(1+i,j,k+1,3)+bt(i,j,k+1,
     1            3)+bt(i,j+1,k+1,3))
               bz2 = 0.25*(bt(1+i,j+1,k+2,3)+bt(1+i,j,k+2,3)+bt(i,j,k+2,
     1            3)+bt(i,j+1,k+2,3))
 
               btmf(1+i,j+1,k+1,1) = 0.5*(bx1 + bx2)
               btmf(1+i,j+1,k+1,2) = 0.5*(by1 + by2)
               btmf(1+i,j+1,k+1,3) = 0.5*(bz1 + bz2)
            end do
 
         end do
      end do
      return 
      end 
c----------------------------------------------------------------------
 
 
c----------------------------------------------------------------------
      SUBROUTINE curlB(b1,nf,np,aj)
c Calculates curl B / n*alpha.  The resulting "current" is called aj
c which is used in several other places in the code.  This curl is
c performed on the main cell where B is covarient.  The resulting
c current is main cell contravarient.  Note that dz_cell is used for
c the cell dimensions since dz_grid is not equal to dz_cell on non-
c uniform grid.
c----------------------------------------------------------------------
      PARAMETER (nx = 79, ny = 19, nz = 37)
 
c grid parameters
      PARAMETER (dx = 10000.0, dy = 10000.0)   !units in km
      PARAMETER (delz = 10000.0)          !dz at release coordinates
                                       !this is not the final dz
                                       !d_lamba determines dz
c time stepping parameters
      PARAMETER (tstart = 0.2)   !starting time
      PARAMETER (dt = 1.)     !main time step
      PARAMETER (nt = 5.0)        !number of time steps
c      PARAMETER (dtsub = 0.01) !subcycle time step
      PARAMETER (dtsub_init = 0.1) !subcycle time step
      PARAMETER (ntsub = 10.0)   !number of subcycle time steps
      PARAMETER (nout = 50)      !number of steps to diagnosic output
 
c logical variable for restart
C...Translated by Crescent Bay Software VAST-F    7.0K8 15:17:00   4/ 1/04    
      logical restart
      PARAMETER (restart = .false.)
      PARAMETER (mrestart = 500)      ! use -1 for no save
      PARAMETER (mbegin = 0)      !mstart
 
c inital ambient density at the release point.
      real nf_init,b0_init,nn_coef
      PARAMETER (nf_init = 4e15)   !was 0.01
      PARAMETER (b0_init = 2e-9)    !was 0.2
      PARAMETER (nn_coef = 1e5)
 
c electron ion collision frequency
      real nu_init, eta_init
      PARAMETER (nu_init = 0.001)
      PARAMETER (eta_init = 0.0)
 
c temporary directory for output
      character*9 temp_dir
      PARAMETER (temp_dir = '/tmp/pad/')
 
c total number of neutrals released
      real No
      PARAMETER (No=1.0e28)
 
c neutral cloud expansion characteristics
      PARAMETER(vo = 20.0)
      PARAMETER(vth = 19.3)
      PARAMETER(vsat = 0.0)
      PARAMETER(vsw = 340.0)
 
c total kinetic energy of released neutrals
      PARAMETER(Etot = 2.696e7)       !units of kg m^2 / s^2
 
c density scaling parameter, alpha, and ion particle array dims
 
      real alpha
 
c      PARAMETER (beta =  8.5227271e-18)
      PARAMETER (alpha = 1.9263418e-20) !mH
 
c max number of ion particles to be produced.  This parameter
c is used to dimension the particle arrays.
      integer*4 Ni_max
      PARAMETER (Ni_max = 10000)
 
c earth's magnetic dipole field geometry
c earth radius and equatorial B alt (defining L shell)
      PARAMETER (r_earth=6.38e3)                !earth radius
      PARAMETER (B_max_alt=2.0e3+r_earth)       !ro defining L shell
 
c misc constants
      real mu0,epsilon,pi,rtod,q,mO,mBa,O_to_Ba,km_to_m,kboltz,melec
      real mproton,tempf0
      PARAMETER (pi = 3.14159)
      PARAMETER (rtod = 180.0/pi)  !radians to degreesc
      PARAMETER (mu0 = pi*4.0e-7)  !magnetic permeability of free space
      PARAMETER (epsilon = 8.85e-12) !dielectric constant
      PARAMETER (q = 1.6e-19)        !electron charge (Coulombs)
      PARAMETER (melec = 9.1e-31)  !mass of electron (kg)
      PARAMETER (mproton = 1.67e-27)
      PARAMETER (mO = 1.67e-27)    !mass of H (kg)
c      PARAMETER (mO = 2.3e-25)    !mass of Ba (kg)
      PARAMETER (mBa = 17.0*mO)    !mass of Ba (kg)
      PARAMETER (O_to_Ba = mO/mBa) !convert E and B for particle move
      PARAMETER (km_to_m = 1e3)    !km to meters conversion
      PARAMETER (kboltz = 1.38e-29)   !kg km^2 / s^2 / K
      PARAMETER (tempf0 = 50*11600.)     !K
 
      real Qo, vrad, N_o, Rp, tau_photo
      integer S_radius
      PARAMETER (Qo = 1.0e28)  !neutral source rate
      PARAMETER (vrad = 1.0)   !escape velocity
c      PARAMETER (N_o = 5e34)   !Steady state neutral particle constant
      PARAMETER (N_o = 1e33)   !Steady state neutral particle constant
      PARAMETER (Rp = 20000.0)  !Pluto radius
      PARAMETER (tau_photo = 1.2e6)
c      PARAMETER (dNi = 2500)
      PARAMETER (S_radius = 3) !units of dx
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
c raw grid coordinate data
 
      real qx(nx), qy(ny), qz(nz)
                            !cartisian coordinates of grid points
      real lambda(nz)       !latitude of given grid point
      integer ri, rj, rk    !release coordinates
      real dz_grid(nz)      !dz for grid, variable in z direction
      real dz_cell(nz)      !dz for cell, variable in z direction
 
      common /coords/ qx,qy,qz,lambda,ri,rj,rk,dz_grid,dz_cell
 
 
c Total number of ions produced at a given time during the simulation
      integer*4 Ni_tot, Nn_tot
 
c Location (indices) of particles with the grid
      integer ijkp(Ni_max,3)
      logical ionized(Ni_max)
      common /ion_info/ Ni_tot, Nn_tot, ijkp, ionized
 
      real seed
      common /rndseed/ seed
 
c Indices for boundary cells
      integer ip,im,jp,jm,kp,km
 
c Total input energy (ions)
      real input_E,input_EeP,prev_Etot, bndry_Eflux
      real input_chex, input_bill
      common /energy/ input_E, input_EeP, prev_Etot, bndry_Eflux,
     x                input_chex, input_bill
 
c Momentum conservation
      real surf_tot(3),graduu_tot(3),ugradu_tot(3)
      common /momentum/ surf_tot,graduu_tot,ugradu_tot
 
c Dummy variables for doing vector operations and particle flux
c calculation....etc.
 
      real a(nx,ny,nz,3),c(nx,ny,nz,3),ct(nx,ny,nz,3)
      common /dummy/ a,c,ct
 
c Variable time step for fluid velocity update used in both
c particle update and fluid update
 
      real dtf,ntf,dtsub
      common /timestep/ ntf,dtsub
 
c Weight variables for trilinear interpolation
 
      real wght(Ni_max,8), wquad(Ni_max,3)
      common /weights/ wght, wquad
 
c variable for anomlous resistivity
 
      real etar(nx,ny,nz,3)
      real nuei
      common /resist/ etar,nuei
 
c variable for particle scale
 
      real beta, dNi
      common /scale/ beta, dNi
 
 
c mass array for multi-ion species
      real mrat,m_arr
      common /mass/ mrat(Ni_max),m_arr(Ni_max)
 
 
 
      real b1(nx,ny,nz,3),
     x     nf(nx,ny,nz),
     x     np(nx,ny,nz),
     x     aj(nx,ny,nz,3)
 
      real curl_B(3)      !dummy for holding curl vector
      real ntot(3)        !total density, np + nf
 
c      call periodic_scalar(np)
      REAL R1, R2, R3, R4, R5, R6, R7, R8, R9
      call periodic_scalar(nf)
      call periodic(b1)
      call fix_normal_b(b1)
 
      DO 10 K = 2, 36
         DO 10 J = 2, 18
            DO 10 I = 2, 78
 
               IP = I + 1
               JP = J + 1
               KP = K + 1
 
c               if (ip .gt. nx) then ip = nx
c               if (jp .gt. ny) then jp = ny
c               if (kp .gt. nz) then kp = nz
 
               NTOT(1) = 0.5*(NF(I,J,K)+NF(IP,J,K)) + 0.5*(NP(I,J,K)+NP(
     1            IP,J,K))
               NTOT(2) = 0.5*(NF(I,J,K)+NF(I,JP,K)) + 0.5*(NP(I,J,K)+NP(
     1            I,JP,K))
               NTOT(3) = 0.5*(NF(I,J,K)+NF(I,J,KP)) + 0.5*(NP(I,J,K)+NP(
     1            I,J,KP))
 
               CURL_B(1) = B1(I,J,K,3)/10000.0 - B1(I,J-1,K,3)/10000.0
     1             + B1(I,J,K-1,2)/DZ_CELL(K) - B1(I,J,K,2)/DZ_CELL(K)
               CURL_B(2) = B1(I,J,K,1)/DZ_CELL(K) - B1(I,J,K-1,1)/
     1            DZ_CELL(K) - B1(I,J,K,3)/10000.0 + B1(I-1,J,K,3)/
     2            10000.0
               CURL_B(3) = B1(I,J,K,2)/10000.0 - B1(I-1,J,K,2)/10000.0
     1             + B1(I,J-1,K,1)/10000.0 - B1(I,J,K,1)/10000.0
 
               R7 = CURL_B(1)/(NTOT(1)*1.9263418E-20)
               R8 = CURL_B(2)/(NTOT(2)*1.9263418E-20)
               R9 = CURL_B(3)/(NTOT(3)*1.9263418E-20)
               AJ(I,J,K,1) = R7
               AJ(I,J,K,2) = R8
               AJ(I,J,K,3) = R9
   10 CONTINUE
 
      call periodic(aj)
 
      return
      end
c----------------------------------------------------------------------
 
 
c----------------------------------------------------------------------
      SUBROUTINE curlE(E,curl_E)
c E is dual cell covarient, and curl_E will be returned as main
c cell covarient...as all magnetic fields are.  All i,j,k exclude
c boundaries.  Boundaries are taken care of in main fluid code.
c----------------------------------------------------------------------
      PARAMETER (nx = 79, ny = 19, nz = 37)
 
c grid parameters
      PARAMETER (dx = 10000.0, dy = 10000.0)   !units in km
      PARAMETER (delz = 10000.0)          !dz at release coordinates
                                       !this is not the final dz
                                       !d_lamba determines dz
c time stepping parameters
      PARAMETER (tstart = 0.2)   !starting time
      PARAMETER (dt = 1.)     !main time step
      PARAMETER (nt = 5.0)        !number of time steps
c      PARAMETER (dtsub = 0.01) !subcycle time step
      PARAMETER (dtsub_init = 0.1) !subcycle time step
      PARAMETER (ntsub = 10.0)   !number of subcycle time steps
      PARAMETER (nout = 50)      !number of steps to diagnosic output
 
c logical variable for restart
C...Translated by Crescent Bay Software VAST-F    7.0K8 15:17:00   4/ 1/04    
      logical restart
      PARAMETER (restart = .false.)
      PARAMETER (mrestart = 500)      ! use -1 for no save
      PARAMETER (mbegin = 0)      !mstart
 
c inital ambient density at the release point.
      real nf_init,b0_init,nn_coef
      PARAMETER (nf_init = 4e15)   !was 0.01
      PARAMETER (b0_init = 2e-9)    !was 0.2
      PARAMETER (nn_coef = 1e5)
 
c electron ion collision frequency
      real nu_init, eta_init
      PARAMETER (nu_init = 0.001)
      PARAMETER (eta_init = 0.0)
 
c temporary directory for output
      character*9 temp_dir
      PARAMETER (temp_dir = '/tmp/pad/')
 
c total number of neutrals released
      real No
      PARAMETER (No=1.0e28)
 
c neutral cloud expansion characteristics
      PARAMETER(vo = 20.0)
      PARAMETER(vth = 19.3)
      PARAMETER(vsat = 0.0)
      PARAMETER(vsw = 340.0)
 
c total kinetic energy of released neutrals
      PARAMETER(Etot = 2.696e7)       !units of kg m^2 / s^2
 
c density scaling parameter, alpha, and ion particle array dims
 
      real alpha
 
c      PARAMETER (beta =  8.5227271e-18)
      PARAMETER (alpha = 1.9263418e-20) !mH
 
c max number of ion particles to be produced.  This parameter
c is used to dimension the particle arrays.
      integer*4 Ni_max
      PARAMETER (Ni_max = 10000)
 
c earth's magnetic dipole field geometry
c earth radius and equatorial B alt (defining L shell)
      PARAMETER (r_earth=6.38e3)                !earth radius
      PARAMETER (B_max_alt=2.0e3+r_earth)       !ro defining L shell
 
c misc constants
      real mu0,epsilon,pi,rtod,q,mO,mBa,O_to_Ba,km_to_m,kboltz,melec
      real mproton,tempf0
      PARAMETER (pi = 3.14159)
      PARAMETER (rtod = 180.0/pi)  !radians to degreesc
      PARAMETER (mu0 = pi*4.0e-7)  !magnetic permeability of free space
      PARAMETER (epsilon = 8.85e-12) !dielectric constant
      PARAMETER (q = 1.6e-19)        !electron charge (Coulombs)
      PARAMETER (melec = 9.1e-31)  !mass of electron (kg)
      PARAMETER (mproton = 1.67e-27)
      PARAMETER (mO = 1.67e-27)    !mass of H (kg)
c      PARAMETER (mO = 2.3e-25)    !mass of Ba (kg)
      PARAMETER (mBa = 17.0*mO)    !mass of Ba (kg)
      PARAMETER (O_to_Ba = mO/mBa) !convert E and B for particle move
      PARAMETER (km_to_m = 1e3)    !km to meters conversion
      PARAMETER (kboltz = 1.38e-29)   !kg km^2 / s^2 / K
      PARAMETER (tempf0 = 50*11600.)     !K
 
      real Qo, vrad, N_o, Rp, tau_photo
      integer S_radius
      PARAMETER (Qo = 1.0e28)  !neutral source rate
      PARAMETER (vrad = 1.0)   !escape velocity
c      PARAMETER (N_o = 5e34)   !Steady state neutral particle constant
      PARAMETER (N_o = 1e33)   !Steady state neutral particle constant
      PARAMETER (Rp = 20000.0)  !Pluto radius
      PARAMETER (tau_photo = 1.2e6)
c      PARAMETER (dNi = 2500)
      PARAMETER (S_radius = 3) !units of dx
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
c raw grid coordinate data
 
      real qx(nx), qy(ny), qz(nz)
                            !cartisian coordinates of grid points
      real lambda(nz)       !latitude of given grid point
      integer ri, rj, rk    !release coordinates
      real dz_grid(nz)      !dz for grid, variable in z direction
      real dz_cell(nz)      !dz for cell, variable in z direction
 
      common /coords/ qx,qy,qz,lambda,ri,rj,rk,dz_grid,dz_cell
 
 
c Total number of ions produced at a given time during the simulation
      integer*4 Ni_tot, Nn_tot
 
c Location (indices) of particles with the grid
      integer ijkp(Ni_max,3)
      logical ionized(Ni_max)
      common /ion_info/ Ni_tot, Nn_tot, ijkp, ionized
 
      real seed
      common /rndseed/ seed
 
c Indices for boundary cells
      integer ip,im,jp,jm,kp,km
 
c Total input energy (ions)
      real input_E,input_EeP,prev_Etot, bndry_Eflux
      real input_chex, input_bill
      common /energy/ input_E, input_EeP, prev_Etot, bndry_Eflux,
     x                input_chex, input_bill
 
c Momentum conservation
      real surf_tot(3),graduu_tot(3),ugradu_tot(3)
      common /momentum/ surf_tot,graduu_tot,ugradu_tot
 
c Dummy variables for doing vector operations and particle flux
c calculation....etc.
 
      real a(nx,ny,nz,3),c(nx,ny,nz,3),ct(nx,ny,nz,3)
      common /dummy/ a,c,ct
 
c Variable time step for fluid velocity update used in both
c particle update and fluid update
 
      real dtf,ntf,dtsub
      common /timestep/ ntf,dtsub
 
c Weight variables for trilinear interpolation
 
      real wght(Ni_max,8), wquad(Ni_max,3)
      common /weights/ wght, wquad
 
c variable for anomlous resistivity
 
      real etar(nx,ny,nz,3)
      real nuei
      common /resist/ etar,nuei
 
c variable for particle scale
 
      real beta, dNi
      common /scale/ beta, dNi
 
 
c mass array for multi-ion species
      real mrat,m_arr
      common /mass/ mrat(Ni_max),m_arr(Ni_max)
 
 
 
      real E(nx,ny,nz,3)      !E field, main cell contravarient
      real curl_E(nx,ny,nz,3) !curl of E, main cell covarient
      real lx, ly, lz         !lengths of dual cell edges
 
      EXTERNAL VP_CURLE1
      INTEGER J3, J4, J5
      call periodic(E)
 
      CALL VP_PCALL (VP_CURLE1, 2, 2, E, CURL_E)
 
      call periodic(curl_E)
 
      return
      end
      SUBROUTINE VP_CURLE1(VP_MYTHREAD, VP_NMTHREADS, E, CURL_E)
      COMMON /COORDS/QX, QY, QZ, LAMBDA, RI, RJ, RK, DZ_GRID, DZ_CELL
      INTEGER RI, RJ, RK
      REAL QX(79), QY(19), QZ(37), E(79,19,37,3), CURL_E(79,19,37,3), 
     1   LAMBDA(37), DZ_GRID(37), DZ_CELL(37)
      REAL LX, LY, LZ
      INTEGER J, K, I, J3, J4, J5, VP_MYTHREAD, VP_NMTHREADS, J1, J2
      J3 = (17 + VP_NMTHREADS - 1)/VP_NMTHREADS
      J4 = VP_MYTHREAD*J3 + 1
      J5 = MIN(17 - J4 + 1,J3)
      DO J = J4, J5 + J4 - 1
         DO K = 1, 35
            DO I = 1, 77
 
               LX = QX(2+I) - QX(1+I)
               LY = QY(J+2) - QY(J+1)
               LZ = QZ(K+2) - QZ(K+1)
 
               CURL_E(1+I,J+1,K+1,1) = E(1+I,J+2,K+1,3)/LY - E(1+I,J+1,K
     1            +1,3)/LY + E(1+I,J+1,K+1,2)/LZ - E(1+I,J+1,K+2,2)/LZ
               CURL_E(1+I,J+1,K+1,2) = E(1+I,J+1,K+1,3)/LX - E(2+I,J+1,K
     1            +1,3)/LX + E(1+I,J+1,K+2,1)/LZ - E(1+I,J+1,K+1,1)/LZ
               CURL_E(1+I,J+1,K+1,3) = E(1+I,J+1,K+1,1)/LY - E(1+I,J+2,K
     1            +1,1)/LY + E(2+I,J+1,K+1,2)/LX - E(1+I,J+1,K+1,2)/LX
            END DO
         END DO
      END DO
      RETURN 
      END 
c----------------------------------------------------------------------
 
 
c----------------------------------------------------------------------
      SUBROUTINE get_um_dot_BB(u,b,cc)
c uf and btmf are gathered at main cell center and uf.B*B
c calculated.  Result returned to main cell contravarient
c postion.
c----------------------------------------------------------------------
      PARAMETER (nx = 79, ny = 19, nz = 37)
 
c grid parameters
      PARAMETER (dx = 10000.0, dy = 10000.0)   !units in km
      PARAMETER (delz = 10000.0)          !dz at release coordinates
                                       !this is not the final dz
                                       !d_lamba determines dz
c time stepping parameters
      PARAMETER (tstart = 0.2)   !starting time
      PARAMETER (dt = 1.)     !main time step
      PARAMETER (nt = 5.0)        !number of time steps
c      PARAMETER (dtsub = 0.01) !subcycle time step
      PARAMETER (dtsub_init = 0.1) !subcycle time step
      PARAMETER (ntsub = 10.0)   !number of subcycle time steps
      PARAMETER (nout = 50)      !number of steps to diagnosic output
 
c logical variable for restart
C...Translated by Crescent Bay Software VAST-F    7.0K8 15:17:00   4/ 1/04    
      logical restart
      PARAMETER (restart = .false.)
      PARAMETER (mrestart = 500)      ! use -1 for no save
      PARAMETER (mbegin = 0)      !mstart
 
c inital ambient density at the release point.
      real nf_init,b0_init,nn_coef
      PARAMETER (nf_init = 4e15)   !was 0.01
      PARAMETER (b0_init = 2e-9)    !was 0.2
      PARAMETER (nn_coef = 1e5)
 
c electron ion collision frequency
      real nu_init, eta_init
      PARAMETER (nu_init = 0.001)
      PARAMETER (eta_init = 0.0)
 
c temporary directory for output
      character*9 temp_dir
      PARAMETER (temp_dir = '/tmp/pad/')
 
c total number of neutrals released
      real No
      PARAMETER (No=1.0e28)
 
c neutral cloud expansion characteristics
      PARAMETER(vo = 20.0)
      PARAMETER(vth = 19.3)
      PARAMETER(vsat = 0.0)
      PARAMETER(vsw = 340.0)
 
c total kinetic energy of released neutrals
      PARAMETER(Etot = 2.696e7)       !units of kg m^2 / s^2
 
c density scaling parameter, alpha, and ion particle array dims
 
      real alpha
 
c      PARAMETER (beta =  8.5227271e-18)
      PARAMETER (alpha = 1.9263418e-20) !mH
 
c max number of ion particles to be produced.  This parameter
c is used to dimension the particle arrays.
      integer*4 Ni_max
      PARAMETER (Ni_max = 10000)
 
c earth's magnetic dipole field geometry
c earth radius and equatorial B alt (defining L shell)
      PARAMETER (r_earth=6.38e3)                !earth radius
      PARAMETER (B_max_alt=2.0e3+r_earth)       !ro defining L shell
 
c misc constants
      real mu0,epsilon,pi,rtod,q,mO,mBa,O_to_Ba,km_to_m,kboltz,melec
      real mproton,tempf0
      PARAMETER (pi = 3.14159)
      PARAMETER (rtod = 180.0/pi)  !radians to degreesc
      PARAMETER (mu0 = pi*4.0e-7)  !magnetic permeability of free space
      PARAMETER (epsilon = 8.85e-12) !dielectric constant
      PARAMETER (q = 1.6e-19)        !electron charge (Coulombs)
      PARAMETER (melec = 9.1e-31)  !mass of electron (kg)
      PARAMETER (mproton = 1.67e-27)
      PARAMETER (mO = 1.67e-27)    !mass of H (kg)
c      PARAMETER (mO = 2.3e-25)    !mass of Ba (kg)
      PARAMETER (mBa = 17.0*mO)    !mass of Ba (kg)
      PARAMETER (O_to_Ba = mO/mBa) !convert E and B for particle move
      PARAMETER (km_to_m = 1e3)    !km to meters conversion
      PARAMETER (kboltz = 1.38e-29)   !kg km^2 / s^2 / K
      PARAMETER (tempf0 = 50*11600.)     !K
 
      real Qo, vrad, N_o, Rp, tau_photo
      integer S_radius
      PARAMETER (Qo = 1.0e28)  !neutral source rate
      PARAMETER (vrad = 1.0)   !escape velocity
c      PARAMETER (N_o = 5e34)   !Steady state neutral particle constant
      PARAMETER (N_o = 1e33)   !Steady state neutral particle constant
      PARAMETER (Rp = 20000.0)  !Pluto radius
      PARAMETER (tau_photo = 1.2e6)
c      PARAMETER (dNi = 2500)
      PARAMETER (S_radius = 3) !units of dx
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
c raw grid coordinate data
 
      real qx(nx), qy(ny), qz(nz)
                            !cartisian coordinates of grid points
      real lambda(nz)       !latitude of given grid point
      integer ri, rj, rk    !release coordinates
      real dz_grid(nz)      !dz for grid, variable in z direction
      real dz_cell(nz)      !dz for cell, variable in z direction
 
      common /coords/ qx,qy,qz,lambda,ri,rj,rk,dz_grid,dz_cell
 
 
c Total number of ions produced at a given time during the simulation
      integer*4 Ni_tot, Nn_tot
 
c Location (indices) of particles with the grid
      integer ijkp(Ni_max,3)
      logical ionized(Ni_max)
      common /ion_info/ Ni_tot, Nn_tot, ijkp, ionized
 
      real seed
      common /rndseed/ seed
 
c Indices for boundary cells
      integer ip,im,jp,jm,kp,km
 
c Total input energy (ions)
      real input_E,input_EeP,prev_Etot, bndry_Eflux
      real input_chex, input_bill
      common /energy/ input_E, input_EeP, prev_Etot, bndry_Eflux,
     x                input_chex, input_bill
 
c Momentum conservation
      real surf_tot(3),graduu_tot(3),ugradu_tot(3)
      common /momentum/ surf_tot,graduu_tot,ugradu_tot
 
c Dummy variables for doing vector operations and particle flux
c calculation....etc.
 
      real a(nx,ny,nz,3),c(nx,ny,nz,3),ct(nx,ny,nz,3)
      common /dummy/ a,c,ct
 
c Variable time step for fluid velocity update used in both
c particle update and fluid update
 
      real dtf,ntf,dtsub
      common /timestep/ ntf,dtsub
 
c Weight variables for trilinear interpolation
 
      real wght(Ni_max,8), wquad(Ni_max,3)
      common /weights/ wght, wquad
 
c variable for anomlous resistivity
 
      real etar(nx,ny,nz,3)
      real nuei
      common /resist/ etar,nuei
 
c variable for particle scale
 
      real beta, dNi
      common /scale/ beta, dNi
 
 
c mass array for multi-ion species
      real mrat,m_arr
      common /mass/ mrat(Ni_max),m_arr(Ni_max)
 
 
 
      real u(nx,ny,nz,3),  !main cell contravarient fluid velocity (uf)
     x     b(nx,ny,nz,3),   !main cell contravarient bt (btmf)
     x     cc(nx,ny,nz,3)   !(uf.B)*B
 
      real ux,uy,uz,bx,by,bz   !temp vars for u and b at cell center
      real temp                !used to vectorize loop
c      real ct(nx,ny,nz,3)      !result are main cell center
      real udotb               !u dot b
      real zfrc(nz)            !0.5*dz_grid(k)/dz_cell(k)
 
! first gather everything at center
 
      external vp_get_um_dot_bb1
      external vp_get_um_dot_bb2
      integer j3, j4, j5, j6, j7, j8
      real r1, r2, r3, r4, r5, r6
      call periodic(u)
      call periodic(b)
      call fix_normal_b(b)
 
      do 5 k=1,nz
         zfrc(k) = 0.5*dz_grid(k)/dz_cell(k)
 5       continue
 
      call vp_pcall (vp_get_um_dot_bb1, 2, 3, b, zfrc, u)
 
      call periodic(ct)
 
 
c extrapolate back to main cell contravarient positions.
c ...just average across cells.
 
      call vp_pcall (vp_get_um_dot_bb2, 2, 1, cc)
 
      call periodic(cc)
 
      return
      end
      subroutine vp_get_um_dot_bb1(vp_mythread,vp_nmthreads,b,zfrc,u)
      common /dummy/a, c, ct
      REAL b(79,19,37,3), zfrc(37), u(79,19,37,3), ct(79,19,37,3), a(79,
     1   19,37,3), c(79,19,37,3)
      REAL ux, bx, uy, by, uz, bz, udotb
      INTEGER k,km,j,jm,i,im,j3,j4,j5,vp_mythread,vp_nmthreads,j1,j2
 
      j3 = (35 + vp_nmthreads - 1)/vp_nmthreads
      j4 = vp_mythread*j3 + 1
      j5 = min(35 - j4 + 1,j3)
      do k = j4, j5 + j4 - 1
         do j = 1, 17
            do i = 1, 77
               bx = 0.5*(b(1+i,j+1,k+1,1)+b(i,j+1,k+1,1))
               by = 0.5*(b(1+i,j+1,k+1,2)+b(1+i,j,k+1,2))
 
               uz = zfrc(k+1)*(u(1+i,j+1,k+1,3)-u(1+i,j+1,k,3)) + u(1+i,
     1            j+1,k,3)
               bz = zfrc(k+1)*(b(1+i,j+1,k+1,3)-b(1+i,j+1,k,3)) + b(1+i,
     1            j+1,k,3)
 
               udotb = 0.5*(u(1+i,j+1,k+1,1)+u(i,j+1,k+1,1))*bx + 0.5*(u
     1            (1+i,j+1,k+1,2)+u(1+i,j,k+1,2))*by + uz*bz
 
               ct(1+i,j+1,k+1,1) = udotb*bx
               ct(1+i,j+1,k+1,2) = udotb*by
               ct(1+i,j+1,k+1,3) = udotb*bz
            end do
 
         end do
      end do
      return 
      end 
      subroutine vp_get_um_dot_bb2(vp_mythread, vp_nmthreads, cc)
      common /dummy/a, c, ct
      REAL ct(79,19,37,3), cc(79,19,37,3), a(79,19,37,3), c(79,19,37,3)
      REAL r1, r2, r3, r4, r5, r6
      INTEGER k,kp,j,jp,i,ip,j6,j7,j8,vp_mythread,vp_nmthreads,j1,j2
 
 
c extrapolate back to main cell contravarient positions.
c ...just average across cells.
 
      j6 = (35 + vp_nmthreads - 1)/vp_nmthreads
      j7 = vp_mythread*j6 + 1
      j8 = min(35 - j7 + 1,j6)
      do k = j7, j8 + j7 - 1
         j = 1
         do i = 1, 77
 
            cc(1+i,2,k+1,1) = 0.5*(ct(1+i,2,k+1,1)+ct(2+i,2,k+1,1))
            cc(1+i,2,k+1,2) = 0.5*(ct(1+i,2,k+1,2)+ct(1+i,3,k+1,2))
            cc(1+i,2,k+1,3) = 0.5*(ct(1+i,2,k+1,3)+ct(1+i,2,k+2,3))
         end do
 
         do j = 1, 8
            do i = 1, 77
               r1 = 0.5*(ct(1+i,j*2+1,k+1,1)+ct(2+i,j*2+1,k+1,1))
               r2 = 0.5*(ct(1+i,(j+1)*2,k+1,1)+ct(2+i,(j+1)*2,k+1,1))
 
               cc(1+i,j*2+1,k+1,1) = r1
               cc(1+i,(j+1)*2,k+1,1) = r2
               r3 = 0.5*(ct(1+i,j*2+1,k+1,2)+ct(1+i,(j+1)*2,k+1,2))
               r4 = 0.5*(ct(1+i,(j+1)*2,k+1,2)+ct(1+i,j*2+3,k+1,2))
               cc(1+i,j*2+1,k+1,2) = r3
               cc(1+i,(j+1)*2,k+1,2) = r4
               r5 = 0.5*(ct(1+i,j*2+1,k+1,3)+ct(1+i,j*2+1,k+2,3))
               r6 = 0.5*(ct(1+i,(j+1)*2,k+1,3)+ct(1+i,(j+1)*2,k+2,3))
               cc(1+i,j*2+1,k+1,3) = r5
               cc(1+i,(j+1)*2,k+1,3) = r6
            end do
 
         end do
      end do
      return 
      end 
c----------------------------------------------------------------------
 
 
c----------------------------------------------------------------------
      SUBROUTINE get_um_dot_BB_old(u,b,cc)
c uf and btmf are gathered at main cell center and uf.B*B
c calculated.  Result returned to main cell contravarient
c postion.
c----------------------------------------------------------------------
      PARAMETER (nx = 79, ny = 19, nz = 37)
 
c grid parameters
      PARAMETER (dx = 10000.0, dy = 10000.0)   !units in km
      PARAMETER (delz = 10000.0)          !dz at release coordinates
                                       !this is not the final dz
                                       !d_lamba determines dz
c time stepping parameters
      PARAMETER (tstart = 0.2)   !starting time
      PARAMETER (dt = 1.)     !main time step
      PARAMETER (nt = 5.0)        !number of time steps
c      PARAMETER (dtsub = 0.01) !subcycle time step
      PARAMETER (dtsub_init = 0.1) !subcycle time step
      PARAMETER (ntsub = 10.0)   !number of subcycle time steps
      PARAMETER (nout = 50)      !number of steps to diagnosic output
 
c logical variable for restart
C...Translated by Crescent Bay Software VAST-F    7.0K8 15:17:00   4/ 1/04    
      integer im1u,im2u
      logical restart
      PARAMETER (restart = .false.)
      PARAMETER (mrestart = 500)      ! use -1 for no save
      PARAMETER (mbegin = 0)      !mstart
 
c inital ambient density at the release point.
      real nf_init,b0_init,nn_coef
      PARAMETER (nf_init = 4e15)   !was 0.01
      PARAMETER (b0_init = 2e-9)    !was 0.2
      PARAMETER (nn_coef = 1e5)
 
c electron ion collision frequency
      real nu_init, eta_init
      PARAMETER (nu_init = 0.001)
      PARAMETER (eta_init = 0.0)
 
c temporary directory for output
      character*9 temp_dir
      PARAMETER (temp_dir = '/tmp/pad/')
 
c total number of neutrals released
      real No
      PARAMETER (No=1.0e28)
 
c neutral cloud expansion characteristics
      PARAMETER(vo = 20.0)
      PARAMETER(vth = 19.3)
      PARAMETER(vsat = 0.0)
      PARAMETER(vsw = 340.0)
 
c total kinetic energy of released neutrals
      PARAMETER(Etot = 2.696e7)       !units of kg m^2 / s^2
 
c density scaling parameter, alpha, and ion particle array dims
 
      real alpha
 
c      PARAMETER (beta =  8.5227271e-18)
      PARAMETER (alpha = 1.9263418e-20) !mH
 
c max number of ion particles to be produced.  This parameter
c is used to dimension the particle arrays.
      integer*4 Ni_max
      PARAMETER (Ni_max = 10000)
 
c earth's magnetic dipole field geometry
c earth radius and equatorial B alt (defining L shell)
      PARAMETER (r_earth=6.38e3)                !earth radius
      PARAMETER (B_max_alt=2.0e3+r_earth)       !ro defining L shell
 
c misc constants
      real mu0,epsilon,pi,rtod,q,mO,mBa,O_to_Ba,km_to_m,kboltz,melec
      real mproton,tempf0
      PARAMETER (pi = 3.14159)
      PARAMETER (rtod = 180.0/pi)  !radians to degreesc
      PARAMETER (mu0 = pi*4.0e-7)  !magnetic permeability of free space
      PARAMETER (epsilon = 8.85e-12) !dielectric constant
      PARAMETER (q = 1.6e-19)        !electron charge (Coulombs)
      PARAMETER (melec = 9.1e-31)  !mass of electron (kg)
      PARAMETER (mproton = 1.67e-27)
      PARAMETER (mO = 1.67e-27)    !mass of H (kg)
c      PARAMETER (mO = 2.3e-25)    !mass of Ba (kg)
      PARAMETER (mBa = 17.0*mO)    !mass of Ba (kg)
      PARAMETER (O_to_Ba = mO/mBa) !convert E and B for particle move
      PARAMETER (km_to_m = 1e3)    !km to meters conversion
      PARAMETER (kboltz = 1.38e-29)   !kg km^2 / s^2 / K
      PARAMETER (tempf0 = 50*11600.)     !K
 
      real Qo, vrad, N_o, Rp, tau_photo
      integer S_radius
      PARAMETER (Qo = 1.0e28)  !neutral source rate
      PARAMETER (vrad = 1.0)   !escape velocity
c      PARAMETER (N_o = 5e34)   !Steady state neutral particle constant
      PARAMETER (N_o = 1e33)   !Steady state neutral particle constant
      PARAMETER (Rp = 20000.0)  !Pluto radius
      PARAMETER (tau_photo = 1.2e6)
c      PARAMETER (dNi = 2500)
      PARAMETER (S_radius = 3) !units of dx
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
c raw grid coordinate data
 
      real qx(nx), qy(ny), qz(nz)
                            !cartisian coordinates of grid points
      real lambda(nz)       !latitude of given grid point
      integer ri, rj, rk    !release coordinates
      real dz_grid(nz)      !dz for grid, variable in z direction
      real dz_cell(nz)      !dz for cell, variable in z direction
 
      common /coords/ qx,qy,qz,lambda,ri,rj,rk,dz_grid,dz_cell
 
 
c Total number of ions produced at a given time during the simulation
      integer*4 Ni_tot, Nn_tot
 
c Location (indices) of particles with the grid
      integer ijkp(Ni_max,3)
      logical ionized(Ni_max)
      common /ion_info/ Ni_tot, Nn_tot, ijkp, ionized
 
      real seed
      common /rndseed/ seed
 
c Indices for boundary cells
      integer ip,im,jp,jm,kp,km
 
c Total input energy (ions)
      real input_E,input_EeP,prev_Etot, bndry_Eflux
      real input_chex, input_bill
      common /energy/ input_E, input_EeP, prev_Etot, bndry_Eflux,
     x                input_chex, input_bill
 
c Momentum conservation
      real surf_tot(3),graduu_tot(3),ugradu_tot(3)
      common /momentum/ surf_tot,graduu_tot,ugradu_tot
 
c Dummy variables for doing vector operations and particle flux
c calculation....etc.
 
      real a(nx,ny,nz,3),c(nx,ny,nz,3),ct(nx,ny,nz,3)
      common /dummy/ a,c,ct
 
c Variable time step for fluid velocity update used in both
c particle update and fluid update
 
      real dtf,ntf,dtsub
      common /timestep/ ntf,dtsub
 
c Weight variables for trilinear interpolation
 
      real wght(Ni_max,8), wquad(Ni_max,3)
      common /weights/ wght, wquad
 
c variable for anomlous resistivity
 
      real etar(nx,ny,nz,3)
      real nuei
      common /resist/ etar,nuei
 
c variable for particle scale
 
      real beta, dNi
      common /scale/ beta, dNi
 
 
c mass array for multi-ion species
      real mrat,m_arr
      common /mass/ mrat(Ni_max),m_arr(Ni_max)
 
 
 
      real u(nx,ny,nz,3),  !main cell contravarient fluid velocity (uf)
     x     b(nx,ny,nz,3),   !main cell contravarient bt (btmf)
     x     cc(nx,ny,nz,3)   !(uf.B)*B
 
      real ux,uy,uz,bx,by,bz   !temp vars for u and b at cell center
      real temp                !used to vectorize loop
c      real ct(nx,ny,nz,3)      !result are main cell center
      real udotb               !u dot b
      real zfrc(nz)            !0.5*dz_grid(k)/dz_cell(k)
 
! first gather everything at center
 
      external vp_get_um_dot_bb_old1
      integer j3, j4, j5
      call periodic(u)
      call periodic(b)
      call fix_normal_b(b)
 
      do 5 k=1,nz
         zfrc(k) = 0.5*dz_grid(k)/dz_cell(k)
 5       continue
 
      call vp_pcall (vp_get_um_dot_bb_old1, 2, 3, u, b, zfrc)
 
      call periodic(ct)
 
c extrapolate back to main cell contravarient positions.
c ...just average across cells.
 
      do 60 i = 2, 79
         if (i .eq. 79-1) then
            do j = 2, 19
               if (j .eq. 19-1) then
                  do k = 1, 36
                     cc(79,j,1+k,1) = 2.0*ct(79,j,1+k,1) - 0.5*(ct(79,j,
     1                  1+k,1)+ct(78,j,1+k,1))
                     cc(i,19,1+k,2) = 2.0*ct(i,19,1+k,2) - 0.5*(ct(i,19,
     1                  1+k,2)+ct(i,18,1+k,2))
 
                     if (1 + k .eq. 36) then
c                  temp = 2.0*ct(i,j,nz,3) -
c     x                           0.5*(ct(i,j,nz,3) + ct(i,j,nz-1,3))
                        temp = (0.5*(ct(i,j,37,3)+ct(i,j,36,3))) + (2.0*
     1                     dz_cell(37)/dz_grid(37))*(ct(i,j,37,3)-(0.5*(
     2                     ct(i,j,37,3)+ct(i,j,36,3))))
                     else
                        temp = 0.5*(ct(i,j,1+k,3)+ct(i,j,2+k,3))
                     endif
 
                     cc(i,j,1+k,3) = temp
                  end do
 
                  ip = i + 1
                  jp = j + 1
               else
                  jp = j + 1
                  do k = 1, 36
                     cc(79,j,1+k,1) = 2.0*ct(79,j,1+k,1) - 0.5*(ct(79,j,
     1                  1+k,1)+ct(78,j,1+k,1))
                     cc(i,j,1+k,2) = 0.5*(ct(i,j,1+k,2)+ct(i,jp,1+k,2))
 
                     if (1 + k .eq. 36) then
c                  temp = 2.0*ct(i,j,nz,3) -
c     x                           0.5*(ct(i,j,nz,3) + ct(i,j,nz-1,3))
                        temp = (0.5*(ct(i,j,37,3)+ct(i,j,36,3))) + (2.0*
     1                     dz_cell(37)/dz_grid(37))*(ct(i,j,37,3)-(0.5*(
     2                     ct(i,j,37,3)+ct(i,j,36,3))))
                     else
                        temp = 0.5*(ct(i,j,1+k,3)+ct(i,j,2+k,3))
                     endif
 
                     cc(i,j,1+k,3) = temp
                  end do
 
                  ip = i + 1
               endif
            end do
         else
            do j = 2, 19
               if (j .eq. 19-1) then
 
                  ip = i + 1
                  do k = 1, 36
                     cc(i,j,1+k,1) = 0.5*(ct(i,j,1+k,1)+ct(ip,j,1+k,1))
                     cc(i,19,1+k,2) = 2.0*ct(i,19,1+k,2) - 0.5*(ct(i,19,
     1                  1+k,2)+ct(i,18,1+k,2))
 
                     if (1 + k .eq. 36) then
c                  temp = 2.0*ct(i,j,nz,3) -
c     x                           0.5*(ct(i,j,nz,3) + ct(i,j,nz-1,3))
                        temp = (0.5*(ct(i,j,37,3)+ct(i,j,36,3))) + (2.0*
     1                     dz_cell(37)/dz_grid(37))*(ct(i,j,37,3)-(0.5*(
     2                     ct(i,j,37,3)+ct(i,j,36,3))))
                     else
                        temp = 0.5*(ct(i,j,1+k,3)+ct(i,j,2+k,3))
                     endif
 
                     cc(i,j,1+k,3) = temp
                  end do
 
                  jp = j + 1
               else
 
                  ip = i + 1
                  jp = j + 1
                  do k = 1, 36
                     cc(i,j,1+k,1) = 0.5*(ct(i,j,1+k,1)+ct(ip,j,1+k,1))
                     cc(i,j,1+k,2) = 0.5*(ct(i,j,1+k,2)+ct(i,jp,1+k,2))
 
                     if (1 + k .eq. 36) then
c                  temp = 2.0*ct(i,j,nz,3) -
c     x                           0.5*(ct(i,j,nz,3) + ct(i,j,nz-1,3))
                        temp = (0.5*(ct(i,j,37,3)+ct(i,j,36,3))) + (2.0*
     1                     dz_cell(37)/dz_grid(37))*(ct(i,j,37,3)-(0.5*(
     2                     ct(i,j,37,3)+ct(i,j,36,3))))
                     else
                        temp = 0.5*(ct(i,j,1+k,3)+ct(i,j,2+k,3))
                     endif
 
                     cc(i,j,1+k,3) = temp
                  end do
 
               endif
            end do
         endif
   60 continue
 
      call periodic(cc)
 
      return
      end
      subroutine vp_get_um_dot_bb_old1(vp_mythread, vp_nmthreads, u, b, 
     1   zfrc)
      common /dummy/a, c, ct
      REAL u(79,19,37,3), b(79,19,37,3), zfrc(37), ct(79,19,37,3), a(79,
     1   19,37,3), c(79,19,37,3)
      REAL ux, bx, uy, by, uz, bz, udotb
      INTEGER i,j,im,jm,k,km,j3,j4,j5,vp_mythread,vp_nmthreads,j1,j2
C...Translated by Crescent Bay Software VAST-F    7.0K8 15:17:00   4/ 1/04    
      integer im1u,im2u
 
      j3 = (77 + vp_nmthreads - 1)/vp_nmthreads
      j4 = vp_mythread*j3 + 1
      j5 = min(77 - j4 + 1,j3)
      do i = j4, j5 + j4 - 1
         if (1 + i .eq. 2) then
            do j = 1, 17
               if (1 + j .eq. 2) then
                  do k = 1, 35
                     ux = 2.0*u(2,j+1,1+k,1) - 0.5*(u(3,j+1,1+k,1)+u(2,j
     1                  +1,1+k,1))
                     bx = 2.0*b(2,j+1,1+k,1) - 0.5*(b(3,j+1,1+k,1)+b(2,j
     1                  +1,1+k,1))
                     uy = 2.0*u(i+1,2,1+k,2) - 0.5*(u(i+1,3,1+k,2)+u(i+1
     1                  ,2,1+k,2))
                     by = 2.0*b(i+1,2,1+k,2) - 0.5*(b(i+1,3,1+k,2)+b(i+1
     1                  ,2,1+k,2))
 
                     if (1 + k .eq. 2) then
                        uz = 2.0*u(i+1,j+1,2,3) - zfrc(1+k)*(u(i+1,j+1,3
     1                     ,3)-u(i+1,j+1,2,3)) + u(i+1,j+1,2,3)
c                  uz = 2.0*u(i,j,2,3) -
c     x                 0.5*(u(i,j,3,3) + u(i,j,2,3))
                        bz = b(i+1,j+1,2,3)
c                  bz = 2.0*b(i,j,2,3) -
c     x                 0.5*(b(i,j,3,3) + b(i,j,2,3))
                     else
                        uz = zfrc(1+k)*(u(i+1,j+1,1+k,3)-u(i+1,j+1,k,3))
     1                      + u(i+1,j+1,k,3)
                        bz = zfrc(1+k)*(b(i+1,j+1,1+k,3)-b(i+1,j+1,k,3))
     1                      + b(i+1,j+1,k,3)
c                  uz = 0.5*(u(i,j,k,3) + u(i,j,km,3))
c                  bz = 0.5*(b(i,j,k,3) + b(i,j,km,3))
                     endif
 
                     udotb = ux*bx + uy*by + uz*bz
 
                     ct(i+1,j+1,1+k,1) = udotb*bx
                     ct(i+1,j+1,1+k,2) = udotb*by
                     ct(i+1,j+1,1+k,3) = udotb*bz
                  end do
 
                  im1u = i
               else
                  do k = 1, 35
                     ux = 2.0*u(2,j+1,1+k,1) - 0.5*(u(3,j+1,1+k,1)+u(2,j
     1                  +1,1+k,1))
                     bx = 2.0*b(2,j+1,1+k,1) - 0.5*(b(3,j+1,1+k,1)+b(2,j
     1                  +1,1+k,1))
                     uy = 0.5*(u(i+1,j+1,1+k,2)+u(i+1,j,1+k,2))
                     by = 0.5*(b(i+1,j+1,1+k,2)+b(i+1,j,1+k,2))
 
                     if (1 + k .eq. 2) then
                        uz = 2.0*u(i+1,j+1,2,3) - zfrc(1+k)*(u(i+1,j+1,3
     1                     ,3)-u(i+1,j+1,2,3)) + u(i+1,j+1,2,3)
c                  uz = 2.0*u(i,j,2,3) -
c     x                 0.5*(u(i,j,3,3) + u(i,j,2,3))
                        bz = b(i+1,j+1,2,3)
c                  bz = 2.0*b(i,j,2,3) -
c     x                 0.5*(b(i,j,3,3) + b(i,j,2,3))
                     else
                        uz = zfrc(1+k)*(u(i+1,j+1,1+k,3)-u(i+1,j+1,k,3))
     1                      + u(i+1,j+1,k,3)
                        bz = zfrc(1+k)*(b(i+1,j+1,1+k,3)-b(i+1,j+1,k,3))
     1                      + b(i+1,j+1,k,3)
c                  uz = 0.5*(u(i,j,k,3) + u(i,j,km,3))
c                  bz = 0.5*(b(i,j,k,3) + b(i,j,km,3))
                     endif
 
                     udotb = ux*bx + uy*by + uz*bz
 
                     ct(i+1,j+1,1+k,1) = udotb*bx
                     ct(i+1,j+1,1+k,2) = udotb*by
                     ct(i+1,j+1,1+k,3) = udotb*bz
                  end do
 
                  im1u = i
               endif
            end do
         else
            do j = 1, 17
               if (1 + j .eq. 2) then
 
                  im2u = i
                  do k = 1, 35
                     ux = 0.5*(u(i+1,j+1,1+k,1)+u(im2u,1+j,1+k,1))
                     bx = 0.5*(b(i+1,j+1,1+k,1)+b(im2u,1+j,1+k,1))
                     uy = 2.0*u(i+1,2,1+k,2) - 0.5*(u(i+1,3,1+k,2)+u(i+1
     1                  ,2,1+k,2))
                     by = 2.0*b(i+1,2,1+k,2) - 0.5*(b(i+1,3,1+k,2)+b(i+1
     1                  ,2,1+k,2))
 
                     if (1 + k .eq. 2) then
                        uz = 2.0*u(i+1,j+1,2,3) - zfrc(1+k)*(u(i+1,j+1,3
     1                     ,3)-u(i+1,j+1,2,3)) + u(i+1,j+1,2,3)
c                  uz = 2.0*u(i,j,2,3) -
c     x                 0.5*(u(i,j,3,3) + u(i,j,2,3))
                        bz = b(i+1,j+1,2,3)
c                  bz = 2.0*b(i,j,2,3) -
c     x                 0.5*(b(i,j,3,3) + b(i,j,2,3))
                     else
                        uz = zfrc(1+k)*(u(i+1,j+1,1+k,3)-u(i+1,j+1,k,3))
     1                      + u(i+1,j+1,k,3)
                        bz = zfrc(1+k)*(b(i+1,j+1,1+k,3)-b(i+1,j+1,k,3))
     1                      + b(i+1,j+1,k,3)
c                  uz = 0.5*(u(i,j,k,3) + u(i,j,km,3))
c                  bz = 0.5*(b(i,j,k,3) + b(i,j,km,3))
                     endif
 
                     udotb = ux*bx + uy*by + uz*bz
 
                     ct(i+1,j+1,1+k,1) = udotb*bx
                     ct(i+1,j+1,1+k,2) = udotb*by
                     ct(i+1,j+1,1+k,3) = udotb*bz
                  end do
 
               else
 
                  im2u = i
                  do k = 1, 35
                     ux = 0.5*(u(i+1,j+1,1+k,1)+u(im2u,1+j,1+k,1))
                     bx = 0.5*(b(i+1,j+1,1+k,1)+b(im2u,1+j,1+k,1))
                     uy = 0.5*(u(i+1,j+1,1+k,2)+u(i+1,j,1+k,2))
                     by = 0.5*(b(i+1,j+1,1+k,2)+b(i+1,j,1+k,2))
 
                     if (1 + k .eq. 2) then
                        uz = 2.0*u(i+1,j+1,2,3) - zfrc(1+k)*(u(i+1,j+1,3
     1                     ,3)-u(i+1,j+1,2,3)) + u(i+1,j+1,2,3)
c                  uz = 2.0*u(i,j,2,3) -
c     x                 0.5*(u(i,j,3,3) + u(i,j,2,3))
                        bz = b(i+1,j+1,2,3)
c                  bz = 2.0*b(i,j,2,3) -
c     x                 0.5*(b(i,j,3,3) + b(i,j,2,3))
                     else
                        uz = zfrc(1+k)*(u(i+1,j+1,1+k,3)-u(i+1,j+1,k,3))
     1                      + u(i+1,j+1,k,3)
                        bz = zfrc(1+k)*(b(i+1,j+1,1+k,3)-b(i+1,j+1,k,3))
     1                      + b(i+1,j+1,k,3)
c                  uz = 0.5*(u(i,j,k,3) + u(i,j,km,3))
c                  bz = 0.5*(b(i,j,k,3) + b(i,j,km,3))
                     endif
 
                     udotb = ux*bx + uy*by + uz*bz
 
                     ct(i+1,j+1,1+k,1) = udotb*bx
                     ct(i+1,j+1,1+k,2) = udotb*by
                     ct(i+1,j+1,1+k,3) = udotb*bz
                  end do
 
               endif
            end do
         endif
      end do
      return 
      end 
c----------------------------------------------------------------------
 
 
c----------------------------------------------------------------------
      SUBROUTINE get_ugradu_Lax(uf,ugradu,delta_t)
c----------------------------------------------------------------------
      PARAMETER (nx = 79, ny = 19, nz = 37)
 
c grid parameters
      PARAMETER (dx = 10000.0, dy = 10000.0)   !units in km
      PARAMETER (delz = 10000.0)          !dz at release coordinates
                                       !this is not the final dz
                                       !d_lamba determines dz
c time stepping parameters
      PARAMETER (tstart = 0.2)   !starting time
      PARAMETER (dt = 1.)     !main time step
      PARAMETER (nt = 5.0)        !number of time steps
c      PARAMETER (dtsub = 0.01) !subcycle time step
      PARAMETER (dtsub_init = 0.1) !subcycle time step
      PARAMETER (ntsub = 10.0)   !number of subcycle time steps
      PARAMETER (nout = 50)      !number of steps to diagnosic output
 
c logical variable for restart
C...Translated by Crescent Bay Software VAST-F    7.0K8 15:17:00   4/ 1/04    
      logical restart
      PARAMETER (restart = .false.)
      PARAMETER (mrestart = 500)      ! use -1 for no save
      PARAMETER (mbegin = 0)      !mstart
 
c inital ambient density at the release point.
      real nf_init,b0_init,nn_coef
      PARAMETER (nf_init = 4e15)   !was 0.01
      PARAMETER (b0_init = 2e-9)    !was 0.2
      PARAMETER (nn_coef = 1e5)
 
c electron ion collision frequency
      real nu_init, eta_init
      PARAMETER (nu_init = 0.001)
      PARAMETER (eta_init = 0.0)
 
c temporary directory for output
      character*9 temp_dir
      PARAMETER (temp_dir = '/tmp/pad/')
 
c total number of neutrals released
      real No
      PARAMETER (No=1.0e28)
 
c neutral cloud expansion characteristics
      PARAMETER(vo = 20.0)
      PARAMETER(vth = 19.3)
      PARAMETER(vsat = 0.0)
      PARAMETER(vsw = 340.0)
 
c total kinetic energy of released neutrals
      PARAMETER(Etot = 2.696e7)       !units of kg m^2 / s^2
 
c density scaling parameter, alpha, and ion particle array dims
 
      real alpha
 
c      PARAMETER (beta =  8.5227271e-18)
      PARAMETER (alpha = 1.9263418e-20) !mH
 
c max number of ion particles to be produced.  This parameter
c is used to dimension the particle arrays.
      integer*4 Ni_max
      PARAMETER (Ni_max = 10000)
 
c earth's magnetic dipole field geometry
c earth radius and equatorial B alt (defining L shell)
      PARAMETER (r_earth=6.38e3)                !earth radius
      PARAMETER (B_max_alt=2.0e3+r_earth)       !ro defining L shell
 
c misc constants
      real mu0,epsilon,pi,rtod,q,mO,mBa,O_to_Ba,km_to_m,kboltz,melec
      real mproton,tempf0
      PARAMETER (pi = 3.14159)
      PARAMETER (rtod = 180.0/pi)  !radians to degreesc
      PARAMETER (mu0 = pi*4.0e-7)  !magnetic permeability of free space
      PARAMETER (epsilon = 8.85e-12) !dielectric constant
      PARAMETER (q = 1.6e-19)        !electron charge (Coulombs)
      PARAMETER (melec = 9.1e-31)  !mass of electron (kg)
      PARAMETER (mproton = 1.67e-27)
      PARAMETER (mO = 1.67e-27)    !mass of H (kg)
c      PARAMETER (mO = 2.3e-25)    !mass of Ba (kg)
      PARAMETER (mBa = 17.0*mO)    !mass of Ba (kg)
      PARAMETER (O_to_Ba = mO/mBa) !convert E and B for particle move
      PARAMETER (km_to_m = 1e3)    !km to meters conversion
      PARAMETER (kboltz = 1.38e-29)   !kg km^2 / s^2 / K
      PARAMETER (tempf0 = 50*11600.)     !K
 
      real Qo, vrad, N_o, Rp, tau_photo
      integer S_radius
      PARAMETER (Qo = 1.0e28)  !neutral source rate
      PARAMETER (vrad = 1.0)   !escape velocity
c      PARAMETER (N_o = 5e34)   !Steady state neutral particle constant
      PARAMETER (N_o = 1e33)   !Steady state neutral particle constant
      PARAMETER (Rp = 20000.0)  !Pluto radius
      PARAMETER (tau_photo = 1.2e6)
c      PARAMETER (dNi = 2500)
      PARAMETER (S_radius = 3) !units of dx
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
c raw grid coordinate data
 
      real qx(nx), qy(ny), qz(nz)
                            !cartisian coordinates of grid points
      real lambda(nz)       !latitude of given grid point
      integer ri, rj, rk    !release coordinates
      real dz_grid(nz)      !dz for grid, variable in z direction
      real dz_cell(nz)      !dz for cell, variable in z direction
 
      common /coords/ qx,qy,qz,lambda,ri,rj,rk,dz_grid,dz_cell
 
 
c Total number of ions produced at a given time during the simulation
      integer*4 Ni_tot, Nn_tot
 
c Location (indices) of particles with the grid
      integer ijkp(Ni_max,3)
      logical ionized(Ni_max)
      common /ion_info/ Ni_tot, Nn_tot, ijkp, ionized
 
      real seed
      common /rndseed/ seed
 
c Indices for boundary cells
      integer ip,im,jp,jm,kp,km
 
c Total input energy (ions)
      real input_E,input_EeP,prev_Etot, bndry_Eflux
      real input_chex, input_bill
      common /energy/ input_E, input_EeP, prev_Etot, bndry_Eflux,
     x                input_chex, input_bill
 
c Momentum conservation
      real surf_tot(3),graduu_tot(3),ugradu_tot(3)
      common /momentum/ surf_tot,graduu_tot,ugradu_tot
 
c Dummy variables for doing vector operations and particle flux
c calculation....etc.
 
      real a(nx,ny,nz,3),c(nx,ny,nz,3),ct(nx,ny,nz,3)
      common /dummy/ a,c,ct
 
c Variable time step for fluid velocity update used in both
c particle update and fluid update
 
      real dtf,ntf,dtsub
      common /timestep/ ntf,dtsub
 
c Weight variables for trilinear interpolation
 
      real wght(Ni_max,8), wquad(Ni_max,3)
      common /weights/ wght, wquad
 
c variable for anomlous resistivity
 
      real etar(nx,ny,nz,3)
      real nuei
      common /resist/ etar,nuei
 
c variable for particle scale
 
      real beta, dNi
      common /scale/ beta, dNi
 
 
c mass array for multi-ion species
      real mrat,m_arr
      common /mass/ mrat(Ni_max),m_arr(Ni_max)
 
 
 
      real uf(nx,ny,nz,3),
     x     ugradu(nx,ny,nz,3)
 
      real ufc(nx,ny,nz,3)
      real ax1,ax2,ay1,ay2,az1,az2       !misc const
      real u1,u2,u3                      !temp vars
 
      parameter(ad = 0.01)                 !coefficient to add extra
                                         !diffusion
      external vp_get_ugradu_lax1
      external vp_get_ugradu_lax2
      integer j3, j4, j5, j6, j7, j8
      real r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, 
     .r15
      call periodic(uf)
 
      call face_to_center(uf,ufc)
 
      call vp_pcall (vp_get_ugradu_lax1, 2, 1, ufc)
 
      call periodic(ct)
 
c interpolate back to contravarient positions.
 
      call vp_pcall (vp_get_ugradu_lax2, 2, 1, ugradu)
 
       call periodic(ugradu)
 
      return
      end
      subroutine vp_get_ugradu_lax1(vp_mythread, vp_nmthreads, ufc)
      common /coords/qx, qy, qz, lambda, ri, rj, rk, dz_grid, dz_cell
      common /dummy/a, c, ct
      INTEGER ri, rj, rk
      REAL ufc(79,19,37,3), dz_grid(37), ct(79,19,37,3), qx(79), qy(19)
     1   , qz(37), lambda(37), dz_cell(37), a(79,19,37,3), c(79,19,37,3)
      REAL ax1,ax2,ay1,ay2,u1,u2,az1,az2,u3,r1,r2,r3,r4,r5,r6,r7,r8,r9
      INTEGER k, kp, km, j, jp, jm, i, ip, im, j3, j4, j5, vp_mythread, 
     1   vp_nmthreads, j1, j2
 
      j3 = (35 + vp_nmthreads - 1)/vp_nmthreads
      j4 = vp_mythread*j3 + 1
      j5 = min(35 - j4 + 1,j3)
      do k = j4, j5 + j4 - 1
         do j = 1, 17
            r1 = 1./10000.0
            r2 = 1./10000.0
            r3 = 1./(dz_grid(k+1)+dz_grid(k+2))
            r4 = 1./10000.0
            r5 = 1./10000.0
            r6 = 1./(dz_grid(k+1)+dz_grid(k+2))
            r7 = 1./10000.0
            r8 = 1./10000.0
            r9 = 1./(dz_grid(k+1)+dz_grid(k+2))
            do i = 1, 77
               u1 = ((0.5*ufc(1+i,j-1+2,k-1+2,1))*r1)*(ufc(2+i,j-1+2,k-1
     1            +2,1)-ufc(i,j-1+2,k-1+2,1)) - 0.01*abs((ufc(2+i,j-1+2,
     2            k-1+2,1)-ufc(i,j-1+2,k-1+2,1)))*(ufc(i,j-1+2,k-1+2,1)-
     3            (2.0*ufc(1+i,j-1+2,k-1+2,1))+ufc(2+i,j-1+2,k-1+2,1))
               u2 = ((0.5*ufc(1+i,j-1+2,k-1+2,2))*r2)*(ufc(1+i,j-1+3,k-1
     1            +2,1)-ufc(1+i,j,k-1+2,1)) - 0.01*abs((ufc(2+i,j-1+2,k-
     2            1+2,2)-ufc(i,j-1+2,k-1+2,2)))*(ufc(1+i,j,k-1+2,1)-(2.0
     3            *ufc(1+i,j-1+2,k-1+2,1))+ufc(1+i,j-1+3,k-1+2,1))
               u3 = (ufc(1+i,j-1+2,k-1+2,3)*r3)*(ufc(1+i,j-1+2,k-1+3,1)-
     1            ufc(1+i,j-1+2,k,1)) - 0.01*abs((ufc(2+i,j-1+2,k-1+2,3)
     2            -ufc(i,j-1+2,k-1+2,3)))*(ufc(1+i,j-1+2,k,1)-(2.0*ufc(1
     3            +i,j-1+2,k-1+2,1))+ufc(1+i,j-1+2,k-1+3,1))
               ct(1+i,j+1,k+1,1) = u1 + u2 + u3
               u1 = ((0.5*ufc(1+i,j-1+2,k-1+2,1))*r4)*(ufc(2+i,j-1+2,k-1
     1            +2,2)-ufc(i,j-1+2,k-1+2,2)) - 0.01*abs((ufc(1+i,j-1+3,
     2            k-1+2,1)-ufc(1+i,j,k-1+2,1)))*(ufc(i,j-1+2,k-1+2,2)-(
     3            2.0*ufc(1+i,j-1+2,k-1+2,2))+ufc(2+i,j-1+2,k-1+2,2))
               u2 = ((0.5*ufc(1+i,j-1+2,k-1+2,2))*r5)*(ufc(1+i,j-1+3,k-1
     1            +2,2)-ufc(1+i,j,k-1+2,2)) - 0.01*abs((ufc(1+i,j-1+3,k-
     2            1+2,2)-ufc(1+i,j,k-1+2,2)))*(ufc(1+i,j,k-1+2,2)-(2.0*
     3            ufc(1+i,j-1+2,k-1+2,2))+ufc(1+i,j-1+3,k-1+2,2))
               u3 = (ufc(1+i,j-1+2,k-1+2,3)*r6)*(ufc(1+i,j-1+2,k-1+3,2)-
     1            ufc(1+i,j-1+2,k,2)) - 0.01*abs((ufc(1+i,j-1+3,k-1+2,3)
     2            -ufc(1+i,j,k-1+2,3)))*(ufc(1+i,j-1+2,k,2)-(2.0*ufc(1+i
     3            ,j-1+2,k-1+2,2))+ufc(1+i,j-1+2,k-1+3,2))
               ct(1+i,j+1,k+1,2) = u1 + u2 + u3
               u1 = ((0.5*ufc(1+i,j-1+2,k-1+2,1))*r7)*(ufc(2+i,j-1+2,k-1
     1            +2,3)-ufc(i,j-1+2,k-1+2,3)) - 0.01*abs((ufc(1+i,j-1+2,
     2            k-1+3,1)-ufc(1+i,j-1+2,k,1)))*(ufc(i,j-1+2,k-1+2,3)-(
     3            2.0*ufc(1+i,j-1+2,k-1+2,3))+ufc(2+i,j-1+2,k-1+2,3))
               u2 = ((0.5*ufc(1+i,j-1+2,k-1+2,2))*r8)*(ufc(1+i,j-1+3,k-1
     1            +2,3)-ufc(1+i,j,k-1+2,3)) - 0.01*abs((ufc(1+i,j-1+2,k-
     2            1+3,2)-ufc(1+i,j-1+2,k,2)))*(ufc(1+i,j,k-1+2,3)-(2.0*
     3            ufc(1+i,j-1+2,k-1+2,3))+ufc(1+i,j-1+3,k-1+2,3))
               u3 = (ufc(1+i,j-1+2,k-1+2,3)*r9)*(ufc(1+i,j-1+2,k-1+3,3)-
     1            ufc(1+i,j-1+2,k,3)) - 0.01*abs((ufc(1+i,j-1+2,k-1+3,3)
     2            -ufc(1+i,j-1+2,k,3)))*(ufc(1+i,j-1+2,k,3)-(2.0*ufc(1+i
     3            ,j-1+2,k-1+2,3))+ufc(1+i,j-1+2,k-1+3,3))
               ct(1+i,j+1,k+1,3) = u1 + u2 + u3
            end do
 
         end do
      end do
      return 
      end 
      subroutine vp_get_ugradu_lax2(vp_mythread, vp_nmthreads, ugradu)
      common /dummy/a, c, ct
      REAL ct(79,19,37,3), ugradu(79,19,37,3), a(79,19,37,3), c(79,19,37
     1   ,3)
      REAL r10, r11, r12, r13, r14, r15
      INTEGER k,kp,j,jp,i,ip,j6,j7,j8,vp_mythread,vp_nmthreads,j1,j2
 
c interpolate back to contravarient positions.
 
      j6 = (35 + vp_nmthreads - 1)/vp_nmthreads
      j7 = vp_mythread*j6 + 1
      j8 = min(35 - j7 + 1,j6)
      do k = j7, j8 + j7 - 1
         j = 1
         do i = 1, 77
            ugradu(1+i,2,k+1,1) = 0.5*(ct(1+i,2,k+1,1)+ct(2+i,2,k+1,1))
            ugradu(1+i,2,k+1,2) = 0.5*(ct(1+i,2,k+1,2)+ct(1+i,3,k+1,2))
            ugradu(1+i,2,k+1,3) = 0.5*(ct(1+i,2,k+1,3)+ct(1+i,2,k+2,3))
         end do
         do j = 1, 8
            do i = 1, 77
               r10 = 0.5*(ct(1+i,j*2+1,k+1,1)+ct(2+i,j*2+1,k+1,1))
               r11 = 0.5*(ct(1+i,(j+1)*2,k+1,1)+ct(2+i,(j+1)*2,k+1,1))
               ugradu(1+i,j*2+1,k+1,1) = r10
               ugradu(1+i,(j+1)*2,k+1,1) = r11
               r12 = 0.5*(ct(1+i,j*2+1,k+1,2)+ct(1+i,(j+1)*2,k+1,2))
               r13 = 0.5*(ct(1+i,(j+1)*2,k+1,2)+ct(1+i,j*2+3,k+1,2))
               ugradu(1+i,j*2+1,k+1,2) = r12
               ugradu(1+i,(j+1)*2,k+1,2) = r13
               r14 = 0.5*(ct(1+i,j*2+1,k+1,3)+ct(1+i,j*2+1,k+2,3))
               r15 = 0.5*(ct(1+i,(j+1)*2,k+1,3)+ct(1+i,(j+1)*2,k+2,3))
               ugradu(1+i,j*2+1,k+1,3) = r14
               ugradu(1+i,(j+1)*2,k+1,3) = r15
            end do
         end do
      end do
      return 
      end 
c----------------------------------------------------------------------
 
 
c----------------------------------------------------------------------
      SUBROUTINE get_Ef(Ef,aj,np,nf,up,uf,btmf,nu,ugradu,delta_t,
     x                  gradPf)
c Need to treat boundaries separately!!
c----------------------------------------------------------------------
      PARAMETER (nx = 79, ny = 19, nz = 37)
 
c grid parameters
      PARAMETER (dx = 10000.0, dy = 10000.0)   !units in km
      PARAMETER (delz = 10000.0)          !dz at release coordinates
                                       !this is not the final dz
                                       !d_lamba determines dz
c time stepping parameters
      PARAMETER (tstart = 0.2)   !starting time
      PARAMETER (dt = 1.)     !main time step
      PARAMETER (nt = 5.0)        !number of time steps
c      PARAMETER (dtsub = 0.01) !subcycle time step
      PARAMETER (dtsub_init = 0.1) !subcycle time step
      PARAMETER (ntsub = 10.0)   !number of subcycle time steps
      PARAMETER (nout = 50)      !number of steps to diagnosic output
 
c logical variable for restart
C...Translated by Crescent Bay Software VAST-F    7.0K8 15:17:00   4/ 1/04    
      logical restart
      PARAMETER (restart = .false.)
      PARAMETER (mrestart = 500)      ! use -1 for no save
      PARAMETER (mbegin = 0)      !mstart
 
c inital ambient density at the release point.
      real nf_init,b0_init,nn_coef
      PARAMETER (nf_init = 4e15)   !was 0.01
      PARAMETER (b0_init = 2e-9)    !was 0.2
      PARAMETER (nn_coef = 1e5)
 
c electron ion collision frequency
      real nu_init, eta_init
      PARAMETER (nu_init = 0.001)
      PARAMETER (eta_init = 0.0)
 
c temporary directory for output
      character*9 temp_dir
      PARAMETER (temp_dir = '/tmp/pad/')
 
c total number of neutrals released
      real No
      PARAMETER (No=1.0e28)
 
c neutral cloud expansion characteristics
      PARAMETER(vo = 20.0)
      PARAMETER(vth = 19.3)
      PARAMETER(vsat = 0.0)
      PARAMETER(vsw = 340.0)
 
c total kinetic energy of released neutrals
      PARAMETER(Etot = 2.696e7)       !units of kg m^2 / s^2
 
c density scaling parameter, alpha, and ion particle array dims
 
      real alpha
 
c      PARAMETER (beta =  8.5227271e-18)
      PARAMETER (alpha = 1.9263418e-20) !mH
 
c max number of ion particles to be produced.  This parameter
c is used to dimension the particle arrays.
      integer*4 Ni_max
      PARAMETER (Ni_max = 10000)
 
c earth's magnetic dipole field geometry
c earth radius and equatorial B alt (defining L shell)
      PARAMETER (r_earth=6.38e3)                !earth radius
      PARAMETER (B_max_alt=2.0e3+r_earth)       !ro defining L shell
 
c misc constants
      real mu0,epsilon,pi,rtod,q,mO,mBa,O_to_Ba,km_to_m,kboltz,melec
      real mproton,tempf0
      PARAMETER (pi = 3.14159)
      PARAMETER (rtod = 180.0/pi)  !radians to degreesc
      PARAMETER (mu0 = pi*4.0e-7)  !magnetic permeability of free space
      PARAMETER (epsilon = 8.85e-12) !dielectric constant
      PARAMETER (q = 1.6e-19)        !electron charge (Coulombs)
      PARAMETER (melec = 9.1e-31)  !mass of electron (kg)
      PARAMETER (mproton = 1.67e-27)
      PARAMETER (mO = 1.67e-27)    !mass of H (kg)
c      PARAMETER (mO = 2.3e-25)    !mass of Ba (kg)
      PARAMETER (mBa = 17.0*mO)    !mass of Ba (kg)
      PARAMETER (O_to_Ba = mO/mBa) !convert E and B for particle move
      PARAMETER (km_to_m = 1e3)    !km to meters conversion
      PARAMETER (kboltz = 1.38e-29)   !kg km^2 / s^2 / K
      PARAMETER (tempf0 = 50*11600.)     !K
 
      real Qo, vrad, N_o, Rp, tau_photo
      integer S_radius
      PARAMETER (Qo = 1.0e28)  !neutral source rate
      PARAMETER (vrad = 1.0)   !escape velocity
c      PARAMETER (N_o = 5e34)   !Steady state neutral particle constant
      PARAMETER (N_o = 1e33)   !Steady state neutral particle constant
      PARAMETER (Rp = 20000.0)  !Pluto radius
      PARAMETER (tau_photo = 1.2e6)
c      PARAMETER (dNi = 2500)
      PARAMETER (S_radius = 3) !units of dx
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
c raw grid coordinate data
 
      real qx(nx), qy(ny), qz(nz)
                            !cartisian coordinates of grid points
      real lambda(nz)       !latitude of given grid point
      integer ri, rj, rk    !release coordinates
      real dz_grid(nz)      !dz for grid, variable in z direction
      real dz_cell(nz)      !dz for cell, variable in z direction
 
      common /coords/ qx,qy,qz,lambda,ri,rj,rk,dz_grid,dz_cell
 
 
c Total number of ions produced at a given time during the simulation
      integer*4 Ni_tot, Nn_tot
 
c Location (indices) of particles with the grid
      integer ijkp(Ni_max,3)
      logical ionized(Ni_max)
      common /ion_info/ Ni_tot, Nn_tot, ijkp, ionized
 
      real seed
      common /rndseed/ seed
 
c Indices for boundary cells
      integer ip,im,jp,jm,kp,km
 
c Total input energy (ions)
      real input_E,input_EeP,prev_Etot, bndry_Eflux
      real input_chex, input_bill
      common /energy/ input_E, input_EeP, prev_Etot, bndry_Eflux,
     x                input_chex, input_bill
 
c Momentum conservation
      real surf_tot(3),graduu_tot(3),ugradu_tot(3)
      common /momentum/ surf_tot,graduu_tot,ugradu_tot
 
c Dummy variables for doing vector operations and particle flux
c calculation....etc.
 
      real a(nx,ny,nz,3),c(nx,ny,nz,3),ct(nx,ny,nz,3)
      common /dummy/ a,c,ct
 
c Variable time step for fluid velocity update used in both
c particle update and fluid update
 
      real dtf,ntf,dtsub
      common /timestep/ ntf,dtsub
 
c Weight variables for trilinear interpolation
 
      real wght(Ni_max,8), wquad(Ni_max,3)
      common /weights/ wght, wquad
 
c variable for anomlous resistivity
 
      real etar(nx,ny,nz,3)
      real nuei
      common /resist/ etar,nuei
 
c variable for particle scale
 
      real beta, dNi
      common /scale/ beta, dNi
 
 
c mass array for multi-ion species
      real mrat,m_arr
      common /mass/ mrat(Ni_max),m_arr(Ni_max)
 
 
 
      real Ef(nx,ny,nz,3),
     x     aj(nx,ny,nz,3),
     x     np(nx,ny,nz),
     x     nf(nx,ny,nz),
     x     up(nx,ny,nz,3),
     x     uf(nx,ny,nz,3),
     x     btmf(nx,ny,nz,3),
     x     nu(nx,ny,nz),
     x     ugradu(nx,ny,nz,3),
     x     gradPf(nx,ny,nz,3)
 
      real ntot(3)                !total plasma density
      real fnp(3)                 !fraction, np/n
      real aac(3),bbc(3),ccc(3)
      real cc(nx,ny,nz,3)
 
      real r1, r2, r3, r4, r5, r6
      call periodic_scalar(np)
      call periodic_scalar(nf)
 
      do 10 k = 2, 36
         do 10 j = 2, 18
            do 10 i = 2, 78
 
               ip = i + 1
               jp = j + 1
               kp = k + 1
 
c               if (ip .gt. nx) then ip = nx
c               if (jp .gt. ny) then jp = ny
c               if (kp .gt. nz) then kp = nz
 
               ntot(1) = 0.5*(nf(i,j,k)+nf(ip,j,k)) + 0.5*(np(i,j,k)+np(
     1            ip,j,k))
               ntot(2) = 0.5*(nf(i,j,k)+nf(i,jp,k)) + 0.5*(np(i,j,k)+np(
     1            i,jp,k))
               ntot(3) = 0.5*(nf(i,j,k)+nf(i,j,kp)) + 0.5*(np(i,j,k)+np(
     1            i,j,kp))
 
               fnp(1) = 0.5*(np(i,j,k)+np(ip,j,k))/ntot(1)
               fnp(2) = 0.5*(np(i,j,k)+np(i,jp,k))/ntot(2)
               fnp(3) = 0.5*(np(i,j,k)+np(i,j,kp))/ntot(3)
 
               r4 = aj(i,j,k,1) - fnp(1)*up(i,j,k,1)
               r5 = aj(i,j,k,2) - fnp(2)*up(i,j,k,2)
               r6 = aj(i,j,k,3) - fnp(3)*up(i,j,k,3)
               a(i,j,k,1) = r4
               a(i,j,k,2) = r5
               a(i,j,k,3) = r6
   10 continue
 
      call crossf(a,btmf,c)
 
      call get_ugradu_Lax(uf,ugradu,delta_t)
 
      do 20 k = 2, 36
         do 20 j = 2, 18
            do 20 i = 2, 78
 
               ip = i + 1
               jp = j + 1
               kp = k + 1
 
c               if (ip .gt. nx) then ip = nx
c               if (jp .gt. ny) then jp = ny
c               if (kp .gt. nz) then kp = nz
 
               ntot(1) = 0.5*(nf(i,j,k)+nf(ip,j,k)) + 0.5*(np(i,j,k)+np(
     1            ip,j,k))
               ntot(2) = 0.5*(nf(i,j,k)+nf(i,jp,k)) + 0.5*(np(i,j,k)+np(
     1            i,jp,k))
               ntot(3) = 0.5*(nf(i,j,k)+nf(i,j,kp)) + 0.5*(np(i,j,k)+np(
     1            i,j,kp))
 
               fnp(1) = 0.5*(np(i,j,k)+np(ip,j,k))/ntot(1)
               fnp(2) = 0.5*(np(i,j,k)+np(i,jp,k))/ntot(2)
               fnp(3) = 0.5*(np(i,j,k)+np(i,j,kp))/ntot(3)
 
               r1 = c(i,j,k,1) - ugradu(i,j,k,1) + nu(i,j,k)*fnp(1)*up(i
     1            ,j,k,1) + nuei*aj(i,j,k,1) - gradpf(i,j,k,1)
               r2 = c(i,j,k,2) - ugradu(i,j,k,2) + nu(i,j,k)*fnp(2)*up(i
     1            ,j,k,2) + nuei*aj(i,j,k,2) - gradpf(i,j,k,2)
               r3 = c(i,j,k,3) - ugradu(i,j,k,3) + nu(i,j,k)*fnp(3)*up(i
     1            ,j,k,3) + nuei*aj(i,j,k,3) - gradpf(i,j,k,3)
               ef(i,j,k,1) = r1
               ef(i,j,k,2) = r2
               ef(i,j,k,3) = r3
c     x                          + etar(i,j,k,m)*aj(i,j,k,m)
c                  Ef(i,j,k,m) = c(i,j,k,m) +
c     x                          nu(i,j,k)*fnp(m)*up(i,j,k,m)
   20 continue
 
      call periodic(Ef)
c      call fix_tangential_E(Ef)
 
      return
      end
c----------------------------------------------------------------------
 
 
 
c----------------------------------------------------------------------
      SUBROUTINE get_uplus_uminus(Ef,btmf,uf,nu,np,nf,uplus,uminus,
     x                            delta_t)
c This is the heart of the fluid velocity update.  It solves eqn. 18
c (Dan's paper) for uf+
c----------------------------------------------------------------------
      PARAMETER (nx = 79, ny = 19, nz = 37)
 
c grid parameters
      PARAMETER (dx = 10000.0, dy = 10000.0)   !units in km
      PARAMETER (delz = 10000.0)          !dz at release coordinates
                                       !this is not the final dz
                                       !d_lamba determines dz
c time stepping parameters
      PARAMETER (tstart = 0.2)   !starting time
      PARAMETER (dt = 1.)     !main time step
      PARAMETER (nt = 5.0)        !number of time steps
c      PARAMETER (dtsub = 0.01) !subcycle time step
      PARAMETER (dtsub_init = 0.1) !subcycle time step
      PARAMETER (ntsub = 10.0)   !number of subcycle time steps
      PARAMETER (nout = 50)      !number of steps to diagnosic output
 
c logical variable for restart
C...Translated by Crescent Bay Software VAST-F    7.0K8 15:17:00   4/ 1/04    
      logical restart
      PARAMETER (restart = .false.)
      PARAMETER (mrestart = 500)      ! use -1 for no save
      PARAMETER (mbegin = 0)      !mstart
 
c inital ambient density at the release point.
      real nf_init,b0_init,nn_coef
      PARAMETER (nf_init = 4e15)   !was 0.01
      PARAMETER (b0_init = 2e-9)    !was 0.2
      PARAMETER (nn_coef = 1e5)
 
c electron ion collision frequency
      real nu_init, eta_init
      PARAMETER (nu_init = 0.001)
      PARAMETER (eta_init = 0.0)
 
c temporary directory for output
      character*9 temp_dir
      PARAMETER (temp_dir = '/tmp/pad/')
 
c total number of neutrals released
      real No
      PARAMETER (No=1.0e28)
 
c neutral cloud expansion characteristics
      PARAMETER(vo = 20.0)
      PARAMETER(vth = 19.3)
      PARAMETER(vsat = 0.0)
      PARAMETER(vsw = 340.0)
 
c total kinetic energy of released neutrals
      PARAMETER(Etot = 2.696e7)       !units of kg m^2 / s^2
 
c density scaling parameter, alpha, and ion particle array dims
 
      real alpha
 
c      PARAMETER (beta =  8.5227271e-18)
      PARAMETER (alpha = 1.9263418e-20) !mH
 
c max number of ion particles to be produced.  This parameter
c is used to dimension the particle arrays.
      integer*4 Ni_max
      PARAMETER (Ni_max = 10000)
 
c earth's magnetic dipole field geometry
c earth radius and equatorial B alt (defining L shell)
      PARAMETER (r_earth=6.38e3)                !earth radius
      PARAMETER (B_max_alt=2.0e3+r_earth)       !ro defining L shell
 
c misc constants
      real mu0,epsilon,pi,rtod,q,mO,mBa,O_to_Ba,km_to_m,kboltz,melec
      real mproton,tempf0
      PARAMETER (pi = 3.14159)
      PARAMETER (rtod = 180.0/pi)  !radians to degreesc
      PARAMETER (mu0 = pi*4.0e-7)  !magnetic permeability of free space
      PARAMETER (epsilon = 8.85e-12) !dielectric constant
      PARAMETER (q = 1.6e-19)        !electron charge (Coulombs)
      PARAMETER (melec = 9.1e-31)  !mass of electron (kg)
      PARAMETER (mproton = 1.67e-27)
      PARAMETER (mO = 1.67e-27)    !mass of H (kg)
c      PARAMETER (mO = 2.3e-25)    !mass of Ba (kg)
      PARAMETER (mBa = 17.0*mO)    !mass of Ba (kg)
      PARAMETER (O_to_Ba = mO/mBa) !convert E and B for particle move
      PARAMETER (km_to_m = 1e3)    !km to meters conversion
      PARAMETER (kboltz = 1.38e-29)   !kg km^2 / s^2 / K
      PARAMETER (tempf0 = 50*11600.)     !K
 
      real Qo, vrad, N_o, Rp, tau_photo
      integer S_radius
      PARAMETER (Qo = 1.0e28)  !neutral source rate
      PARAMETER (vrad = 1.0)   !escape velocity
c      PARAMETER (N_o = 5e34)   !Steady state neutral particle constant
      PARAMETER (N_o = 1e33)   !Steady state neutral particle constant
      PARAMETER (Rp = 20000.0)  !Pluto radius
      PARAMETER (tau_photo = 1.2e6)
c      PARAMETER (dNi = 2500)
      PARAMETER (S_radius = 3) !units of dx
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
c raw grid coordinate data
 
      real qx(nx), qy(ny), qz(nz)
                            !cartisian coordinates of grid points
      real lambda(nz)       !latitude of given grid point
      integer ri, rj, rk    !release coordinates
      real dz_grid(nz)      !dz for grid, variable in z direction
      real dz_cell(nz)      !dz for cell, variable in z direction
 
      common /coords/ qx,qy,qz,lambda,ri,rj,rk,dz_grid,dz_cell
 
 
c Total number of ions produced at a given time during the simulation
      integer*4 Ni_tot, Nn_tot
 
c Location (indices) of particles with the grid
      integer ijkp(Ni_max,3)
      logical ionized(Ni_max)
      common /ion_info/ Ni_tot, Nn_tot, ijkp, ionized
 
      real seed
      common /rndseed/ seed
 
c Indices for boundary cells
      integer ip,im,jp,jm,kp,km
 
c Total input energy (ions)
      real input_E,input_EeP,prev_Etot, bndry_Eflux
      real input_chex, input_bill
      common /energy/ input_E, input_EeP, prev_Etot, bndry_Eflux,
     x                input_chex, input_bill
 
c Momentum conservation
      real surf_tot(3),graduu_tot(3),ugradu_tot(3)
      common /momentum/ surf_tot,graduu_tot,ugradu_tot
 
c Dummy variables for doing vector operations and particle flux
c calculation....etc.
 
      real a(nx,ny,nz,3),c(nx,ny,nz,3),ct(nx,ny,nz,3)
      common /dummy/ a,c,ct
 
c Variable time step for fluid velocity update used in both
c particle update and fluid update
 
      real dtf,ntf,dtsub
      common /timestep/ ntf,dtsub
 
c Weight variables for trilinear interpolation
 
      real wght(Ni_max,8), wquad(Ni_max,3)
      common /weights/ wght, wquad
 
c variable for anomlous resistivity
 
      real etar(nx,ny,nz,3)
      real nuei
      common /resist/ etar,nuei
 
c variable for particle scale
 
      real beta, dNi
      common /scale/ beta, dNi
 
 
c mass array for multi-ion species
      real mrat,m_arr
      common /mass/ mrat(Ni_max),m_arr(Ni_max)
 
 
 
      real Ef(nx,ny,nz,3),
     x     btmf(nx,ny,nz,3),
     x     uf(nx,ny,nz,3), !using particle velocities at t level n-1/2
     x     nu(nx,ny,nz),
     x     np(nx,ny,nz),
     x     nf(nx,ny,nz),
     x     uplus(nx,ny,nz,3),
     x     uminus(nx,ny,nz,3)
 
      real a1,b1,b2,b3      !4 coefficients (see notebook solution)
      real PP,QQ            !intermediate variables for a1,b1,b2,b3
      real eta              !intermediate variable for PP,QQ
c      real B(3)             !B for cross product call
c      real Bsqrd                     !B*B
      real um_x_B(nx,ny,nz,3)        !uf- X B
      real um_dot_BB(nx,ny,nz,3)     !uf- . B
      real ntot(3)                   !total density np + nf
      real npave(3)
      real btc(nx,ny,nz,3)
      real bsqrd(nx,ny,nz),bsq(3)
 
      external vp_get_uplus_uminus1
      external vp_get_uplus_uminus2
      integer j3, j4, j5, j6, j7, j8
      real r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, 
     .r15, r16
      call vp_pcall(vp_get_uplus_uminus1,2,4,delta_t,uf,ef,uminus)
 
      call crossf(uminus, btmf, um_x_B)
      call get_um_dot_BB(uminus , btmf, um_dot_BB)
 
      call face_to_center(btmf,btc)
 
      call vp_pcall (vp_get_uplus_uminus2, 2, 2, btc, bsqrd)
 
      call periodic_scalar(bsqrd)
 
      do 20 k = 2, 36
         do 20 j = 2, 18
            do 20 i = 2, 78
 
               ip = i + 1
               jp = j + 1
               kp = k + 1
 
c               if (ip .gt. nx) then ip = nx
c               if (jp .gt. ny) then jp = ny
c               if (kp .gt. nz) then kp = nz
 
               bsq(1) = 0.5*(bsqrd(i,j,k)+bsqrd(ip,j,k))
               bsq(2) = 0.5*(bsqrd(i,j,k)+bsqrd(i,jp,k))
               bsq(3) = 0.5*(bsqrd(i,j,k)+bsqrd(i,j,kp))
 
               npave(1) = 0.5*(np(i,j,k)+np(ip,j,k))
               npave(2) = 0.5*(np(i,j,k)+np(i,jp,k))
               npave(3) = 0.5*(np(i,j,k)+np(i,j,kp))
 
               ntot(1) = 0.5*(nf(i,j,k)+nf(ip,j,k)) + npave(1)
               ntot(2) = 0.5*(nf(i,j,k)+nf(i,jp,k)) + npave(2)
               ntot(3) = 0.5*(nf(i,j,k)+nf(i,j,kp)) + npave(3)
 
               eta = npave(1)*delta_t/(2.0*ntot(1))
               qq = eta/(1.0 + (nu(i,j,k)*eta))
               pp = (1.0 - (nu(i,j,k)*eta))/(1.0 + (nu(i,j,k)*eta))
               a1 = 1.0/(1.0 + ((qq*qq)*bsq(1)))
               b1 = pp - ((qq*qq)*bsq(1))
               b3 = (qq*qq)*pp + (qq*qq)
               r14 = a1*(b1*uminus(i,j,k,1)+(qq*pp+qq)*um_x_b(i,j,k,1)+
     1            b3*um_dot_bb(i,j,k,1))
               eta = npave(2)*delta_t/(2.0*ntot(2))
               qq = eta/(1.0 + (nu(i,j,k)*eta))
               pp = (1.0 - (nu(i,j,k)*eta))/(1.0 + (nu(i,j,k)*eta))
               a1 = 1.0/(1.0 + ((qq*qq)*bsq(2)))
               b1 = pp - ((qq*qq)*bsq(2))
               b3 = (qq*qq)*pp + (qq*qq)
               r15 = a1*(b1*uminus(i,j,k,2)+(qq*pp+qq)*um_x_b(i,j,k,2)+
     1            b3*um_dot_bb(i,j,k,2))
               eta = npave(3)*delta_t/(2.0*ntot(3))
               qq = eta/(1.0 + (nu(i,j,k)*eta))
               pp = (1.0 - (nu(i,j,k)*eta))/(1.0 + (nu(i,j,k)*eta))
               a1 = 1.0/(1.0 + ((qq*qq)*bsq(3)))
               b1 = pp - ((qq*qq)*bsq(3))
               b3 = (qq*qq)*pp + (qq*qq)
               r16 = a1*(b1*uminus(i,j,k,3)+(qq*pp+qq)*um_x_b(i,j,k,3)+
     1            b3*um_dot_bb(i,j,k,3))
               uplus(i,j,k,1) = r14
               uplus(i,j,k,2) = r15
               uplus(i,j,k,3) = r16
 
   20 continue
 
      return
      end
      subroutine vp_get_uplus_uminus1(vp_mythread, vp_nmthreads, delta_t
     1   , uf, ef, uminus)
      REAL delta_t
      REAL uf(79,19,37,3), ef(79,19,37,3), uminus(79,19,37,3)
      INTEGER j, k, m, i, j3, j4, j5, vp_mythread, vp_nmthreads, j1, j2
      j3 = (17 + vp_nmthreads - 1)/vp_nmthreads
      j4 = vp_mythread*j3 + 1
      j5 = min(17 - j4 + 1,j3)
      do j = j4, j5 + j4 - 1
         do k = 1, 35
            do i = 1, 77
               uminus(1+i,j+1,k+1,1) = uf(1+i,j+1,k+1,1) + (0.5*delta_t)
     1            *ef(1+i,j+1,k+1,1)
               uminus(1+i,j+1,k+1,2) = uf(1+i,j+1,k+1,2) + (0.5*delta_t)
     1            *ef(1+i,j+1,k+1,2)
               uminus(1+i,j+1,k+1,3) = uf(1+i,j+1,k+1,3) + (0.5*delta_t)
     1            *ef(1+i,j+1,k+1,3)
            end do
         end do
      end do
      return 
      end 
      subroutine vp_get_uplus_uminus2(vp_mythread, vp_nmthreads, btc, 
     1   bsqrd)
      REAL btc(79,19,37,3), bsqrd(79,19,37)
      REAL r1, r2, r3, r4
      INTEGER k, j, i, j6, j7, j8, vp_mythread, vp_nmthreads, j1, j2
 
      j6 = (35 + vp_nmthreads - 1)/vp_nmthreads
      j7 = vp_mythread*j6 + 1
      j8 = min(35 - j7 + 1,j6)
      do k = j7, j8 + j7 - 1
         j = 1
         do i = 1, 77
            bsqrd(1+i,2,k+1) = btc(1+i,2,k+1,1)**2 + btc(1+i,2,k+1,2)**2
     1          + btc(1+i,2,k+1,3)**2
         end do
         do j = 1, 4
            do i = 1, 77
               r1 = btc(1+i,j*4-1,k+1,1)**2 + btc(1+i,j*4-1,k+1,2)**2 + 
     1            btc(1+i,j*4-1,k+1,3)**2
               r2 = btc(1+i,j*4,k+1,1)**2 + btc(1+i,j*4,k+1,2)**2 + btc(
     1            1+i,j*4,k+1,3)**2
               r3 = btc(1+i,j*4+1,k+1,1)**2 + btc(1+i,j*4+1,k+1,2)**2 + 
     1            btc(1+i,j*4+1,k+1,3)**2
               r4 = btc(1+i,j*4+2,k+1,1)**2 + btc(1+i,j*4+2,k+1,2)**2 + 
     1            btc(1+i,j*4+2,k+1,3)**2
               bsqrd(1+i,j*4-1,k+1) = r1
               bsqrd(1+i,j*4,k+1) = r2
               bsqrd(1+i,j*4+1,k+1) = r3
               bsqrd(1+i,j*4+2,k+1) = r4
            end do
         end do
      end do
      return 
      end 
c----------------------------------------------------------------------
 
 
c----------------------------------------------------------------------
      SUBROUTINE predict_uf(Ef,b0,b1,b12,uf,uf2,ufp2,nu,np,nf,uplus,
     x                      uminus,ugradu,up,gradP,nuin,bdp,pf1)
c Calculate the fluid velocity, uf,  at the new time step and replace
c uf1 with the new value, uf, in preparation for the next time step.
c----------------------------------------------------------------------
      PARAMETER (nx = 79, ny = 19, nz = 37)
 
c grid parameters
      PARAMETER (dx = 10000.0, dy = 10000.0)   !units in km
      PARAMETER (delz = 10000.0)          !dz at release coordinates
                                       !this is not the final dz
                                       !d_lamba determines dz
c time stepping parameters
      PARAMETER (tstart = 0.2)   !starting time
      PARAMETER (dt = 1.)     !main time step
      PARAMETER (nt = 5.0)        !number of time steps
c      PARAMETER (dtsub = 0.01) !subcycle time step
      PARAMETER (dtsub_init = 0.1) !subcycle time step
      PARAMETER (ntsub = 10.0)   !number of subcycle time steps
      PARAMETER (nout = 50)      !number of steps to diagnosic output
 
c logical variable for restart
C...Translated by Crescent Bay Software VAST-F    7.0K8 15:17:00   4/ 1/04    
      integer j1, j2, j3
      logical restart
      PARAMETER (restart = .false.)
      PARAMETER (mrestart = 500)      ! use -1 for no save
      PARAMETER (mbegin = 0)      !mstart
 
c inital ambient density at the release point.
      real nf_init,b0_init,nn_coef
      PARAMETER (nf_init = 4e15)   !was 0.01
      PARAMETER (b0_init = 2e-9)    !was 0.2
      PARAMETER (nn_coef = 1e5)
 
c electron ion collision frequency
      real nu_init, eta_init
      PARAMETER (nu_init = 0.001)
      PARAMETER (eta_init = 0.0)
 
c temporary directory for output
      character*9 temp_dir
      PARAMETER (temp_dir = '/tmp/pad/')
 
c total number of neutrals released
      real No
      PARAMETER (No=1.0e28)
 
c neutral cloud expansion characteristics
      PARAMETER(vo = 20.0)
      PARAMETER(vth = 19.3)
      PARAMETER(vsat = 0.0)
      PARAMETER(vsw = 340.0)
 
c total kinetic energy of released neutrals
      PARAMETER(Etot = 2.696e7)       !units of kg m^2 / s^2
 
c density scaling parameter, alpha, and ion particle array dims
 
      real alpha
 
c      PARAMETER (beta =  8.5227271e-18)
      PARAMETER (alpha = 1.9263418e-20) !mH
 
c max number of ion particles to be produced.  This parameter
c is used to dimension the particle arrays.
      integer*4 Ni_max
      PARAMETER (Ni_max = 10000)
 
c earth's magnetic dipole field geometry
c earth radius and equatorial B alt (defining L shell)
      PARAMETER (r_earth=6.38e3)                !earth radius
      PARAMETER (B_max_alt=2.0e3+r_earth)       !ro defining L shell
 
c misc constants
      real mu0,epsilon,pi,rtod,q,mO,mBa,O_to_Ba,km_to_m,kboltz,melec
      real mproton,tempf0
      PARAMETER (pi = 3.14159)
      PARAMETER (rtod = 180.0/pi)  !radians to degreesc
      PARAMETER (mu0 = pi*4.0e-7)  !magnetic permeability of free space
      PARAMETER (epsilon = 8.85e-12) !dielectric constant
      PARAMETER (q = 1.6e-19)        !electron charge (Coulombs)
      PARAMETER (melec = 9.1e-31)  !mass of electron (kg)
      PARAMETER (mproton = 1.67e-27)
      PARAMETER (mO = 1.67e-27)    !mass of H (kg)
c      PARAMETER (mO = 2.3e-25)    !mass of Ba (kg)
      PARAMETER (mBa = 17.0*mO)    !mass of Ba (kg)
      PARAMETER (O_to_Ba = mO/mBa) !convert E and B for particle move
      PARAMETER (km_to_m = 1e3)    !km to meters conversion
      PARAMETER (kboltz = 1.38e-29)   !kg km^2 / s^2 / K
      PARAMETER (tempf0 = 50*11600.)     !K
 
      real Qo, vrad, N_o, Rp, tau_photo
      integer S_radius
      PARAMETER (Qo = 1.0e28)  !neutral source rate
      PARAMETER (vrad = 1.0)   !escape velocity
c      PARAMETER (N_o = 5e34)   !Steady state neutral particle constant
      PARAMETER (N_o = 1e33)   !Steady state neutral particle constant
      PARAMETER (Rp = 20000.0)  !Pluto radius
      PARAMETER (tau_photo = 1.2e6)
c      PARAMETER (dNi = 2500)
      PARAMETER (S_radius = 3) !units of dx
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
c raw grid coordinate data
 
      real qx(nx), qy(ny), qz(nz)
                            !cartisian coordinates of grid points
      real lambda(nz)       !latitude of given grid point
      integer ri, rj, rk    !release coordinates
      real dz_grid(nz)      !dz for grid, variable in z direction
      real dz_cell(nz)      !dz for cell, variable in z direction
 
      common /coords/ qx,qy,qz,lambda,ri,rj,rk,dz_grid,dz_cell
 
 
c Total number of ions produced at a given time during the simulation
      integer*4 Ni_tot, Nn_tot
 
c Location (indices) of particles with the grid
      integer ijkp(Ni_max,3)
      logical ionized(Ni_max)
      common /ion_info/ Ni_tot, Nn_tot, ijkp, ionized
 
      real seed
      common /rndseed/ seed
 
c Indices for boundary cells
      integer ip,im,jp,jm,kp,km
 
c Total input energy (ions)
      real input_E,input_EeP,prev_Etot, bndry_Eflux
      real input_chex, input_bill
      common /energy/ input_E, input_EeP, prev_Etot, bndry_Eflux,
     x                input_chex, input_bill
 
c Momentum conservation
      real surf_tot(3),graduu_tot(3),ugradu_tot(3)
      common /momentum/ surf_tot,graduu_tot,ugradu_tot
 
c Dummy variables for doing vector operations and particle flux
c calculation....etc.
 
      real a(nx,ny,nz,3),c(nx,ny,nz,3),ct(nx,ny,nz,3)
      common /dummy/ a,c,ct
 
c Variable time step for fluid velocity update used in both
c particle update and fluid update
 
      real dtf,ntf,dtsub
      common /timestep/ ntf,dtsub
 
c Weight variables for trilinear interpolation
 
      real wght(Ni_max,8), wquad(Ni_max,3)
      common /weights/ wght, wquad
 
c variable for anomlous resistivity
 
      real etar(nx,ny,nz,3)
      real nuei
      common /resist/ etar,nuei
 
c variable for particle scale
 
      real beta, dNi
      common /scale/ beta, dNi
 
 
c mass array for multi-ion species
      real mrat,m_arr
      common /mass/ mrat(Ni_max),m_arr(Ni_max)
 
 
 
      real Ef(nx,ny,nz,3),
     x     b0(ny),
     x     b1(nx,ny,nz,3),
     x     b12(nx,ny,nz,3),
     x     uf(nx,ny,nz,3),
     x     uf2(nx,ny,nz,3),
     x     ufp2(nx,ny,nz,3),
     x     nu(nx,ny,nz),
     x     np(nx,ny,nz),
     x     nf(nx,ny,nz),
     x     uplus(nx,ny,nz,3),
     x     uminus(nx,ny,nz,3),
     x     ugradu(nx,ny,nz,3),
     x     up(nx,ny,nz,3),
c     x     gradP(nx,ny,nz,3),
     x     nuin(nx,ny,nz),
     x     bdp(nx,ny,nz,3),
     x     pf1(nx,ny,nz)
 
      real b1h(nx,ny,nz,3)
      real bth(nx,ny,nz,3)
      real btmfh(nx,ny,nz,3)
      real ajh(nx,ny,nz,3)
      real gradPf(nx,ny,nz,3)
 
      real m_den
      real delta_t
 
      external vp_predict_uf1
      external vp_predict_uf2
      external vp_predict_uf3
      integer j6, j7, j8, j9, j10, j11, j12, j13, j14, j15, j16, j17
      real r1
      delta_t = 2.0*dtsub
 
      call vp_pcall (vp_predict_uf1, 2, 9, b1, b12, b1h, bdp, bth, b0, 
     1   nf, pf1, gradpf)
 
      call cov_to_contra(bth,btmfh)
      call curlB(b1h,nf,np,ajh)
 
      call get_Ef(Ef,ajh,np,nf,up,uf,btmfh,nu,ugradu,delta_t,gradPf)
      call get_uplus_uminus(Ef,btmfh,uf2,nu,np,nf,uplus,uminus,
     x                      delta_t)
 
      call vp_pcall(vp_predict_uf2,2,5,delta_t,uplus,ef,nuin,ufp2)
 
      call vp_pcall (vp_predict_uf3, 2, 1, ufp2)
 
      return
      end
      subroutine vp_predict_uf1(vp_mythread, vp_nmthreads, b1, b12, b1h
     1   , bdp, bth, b0, nf, pf1, gradpf)
      common /coords/qx, qy, qz, lambda, ri, rj, rk, dz_grid, dz_cell
      INTEGER ri, rj, rk
      REAL b1(79,19,37,3), b12(79,19,37,3), b1h(79,19,37,3), bdp(79,19,
     1   37,3), bth(79,19,37,3), b0(19), nf(79,19,37), pf1(79,19,37), 
     2   gradpf(79,19,37,3), dz_grid(37), qx(79), qy(19), qz(37), lambda
     3   (37), dz_cell(37)
      REAL m_den
      INTEGER k, j, i, j6, j7, j8, j9, j10, j11, vp_mythread, 
     1   vp_nmthreads, j4, j5
 
      j6 = (35 + vp_nmthreads - 1)/vp_nmthreads
      j7 = vp_mythread*j6 + 1
      j8 = min(35 - j7 + 1,j6)
      do k = j7, j8 + j7 - 1
         do j = 1, 17
            do i = 1, 77
               m_den = (1.67E-27*0.5)*(nf(2+i,j+1,k+1)+nf(1+i,j+1,k+1))
               gradpf(1+i,j+1,k+1,1) = (pf1(2+i,j+1,k+1)-pf1(1+i,j+1,k+1
     1            ))/(10000.0*m_den)
               m_den = (1.67E-27*0.5)*(nf(1+i,j+2,k+1)+nf(1+i,j+1,k+1))
               gradpf(1+i,j+1,k+1,2) = (pf1(1+i,j+2,k+1)-pf1(1+i,j+1,k+1
     1            ))/(10000.0*m_den)
               m_den = (1.67E-27*0.5)*(nf(1+i,j+1,k+2)+nf(1+i,j+1,k+1))
               gradpf(1+i,j+1,k+1,3) = (pf1(1+i,j+1,k+2)-pf1(1+i,j+1,k+1
     1            ))/(dz_grid(k+1)*m_den)
            end do
         end do
      end do
 
 
      j9 = (37 + vp_nmthreads - 1)/vp_nmthreads
      j10 = vp_mythread*j9 + 1
      j11 = min(37 - j10 + 1,j9)
      do k = j10, j11 + j10 - 1
         do j = 1, 19
            do i = 1, 79
               b1h(i,j,k,1) = 0.5*(b1(i,j,k,1)+b12(i,j,k,1))
               b1h(i,j,k,2) = 0.5*(b1(i,j,k,2)+b12(i,j,k,2))
               b1h(i,j,k,3) = 0.5*(b1(i,j,k,3)+b12(i,j,k,3))
               bth(i,j,k,1) = b1h(i,j,k,1) + bdp(i,j,k,1)
               bth(i,j,k,2) = b1h(i,j,k,2) + b0(j) + bdp(i,j,k,2)
               bth(i,j,k,3) = b1h(i,j,k,3) + bdp(i,j,k,3)
            end do
         end do
      end do
      return 
      end 
      subroutine vp_predict_uf2(vp_mythread, vp_nmthreads, delta_t, 
     1   uplus, ef, nuin, ufp2)
      REAL delta_t
      REAL uplus(79,19,37,3), ef(79,19,37,3), nuin(79,19,37), ufp2(79,19
     1   ,37,3)
      REAL r1
      INTEGER k,j,i,m,j12,j13,j14,vp_mythread,vp_nmthreads,j4,j5
 
      j12 = (35 + vp_nmthreads - 1)/vp_nmthreads
      j13 = vp_mythread*j12 + 1
      j14 = min(35 - j13 + 1,j12)
      do k = j13, j14 + j13 - 1
         do j = 1, 17
            do i = 1, 77
               r1 = nuin(1+i,j+1,k+1)
               ufp2(1+i,j+1,k+1,1) = uplus(1+i,j+1,k+1,1) + (0.5*delta_t
     1            )*ef(1+i,j+1,k+1,1) - ((0.5*delta_t)*r1)*uplus(1+i,j+1
     2            ,k+1,1)
               ufp2(1+i,j+1,k+1,2) = uplus(1+i,j+1,k+1,2) + (0.5*delta_t
     1            )*ef(1+i,j+1,k+1,2) - ((0.5*delta_t)*r1)*uplus(1+i,j+1
     2            ,k+1,2)
               ufp2(1+i,j+1,k+1,3) = uplus(1+i,j+1,k+1,3) + (0.5*delta_t
     1            )*ef(1+i,j+1,k+1,3) - ((0.5*delta_t)*r1)*uplus(1+i,j+1
     2            ,k+1,3)
            end do
         end do
      end do
      return 
      end 
      subroutine vp_predict_uf3(vp_mythread, vp_nmthreads, ufp2)
      REAL ufp2(79,19,37,3)
      INTEGER j3, j2, j15, j16, j17, vp_mythread, vp_nmthreads, j4, j5
 
      j15 = (703 + vp_nmthreads - 1)/vp_nmthreads
      j16 = vp_mythread*j15 + 1
      j17 = min(703 - j16 + 1,j15)
      do j3 = j16, j17 + j16 - 1
         ufp2(78,j3,1,1) = ((-340.0))
         ufp2(79,j3,1,1) = ((-340.0))
      end do
      return 
      end 
c----------------------------------------------------------------------
 
 
c----------------------------------------------------------------------
      SUBROUTINE correct_uf(Ef,btmf,uf,uf2,ufp2,nu,np,nf,uplus,uminus,
     x                      ugradu,aj,up,ufp1,gradP,nuin,pf)
c Calculate the fluid velocity, uf,  at the new time step and replace
c uf1 with the new value, uf, in preparation for the next time step.
c----------------------------------------------------------------------
      PARAMETER (nx = 79, ny = 19, nz = 37)
 
c grid parameters
      PARAMETER (dx = 10000.0, dy = 10000.0)   !units in km
      PARAMETER (delz = 10000.0)          !dz at release coordinates
                                       !this is not the final dz
                                       !d_lamba determines dz
c time stepping parameters
      PARAMETER (tstart = 0.2)   !starting time
      PARAMETER (dt = 1.)     !main time step
      PARAMETER (nt = 5.0)        !number of time steps
c      PARAMETER (dtsub = 0.01) !subcycle time step
      PARAMETER (dtsub_init = 0.1) !subcycle time step
      PARAMETER (ntsub = 10.0)   !number of subcycle time steps
      PARAMETER (nout = 50)      !number of steps to diagnosic output
 
c logical variable for restart
C...Translated by Crescent Bay Software VAST-F    7.0K8 15:17:00   4/ 1/04    
      integer j1, j2, j3
      logical restart
      PARAMETER (restart = .false.)
      PARAMETER (mrestart = 500)      ! use -1 for no save
      PARAMETER (mbegin = 0)      !mstart
 
c inital ambient density at the release point.
      real nf_init,b0_init,nn_coef
      PARAMETER (nf_init = 4e15)   !was 0.01
      PARAMETER (b0_init = 2e-9)    !was 0.2
      PARAMETER (nn_coef = 1e5)
 
c electron ion collision frequency
      real nu_init, eta_init
      PARAMETER (nu_init = 0.001)
      PARAMETER (eta_init = 0.0)
 
c temporary directory for output
      character*9 temp_dir
      PARAMETER (temp_dir = '/tmp/pad/')
 
c total number of neutrals released
      real No
      PARAMETER (No=1.0e28)
 
c neutral cloud expansion characteristics
      PARAMETER(vo = 20.0)
      PARAMETER(vth = 19.3)
      PARAMETER(vsat = 0.0)
      PARAMETER(vsw = 340.0)
 
c total kinetic energy of released neutrals
      PARAMETER(Etot = 2.696e7)       !units of kg m^2 / s^2
 
c density scaling parameter, alpha, and ion particle array dims
 
      real alpha
 
c      PARAMETER (beta =  8.5227271e-18)
      PARAMETER (alpha = 1.9263418e-20) !mH
 
c max number of ion particles to be produced.  This parameter
c is used to dimension the particle arrays.
      integer*4 Ni_max
      PARAMETER (Ni_max = 10000)
 
c earth's magnetic dipole field geometry
c earth radius and equatorial B alt (defining L shell)
      PARAMETER (r_earth=6.38e3)                !earth radius
      PARAMETER (B_max_alt=2.0e3+r_earth)       !ro defining L shell
 
c misc constants
      real mu0,epsilon,pi,rtod,q,mO,mBa,O_to_Ba,km_to_m,kboltz,melec
      real mproton,tempf0
      PARAMETER (pi = 3.14159)
      PARAMETER (rtod = 180.0/pi)  !radians to degreesc
      PARAMETER (mu0 = pi*4.0e-7)  !magnetic permeability of free space
      PARAMETER (epsilon = 8.85e-12) !dielectric constant
      PARAMETER (q = 1.6e-19)        !electron charge (Coulombs)
      PARAMETER (melec = 9.1e-31)  !mass of electron (kg)
      PARAMETER (mproton = 1.67e-27)
      PARAMETER (mO = 1.67e-27)    !mass of H (kg)
c      PARAMETER (mO = 2.3e-25)    !mass of Ba (kg)
      PARAMETER (mBa = 17.0*mO)    !mass of Ba (kg)
      PARAMETER (O_to_Ba = mO/mBa) !convert E and B for particle move
      PARAMETER (km_to_m = 1e3)    !km to meters conversion
      PARAMETER (kboltz = 1.38e-29)   !kg km^2 / s^2 / K
      PARAMETER (tempf0 = 50*11600.)     !K
 
      real Qo, vrad, N_o, Rp, tau_photo
      integer S_radius
      PARAMETER (Qo = 1.0e28)  !neutral source rate
      PARAMETER (vrad = 1.0)   !escape velocity
c      PARAMETER (N_o = 5e34)   !Steady state neutral particle constant
      PARAMETER (N_o = 1e33)   !Steady state neutral particle constant
      PARAMETER (Rp = 20000.0)  !Pluto radius
      PARAMETER (tau_photo = 1.2e6)
c      PARAMETER (dNi = 2500)
      PARAMETER (S_radius = 3) !units of dx
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
c raw grid coordinate data
 
      real qx(nx), qy(ny), qz(nz)
                            !cartisian coordinates of grid points
      real lambda(nz)       !latitude of given grid point
      integer ri, rj, rk    !release coordinates
      real dz_grid(nz)      !dz for grid, variable in z direction
      real dz_cell(nz)      !dz for cell, variable in z direction
 
      common /coords/ qx,qy,qz,lambda,ri,rj,rk,dz_grid,dz_cell
 
 
c Total number of ions produced at a given time during the simulation
      integer*4 Ni_tot, Nn_tot
 
c Location (indices) of particles with the grid
      integer ijkp(Ni_max,3)
      logical ionized(Ni_max)
      common /ion_info/ Ni_tot, Nn_tot, ijkp, ionized
 
      real seed
      common /rndseed/ seed
 
c Indices for boundary cells
      integer ip,im,jp,jm,kp,km
 
c Total input energy (ions)
      real input_E,input_EeP,prev_Etot, bndry_Eflux
      real input_chex, input_bill
      common /energy/ input_E, input_EeP, prev_Etot, bndry_Eflux,
     x                input_chex, input_bill
 
c Momentum conservation
      real surf_tot(3),graduu_tot(3),ugradu_tot(3)
      common /momentum/ surf_tot,graduu_tot,ugradu_tot
 
c Dummy variables for doing vector operations and particle flux
c calculation....etc.
 
      real a(nx,ny,nz,3),c(nx,ny,nz,3),ct(nx,ny,nz,3)
      common /dummy/ a,c,ct
 
c Variable time step for fluid velocity update used in both
c particle update and fluid update
 
      real dtf,ntf,dtsub
      common /timestep/ ntf,dtsub
 
c Weight variables for trilinear interpolation
 
      real wght(Ni_max,8), wquad(Ni_max,3)
      common /weights/ wght, wquad
 
c variable for anomlous resistivity
 
      real etar(nx,ny,nz,3)
      real nuei
      common /resist/ etar,nuei
 
c variable for particle scale
 
      real beta, dNi
      common /scale/ beta, dNi
 
 
c mass array for multi-ion species
      real mrat,m_arr
      common /mass/ mrat(Ni_max),m_arr(Ni_max)
 
 
 
      real Ef(nx,ny,nz,3),
     x     btmf(nx,ny,nz,3),
     x     uf(nx,ny,nz,3),
     x     uf2(nx,ny,nz,3),
     x     ufp2(nx,ny,nz,3),
     x     nu(nx,ny,nz),
     x     np(nx,ny,nz),
     x     nf(nx,ny,nz),
     x     uplus(nx,ny,nz,3),
     x     uminus(nx,ny,nz,3),
     x     ugradu(nx,ny,nz,3),
     x     aj(nx,ny,nz,3),
     x     up(nx,ny,nz,3),
     x     ufp1(nx,ny,nz,3),
c     x     gradP(nx,ny,nz,3),
     x     nuin(nx,ny,nz),
     x     pf(nx,ny,nz)
 
      real m_den
      real delta_t
      real gradPf(nx,ny,nz,3)
 
      external vp_correct_uf1
      external vp_correct_uf2
      external vp_correct_uf3
      integer j6, j7, j8, j9, j10, j11, j12, j13, j14, j15, j16, j17
      real r1
      delta_t = dtsub
 
      call vp_pcall(vp_correct_uf1,2,6,uf,ufp2,ufp1,nf,pf,gradpf)
 
      call get_Ef(Ef,aj,np,nf,up,ufp1,btmf,nu,ugradu,delta_t,gradPf)
      call get_uplus_uminus(Ef,btmf,uf,nu,np,nf,uplus,uminus,
     x                      delta_t)
 
      call vp_pcall (vp_correct_uf2, 2, 5, uf, uf2, uplus, ef, nuin)
 
      call vp_pcall (vp_correct_uf3, 2, 1, uf)
 
      return
      end
      subroutine vp_correct_uf1(vp_mythread, vp_nmthreads, uf, ufp2, 
     1   ufp1, nf, pf, gradpf)
      common /coords/qx, qy, qz, lambda, ri, rj, rk, dz_grid, dz_cell
      INTEGER ri, rj, rk
      REAL uf(79,19,37,3), ufp2(79,19,37,3), ufp1(79,19,37,3), nf(79,19,
     1   37), pf(79,19,37), gradpf(79,19,37,3), dz_grid(37), qx(79), qy(
     2   19), qz(37), lambda(37), dz_cell(37)
      REAL m_den
      INTEGER k, j, i, m, j6, j7, j8, j9, j10, j11, vp_mythread, 
     1   vp_nmthreads, j4, j5
 
      j6 = (35 + vp_nmthreads - 1)/vp_nmthreads
      j7 = vp_mythread*j6 + 1
      j8 = min(35 - j7 + 1,j6)
      do k = j7, j8 + j7 - 1
         do j = 1, 17
            do i = 1, 77
               m_den = (1.67E-27*0.5)*(nf(2+i,j+1,k+1)+nf(1+i,j+1,k+1))
               gradpf(1+i,j+1,k+1,1) = (pf(2+i,j+1,k+1)-pf(1+i,j+1,k+1))
     1            /(10000.0*m_den)
               m_den = (1.67E-27*0.5)*(nf(1+i,j+2,k+1)+nf(1+i,j+1,k+1))
               gradpf(1+i,j+1,k+1,2) = (pf(1+i,j+2,k+1)-pf(1+i,j+1,k+1))
     1            /(10000.0*m_den)
               m_den = (1.67E-27*0.5)*(nf(1+i,j+1,k+2)+nf(1+i,j+1,k+1))
               gradpf(1+i,j+1,k+1,3) = (pf(1+i,j+1,k+2)-pf(1+i,j+1,k+1))
     1            /(dz_grid(k+1)*m_den)
            end do
         end do
      end do
 
 
      j9 = (35 + vp_nmthreads - 1)/vp_nmthreads
      j10 = vp_mythread*j9 + 1
      j11 = min(35 - j10 + 1,j9)
      do k = j10, j11 + j10 - 1
         do j = 1, 17
            do i = 1, 77
               ufp1(1+i,j+1,k+1,1) = 0.5*(uf(1+i,j+1,k+1,1)+ufp2(1+i,j+1
     1            ,k+1,1))
               ufp1(1+i,j+1,k+1,2) = 0.5*(uf(1+i,j+1,k+1,2)+ufp2(1+i,j+1
     1            ,k+1,2))
               ufp1(1+i,j+1,k+1,3) = 0.5*(uf(1+i,j+1,k+1,3)+ufp2(1+i,j+1
     1            ,k+1,3))
            end do
         end do
      end do
      return 
      end 
      subroutine vp_correct_uf2(vp_mythread, vp_nmthreads, uf, uf2, 
     1   uplus, ef, nuin)
      common /timestep/ntf, dtsub
      REAL dtsub, ntf
      REAL uf(79,19,37,3), uf2(79,19,37,3), uplus(79,19,37,3), ef(79,19,
     1   37,3), nuin(79,19,37)
      REAL r1
      INTEGER k,j,i,m,j12,j13,j14,vp_mythread,vp_nmthreads,j4,j5
 
      j12 = (35 + vp_nmthreads - 1)/vp_nmthreads
      j13 = vp_mythread*j12 + 1
      j14 = min(35 - j13 + 1,j12)
      do k = j13, j14 + j13 - 1
         do j = 1, 17
            do i = 1, 77
               uf2(1+i,j+1,k+1,1) = uf(1+i,j+1,k+1,1)
               uf2(1+i,j+1,k+1,2) = uf(1+i,j+1,k+1,2)
               uf2(1+i,j+1,k+1,3) = uf(1+i,j+1,k+1,3)
               r1 = nuin(1+i,j+1,k+1)
               uf(1+i,j+1,k+1,1) = uplus(1+i,j+1,k+1,1) + (0.5*dtsub)*ef
     1            (1+i,j+1,k+1,1) - ((0.5*dtsub)*r1)*uplus(1+i,j+1,k+1,1
     2            )
               uf(1+i,j+1,k+1,2) = uplus(1+i,j+1,k+1,2) + (0.5*dtsub)*ef
     1            (1+i,j+1,k+1,2) - ((0.5*dtsub)*r1)*uplus(1+i,j+1,k+1,2
     2            )
               uf(1+i,j+1,k+1,3) = uplus(1+i,j+1,k+1,3) + (0.5*dtsub)*ef
     1            (1+i,j+1,k+1,3) - ((0.5*dtsub)*r1)*uplus(1+i,j+1,k+1,3
     2            )
            end do
         end do
      end do
      return 
      end 
      subroutine vp_correct_uf3(vp_mythread, vp_nmthreads, uf)
      REAL uf(79,19,37,3)
      INTEGER j3, j2, j15, j16, j17, vp_mythread, vp_nmthreads, j4, j5
 
      j15 = (703 + vp_nmthreads - 1)/vp_nmthreads
      j16 = vp_mythread*j15 + 1
      j17 = min(703 - j16 + 1,j15)
      do j3 = j16, j17 + j16 - 1
         uf(78,j3,1,1) = ((-340.0))
         uf(79,j3,1,1) = ((-340.0))
      end do
      return 
      end 
c----------------------------------------------------------------------
 
 
c----------------------------------------------------------------------
      SUBROUTINE get_E(E,bt,btmf,aj,up,uf,uf2,np,nf,nu,gradP)
c E must be at time level m. We have uf at levels m-1/2 and m+1/2, so
c the average value is used for uf in the calculation of ui.
c----------------------------------------------------------------------
      PARAMETER (nx = 79, ny = 19, nz = 37)
 
c grid parameters
      PARAMETER (dx = 10000.0, dy = 10000.0)   !units in km
      PARAMETER (delz = 10000.0)          !dz at release coordinates
                                       !this is not the final dz
                                       !d_lamba determines dz
c time stepping parameters
      PARAMETER (tstart = 0.2)   !starting time
      PARAMETER (dt = 1.)     !main time step
      PARAMETER (nt = 5.0)        !number of time steps
c      PARAMETER (dtsub = 0.01) !subcycle time step
      PARAMETER (dtsub_init = 0.1) !subcycle time step
      PARAMETER (ntsub = 10.0)   !number of subcycle time steps
      PARAMETER (nout = 50)      !number of steps to diagnosic output
 
c logical variable for restart
C...Translated by Crescent Bay Software VAST-F    7.0K8 15:17:00   4/ 1/04    
      logical restart
      PARAMETER (restart = .false.)
      PARAMETER (mrestart = 500)      ! use -1 for no save
      PARAMETER (mbegin = 0)      !mstart
 
c inital ambient density at the release point.
      real nf_init,b0_init,nn_coef
      PARAMETER (nf_init = 4e15)   !was 0.01
      PARAMETER (b0_init = 2e-9)    !was 0.2
      PARAMETER (nn_coef = 1e5)
 
c electron ion collision frequency
      real nu_init, eta_init
      PARAMETER (nu_init = 0.001)
      PARAMETER (eta_init = 0.0)
 
c temporary directory for output
      character*9 temp_dir
      PARAMETER (temp_dir = '/tmp/pad/')
 
c total number of neutrals released
      real No
      PARAMETER (No=1.0e28)
 
c neutral cloud expansion characteristics
      PARAMETER(vo = 20.0)
      PARAMETER(vth = 19.3)
      PARAMETER(vsat = 0.0)
      PARAMETER(vsw = 340.0)
 
c total kinetic energy of released neutrals
      PARAMETER(Etot = 2.696e7)       !units of kg m^2 / s^2
 
c density scaling parameter, alpha, and ion particle array dims
 
      real alpha
 
c      PARAMETER (beta =  8.5227271e-18)
      PARAMETER (alpha = 1.9263418e-20) !mH
 
c max number of ion particles to be produced.  This parameter
c is used to dimension the particle arrays.
      integer*4 Ni_max
      PARAMETER (Ni_max = 10000)
 
c earth's magnetic dipole field geometry
c earth radius and equatorial B alt (defining L shell)
      PARAMETER (r_earth=6.38e3)                !earth radius
      PARAMETER (B_max_alt=2.0e3+r_earth)       !ro defining L shell
 
c misc constants
      real mu0,epsilon,pi,rtod,q,mO,mBa,O_to_Ba,km_to_m,kboltz,melec
      real mproton,tempf0
      PARAMETER (pi = 3.14159)
      PARAMETER (rtod = 180.0/pi)  !radians to degreesc
      PARAMETER (mu0 = pi*4.0e-7)  !magnetic permeability of free space
      PARAMETER (epsilon = 8.85e-12) !dielectric constant
      PARAMETER (q = 1.6e-19)        !electron charge (Coulombs)
      PARAMETER (melec = 9.1e-31)  !mass of electron (kg)
      PARAMETER (mproton = 1.67e-27)
      PARAMETER (mO = 1.67e-27)    !mass of H (kg)
c      PARAMETER (mO = 2.3e-25)    !mass of Ba (kg)
      PARAMETER (mBa = 17.0*mO)    !mass of Ba (kg)
      PARAMETER (O_to_Ba = mO/mBa) !convert E and B for particle move
      PARAMETER (km_to_m = 1e3)    !km to meters conversion
      PARAMETER (kboltz = 1.38e-29)   !kg km^2 / s^2 / K
      PARAMETER (tempf0 = 50*11600.)     !K
 
      real Qo, vrad, N_o, Rp, tau_photo
      integer S_radius
      PARAMETER (Qo = 1.0e28)  !neutral source rate
      PARAMETER (vrad = 1.0)   !escape velocity
c      PARAMETER (N_o = 5e34)   !Steady state neutral particle constant
      PARAMETER (N_o = 1e33)   !Steady state neutral particle constant
      PARAMETER (Rp = 20000.0)  !Pluto radius
      PARAMETER (tau_photo = 1.2e6)
c      PARAMETER (dNi = 2500)
      PARAMETER (S_radius = 3) !units of dx
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
c raw grid coordinate data
 
      real qx(nx), qy(ny), qz(nz)
                            !cartisian coordinates of grid points
      real lambda(nz)       !latitude of given grid point
      integer ri, rj, rk    !release coordinates
      real dz_grid(nz)      !dz for grid, variable in z direction
      real dz_cell(nz)      !dz for cell, variable in z direction
 
      common /coords/ qx,qy,qz,lambda,ri,rj,rk,dz_grid,dz_cell
 
 
c Total number of ions produced at a given time during the simulation
      integer*4 Ni_tot, Nn_tot
 
c Location (indices) of particles with the grid
      integer ijkp(Ni_max,3)
      logical ionized(Ni_max)
      common /ion_info/ Ni_tot, Nn_tot, ijkp, ionized
 
      real seed
      common /rndseed/ seed
 
c Indices for boundary cells
      integer ip,im,jp,jm,kp,km
 
c Total input energy (ions)
      real input_E,input_EeP,prev_Etot, bndry_Eflux
      real input_chex, input_bill
      common /energy/ input_E, input_EeP, prev_Etot, bndry_Eflux,
     x                input_chex, input_bill
 
c Momentum conservation
      real surf_tot(3),graduu_tot(3),ugradu_tot(3)
      common /momentum/ surf_tot,graduu_tot,ugradu_tot
 
c Dummy variables for doing vector operations and particle flux
c calculation....etc.
 
      real a(nx,ny,nz,3),c(nx,ny,nz,3),ct(nx,ny,nz,3)
      common /dummy/ a,c,ct
 
c Variable time step for fluid velocity update used in both
c particle update and fluid update
 
      real dtf,ntf,dtsub
      common /timestep/ ntf,dtsub
 
c Weight variables for trilinear interpolation
 
      real wght(Ni_max,8), wquad(Ni_max,3)
      common /weights/ wght, wquad
 
c variable for anomlous resistivity
 
      real etar(nx,ny,nz,3)
      real nuei
      common /resist/ etar,nuei
 
c variable for particle scale
 
      real beta, dNi
      common /scale/ beta, dNi
 
 
c mass array for multi-ion species
      real mrat,m_arr
      common /mass/ mrat(Ni_max),m_arr(Ni_max)
 
 
 
      real E(nx,ny,nz,3),
     x     bt(nx,ny,nz,3),
     x     btmf(nx,ny,nz,3),
     x     aj(nx,ny,nz,3),
     x     up(nx,ny,nz,3),
     x     uf(nx,ny,nz,3),
     x     uf2(nx,ny,nz,3),
     x     np(nx,ny,nz),
     x     nf(nx,ny,nz),
     x     nu(nx,ny,nz)
c     x     gradP(nx,ny,nz,3)
 
      real ntot(3)         !total density np + nf
      real fnp(3),fnf(3)   !fraction np and nf of n
      real npave(3)
 
c      real a(nx,ny,nz,3),
c     x     c(nx,ny,nz,3)  !dummy vars for doing cross product
 
      external vp_get_e1
      integer j3, j4, j5
      real r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, 
     .r15, r16
      call periodic_scalar(np)
      call periodic_scalar(nf)
 
      do 10 k = 2, 36
         do 10 j = 2, 18
            do 10 i = 2, 78
 
               ip = i + 1
               jp = j + 1
               kp = k + 1
 
c               if (ip .gt. nx) then ip = nx
c               if (jp .gt. ny) then jp = ny
c               if (kp .gt. nz) then kp = nz
 
               npave(1) = 0.5*(np(i,j,k)+np(ip,j,k))
               npave(2) = 0.5*(np(i,j,k)+np(i,jp,k))
               npave(3) = 0.5*(np(i,j,k)+np(i,j,kp))
 
               ntot(1) = npave(1) + 0.5*(nf(i,j,k)+nf(ip,j,k))
               ntot(2) = npave(2) + 0.5*(nf(i,j,k)+nf(i,jp,k))
               ntot(3) = npave(3) + 0.5*(nf(i,j,k)+nf(i,j,kp))
 
               fnp(1) = npave(1)/ntot(1)
               fnp(2) = npave(2)/ntot(2)
               fnp(3) = npave(3)/ntot(3)
 
               fnf(1) = 0.5*(nf(i,j,k)+nf(ip,j,k))/ntot(1)
               fnf(2) = 0.5*(nf(i,j,k)+nf(i,jp,k))/ntot(2)
               fnf(3) = 0.5*(nf(i,j,k)+nf(i,j,kp))/ntot(3)
 
c               ntot = np(i,j,k) + nf(i,j,k)
c               fnp = np(i,j,k)/ntot
c               fnf = nf(i,j,k)/ntot
 
               r14 = aj(i,j,k,1) - fnp(1)*up(i,j,k,1) - fnf(1)*0.5*(uf2(
     1            i,j,k,1)+uf(i,j,k,1))
               r15 = aj(i,j,k,2) - fnp(2)*up(i,j,k,2) - fnf(2)*0.5*(uf2(
     1            i,j,k,2)+uf(i,j,k,2))
               r16 = aj(i,j,k,3) - fnp(3)*up(i,j,k,3) - fnf(3)*0.5*(uf2(
     1            i,j,k,3)+uf(i,j,k,3))
               a(i,j,k,1) = r14
               a(i,j,k,2) = r15
               a(i,j,k,3) = r16
c                  a(i,j,k,m) = - fnp(m)*up(i,j,k,m) -
c     x                         fnf(m)*0.5*(uf2(i,j,k,m)+uf(i,j,k,m))
   10 continue
 
      call crossf(a,btmf,c)
 
 
      call vp_pcall (vp_get_e1, 2, 3, nu, aj, e)
 
 
c      call fix_tangential_E(E)
      call periodic(E)
c      call fix_tangential_E(E)
 
      return
      end
      subroutine vp_get_e1(vp_mythread, vp_nmthreads, nu, aj, e)
      common /dummy/a, c, ct
      common /resist/etar, nuei
      REAL nuei
      REAL c(79,19,37,3), nu(79,19,37), aj(79,19,37,3), e(79,19,37,3), a
     1   (79,19,37,3), ct(79,19,37,3), etar(79,19,37,3)
      REAL r13
      INTEGER k, j, i, m, j3, j4, j5, vp_mythread, vp_nmthreads, j1, j2
 
 
      j3 = (35 + vp_nmthreads - 1)/vp_nmthreads
      j4 = vp_mythread*j3 + 1
      j5 = min(35 - j4 + 1,j3)
      do k = j4, j5 + j4 - 1
         do j = 1, 17
            do i = 1, 77
               r13 = nu(1+i,j+1,k+1)
               e(1+i,j+1,k+1,1) = c(1+i,j+1,k+1,1) + r13*aj(1+i,j+1,k+1,
     1            1) + nuei*aj(1+i,j+1,k+1,1)    !- gradP(i,j,k,m)
               e(1+i,j+1,k+1,2) = c(1+i,j+1,k+1,2) + r13*aj(1+i,j+1,k+1,
     1            2) + nuei*aj(1+i,j+1,k+1,2)
               e(1+i,j+1,k+1,3) = c(1+i,j+1,k+1,3) + r13*aj(1+i,j+1,k+1,
     1            3) + nuei*aj(1+i,j+1,k+1,3)
            end do
         end do
      end do
      return 
      end 
c----------------------------------------------------------------------
 
 
c----------------------------------------------------------------------
      SUBROUTINE predict_B(b1,b12,b1p2,bt,btmf,E,aj,up,uf,uf2,np,nf,nu,
     x                     gradP)
c Predictor step in magnetic field update.
c----------------------------------------------------------------------
      PARAMETER (nx = 79, ny = 19, nz = 37)
 
c grid parameters
      PARAMETER (dx = 10000.0, dy = 10000.0)   !units in km
      PARAMETER (delz = 10000.0)          !dz at release coordinates
                                       !this is not the final dz
                                       !d_lamba determines dz
c time stepping parameters
      PARAMETER (tstart = 0.2)   !starting time
      PARAMETER (dt = 1.)     !main time step
      PARAMETER (nt = 5.0)        !number of time steps
c      PARAMETER (dtsub = 0.01) !subcycle time step
      PARAMETER (dtsub_init = 0.1) !subcycle time step
      PARAMETER (ntsub = 10.0)   !number of subcycle time steps
      PARAMETER (nout = 50)      !number of steps to diagnosic output
 
c logical variable for restart
C...Translated by Crescent Bay Software VAST-F    7.0K8 15:17:00   4/ 1/04    
      logical restart
      PARAMETER (restart = .false.)
      PARAMETER (mrestart = 500)      ! use -1 for no save
      PARAMETER (mbegin = 0)      !mstart
 
c inital ambient density at the release point.
      real nf_init,b0_init,nn_coef
      PARAMETER (nf_init = 4e15)   !was 0.01
      PARAMETER (b0_init = 2e-9)    !was 0.2
      PARAMETER (nn_coef = 1e5)
 
c electron ion collision frequency
      real nu_init, eta_init
      PARAMETER (nu_init = 0.001)
      PARAMETER (eta_init = 0.0)
 
c temporary directory for output
      character*9 temp_dir
      PARAMETER (temp_dir = '/tmp/pad/')
 
c total number of neutrals released
      real No
      PARAMETER (No=1.0e28)
 
c neutral cloud expansion characteristics
      PARAMETER(vo = 20.0)
      PARAMETER(vth = 19.3)
      PARAMETER(vsat = 0.0)
      PARAMETER(vsw = 340.0)
 
c total kinetic energy of released neutrals
      PARAMETER(Etot = 2.696e7)       !units of kg m^2 / s^2
 
c density scaling parameter, alpha, and ion particle array dims
 
      real alpha
 
c      PARAMETER (beta =  8.5227271e-18)
      PARAMETER (alpha = 1.9263418e-20) !mH
 
c max number of ion particles to be produced.  This parameter
c is used to dimension the particle arrays.
      integer*4 Ni_max
      PARAMETER (Ni_max = 10000)
 
c earth's magnetic dipole field geometry
c earth radius and equatorial B alt (defining L shell)
      PARAMETER (r_earth=6.38e3)                !earth radius
      PARAMETER (B_max_alt=2.0e3+r_earth)       !ro defining L shell
 
c misc constants
      real mu0,epsilon,pi,rtod,q,mO,mBa,O_to_Ba,km_to_m,kboltz,melec
      real mproton,tempf0
      PARAMETER (pi = 3.14159)
      PARAMETER (rtod = 180.0/pi)  !radians to degreesc
      PARAMETER (mu0 = pi*4.0e-7)  !magnetic permeability of free space
      PARAMETER (epsilon = 8.85e-12) !dielectric constant
      PARAMETER (q = 1.6e-19)        !electron charge (Coulombs)
      PARAMETER (melec = 9.1e-31)  !mass of electron (kg)
      PARAMETER (mproton = 1.67e-27)
      PARAMETER (mO = 1.67e-27)    !mass of H (kg)
c      PARAMETER (mO = 2.3e-25)    !mass of Ba (kg)
      PARAMETER (mBa = 17.0*mO)    !mass of Ba (kg)
      PARAMETER (O_to_Ba = mO/mBa) !convert E and B for particle move
      PARAMETER (km_to_m = 1e3)    !km to meters conversion
      PARAMETER (kboltz = 1.38e-29)   !kg km^2 / s^2 / K
      PARAMETER (tempf0 = 50*11600.)     !K
 
      real Qo, vrad, N_o, Rp, tau_photo
      integer S_radius
      PARAMETER (Qo = 1.0e28)  !neutral source rate
      PARAMETER (vrad = 1.0)   !escape velocity
c      PARAMETER (N_o = 5e34)   !Steady state neutral particle constant
      PARAMETER (N_o = 1e33)   !Steady state neutral particle constant
      PARAMETER (Rp = 20000.0)  !Pluto radius
      PARAMETER (tau_photo = 1.2e6)
c      PARAMETER (dNi = 2500)
      PARAMETER (S_radius = 3) !units of dx
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
c raw grid coordinate data
 
      real qx(nx), qy(ny), qz(nz)
                            !cartisian coordinates of grid points
      real lambda(nz)       !latitude of given grid point
      integer ri, rj, rk    !release coordinates
      real dz_grid(nz)      !dz for grid, variable in z direction
      real dz_cell(nz)      !dz for cell, variable in z direction
 
      common /coords/ qx,qy,qz,lambda,ri,rj,rk,dz_grid,dz_cell
 
 
c Total number of ions produced at a given time during the simulation
      integer*4 Ni_tot, Nn_tot
 
c Location (indices) of particles with the grid
      integer ijkp(Ni_max,3)
      logical ionized(Ni_max)
      common /ion_info/ Ni_tot, Nn_tot, ijkp, ionized
 
      real seed
      common /rndseed/ seed
 
c Indices for boundary cells
      integer ip,im,jp,jm,kp,km
 
c Total input energy (ions)
      real input_E,input_EeP,prev_Etot, bndry_Eflux
      real input_chex, input_bill
      common /energy/ input_E, input_EeP, prev_Etot, bndry_Eflux,
     x                input_chex, input_bill
 
c Momentum conservation
      real surf_tot(3),graduu_tot(3),ugradu_tot(3)
      common /momentum/ surf_tot,graduu_tot,ugradu_tot
 
c Dummy variables for doing vector operations and particle flux
c calculation....etc.
 
      real a(nx,ny,nz,3),c(nx,ny,nz,3),ct(nx,ny,nz,3)
      common /dummy/ a,c,ct
 
c Variable time step for fluid velocity update used in both
c particle update and fluid update
 
      real dtf,ntf,dtsub
      common /timestep/ ntf,dtsub
 
c Weight variables for trilinear interpolation
 
      real wght(Ni_max,8), wquad(Ni_max,3)
      common /weights/ wght, wquad
 
c variable for anomlous resistivity
 
      real etar(nx,ny,nz,3)
      real nuei
      common /resist/ etar,nuei
 
c variable for particle scale
 
      real beta, dNi
      common /scale/ beta, dNi
 
 
c mass array for multi-ion species
      real mrat,m_arr
      common /mass/ mrat(Ni_max),m_arr(Ni_max)
 
 
 
      real b1(nx,ny,nz,3),
     x     b12(nx,ny,nz,3),
     x     b1p2(nx,ny,nz,3),
     x     bt(nx,ny,nz,3),
     x     btmf(nx,ny,nz,3),
     x     E(nx,ny,nz,3),
     x     aj(nx,ny,nz,3),
     x     up(nx,ny,nz,3),
     x     uf(nx,ny,nz,3),
     x     uf2(nx,ny,nz,3),
     x     np(nx,ny,nz),
     x     nf(nx,ny,nz),
     x     nu(nx,ny,nz)
c     x     gradP(nx,ny,nz,3)
 
      real curl_E(nx,ny,nz,3)   !curl of E
 
      EXTERNAL VP_PREDICT_B1
      INTEGER J3, J4, J5
      call get_E(E,bt,btmf,aj,up,uf,uf2,np,nf,nu,gradP)  !E at time level m
 
 
      call curlE(E,curl_E)
c      call fix_tangential_E(E)
      call periodic(E)
c      call fix_tangential_E(E)
 
      CALL VP_PCALL (VP_PREDICT_B1, 2, 3, B12, CURL_E, B1P2)
 
c      call boundaries(b1p2)
c      call damp(b1p2)
      call periodic(b1p2)
      call fix_normal_b(b1p2)
 
      return
      end
      SUBROUTINE VP_PREDICT_B1(VP_MYTHREAD, VP_NMTHREADS, B12, CURL_E, 
     1   B1P2)
      COMMON /TIMESTEP/NTF, DTSUB
      REAL DTSUB, NTF
      REAL B12(79,19,37,3), CURL_E(79,19,37,3), B1P2(79,19,37,3)
      INTEGER K, J, I, M, J3, J4, J5, VP_MYTHREAD, VP_NMTHREADS, J1, J2
c      call fix_tangential_E(E)
 
      J3 = (35 + VP_NMTHREADS - 1)/VP_NMTHREADS
      J4 = VP_MYTHREAD*J3 + 1
      J5 = MIN(35 - J4 + 1,J3)
      DO K = J4, J5 + J4 - 1
         DO J = 1, 17
            DO I = 1, 77
               B1P2(1+I,J+1,K+1,1) = B12(1+I,J+1,K+1,1) - (2.0*DTSUB)*
     1            CURL_E(1+I,J+1,K+1,1)
               B1P2(1+I,J+1,K+1,2) = B12(1+I,J+1,K+1,2) - (2.0*DTSUB)*
     1            CURL_E(1+I,J+1,K+1,2)
               B1P2(1+I,J+1,K+1,3) = B12(1+I,J+1,K+1,3) - (2.0*DTSUB)*
     1            CURL_E(1+I,J+1,K+1,3)
            END DO
         END DO
      END DO
      RETURN 
      END 
c----------------------------------------------------------------------
 
 
c----------------------------------------------------------------------
      SUBROUTINE get_Ep1(E,b0,b1,b1p2,aj,up,uf,np,nf,nu,gradP,bdp)
c The main feature here is that E must be calculated at time level
c m + 1/2.  That means that we need B at m + 1/2.  So b1p1 is
c calculated as 0.5*(b1 + b1p2).  uf and np are already at time level
c m + 1/2, so they are used as is.
c----------------------------------------------------------------------
      PARAMETER (nx = 79, ny = 19, nz = 37)
 
c grid parameters
      PARAMETER (dx = 10000.0, dy = 10000.0)   !units in km
      PARAMETER (delz = 10000.0)          !dz at release coordinates
                                       !this is not the final dz
                                       !d_lamba determines dz
c time stepping parameters
      PARAMETER (tstart = 0.2)   !starting time
      PARAMETER (dt = 1.)     !main time step
      PARAMETER (nt = 5.0)        !number of time steps
c      PARAMETER (dtsub = 0.01) !subcycle time step
      PARAMETER (dtsub_init = 0.1) !subcycle time step
      PARAMETER (ntsub = 10.0)   !number of subcycle time steps
      PARAMETER (nout = 50)      !number of steps to diagnosic output
 
c logical variable for restart
C...Translated by Crescent Bay Software VAST-F    7.0K8 15:17:00   4/ 1/04    
      logical restart
      PARAMETER (restart = .false.)
      PARAMETER (mrestart = 500)      ! use -1 for no save
      PARAMETER (mbegin = 0)      !mstart
 
c inital ambient density at the release point.
      real nf_init,b0_init,nn_coef
      PARAMETER (nf_init = 4e15)   !was 0.01
      PARAMETER (b0_init = 2e-9)    !was 0.2
      PARAMETER (nn_coef = 1e5)
 
c electron ion collision frequency
      real nu_init, eta_init
      PARAMETER (nu_init = 0.001)
      PARAMETER (eta_init = 0.0)
 
c temporary directory for output
      character*9 temp_dir
      PARAMETER (temp_dir = '/tmp/pad/')
 
c total number of neutrals released
      real No
      PARAMETER (No=1.0e28)
 
c neutral cloud expansion characteristics
      PARAMETER(vo = 20.0)
      PARAMETER(vth = 19.3)
      PARAMETER(vsat = 0.0)
      PARAMETER(vsw = 340.0)
 
c total kinetic energy of released neutrals
      PARAMETER(Etot = 2.696e7)       !units of kg m^2 / s^2
 
c density scaling parameter, alpha, and ion particle array dims
 
      real alpha
 
c      PARAMETER (beta =  8.5227271e-18)
      PARAMETER (alpha = 1.9263418e-20) !mH
 
c max number of ion particles to be produced.  This parameter
c is used to dimension the particle arrays.
      integer*4 Ni_max
      PARAMETER (Ni_max = 10000)
 
c earth's magnetic dipole field geometry
c earth radius and equatorial B alt (defining L shell)
      PARAMETER (r_earth=6.38e3)                !earth radius
      PARAMETER (B_max_alt=2.0e3+r_earth)       !ro defining L shell
 
c misc constants
      real mu0,epsilon,pi,rtod,q,mO,mBa,O_to_Ba,km_to_m,kboltz,melec
      real mproton,tempf0
      PARAMETER (pi = 3.14159)
      PARAMETER (rtod = 180.0/pi)  !radians to degreesc
      PARAMETER (mu0 = pi*4.0e-7)  !magnetic permeability of free space
      PARAMETER (epsilon = 8.85e-12) !dielectric constant
      PARAMETER (q = 1.6e-19)        !electron charge (Coulombs)
      PARAMETER (melec = 9.1e-31)  !mass of electron (kg)
      PARAMETER (mproton = 1.67e-27)
      PARAMETER (mO = 1.67e-27)    !mass of H (kg)
c      PARAMETER (mO = 2.3e-25)    !mass of Ba (kg)
      PARAMETER (mBa = 17.0*mO)    !mass of Ba (kg)
      PARAMETER (O_to_Ba = mO/mBa) !convert E and B for particle move
      PARAMETER (km_to_m = 1e3)    !km to meters conversion
      PARAMETER (kboltz = 1.38e-29)   !kg km^2 / s^2 / K
      PARAMETER (tempf0 = 50*11600.)     !K
 
      real Qo, vrad, N_o, Rp, tau_photo
      integer S_radius
      PARAMETER (Qo = 1.0e28)  !neutral source rate
      PARAMETER (vrad = 1.0)   !escape velocity
c      PARAMETER (N_o = 5e34)   !Steady state neutral particle constant
      PARAMETER (N_o = 1e33)   !Steady state neutral particle constant
      PARAMETER (Rp = 20000.0)  !Pluto radius
      PARAMETER (tau_photo = 1.2e6)
c      PARAMETER (dNi = 2500)
      PARAMETER (S_radius = 3) !units of dx
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
c raw grid coordinate data
 
      real qx(nx), qy(ny), qz(nz)
                            !cartisian coordinates of grid points
      real lambda(nz)       !latitude of given grid point
      integer ri, rj, rk    !release coordinates
      real dz_grid(nz)      !dz for grid, variable in z direction
      real dz_cell(nz)      !dz for cell, variable in z direction
 
      common /coords/ qx,qy,qz,lambda,ri,rj,rk,dz_grid,dz_cell
 
 
c Total number of ions produced at a given time during the simulation
      integer*4 Ni_tot, Nn_tot
 
c Location (indices) of particles with the grid
      integer ijkp(Ni_max,3)
      logical ionized(Ni_max)
      common /ion_info/ Ni_tot, Nn_tot, ijkp, ionized
 
      real seed
      common /rndseed/ seed
 
c Indices for boundary cells
      integer ip,im,jp,jm,kp,km
 
c Total input energy (ions)
      real input_E,input_EeP,prev_Etot, bndry_Eflux
      real input_chex, input_bill
      common /energy/ input_E, input_EeP, prev_Etot, bndry_Eflux,
     x                input_chex, input_bill
 
c Momentum conservation
      real surf_tot(3),graduu_tot(3),ugradu_tot(3)
      common /momentum/ surf_tot,graduu_tot,ugradu_tot
 
c Dummy variables for doing vector operations and particle flux
c calculation....etc.
 
      real a(nx,ny,nz,3),c(nx,ny,nz,3),ct(nx,ny,nz,3)
      common /dummy/ a,c,ct
 
c Variable time step for fluid velocity update used in both
c particle update and fluid update
 
      real dtf,ntf,dtsub
      common /timestep/ ntf,dtsub
 
c Weight variables for trilinear interpolation
 
      real wght(Ni_max,8), wquad(Ni_max,3)
      common /weights/ wght, wquad
 
c variable for anomlous resistivity
 
      real etar(nx,ny,nz,3)
      real nuei
      common /resist/ etar,nuei
 
c variable for particle scale
 
      real beta, dNi
      common /scale/ beta, dNi
 
 
c mass array for multi-ion species
      real mrat,m_arr
      common /mass/ mrat(Ni_max),m_arr(Ni_max)
 
 
 
      real E(nx,ny,nz,3),
     x     b0(ny),
     x     b1(nx,ny,nz,3),
     x     b1p2(nx,ny,nz,3),
     x     aj(nx,ny,nz,3),
     x     up(nx,ny,nz,3),
     x     uf(nx,ny,nz,3),
     x     np(nx,ny,nz),
     x     nf(nx,ny,nz),
     x     nu(nx,ny,nz),
c     x     gradP(nx,ny,nz,3),
     x     bdp(nx,ny,nz,3)
 
      real b1p1(nx,ny,nz,3)   !b1 at time level m + 1/2
      real btp1(nx,ny,nz,3)   !bt at time level m + 1/2
      real btp1mf(nx,ny,nz,3) !btp1 at contravarient position
      real ntot(3)            !total density np + nf
      real fnp(3),fnf(3)      !fraction np and nf of n
      real npave(3)
 
c      real a(nx,ny,nz,3),
c     x     c(nx,ny,nz,3)    !dummy vars for doing cross product
 
 
      external vp_get_ep11
      external vp_get_ep12
      integer j3, j4, j5, j6, j7, j8
      real r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, 
     .r15, r16
      call vp_pcall (vp_get_ep11, 2, 6, b1p2, b1, bdp, btp1, b1p1, b0)
 
      call curlB(b1p1,nf,np,aj)
 
      call periodic_scalar(np)
      call periodic_scalar(nf)
 
      do 10 k = 2, 36
         do 10 j = 2, 18
            do 10 i = 2, 78
 
               ip = i + 1
               jp = j + 1
               kp = k + 1
 
c               if (ip .gt. nx) then ip = nx
c               if (jp .gt. ny) then jp = ny
c               if (kp .gt. nz) then kp = nz
 
               npave(1) = 0.5*(np(i,j,k)+np(ip,j,k))
               npave(2) = 0.5*(np(i,j,k)+np(i,jp,k))
               npave(3) = 0.5*(np(i,j,k)+np(i,j,kp))
 
               ntot(1) = npave(1) + 0.5*(nf(i,j,k)+nf(ip,j,k))
               ntot(2) = npave(2) + 0.5*(nf(i,j,k)+nf(i,jp,k))
               ntot(3) = npave(3) + 0.5*(nf(i,j,k)+nf(i,j,kp))
 
               fnp(1) = npave(1)/ntot(1)
               fnp(2) = npave(2)/ntot(2)
               fnp(3) = npave(3)/ntot(3)
 
               fnf(1) = 0.5*(nf(i,j,k)+nf(ip,j,k))/ntot(1)
               fnf(2) = 0.5*(nf(i,j,k)+nf(i,jp,k))/ntot(2)
               fnf(3) = 0.5*(nf(i,j,k)+nf(i,j,kp))/ntot(3)
 
c               ntot = np(i,j,k) + nf(i,j,k)
c               fnp = np(i,j,k)/ntot
c               fnf = nf(i,j,k)/ntot
 
               r14=aj(i,j,k,1)-fnp(1)*up(i,j,k,1)-fnf(1)*uf(i,j,k,1)
               r15=aj(i,j,k,2)-fnp(2)*up(i,j,k,2)-fnf(2)*uf(i,j,k,2)
               r16=aj(i,j,k,3)-fnp(3)*up(i,j,k,3)-fnf(3)*uf(i,j,k,3)
               a(i,j,k,1) = r14
               a(i,j,k,2) = r15
               a(i,j,k,3) = r16
c                  a(i,j,k,m) = - fnp(m)*up(i,j,k,m) -
c     x                           fnf(m)*uf(i,j,k,m)
   10 continue
 
      call cov_to_contra(btp1,btp1mf)
      call crossf(a,btp1mf,c)
 
 
      call vp_pcall (vp_get_ep12, 2, 3, nu, aj, e)
 
c      call fix_tangential_E(E)
      call periodic(E)
c      call fix_tangential_E(E)
 
      return
      end
      subroutine vp_get_ep11(vp_mythread, vp_nmthreads, b1p2, b1, bdp, 
     1   btp1, b1p1, b0)
      REAL b1p2(79,19,37,3), b1(79,19,37,3), bdp(79,19,37,3), btp1(79,19
     1   ,37,3), b1p1(79,19,37,3), b0(19)
      INTEGER k, j, i, j3, j4, j5, vp_mythread, vp_nmthreads, j1, j2
 
c      real a(nx,ny,nz,3),
c     x     c(nx,ny,nz,3)    !dummy vars for doing cross product
 
 
      j3 = (37 + vp_nmthreads - 1)/vp_nmthreads
      j4 = vp_mythread*j3 + 1
      j5 = min(37 - j4 + 1,j3)
      do k = j4, j5 + j4 - 1
         do j = 1, 19
            do i = 1, 79
               btp1(i,j,k,1) = (0.5*(b1p2(i,j,k,1)+b1(i,j,k,1))) + bdp(i
     1            ,j,k,1)
               b1p1(i,j,k,1) = (0.5*(b1p2(i,j,k,1)+b1(i,j,k,1)))
               btp1(i,j,k,2) = b0(j) + (0.5*(b1p2(i,j,k,2)+b1(i,j,k,2)))
     1             + bdp(i,j,k,2)
               b1p1(i,j,k,2) = (0.5*(b1p2(i,j,k,2)+b1(i,j,k,2)))
               btp1(i,j,k,3) = bdp(i,j,k,3) + (0.5*(b1p2(i,j,k,3)+b1(i,j
     1            ,k,3)))
               b1p1(i,j,k,3) = (0.5*(b1p2(i,j,k,3)+b1(i,j,k,3)))
            end do
         end do
      end do
      return 
      end 
      subroutine vp_get_ep12(vp_mythread, vp_nmthreads, nu, aj, e)
      common /dummy/a, c, ct
      common /resist/etar, nuei
      REAL nuei
      REAL c(79,19,37,3), nu(79,19,37), aj(79,19,37,3), e(79,19,37,3), a
     1   (79,19,37,3), ct(79,19,37,3), etar(79,19,37,3)
      REAL r13
      INTEGER k, j, i, m, j6, j7, j8, vp_mythread, vp_nmthreads, j1, j2
 
 
      j6 = (35 + vp_nmthreads - 1)/vp_nmthreads
      j7 = vp_mythread*j6 + 1
      j8 = min(35 - j7 + 1,j6)
      do k = j7, j8 + j7 - 1
         do j = 1, 17
            do i = 1, 77
               r13 = nu(1+i,j+1,k+1)
               e(1+i,j+1,k+1,1) = c(1+i,j+1,k+1,1) + r13*aj(1+i,j+1,k+1,
     1            1) + nuei*aj(1+i,j+1,k+1,1)    !- gradP(i,j,k,m)
               e(1+i,j+1,k+1,2) = c(1+i,j+1,k+1,2) + r13*aj(1+i,j+1,k+1,
     1            2) + nuei*aj(1+i,j+1,k+1,2)
               e(1+i,j+1,k+1,3) = c(1+i,j+1,k+1,3) + r13*aj(1+i,j+1,k+1,
     1            3) + nuei*aj(1+i,j+1,k+1,3)
            end do
         end do
      end do
      return 
      end 
c----------------------------------------------------------------------
 
 
c----------------------------------------------------------------------
      SUBROUTINE correct_B(b0,b1,b1p2,E,aj,up,uf,np,nf,nu,gradP,bdp)
c Corrector step in magnetic field update.
c----------------------------------------------------------------------
      PARAMETER (nx = 79, ny = 19, nz = 37)
 
c grid parameters
      PARAMETER (dx = 10000.0, dy = 10000.0)   !units in km
      PARAMETER (delz = 10000.0)          !dz at release coordinates
                                       !this is not the final dz
                                       !d_lamba determines dz
c time stepping parameters
      PARAMETER (tstart = 0.2)   !starting time
      PARAMETER (dt = 1.)     !main time step
      PARAMETER (nt = 5.0)        !number of time steps
c      PARAMETER (dtsub = 0.01) !subcycle time step
      PARAMETER (dtsub_init = 0.1) !subcycle time step
      PARAMETER (ntsub = 10.0)   !number of subcycle time steps
      PARAMETER (nout = 50)      !number of steps to diagnosic output
 
c logical variable for restart
C...Translated by Crescent Bay Software VAST-F    7.0K8 15:17:00   4/ 1/04    
      logical restart
      PARAMETER (restart = .false.)
      PARAMETER (mrestart = 500)      ! use -1 for no save
      PARAMETER (mbegin = 0)      !mstart
 
c inital ambient density at the release point.
      real nf_init,b0_init,nn_coef
      PARAMETER (nf_init = 4e15)   !was 0.01
      PARAMETER (b0_init = 2e-9)    !was 0.2
      PARAMETER (nn_coef = 1e5)
 
c electron ion collision frequency
      real nu_init, eta_init
      PARAMETER (nu_init = 0.001)
      PARAMETER (eta_init = 0.0)
 
c temporary directory for output
      character*9 temp_dir
      PARAMETER (temp_dir = '/tmp/pad/')
 
c total number of neutrals released
      real No
      PARAMETER (No=1.0e28)
 
c neutral cloud expansion characteristics
      PARAMETER(vo = 20.0)
      PARAMETER(vth = 19.3)
      PARAMETER(vsat = 0.0)
      PARAMETER(vsw = 340.0)
 
c total kinetic energy of released neutrals
      PARAMETER(Etot = 2.696e7)       !units of kg m^2 / s^2
 
c density scaling parameter, alpha, and ion particle array dims
 
      real alpha
 
c      PARAMETER (beta =  8.5227271e-18)
      PARAMETER (alpha = 1.9263418e-20) !mH
 
c max number of ion particles to be produced.  This parameter
c is used to dimension the particle arrays.
      integer*4 Ni_max
      PARAMETER (Ni_max = 10000)
 
c earth's magnetic dipole field geometry
c earth radius and equatorial B alt (defining L shell)
      PARAMETER (r_earth=6.38e3)                !earth radius
      PARAMETER (B_max_alt=2.0e3+r_earth)       !ro defining L shell
 
c misc constants
      real mu0,epsilon,pi,rtod,q,mO,mBa,O_to_Ba,km_to_m,kboltz,melec
      real mproton,tempf0
      PARAMETER (pi = 3.14159)
      PARAMETER (rtod = 180.0/pi)  !radians to degreesc
      PARAMETER (mu0 = pi*4.0e-7)  !magnetic permeability of free space
      PARAMETER (epsilon = 8.85e-12) !dielectric constant
      PARAMETER (q = 1.6e-19)        !electron charge (Coulombs)
      PARAMETER (melec = 9.1e-31)  !mass of electron (kg)
      PARAMETER (mproton = 1.67e-27)
      PARAMETER (mO = 1.67e-27)    !mass of H (kg)
c      PARAMETER (mO = 2.3e-25)    !mass of Ba (kg)
      PARAMETER (mBa = 17.0*mO)    !mass of Ba (kg)
      PARAMETER (O_to_Ba = mO/mBa) !convert E and B for particle move
      PARAMETER (km_to_m = 1e3)    !km to meters conversion
      PARAMETER (kboltz = 1.38e-29)   !kg km^2 / s^2 / K
      PARAMETER (tempf0 = 50*11600.)     !K
 
      real Qo, vrad, N_o, Rp, tau_photo
      integer S_radius
      PARAMETER (Qo = 1.0e28)  !neutral source rate
      PARAMETER (vrad = 1.0)   !escape velocity
c      PARAMETER (N_o = 5e34)   !Steady state neutral particle constant
      PARAMETER (N_o = 1e33)   !Steady state neutral particle constant
      PARAMETER (Rp = 20000.0)  !Pluto radius
      PARAMETER (tau_photo = 1.2e6)
c      PARAMETER (dNi = 2500)
      PARAMETER (S_radius = 3) !units of dx
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
c raw grid coordinate data
 
      real qx(nx), qy(ny), qz(nz)
                            !cartisian coordinates of grid points
      real lambda(nz)       !latitude of given grid point
      integer ri, rj, rk    !release coordinates
      real dz_grid(nz)      !dz for grid, variable in z direction
      real dz_cell(nz)      !dz for cell, variable in z direction
 
      common /coords/ qx,qy,qz,lambda,ri,rj,rk,dz_grid,dz_cell
 
 
c Total number of ions produced at a given time during the simulation
      integer*4 Ni_tot, Nn_tot
 
c Location (indices) of particles with the grid
      integer ijkp(Ni_max,3)
      logical ionized(Ni_max)
      common /ion_info/ Ni_tot, Nn_tot, ijkp, ionized
 
      real seed
      common /rndseed/ seed
 
c Indices for boundary cells
      integer ip,im,jp,jm,kp,km
 
c Total input energy (ions)
      real input_E,input_EeP,prev_Etot, bndry_Eflux
      real input_chex, input_bill
      common /energy/ input_E, input_EeP, prev_Etot, bndry_Eflux,
     x                input_chex, input_bill
 
c Momentum conservation
      real surf_tot(3),graduu_tot(3),ugradu_tot(3)
      common /momentum/ surf_tot,graduu_tot,ugradu_tot
 
c Dummy variables for doing vector operations and particle flux
c calculation....etc.
 
      real a(nx,ny,nz,3),c(nx,ny,nz,3),ct(nx,ny,nz,3)
      common /dummy/ a,c,ct
 
c Variable time step for fluid velocity update used in both
c particle update and fluid update
 
      real dtf,ntf,dtsub
      common /timestep/ ntf,dtsub
 
c Weight variables for trilinear interpolation
 
      real wght(Ni_max,8), wquad(Ni_max,3)
      common /weights/ wght, wquad
 
c variable for anomlous resistivity
 
      real etar(nx,ny,nz,3)
      real nuei
      common /resist/ etar,nuei
 
c variable for particle scale
 
      real beta, dNi
      common /scale/ beta, dNi
 
 
c mass array for multi-ion species
      real mrat,m_arr
      common /mass/ mrat(Ni_max),m_arr(Ni_max)
 
 
 
      real b0(ny),
     x     b1(nx,ny,nz,3),
     x     b1p2(nx,ny,nz,3),
     x     E(nx,ny,nz,3),
     x     aj(nx,ny,nz,3),
     x     up(nx,ny,nz,3),
     x     uf(nx,ny,nz,3),
     x     np(nx,ny,nz),
     x     nf(nx,ny,nz),
     x     nu(nx,ny,nz),
c     x     gradP(nx,ny,nz,3),
     x     bdp(nx,ny,nz,3)
 
      real curl_E(nx,ny,nz,3)            !curl of E
 
      EXTERNAL VP_CORRECT_B1
      INTEGER J3, J4, J5
      call get_Ep1(E,b0,b1,b1p2,aj,up,uf,np,nf,nu,gradP,bdp)
                                                   !E at time level m
 
      call curlE(E,curl_E)
c      call fix_tangential_E(E)
      call periodic(E)
c      call fix_tangential_E(E)
 
c      write(*,*) 'E cb...',E(23,8,14,1),E(23,8,14,2),E(23,8,14,3)
 
      CALL VP_PCALL (VP_CORRECT_B1, 2, 3, B1, CURL_E, B1P2)
 
c      call boundaries(b1p2)
c      call damp(b1p2)
      call periodic(b1p2)
      call fix_normal_b(b1p2)
 
      return
      end
      SUBROUTINE VP_CORRECT_B1(VP_MYTHREAD, VP_NMTHREADS, B1, CURL_E, 
     1   B1P2)
      COMMON /TIMESTEP/NTF, DTSUB
      REAL DTSUB, NTF
      REAL B1(79,19,37,3), CURL_E(79,19,37,3), B1P2(79,19,37,3)
      INTEGER K, J, I, M, J3, J4, J5, VP_MYTHREAD, VP_NMTHREADS, J1, J2
c      call fix_tangential_E(E)
 
c      write(*,*) 'E cb...',E(23,8,14,1),E(23,8,14,2),E(23,8,14,3)
 
      J3 = (35 + VP_NMTHREADS - 1)/VP_NMTHREADS
      J4 = VP_MYTHREAD*J3 + 1
      J5 = MIN(35 - J4 + 1,J3)
      DO K = J4, J5 + J4 - 1
         DO J = 1, 17
            DO I = 1, 77
               B1P2(1+I,J+1,K+1,1) = B1(1+I,J+1,K+1,1) - DTSUB*CURL_E(1+
     1            I,J+1,K+1,1)
               B1P2(1+I,J+1,K+1,2) = B1(1+I,J+1,K+1,2) - DTSUB*CURL_E(1+
     1            I,J+1,K+1,2)
               B1P2(1+I,J+1,K+1,3) = B1(1+I,J+1,K+1,3) - DTSUB*CURL_E(1+
     1            I,J+1,K+1,3)
            END DO
         END DO
      END DO
      RETURN 
      END 
c----------------------------------------------------------------------
 
 
c----------------------------------------------------------------------
      SUBROUTINE predict_nf(nf,nf1,nf3,nfp1,uf,divu,b1)
c----------------------------------------------------------------------
      PARAMETER (nx = 79, ny = 19, nz = 37)
 
c grid parameters
      PARAMETER (dx = 10000.0, dy = 10000.0)   !units in km
      PARAMETER (delz = 10000.0)          !dz at release coordinates
                                       !this is not the final dz
                                       !d_lamba determines dz
c time stepping parameters
      PARAMETER (tstart = 0.2)   !starting time
      PARAMETER (dt = 1.)     !main time step
      PARAMETER (nt = 5.0)        !number of time steps
c      PARAMETER (dtsub = 0.01) !subcycle time step
      PARAMETER (dtsub_init = 0.1) !subcycle time step
      PARAMETER (ntsub = 10.0)   !number of subcycle time steps
      PARAMETER (nout = 50)      !number of steps to diagnosic output
 
c logical variable for restart
C...Translated by Crescent Bay Software VAST-F    7.0K8 15:17:00   4/ 1/04    
      logical restart
      PARAMETER (restart = .false.)
      PARAMETER (mrestart = 500)      ! use -1 for no save
      PARAMETER (mbegin = 0)      !mstart
 
c inital ambient density at the release point.
      real nf_init,b0_init,nn_coef
      PARAMETER (nf_init = 4e15)   !was 0.01
      PARAMETER (b0_init = 2e-9)    !was 0.2
      PARAMETER (nn_coef = 1e5)
 
c electron ion collision frequency
      real nu_init, eta_init
      PARAMETER (nu_init = 0.001)
      PARAMETER (eta_init = 0.0)
 
c temporary directory for output
      character*9 temp_dir
      PARAMETER (temp_dir = '/tmp/pad/')
 
c total number of neutrals released
      real No
      PARAMETER (No=1.0e28)
 
c neutral cloud expansion characteristics
      PARAMETER(vo = 20.0)
      PARAMETER(vth = 19.3)
      PARAMETER(vsat = 0.0)
      PARAMETER(vsw = 340.0)
 
c total kinetic energy of released neutrals
      PARAMETER(Etot = 2.696e7)       !units of kg m^2 / s^2
 
c density scaling parameter, alpha, and ion particle array dims
 
      real alpha
 
c      PARAMETER (beta =  8.5227271e-18)
      PARAMETER (alpha = 1.9263418e-20) !mH
 
c max number of ion particles to be produced.  This parameter
c is used to dimension the particle arrays.
      integer*4 Ni_max
      PARAMETER (Ni_max = 10000)
 
c earth's magnetic dipole field geometry
c earth radius and equatorial B alt (defining L shell)
      PARAMETER (r_earth=6.38e3)                !earth radius
      PARAMETER (B_max_alt=2.0e3+r_earth)       !ro defining L shell
 
c misc constants
      real mu0,epsilon,pi,rtod,q,mO,mBa,O_to_Ba,km_to_m,kboltz,melec
      real mproton,tempf0
      PARAMETER (pi = 3.14159)
      PARAMETER (rtod = 180.0/pi)  !radians to degreesc
      PARAMETER (mu0 = pi*4.0e-7)  !magnetic permeability of free space
      PARAMETER (epsilon = 8.85e-12) !dielectric constant
      PARAMETER (q = 1.6e-19)        !electron charge (Coulombs)
      PARAMETER (melec = 9.1e-31)  !mass of electron (kg)
      PARAMETER (mproton = 1.67e-27)
      PARAMETER (mO = 1.67e-27)    !mass of H (kg)
c      PARAMETER (mO = 2.3e-25)    !mass of Ba (kg)
      PARAMETER (mBa = 17.0*mO)    !mass of Ba (kg)
      PARAMETER (O_to_Ba = mO/mBa) !convert E and B for particle move
      PARAMETER (km_to_m = 1e3)    !km to meters conversion
      PARAMETER (kboltz = 1.38e-29)   !kg km^2 / s^2 / K
      PARAMETER (tempf0 = 50*11600.)     !K
 
      real Qo, vrad, N_o, Rp, tau_photo
      integer S_radius
      PARAMETER (Qo = 1.0e28)  !neutral source rate
      PARAMETER (vrad = 1.0)   !escape velocity
c      PARAMETER (N_o = 5e34)   !Steady state neutral particle constant
      PARAMETER (N_o = 1e33)   !Steady state neutral particle constant
      PARAMETER (Rp = 20000.0)  !Pluto radius
      PARAMETER (tau_photo = 1.2e6)
c      PARAMETER (dNi = 2500)
      PARAMETER (S_radius = 3) !units of dx
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
c raw grid coordinate data
 
      real qx(nx), qy(ny), qz(nz)
                            !cartisian coordinates of grid points
      real lambda(nz)       !latitude of given grid point
      integer ri, rj, rk    !release coordinates
      real dz_grid(nz)      !dz for grid, variable in z direction
      real dz_cell(nz)      !dz for cell, variable in z direction
 
      common /coords/ qx,qy,qz,lambda,ri,rj,rk,dz_grid,dz_cell
 
 
c Total number of ions produced at a given time during the simulation
      integer*4 Ni_tot, Nn_tot
 
c Location (indices) of particles with the grid
      integer ijkp(Ni_max,3)
      logical ionized(Ni_max)
      common /ion_info/ Ni_tot, Nn_tot, ijkp, ionized
 
      real seed
      common /rndseed/ seed
 
c Indices for boundary cells
      integer ip,im,jp,jm,kp,km
 
c Total input energy (ions)
      real input_E,input_EeP,prev_Etot, bndry_Eflux
      real input_chex, input_bill
      common /energy/ input_E, input_EeP, prev_Etot, bndry_Eflux,
     x                input_chex, input_bill
 
c Momentum conservation
      real surf_tot(3),graduu_tot(3),ugradu_tot(3)
      common /momentum/ surf_tot,graduu_tot,ugradu_tot
 
c Dummy variables for doing vector operations and particle flux
c calculation....etc.
 
      real a(nx,ny,nz,3),c(nx,ny,nz,3),ct(nx,ny,nz,3)
      common /dummy/ a,c,ct
 
c Variable time step for fluid velocity update used in both
c particle update and fluid update
 
      real dtf,ntf,dtsub
      common /timestep/ ntf,dtsub
 
c Weight variables for trilinear interpolation
 
      real wght(Ni_max,8), wquad(Ni_max,3)
      common /weights/ wght, wquad
 
c variable for anomlous resistivity
 
      real etar(nx,ny,nz,3)
      real nuei
      common /resist/ etar,nuei
 
c variable for particle scale
 
      real beta, dNi
      common /scale/ beta, dNi
 
 
c mass array for multi-ion species
      real mrat,m_arr
      common /mass/ mrat(Ni_max),m_arr(Ni_max)
 
 
 
      real nf(nx,ny,nz),
     x     nf1(nx,ny,nz),
     x     nf3(nx,ny,nz),
     x     nfp1(nx,ny,nz),
     x     uf(nx,ny,nz,3),
     x     divu(nx,ny,nz),
     x     b1(nx,ny,nz,3)
 
      real minnf,maxnf
      real flx(nx,ny,nz,3)
      real ufc(nx,ny,nz,3)
      real b1c(nx,ny,nz,3)
 
c      call face_to_center(uf,ufc)
c      call face_to_center(b1,b1c)
 
      external vp_predict_nf1
      external vp_predict_nf2
      integer j3, j4, j5, j6, j7, j8
      real r1, r2, r3, r4, r5, r6, r7
      minnf = 10.0e20
      maxnf = 10.0e20
 
      call periodic_scalar(nf1)
 
      call vp_pcall (vp_predict_nf1, 2, 3, nf1, uf, flx)
 
      call periodic(flx)
 
      call vp_pcall(vp_predict_nf2,2,7,nf3,flx,nfp1,uf,divu,nf1,nf)
 
      return
      end
      subroutine vp_predict_nf1(vp_mythread, vp_nmthreads, nf1, uf, flx)
      REAL nf1(79,19,37), uf(79,19,37,3), flx(79,19,37,3)
      REAL r1, r2, r3, r4, r5, r6
      INTEGER j, k, i, j3, j4, j5, vp_mythread, vp_nmthreads, j1, j2
      j3 = (17 + vp_nmthreads - 1)/vp_nmthreads
      j4 = vp_mythread*j3 + 1
      j5 = min(17 - j4 + 1,j3)
      do j = j4, j5 + j4 - 1
         k = 1
         do i = 1, 77
            flx(1+i,j+1,2,1) = 0.5*(nf1(1+i,j+1,2)+nf1(2+i,j+1,2))*uf(1+
     1         i,j+1,2,1)
            flx(1+i,j+1,2,2) = 0.5*(nf1(1+i,j+1,2)+nf1(1+i,j+2,2))*uf(1+
     1         i,j+1,2,2)
            flx(1+i,j+1,2,3) = 0.5*(nf1(1+i,j+1,2)+nf1(1+i,j+1,3))*uf(1+
     1         i,j+1,2,3)
         end do
         do k = 1, 17
            do i = 1, 77
               r1 = 0.5*(nf1(1+i,j+1,k*2+1)+nf1(2+i,j+1,k*2+1))*uf(1+i,j
     1            +1,k*2+1,1)
               r2 = 0.5*(nf1(1+i,j+1,(k+1)*2)+nf1(2+i,j+1,(k+1)*2))*uf(1
     1            +i,j+1,(k+1)*2,1)
               flx(1+i,j+1,k*2+1,1) = r1
               flx(1+i,j+1,(k+1)*2,1) = r2
               r3 = 0.5*(nf1(1+i,j+1,k*2+1)+nf1(1+i,j+2,k*2+1))*uf(1+i,j
     1            +1,k*2+1,2)
               r4 = 0.5*(nf1(1+i,j+1,(k+1)*2)+nf1(1+i,j+2,(k+1)*2))*uf(1
     1            +i,j+1,(k+1)*2,2)
               flx(1+i,j+1,k*2+1,2) = r3
               flx(1+i,j+1,(k+1)*2,2) = r4
               r5 = 0.5*(nf1(1+i,j+1,k*2+1)+nf1(1+i,j+1,(k+1)*2))*uf(1+i
     1            ,j+1,k*2+1,3)
               r6 = 0.5*(nf1(1+i,j+1,(k+1)*2)+nf1(1+i,j+1,k*2+3))*uf(1+i
     1            ,j+1,(k+1)*2,3)
               flx(1+i,j+1,k*2+1,3) = r5
               flx(1+i,j+1,(k+1)*2,3) = r6
            end do
         end do
      end do
      return 
      end 
      subroutine vp_predict_nf2(vp_mythread, vp_nmthreads, nf3, flx, 
     1   nfp1, uf, divu, nf1, nf)
      common /coords/qx, qy, qz, lambda, ri, rj, rk, dz_grid, dz_cell
      common /timestep/ntf, dtsub
      REAL dtsub, ntf
      INTEGER ri, rj, rk
      REAL nf3(79,19,37), flx(79,19,37,3), dz_cell(37), nfp1(79,19,37), 
     1   uf(79,19,37,3), divu(79,19,37), nf1(79,19,37), nf(79,19,37), qx
     2   (79), qy(19), qz(37), lambda(37), dz_grid(37)
      REAL r7
      INTEGER j, k, i, j6, j7, j8, vp_mythread, vp_nmthreads, j1, j2
      j6 = (17 + vp_nmthreads - 1)/vp_nmthreads
      j7 = vp_mythread*j6 + 1
      j8 = min(17 - j7 + 1,j6)
      do j = j7, j8 + j7 - 1
         do k = 1, 35
            r7 = 1./10000.0
            do i = 1, 77
               nfp1(1+i,j+1,k+1) = nf3(1+i,j-1+2,k-1+2) - ((2.0*dtsub)/
     1            10000.0)*(flx(1+i,j-1+2,k-1+2,1)-flx(i,j-1+2,k-1+2,1))
     2             - ((2.0*dtsub)/10000.0)*(flx(1+i,j-1+2,k-1+2,2)-flx(1
     3            +i,j,k-1+2,2)) - ((2.0*dtsub)/dz_cell(k-1+2))*(flx(1+i
     4            ,j-1+2,k-1+2,3)-flx(1+i,j-1+2,k,3))
               divu(1+i,j+1,k+1) = (uf(1+i,j-1+2,k-1+2,1)-uf(i,j-1+2,k-1
     1            +2,1))*r7 + (uf(1+i,j-1+2,k-1+2,2)-uf(1+i,j,k-1+2,2))*
     2            r7 + (uf(1+i,j-1+2,k-1+2,3)-uf(1+i,j-1+2,k,3))/dz_cell
     3            (k-1+2)
               nf(1+i,j+1,k+1)=0.5*(nfp1(1+i,j+1,k+1)+nf1(1+i,j+1,k+1))
               nf3(1+i,j+1,k+1) = nf1(1+i,j+1,k+1)
            end do
         end do
      end do
      return 
      end 
c----------------------------------------------------------------------
 
 
c----------------------------------------------------------------------
      SUBROUTINE correct_nf(nf,nf1,ufp1)
c----------------------------------------------------------------------
      PARAMETER (nx = 79, ny = 19, nz = 37)
 
c grid parameters
      PARAMETER (dx = 10000.0, dy = 10000.0)   !units in km
      PARAMETER (delz = 10000.0)          !dz at release coordinates
                                       !this is not the final dz
                                       !d_lamba determines dz
c time stepping parameters
      PARAMETER (tstart = 0.2)   !starting time
      PARAMETER (dt = 1.)     !main time step
      PARAMETER (nt = 5.0)        !number of time steps
c      PARAMETER (dtsub = 0.01) !subcycle time step
      PARAMETER (dtsub_init = 0.1) !subcycle time step
      PARAMETER (ntsub = 10.0)   !number of subcycle time steps
      PARAMETER (nout = 50)      !number of steps to diagnosic output
 
c logical variable for restart
C...Translated by Crescent Bay Software VAST-F    7.0K8 15:17:00   4/ 1/04    
      logical restart
      PARAMETER (restart = .false.)
      PARAMETER (mrestart = 500)      ! use -1 for no save
      PARAMETER (mbegin = 0)      !mstart
 
c inital ambient density at the release point.
      real nf_init,b0_init,nn_coef
      PARAMETER (nf_init = 4e15)   !was 0.01
      PARAMETER (b0_init = 2e-9)    !was 0.2
      PARAMETER (nn_coef = 1e5)
 
c electron ion collision frequency
      real nu_init, eta_init
      PARAMETER (nu_init = 0.001)
      PARAMETER (eta_init = 0.0)
 
c temporary directory for output
      character*9 temp_dir
      PARAMETER (temp_dir = '/tmp/pad/')
 
c total number of neutrals released
      real No
      PARAMETER (No=1.0e28)
 
c neutral cloud expansion characteristics
      PARAMETER(vo = 20.0)
      PARAMETER(vth = 19.3)
      PARAMETER(vsat = 0.0)
      PARAMETER(vsw = 340.0)
 
c total kinetic energy of released neutrals
      PARAMETER(Etot = 2.696e7)       !units of kg m^2 / s^2
 
c density scaling parameter, alpha, and ion particle array dims
 
      real alpha
 
c      PARAMETER (beta =  8.5227271e-18)
      PARAMETER (alpha = 1.9263418e-20) !mH
 
c max number of ion particles to be produced.  This parameter
c is used to dimension the particle arrays.
      integer*4 Ni_max
      PARAMETER (Ni_max = 10000)
 
c earth's magnetic dipole field geometry
c earth radius and equatorial B alt (defining L shell)
      PARAMETER (r_earth=6.38e3)                !earth radius
      PARAMETER (B_max_alt=2.0e3+r_earth)       !ro defining L shell
 
c misc constants
      real mu0,epsilon,pi,rtod,q,mO,mBa,O_to_Ba,km_to_m,kboltz,melec
      real mproton,tempf0
      PARAMETER (pi = 3.14159)
      PARAMETER (rtod = 180.0/pi)  !radians to degreesc
      PARAMETER (mu0 = pi*4.0e-7)  !magnetic permeability of free space
      PARAMETER (epsilon = 8.85e-12) !dielectric constant
      PARAMETER (q = 1.6e-19)        !electron charge (Coulombs)
      PARAMETER (melec = 9.1e-31)  !mass of electron (kg)
      PARAMETER (mproton = 1.67e-27)
      PARAMETER (mO = 1.67e-27)    !mass of H (kg)
c      PARAMETER (mO = 2.3e-25)    !mass of Ba (kg)
      PARAMETER (mBa = 17.0*mO)    !mass of Ba (kg)
      PARAMETER (O_to_Ba = mO/mBa) !convert E and B for particle move
      PARAMETER (km_to_m = 1e3)    !km to meters conversion
      PARAMETER (kboltz = 1.38e-29)   !kg km^2 / s^2 / K
      PARAMETER (tempf0 = 50*11600.)     !K
 
      real Qo, vrad, N_o, Rp, tau_photo
      integer S_radius
      PARAMETER (Qo = 1.0e28)  !neutral source rate
      PARAMETER (vrad = 1.0)   !escape velocity
c      PARAMETER (N_o = 5e34)   !Steady state neutral particle constant
      PARAMETER (N_o = 1e33)   !Steady state neutral particle constant
      PARAMETER (Rp = 20000.0)  !Pluto radius
      PARAMETER (tau_photo = 1.2e6)
c      PARAMETER (dNi = 2500)
      PARAMETER (S_radius = 3) !units of dx
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
c raw grid coordinate data
 
      real qx(nx), qy(ny), qz(nz)
                            !cartisian coordinates of grid points
      real lambda(nz)       !latitude of given grid point
      integer ri, rj, rk    !release coordinates
      real dz_grid(nz)      !dz for grid, variable in z direction
      real dz_cell(nz)      !dz for cell, variable in z direction
 
      common /coords/ qx,qy,qz,lambda,ri,rj,rk,dz_grid,dz_cell
 
 
c Total number of ions produced at a given time during the simulation
      integer*4 Ni_tot, Nn_tot
 
c Location (indices) of particles with the grid
      integer ijkp(Ni_max,3)
      logical ionized(Ni_max)
      common /ion_info/ Ni_tot, Nn_tot, ijkp, ionized
 
      real seed
      common /rndseed/ seed
 
c Indices for boundary cells
      integer ip,im,jp,jm,kp,km
 
c Total input energy (ions)
      real input_E,input_EeP,prev_Etot, bndry_Eflux
      real input_chex, input_bill
      common /energy/ input_E, input_EeP, prev_Etot, bndry_Eflux,
     x                input_chex, input_bill
 
c Momentum conservation
      real surf_tot(3),graduu_tot(3),ugradu_tot(3)
      common /momentum/ surf_tot,graduu_tot,ugradu_tot
 
c Dummy variables for doing vector operations and particle flux
c calculation....etc.
 
      real a(nx,ny,nz,3),c(nx,ny,nz,3),ct(nx,ny,nz,3)
      common /dummy/ a,c,ct
 
c Variable time step for fluid velocity update used in both
c particle update and fluid update
 
      real dtf,ntf,dtsub
      common /timestep/ ntf,dtsub
 
c Weight variables for trilinear interpolation
 
      real wght(Ni_max,8), wquad(Ni_max,3)
      common /weights/ wght, wquad
 
c variable for anomlous resistivity
 
      real etar(nx,ny,nz,3)
      real nuei
      common /resist/ etar,nuei
 
c variable for particle scale
 
      real beta, dNi
      common /scale/ beta, dNi
 
 
c mass array for multi-ion species
      real mrat,m_arr
      common /mass/ mrat(Ni_max),m_arr(Ni_max)
 
 
 
      real nf(nx,ny,nz),
     x     nf1(nx,ny,nz),
     x     ufp1(nx,ny,nz,3)
 
      real flx(nx,ny,nz,3)
      real minnf,maxnf
      real ufp1c(nx,ny,nz,3)
 
      external vp_correct_nf1
      external vp_correct_nf2
      integer j3, j4, j5, j6, j7, j8
      real r1, r2, r3, r4, r5, r6
      minnf = 10.0000000000000000e20
      maxnf = 10.0000000000000000e20
 
c      call face_to_center(ufp1,ufp1c)
 
      call vp_pcall (vp_correct_nf1, 2, 3, nf, ufp1, flx)
 
      call periodic(flx)
 
      call vp_pcall (vp_correct_nf2, 2, 5, nf1, flx, minnf, maxnf, nf)
 
 
      call periodic_scalar(nf1)
 
      return
      end
      subroutine vp_correct_nf1(vp_mythread,vp_nmthreads,nf,ufp1,flx)
      REAL nf(79,19,37), ufp1(79,19,37,3), flx(79,19,37,3)
      REAL r1, r2, r3, r4, r5, r6
      INTEGER j, k, i, j3, j4, j5, vp_mythread, vp_nmthreads, j1, j2
      j3 = (17 + vp_nmthreads - 1)/vp_nmthreads
      j4 = vp_mythread*j3 + 1
      j5 = min(17 - j4 + 1,j3)
      do j = j4, j5 + j4 - 1
         k = 1
         do i = 1, 77
            flx(1+i,j+1,2,1) = 0.5*(nf(1+i,j+1,2)+nf(2+i,j+1,2))*ufp1(1+
     1         i,j+1,2,1)
            flx(1+i,j+1,2,2) = 0.5*(nf(1+i,j+1,2)+nf(1+i,j+2,2))*ufp1(1+
     1         i,j+1,2,2)
            flx(1+i,j+1,2,3) = 0.5*(nf(1+i,j+1,2)+nf(1+i,j+1,3))*ufp1(1+
     1         i,j+1,2,3)
         end do
         do k = 1, 17
            do i = 1, 77
               r1 = 0.5*(nf(1+i,j+1,k*2+1)+nf(2+i,j+1,k*2+1))*ufp1(1+i,j
     1            +1,k*2+1,1)
               r2 = 0.5*(nf(1+i,j+1,(k+1)*2)+nf(2+i,j+1,(k+1)*2))*ufp1(1
     1            +i,j+1,(k+1)*2,1)
               flx(1+i,j+1,k*2+1,1) = r1
               flx(1+i,j+1,(k+1)*2,1) = r2
               r3 = 0.5*(nf(1+i,j+1,k*2+1)+nf(1+i,j+2,k*2+1))*ufp1(1+i,j
     1            +1,k*2+1,2)
               r4 = 0.5*(nf(1+i,j+1,(k+1)*2)+nf(1+i,j+2,(k+1)*2))*ufp1(1
     1            +i,j+1,(k+1)*2,2)
               flx(1+i,j+1,k*2+1,2) = r3
               flx(1+i,j+1,(k+1)*2,2) = r4
               r5 = 0.5*(nf(1+i,j+1,k*2+1)+nf(1+i,j+1,(k+1)*2))*ufp1(1+i
     1            ,j+1,k*2+1,3)
               r6 = 0.5*(nf(1+i,j+1,(k+1)*2)+nf(1+i,j+1,k*2+3))*ufp1(1+i
     1            ,j+1,(k+1)*2,3)
               flx(1+i,j+1,k*2+1,3) = r5
               flx(1+i,j+1,(k+1)*2,3) = r6
            end do
         end do
      end do
      return 
      end 
      subroutine vp_correct_nf2(vp_mythread, vp_nmthreads, nf1, flx, 
     1   minnf, maxnf, nf)
      common /coords/qx, qy, qz, lambda, ri, rj, rk, dz_grid, dz_cell
      common /timestep/ntf, dtsub
      REAL dtsub, minnf, maxnf, ntf
      INTEGER ri, rj, rk
      REAL nf1(79,19,37), flx(79,19,37,3), dz_cell(37), nf(79,19,37), qx
     1   (79), qy(19), qz(37), lambda(37), dz_grid(37)
      REAL thenminnf, thenmaxnf
      INTEGER j, k, i, j6, j7, j8, vp_mythread, vp_nmthreads, j1, j2
      j6 = (17 + vp_nmthreads - 1)/vp_nmthreads
      j7 = vp_mythread*j6 + 1
      j8 = min(17 - j7 + 1,j6)
      do j = j7, j8 + j7 - 1
         do k = 1, 35
            do i = 1, 77
               nf1(1+i,j+1,k+1) = nf1(1+i,j-1+2,k-1+2) - (dtsub/10000.0)
     1            *(flx(1+i,j-1+2,k-1+2,1)-flx(i,j-1+2,k-1+2,1)) - (
     2            dtsub/10000.0)*(flx(1+i,j-1+2,k-1+2,2)-flx(1+i,j,k-1+2
     3            ,2)) - (dtsub/dz_cell(k-1+2))*(flx(1+i,j-1+2,k-1+2,3)-
     4            flx(1+i,j-1+2,k,3))
            end do
         end do
      end do
      return 
      end 
c----------------------------------------------------------------------
 
 
c----------------------------------------------------------------------
      SUBROUTINE trans_nf_Lax(nf,nf1,nfp1,uf)
c----------------------------------------------------------------------
      PARAMETER (nx = 79, ny = 19, nz = 37)
 
c grid parameters
      PARAMETER (dx = 10000.0, dy = 10000.0)   !units in km
      PARAMETER (delz = 10000.0)          !dz at release coordinates
                                       !this is not the final dz
                                       !d_lamba determines dz
c time stepping parameters
      PARAMETER (tstart = 0.2)   !starting time
      PARAMETER (dt = 1.)     !main time step
      PARAMETER (nt = 5.0)        !number of time steps
c      PARAMETER (dtsub = 0.01) !subcycle time step
      PARAMETER (dtsub_init = 0.1) !subcycle time step
      PARAMETER (ntsub = 10.0)   !number of subcycle time steps
      PARAMETER (nout = 50)      !number of steps to diagnosic output
 
c logical variable for restart
C...Translated by Crescent Bay Software VAST-F    7.0K8 15:17:00   4/ 1/04    
      logical restart
      PARAMETER (restart = .false.)
      PARAMETER (mrestart = 500)      ! use -1 for no save
      PARAMETER (mbegin = 0)      !mstart
 
c inital ambient density at the release point.
      real nf_init,b0_init,nn_coef
      PARAMETER (nf_init = 4e15)   !was 0.01
      PARAMETER (b0_init = 2e-9)    !was 0.2
      PARAMETER (nn_coef = 1e5)
 
c electron ion collision frequency
      real nu_init, eta_init
      PARAMETER (nu_init = 0.001)
      PARAMETER (eta_init = 0.0)
 
c temporary directory for output
      character*9 temp_dir
      PARAMETER (temp_dir = '/tmp/pad/')
 
c total number of neutrals released
      real No
      PARAMETER (No=1.0e28)
 
c neutral cloud expansion characteristics
      PARAMETER(vo = 20.0)
      PARAMETER(vth = 19.3)
      PARAMETER(vsat = 0.0)
      PARAMETER(vsw = 340.0)
 
c total kinetic energy of released neutrals
      PARAMETER(Etot = 2.696e7)       !units of kg m^2 / s^2
 
c density scaling parameter, alpha, and ion particle array dims
 
      real alpha
 
c      PARAMETER (beta =  8.5227271e-18)
      PARAMETER (alpha = 1.9263418e-20) !mH
 
c max number of ion particles to be produced.  This parameter
c is used to dimension the particle arrays.
      integer*4 Ni_max
      PARAMETER (Ni_max = 10000)
 
c earth's magnetic dipole field geometry
c earth radius and equatorial B alt (defining L shell)
      PARAMETER (r_earth=6.38e3)                !earth radius
      PARAMETER (B_max_alt=2.0e3+r_earth)       !ro defining L shell
 
c misc constants
      real mu0,epsilon,pi,rtod,q,mO,mBa,O_to_Ba,km_to_m,kboltz,melec
      real mproton,tempf0
      PARAMETER (pi = 3.14159)
      PARAMETER (rtod = 180.0/pi)  !radians to degreesc
      PARAMETER (mu0 = pi*4.0e-7)  !magnetic permeability of free space
      PARAMETER (epsilon = 8.85e-12) !dielectric constant
      PARAMETER (q = 1.6e-19)        !electron charge (Coulombs)
      PARAMETER (melec = 9.1e-31)  !mass of electron (kg)
      PARAMETER (mproton = 1.67e-27)
      PARAMETER (mO = 1.67e-27)    !mass of H (kg)
c      PARAMETER (mO = 2.3e-25)    !mass of Ba (kg)
      PARAMETER (mBa = 17.0*mO)    !mass of Ba (kg)
      PARAMETER (O_to_Ba = mO/mBa) !convert E and B for particle move
      PARAMETER (km_to_m = 1e3)    !km to meters conversion
      PARAMETER (kboltz = 1.38e-29)   !kg km^2 / s^2 / K
      PARAMETER (tempf0 = 50*11600.)     !K
 
      real Qo, vrad, N_o, Rp, tau_photo
      integer S_radius
      PARAMETER (Qo = 1.0e28)  !neutral source rate
      PARAMETER (vrad = 1.0)   !escape velocity
c      PARAMETER (N_o = 5e34)   !Steady state neutral particle constant
      PARAMETER (N_o = 1e33)   !Steady state neutral particle constant
      PARAMETER (Rp = 20000.0)  !Pluto radius
      PARAMETER (tau_photo = 1.2e6)
c      PARAMETER (dNi = 2500)
      PARAMETER (S_radius = 3) !units of dx
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
c raw grid coordinate data
 
      real qx(nx), qy(ny), qz(nz)
                            !cartisian coordinates of grid points
      real lambda(nz)       !latitude of given grid point
      integer ri, rj, rk    !release coordinates
      real dz_grid(nz)      !dz for grid, variable in z direction
      real dz_cell(nz)      !dz for cell, variable in z direction
 
      common /coords/ qx,qy,qz,lambda,ri,rj,rk,dz_grid,dz_cell
 
 
c Total number of ions produced at a given time during the simulation
      integer*4 Ni_tot, Nn_tot
 
c Location (indices) of particles with the grid
      integer ijkp(Ni_max,3)
      logical ionized(Ni_max)
      common /ion_info/ Ni_tot, Nn_tot, ijkp, ionized
 
      real seed
      common /rndseed/ seed
 
c Indices for boundary cells
      integer ip,im,jp,jm,kp,km
 
c Total input energy (ions)
      real input_E,input_EeP,prev_Etot, bndry_Eflux
      real input_chex, input_bill
      common /energy/ input_E, input_EeP, prev_Etot, bndry_Eflux,
     x                input_chex, input_bill
 
c Momentum conservation
      real surf_tot(3),graduu_tot(3),ugradu_tot(3)
      common /momentum/ surf_tot,graduu_tot,ugradu_tot
 
c Dummy variables for doing vector operations and particle flux
c calculation....etc.
 
      real a(nx,ny,nz,3),c(nx,ny,nz,3),ct(nx,ny,nz,3)
      common /dummy/ a,c,ct
 
c Variable time step for fluid velocity update used in both
c particle update and fluid update
 
      real dtf,ntf,dtsub
      common /timestep/ ntf,dtsub
 
c Weight variables for trilinear interpolation
 
      real wght(Ni_max,8), wquad(Ni_max,3)
      common /weights/ wght, wquad
 
c variable for anomlous resistivity
 
      real etar(nx,ny,nz,3)
      real nuei
      common /resist/ etar,nuei
 
c variable for particle scale
 
      real beta, dNi
      common /scale/ beta, dNi
 
 
c mass array for multi-ion species
      real mrat,m_arr
      common /mass/ mrat(Ni_max),m_arr(Ni_max)
 
 
 
      real nf(nx,ny,nz),
     x     nf1(nx,ny,nz),
     x     nfp1(nx,ny,nz),
     x     uf(nx,ny,nz,3)
 
      real minnf,maxnf
      real flx(nx,ny,nz,3)
      real ufc(nx,ny,nz,3)
 
c      call face_to_center(uf,ufc)
 
      external vp_trans_nf_lax1
      external vp_trans_nf_lax2
      external vp_trans_nf_lax3
      integer j3, j4, j5, j6, j7, j8, j9, j10, j11
      real r1, r2, r3, r4, r5, r6, r7, r8, r9, r10
      minnf = 10.0e20
      maxnf = 10.0e20
 
      call vp_pcall (vp_trans_nf_lax1, 2, 2, nfp1, nf1)
 
      call periodic_scalar(nf1)
 
      call vp_pcall (vp_trans_nf_lax2, 2, 3, nf1, uf, flx)
 
      call periodic(flx)
 
      call vp_pcall (vp_trans_nf_lax3, 2, 4, nf1, flx, nfp1, nf)
 
      call periodic_scalar(nf)
 
      return
      end
      subroutine vp_trans_nf_lax1(vp_mythread, vp_nmthreads, nfp1, nf1)
      REAL nfp1(79,19,37), nf1(79,19,37)
      REAL r1, r2, r3, r4
      INTEGER j, k, i, j3, j4, j5, vp_mythread, vp_nmthreads, j1, j2
      j3 = (17 + vp_nmthreads - 1)/vp_nmthreads
      j4 = vp_mythread*j3 + 1
      j5 = min(17 - j4 + 1,j3)
      do j = j4, j5 + j4 - 1
         do k = 1, 3
            do i = 1, 77
               nf1(1+i,j+1,k+1) = nfp1(1+i,j+1,k+1)
            end do
         end do
         do k = 1, 8
            do i = 1, 77
               r1 = nfp1(1+i,j+1,k*4+1)
               r2 = nfp1(1+i,j+1,k*4+2)
               r3 = nfp1(1+i,j+1,k*4+3)
               r4 = nfp1(1+i,j+1,(k+1)*4)
               nf1(1+i,j+1,k*4+1) = r1
               nf1(1+i,j+1,k*4+2) = r2
               nf1(1+i,j+1,k*4+3) = r3
               nf1(1+i,j+1,(k+1)*4) = r4
            end do
         end do
      end do
      return 
      end 
      subroutine vp_trans_nf_lax2(vp_mythread, vp_nmthreads, nf1, uf, 
     1   flx)
      REAL nf1(79,19,37), uf(79,19,37,3), flx(79,19,37,3)
      REAL r5, r6, r7, r8, r9, r10
      INTEGER j, k, i, j6, j7, j8, vp_mythread, vp_nmthreads, j1, j2
      j6 = (17 + vp_nmthreads - 1)/vp_nmthreads
      j7 = vp_mythread*j6 + 1
      j8 = min(17 - j7 + 1,j6)
      do j = j7, j8 + j7 - 1
         k = 1
         do i = 1, 77
            flx(1+i,j+1,2,1) = 0.5*(nf1(1+i,j+1,2)+nf1(2+i,j+1,2))*uf(1+
     1         i,j+1,2,1)
            flx(1+i,j+1,2,2) = 0.5*(nf1(1+i,j+1,2)+nf1(1+i,j+2,2))*uf(1+
     1         i,j+1,2,2)
            flx(1+i,j+1,2,3) = 0.5*(nf1(1+i,j+1,2)+nf1(1+i,j+1,3))*uf(1+
     1         i,j+1,2,3)
         end do
         do k = 1, 17
            do i = 1, 77
               r5 = 0.5*(nf1(1+i,j+1,k*2+1)+nf1(2+i,j+1,k*2+1))*uf(1+i,j
     1            +1,k*2+1,1)
               r6 = 0.5*(nf1(1+i,j+1,(k+1)*2)+nf1(2+i,j+1,(k+1)*2))*uf(1
     1            +i,j+1,(k+1)*2,1)
               flx(1+i,j+1,k*2+1,1) = r5
               flx(1+i,j+1,(k+1)*2,1) = r6
               r7 = 0.5*(nf1(1+i,j+1,k*2+1)+nf1(1+i,j+2,k*2+1))*uf(1+i,j
     1            +1,k*2+1,2)
               r8 = 0.5*(nf1(1+i,j+1,(k+1)*2)+nf1(1+i,j+2,(k+1)*2))*uf(1
     1            +i,j+1,(k+1)*2,2)
               flx(1+i,j+1,k*2+1,2) = r7
               flx(1+i,j+1,(k+1)*2,2) = r8
               r9 = 0.5*(nf1(1+i,j+1,k*2+1)+nf1(1+i,j+1,(k+1)*2))*uf(1+i
     1            ,j+1,k*2+1,3)
               r10 = 0.5*(nf1(1+i,j+1,(k+1)*2)+nf1(1+i,j+1,k*2+3))*uf(1+
     1            i,j+1,(k+1)*2,3)
               flx(1+i,j+1,k*2+1,3) = r9
               flx(1+i,j+1,(k+1)*2,3) = r10
            end do
         end do
      end do
      return 
      end 
      subroutine vp_trans_nf_lax3(vp_mythread, vp_nmthreads, nf1, flx, 
     1   nfp1, nf)
      common /coords/qx, qy, qz, lambda, ri, rj, rk, dz_grid, dz_cell
      common /timestep/ntf, dtsub
      REAL dtsub, ntf
      INTEGER ri, rj, rk
      REAL nf1(79,19,37), flx(79,19,37,3), dz_grid(37), nfp1(79,19,37), 
     1   nf(79,19,37), qx(79), qy(19), qz(37), lambda(37), dz_cell(37)
      INTEGER j, k, i, j9, j10, j11, vp_mythread, vp_nmthreads, j1, j2
      j9 = (17 + vp_nmthreads - 1)/vp_nmthreads
      j10 = vp_mythread*j9 + 1
      j11 = min(17 - j10 + 1,j9)
      do j = j10, j11 + j10 - 1
         do k = 1, 35
            do i = 1, 77
               nfp1(1+i,j+1,k+1) = (1.0/6.0)*(nf1(2+i,j-1+2,k-1+2)+nf1(i
     1            ,j-1+2,k-1+2)+nf1(1+i,j-1+3,k-1+2)+nf1(1+i,j,k-1+2)+
     2            nf1(1+i,j-1+2,k-1+3)+nf1(1+i,j-1+2,k)) - 0.5*(dtsub/
     3            10000.0)*(flx(2+i,j-1+2,k-1+2,1)-flx(i,j-1+2,k-1+2,1))
     4             - 0.5*(dtsub/10000.0)*(flx(1+i,j-1+3,k-1+2,2)-flx(1+i
     5            ,j,k-1+2,2)) - (dtsub/(dz_grid(k-1+2)+dz_grid(k-1+3)))
     6            *(flx(1+i,j-1+2,k-1+3,3)-flx(1+i,j-1+2,k,3))
               nf(1+i,j+1,k+1)=0.5*(nf1(1+i,j+1,k+1)+nfp1(1+i,j+1,k+1))
            end do
         end do
      end do
      return 
      end 
c----------------------------------------------------------------------
 
 
c----------------------------------------------------------------------
      SUBROUTINE trans_nf_LaxWend1(nf,nf1,nfp1,uf)
c----------------------------------------------------------------------
      PARAMETER (nx = 79, ny = 19, nz = 37)
 
c grid parameters
      PARAMETER (dx = 10000.0, dy = 10000.0)   !units in km
      PARAMETER (delz = 10000.0)          !dz at release coordinates
                                       !this is not the final dz
                                       !d_lamba determines dz
c time stepping parameters
      PARAMETER (tstart = 0.2)   !starting time
      PARAMETER (dt = 1.)     !main time step
      PARAMETER (nt = 5.0)        !number of time steps
c      PARAMETER (dtsub = 0.01) !subcycle time step
      PARAMETER (dtsub_init = 0.1) !subcycle time step
      PARAMETER (ntsub = 10.0)   !number of subcycle time steps
      PARAMETER (nout = 50)      !number of steps to diagnosic output
 
c logical variable for restart
C...Translated by Crescent Bay Software VAST-F    7.0K8 15:17:00   4/ 1/04    
      integer j1, j2, j3
      logical restart
      PARAMETER (restart = .false.)
      PARAMETER (mrestart = 500)      ! use -1 for no save
      PARAMETER (mbegin = 0)      !mstart
 
c inital ambient density at the release point.
      real nf_init,b0_init,nn_coef
      PARAMETER (nf_init = 4e15)   !was 0.01
      PARAMETER (b0_init = 2e-9)    !was 0.2
      PARAMETER (nn_coef = 1e5)
 
c electron ion collision frequency
      real nu_init, eta_init
      PARAMETER (nu_init = 0.001)
      PARAMETER (eta_init = 0.0)
 
c temporary directory for output
      character*9 temp_dir
      PARAMETER (temp_dir = '/tmp/pad/')
 
c total number of neutrals released
      real No
      PARAMETER (No=1.0e28)
 
c neutral cloud expansion characteristics
      PARAMETER(vo = 20.0)
      PARAMETER(vth = 19.3)
      PARAMETER(vsat = 0.0)
      PARAMETER(vsw = 340.0)
 
c total kinetic energy of released neutrals
      PARAMETER(Etot = 2.696e7)       !units of kg m^2 / s^2
 
c density scaling parameter, alpha, and ion particle array dims
 
      real alpha
 
c      PARAMETER (beta =  8.5227271e-18)
      PARAMETER (alpha = 1.9263418e-20) !mH
 
c max number of ion particles to be produced.  This parameter
c is used to dimension the particle arrays.
      integer*4 Ni_max
      PARAMETER (Ni_max = 10000)
 
c earth's magnetic dipole field geometry
c earth radius and equatorial B alt (defining L shell)
      PARAMETER (r_earth=6.38e3)                !earth radius
      PARAMETER (B_max_alt=2.0e3+r_earth)       !ro defining L shell
 
c misc constants
      real mu0,epsilon,pi,rtod,q,mO,mBa,O_to_Ba,km_to_m,kboltz,melec
      real mproton,tempf0
      PARAMETER (pi = 3.14159)
      PARAMETER (rtod = 180.0/pi)  !radians to degreesc
      PARAMETER (mu0 = pi*4.0e-7)  !magnetic permeability of free space
      PARAMETER (epsilon = 8.85e-12) !dielectric constant
      PARAMETER (q = 1.6e-19)        !electron charge (Coulombs)
      PARAMETER (melec = 9.1e-31)  !mass of electron (kg)
      PARAMETER (mproton = 1.67e-27)
      PARAMETER (mO = 1.67e-27)    !mass of H (kg)
c      PARAMETER (mO = 2.3e-25)    !mass of Ba (kg)
      PARAMETER (mBa = 17.0*mO)    !mass of Ba (kg)
      PARAMETER (O_to_Ba = mO/mBa) !convert E and B for particle move
      PARAMETER (km_to_m = 1e3)    !km to meters conversion
      PARAMETER (kboltz = 1.38e-29)   !kg km^2 / s^2 / K
      PARAMETER (tempf0 = 50*11600.)     !K
 
      real Qo, vrad, N_o, Rp, tau_photo
      integer S_radius
      PARAMETER (Qo = 1.0e28)  !neutral source rate
      PARAMETER (vrad = 1.0)   !escape velocity
c      PARAMETER (N_o = 5e34)   !Steady state neutral particle constant
      PARAMETER (N_o = 1e33)   !Steady state neutral particle constant
      PARAMETER (Rp = 20000.0)  !Pluto radius
      PARAMETER (tau_photo = 1.2e6)
c      PARAMETER (dNi = 2500)
      PARAMETER (S_radius = 3) !units of dx
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
c raw grid coordinate data
 
      real qx(nx), qy(ny), qz(nz)
                            !cartisian coordinates of grid points
      real lambda(nz)       !latitude of given grid point
      integer ri, rj, rk    !release coordinates
      real dz_grid(nz)      !dz for grid, variable in z direction
      real dz_cell(nz)      !dz for cell, variable in z direction
 
      common /coords/ qx,qy,qz,lambda,ri,rj,rk,dz_grid,dz_cell
 
 
c Total number of ions produced at a given time during the simulation
      integer*4 Ni_tot, Nn_tot
 
c Location (indices) of particles with the grid
      integer ijkp(Ni_max,3)
      logical ionized(Ni_max)
      common /ion_info/ Ni_tot, Nn_tot, ijkp, ionized
 
      real seed
      common /rndseed/ seed
 
c Indices for boundary cells
      integer ip,im,jp,jm,kp,km
 
c Total input energy (ions)
      real input_E,input_EeP,prev_Etot, bndry_Eflux
      real input_chex, input_bill
      common /energy/ input_E, input_EeP, prev_Etot, bndry_Eflux,
     x                input_chex, input_bill
 
c Momentum conservation
      real surf_tot(3),graduu_tot(3),ugradu_tot(3)
      common /momentum/ surf_tot,graduu_tot,ugradu_tot
 
c Dummy variables for doing vector operations and particle flux
c calculation....etc.
 
      real a(nx,ny,nz,3),c(nx,ny,nz,3),ct(nx,ny,nz,3)
      common /dummy/ a,c,ct
 
c Variable time step for fluid velocity update used in both
c particle update and fluid update
 
      real dtf,ntf,dtsub
      common /timestep/ ntf,dtsub
 
c Weight variables for trilinear interpolation
 
      real wght(Ni_max,8), wquad(Ni_max,3)
      common /weights/ wght, wquad
 
c variable for anomlous resistivity
 
      real etar(nx,ny,nz,3)
      real nuei
      common /resist/ etar,nuei
 
c variable for particle scale
 
      real beta, dNi
      common /scale/ beta, dNi
 
 
c mass array for multi-ion species
      real mrat,m_arr
      common /mass/ mrat(Ni_max),m_arr(Ni_max)
 
 
 
      real nf(nx,ny,nz),
     x     nf1(nx,ny,nz),
     x     nfp1(nx,ny,nz),
     x     uf(nx,ny,nz,3)
 
      real minnf,maxnf
      real flx(nx,ny,nz,3)
 
      external vp_trans_nf_laxwend11
      external vp_trans_nf_laxwend12
      external vp_trans_nf_laxwend13
      external vp_trans_nf_laxwend14
      integer j6, j7, j8, j9, j10, j11, j12, j13, j14, j15, j16, j17
      real r1, r2, r3, r4, r5, r6, r7, r8, r9, r10
      call vp_pcall (vp_trans_nf_laxwend11, 2, 2, nfp1, nf1)
 
      call periodic_scalar(nf1)
 
      call vp_pcall (vp_trans_nf_laxwend12, 2, 3, nf1, uf, flx)
 
      call periodic(flx)
      call periodic_scalar(nf1)
 
      call vp_pcall (vp_trans_nf_laxwend13, 2, 3, nf1, flx, nf)
 
      call periodic_scalar(nf)
      call vp_pcall (vp_trans_nf_laxwend14, 2, 1, nf)
 
      return
      end
      subroutine vp_trans_nf_laxwend11(vp_mythread, vp_nmthreads, nfp1, 
     1   nf1)
      REAL nfp1(79,19,37), nf1(79,19,37)
      REAL r1, r2, r3, r4
      INTEGER k, j, i, j6, j7, j8, vp_mythread, vp_nmthreads, j4, j5
 
      j6 = (35 + vp_nmthreads - 1)/vp_nmthreads
      j7 = vp_mythread*j6 + 1
      j8 = min(35 - j7 + 1,j6)
      do k = j7, j8 + j7 - 1
         j = 1
         do i = 1, 77
            nf1(1+i,2,k+1) = nfp1(1+i,2,k+1)
         end do
         do j = 1, 4
            do i = 1, 77
               r1 = nfp1(1+i,j*4-1,k+1)
               r2 = nfp1(1+i,j*4,k+1)
               r3 = nfp1(1+i,j*4+1,k+1)
               r4 = nfp1(1+i,j*4+2,k+1)
               nf1(1+i,j*4-1,k+1) = r1
               nf1(1+i,j*4,k+1) = r2
               nf1(1+i,j*4+1,k+1) = r3
               nf1(1+i,j*4+2,k+1) = r4
            end do
         end do
      end do
      return 
      end 
      subroutine vp_trans_nf_laxwend12(vp_mythread, vp_nmthreads, nf1, 
     1   uf, flx)
      REAL nf1(79,19,37), uf(79,19,37,3), flx(79,19,37,3)
      REAL r5, r6, r7, r8, r9, r10
      INTEGER k, j, i, j9, j10, j11, vp_mythread, vp_nmthreads, j4, j5
 
      j9 = (35 + vp_nmthreads - 1)/vp_nmthreads
      j10 = vp_mythread*j9 + 1
      j11 = min(35 - j10 + 1,j9)
      do k = j10, j11 + j10 - 1
         j = 1
         do i = 1, 77
            flx(1+i,2,k+1,1) = 0.5*(nf1(1+i,2,k+1)+nf1(2+i,2,k+1))*uf(1+
     1         i,2,k+1,1)
            flx(1+i,2,k+1,2) = 0.5*(nf1(1+i,2,k+1)+nf1(1+i,3,k+1))*uf(1+
     1         i,2,k+1,2)
            flx(1+i,2,k+1,3) = 0.5*(nf1(1+i,2,k+1)+nf1(1+i,2,k+2))*uf(1+
     1         i,2,k+1,3)
         end do
         do j = 1, 8
            do i = 1, 77
               r5 = 0.5*(nf1(1+i,j*2+1,k+1)+nf1(2+i,j*2+1,k+1))*uf(1+i,j
     1            *2+1,k+1,1)
               r6 = 0.5*(nf1(1+i,(j+1)*2,k+1)+nf1(2+i,(j+1)*2,k+1))*uf(1
     1            +i,(j+1)*2,k+1,1)
               flx(1+i,j*2+1,k+1,1) = r5
               flx(1+i,(j+1)*2,k+1,1) = r6
               r7 = 0.5*(nf1(1+i,j*2+1,k+1)+nf1(1+i,(j+1)*2,k+1))*uf(1+i
     1            ,j*2+1,k+1,2)
               r8 = 0.5*(nf1(1+i,(j+1)*2,k+1)+nf1(1+i,j*2+3,k+1))*uf(1+i
     1            ,(j+1)*2,k+1,2)
               flx(1+i,j*2+1,k+1,2) = r7
               flx(1+i,(j+1)*2,k+1,2) = r8
               r9 = 0.5*(nf1(1+i,j*2+1,k+1)+nf1(1+i,j*2+1,k+2))*uf(1+i,j
     1            *2+1,k+1,3)
               r10 = 0.5*(nf1(1+i,(j+1)*2,k+1)+nf1(1+i,(j+1)*2,k+2))*uf(
     1            1+i,(j+1)*2,k+1,3)
               flx(1+i,j*2+1,k+1,3) = r9
               flx(1+i,(j+1)*2,k+1,3) = r10
            end do
         end do
      end do
      return 
      end 
      subroutine vp_trans_nf_laxwend13(vp_mythread, vp_nmthreads, nf1, 
     1   flx, nf)
      common /coords/qx, qy, qz, lambda, ri, rj, rk, dz_grid, dz_cell
      common /timestep/ntf, dtsub
      REAL dtsub, ntf
      INTEGER ri, rj, rk
      REAL nf1(79,19,37), flx(79,19,37,3), dz_grid(37), nf(79,19,37), qx
     1   (79), qy(19), qz(37), lambda(37), dz_cell(37)
      INTEGER k, j, i, j12, j13, j14, vp_mythread, vp_nmthreads, j4, j5
 
      j12 = (35 + vp_nmthreads - 1)/vp_nmthreads
      j13 = vp_mythread*j12 + 1
      j14 = min(35 - j13 + 1,j12)
      do k = j13, j14 + j13 - 1
         do j = 1, 17
            do i = 1, 77
               nf(1+i,j+1,k+1) = (1.0/12.0)*(nf1(2+i,j-1+2,k-1+2)+nf1(i,
     1            j-1+2,k-1+2)+nf1(1+i,j-1+3,k-1+2)+nf1(1+i,j,k-1+2)+nf1
     2            (1+i,j-1+2,k-1+3)+nf1(1+i,j-1+2,k)+6.0*nf1(1+i,j-1+2,k
     3            -1+2)) - 0.5*(dtsub/10000.0)*(flx(1+i,j-1+2,k-1+2,1)-
     4            flx(i,j-1+2,k-1+2,1)) - 0.5*(dtsub/10000.0)*(flx(1+i,j
     5            -1+2,k-1+2,2)-flx(1+i,j,k-1+2,2)) - 0.5*(dtsub/dz_grid
     6            (k-1+2))*(flx(1+i,j-1+2,k-1+2,3)-flx(1+i,j-1+2,k,3))
            end do
         end do
      end do
      return 
      end 
      subroutine vp_trans_nf_laxwend14(vp_mythread, vp_nmthreads, nf)
      REAL nf(79,19,37)
      INTEGER j3, j2, j15, j16, j17, vp_mythread, vp_nmthreads, j4, j5
      j15 = (703 + vp_nmthreads - 1)/vp_nmthreads
      j16 = vp_mythread*j15 + 1
      j17 = min(703 - j16 + 1,j15)
      do j3 = j16, j17 + j16 - 1
         nf(78,j3,1) = 4E15
         nf(79,j3,1) = 4E15
      end do
      return 
      end 
c----------------------------------------------------------------------
 
 
 
c----------------------------------------------------------------------
      SUBROUTINE trans_nf_LaxWend2(nf,nf1,nfp1,ufp1)
c----------------------------------------------------------------------
      PARAMETER (nx = 79, ny = 19, nz = 37)
 
c grid parameters
      PARAMETER (dx = 10000.0, dy = 10000.0)   !units in km
      PARAMETER (delz = 10000.0)          !dz at release coordinates
                                       !this is not the final dz
                                       !d_lamba determines dz
c time stepping parameters
      PARAMETER (tstart = 0.2)   !starting time
      PARAMETER (dt = 1.)     !main time step
      PARAMETER (nt = 5.0)        !number of time steps
c      PARAMETER (dtsub = 0.01) !subcycle time step
      PARAMETER (dtsub_init = 0.1) !subcycle time step
      PARAMETER (ntsub = 10.0)   !number of subcycle time steps
      PARAMETER (nout = 50)      !number of steps to diagnosic output
 
c logical variable for restart
C...Translated by Crescent Bay Software VAST-F    7.0K8 15:17:00   4/ 1/04    
      integer j1, j2, j3
      logical restart
      PARAMETER (restart = .false.)
      PARAMETER (mrestart = 500)      ! use -1 for no save
      PARAMETER (mbegin = 0)      !mstart
 
c inital ambient density at the release point.
      real nf_init,b0_init,nn_coef
      PARAMETER (nf_init = 4e15)   !was 0.01
      PARAMETER (b0_init = 2e-9)    !was 0.2
      PARAMETER (nn_coef = 1e5)
 
c electron ion collision frequency
      real nu_init, eta_init
      PARAMETER (nu_init = 0.001)
      PARAMETER (eta_init = 0.0)
 
c temporary directory for output
      character*9 temp_dir
      PARAMETER (temp_dir = '/tmp/pad/')
 
c total number of neutrals released
      real No
      PARAMETER (No=1.0e28)
 
c neutral cloud expansion characteristics
      PARAMETER(vo = 20.0)
      PARAMETER(vth = 19.3)
      PARAMETER(vsat = 0.0)
      PARAMETER(vsw = 340.0)
 
c total kinetic energy of released neutrals
      PARAMETER(Etot = 2.696e7)       !units of kg m^2 / s^2
 
c density scaling parameter, alpha, and ion particle array dims
 
      real alpha
 
c      PARAMETER (beta =  8.5227271e-18)
      PARAMETER (alpha = 1.9263418e-20) !mH
 
c max number of ion particles to be produced.  This parameter
c is used to dimension the particle arrays.
      integer*4 Ni_max
      PARAMETER (Ni_max = 10000)
 
c earth's magnetic dipole field geometry
c earth radius and equatorial B alt (defining L shell)
      PARAMETER (r_earth=6.38e3)                !earth radius
      PARAMETER (B_max_alt=2.0e3+r_earth)       !ro defining L shell
 
c misc constants
      real mu0,epsilon,pi,rtod,q,mO,mBa,O_to_Ba,km_to_m,kboltz,melec
      real mproton,tempf0
      PARAMETER (pi = 3.14159)
      PARAMETER (rtod = 180.0/pi)  !radians to degreesc
      PARAMETER (mu0 = pi*4.0e-7)  !magnetic permeability of free space
      PARAMETER (epsilon = 8.85e-12) !dielectric constant
      PARAMETER (q = 1.6e-19)        !electron charge (Coulombs)
      PARAMETER (melec = 9.1e-31)  !mass of electron (kg)
      PARAMETER (mproton = 1.67e-27)
      PARAMETER (mO = 1.67e-27)    !mass of H (kg)
c      PARAMETER (mO = 2.3e-25)    !mass of Ba (kg)
      PARAMETER (mBa = 17.0*mO)    !mass of Ba (kg)
      PARAMETER (O_to_Ba = mO/mBa) !convert E and B for particle move
      PARAMETER (km_to_m = 1e3)    !km to meters conversion
      PARAMETER (kboltz = 1.38e-29)   !kg km^2 / s^2 / K
      PARAMETER (tempf0 = 50*11600.)     !K
 
      real Qo, vrad, N_o, Rp, tau_photo
      integer S_radius
      PARAMETER (Qo = 1.0e28)  !neutral source rate
      PARAMETER (vrad = 1.0)   !escape velocity
c      PARAMETER (N_o = 5e34)   !Steady state neutral particle constant
      PARAMETER (N_o = 1e33)   !Steady state neutral particle constant
      PARAMETER (Rp = 20000.0)  !Pluto radius
      PARAMETER (tau_photo = 1.2e6)
c      PARAMETER (dNi = 2500)
      PARAMETER (S_radius = 3) !units of dx
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
c raw grid coordinate data
 
      real qx(nx), qy(ny), qz(nz)
                            !cartisian coordinates of grid points
      real lambda(nz)       !latitude of given grid point
      integer ri, rj, rk    !release coordinates
      real dz_grid(nz)      !dz for grid, variable in z direction
      real dz_cell(nz)      !dz for cell, variable in z direction
 
      common /coords/ qx,qy,qz,lambda,ri,rj,rk,dz_grid,dz_cell
 
 
c Total number of ions produced at a given time during the simulation
      integer*4 Ni_tot, Nn_tot
 
c Location (indices) of particles with the grid
      integer ijkp(Ni_max,3)
      logical ionized(Ni_max)
      common /ion_info/ Ni_tot, Nn_tot, ijkp, ionized
 
      real seed
      common /rndseed/ seed
 
c Indices for boundary cells
      integer ip,im,jp,jm,kp,km
 
c Total input energy (ions)
      real input_E,input_EeP,prev_Etot, bndry_Eflux
      real input_chex, input_bill
      common /energy/ input_E, input_EeP, prev_Etot, bndry_Eflux,
     x                input_chex, input_bill
 
c Momentum conservation
      real surf_tot(3),graduu_tot(3),ugradu_tot(3)
      common /momentum/ surf_tot,graduu_tot,ugradu_tot
 
c Dummy variables for doing vector operations and particle flux
c calculation....etc.
 
      real a(nx,ny,nz,3),c(nx,ny,nz,3),ct(nx,ny,nz,3)
      common /dummy/ a,c,ct
 
c Variable time step for fluid velocity update used in both
c particle update and fluid update
 
      real dtf,ntf,dtsub
      common /timestep/ ntf,dtsub
 
c Weight variables for trilinear interpolation
 
      real wght(Ni_max,8), wquad(Ni_max,3)
      common /weights/ wght, wquad
 
c variable for anomlous resistivity
 
      real etar(nx,ny,nz,3)
      real nuei
      common /resist/ etar,nuei
 
c variable for particle scale
 
      real beta, dNi
      common /scale/ beta, dNi
 
 
c mass array for multi-ion species
      real mrat,m_arr
      common /mass/ mrat(Ni_max),m_arr(Ni_max)
 
 
 
      real nf(nx,ny,nz),
     x     nf1(nx,ny,nz),
     x     nfp1(nx,ny,nz),
     x     ufp1(nx,ny,nz,3)
 
      real minnf,maxnf
      real flx(nx,ny,nz,3)
 
      external vp_trans_nf_laxwend21
      external vp_trans_nf_laxwend22
      external vp_trans_nf_laxwend23
      integer j6, j7, j8, j9, j10, j11, j12, j13, j14
      real r1, r2, r3, r4, r5, r6
      call periodic_scalar(nf)
 
      call vp_pcall (vp_trans_nf_laxwend21, 2, 3, nf, ufp1, flx)
 
      call periodic(flx)
      call periodic_scalar(nf1)
 
      call vp_pcall (vp_trans_nf_laxwend22, 2, 3, nf1, flx, nfp1)
 
      call periodic_scalar(nfp1)
      call vp_pcall (vp_trans_nf_laxwend23, 2, 1, nfp1)
 
      return
      end
      subroutine vp_trans_nf_laxwend21(vp_mythread, vp_nmthreads, nf, 
     1   ufp1, flx)
      REAL nf(79,19,37), ufp1(79,19,37,3), flx(79,19,37,3)
      REAL r1, r2, r3, r4, r5, r6
      INTEGER k, j, i, j6, j7, j8, vp_mythread, vp_nmthreads, j4, j5
 
      j6 = (35 + vp_nmthreads - 1)/vp_nmthreads
      j7 = vp_mythread*j6 + 1
      j8 = min(35 - j7 + 1,j6)
      do k = j7, j8 + j7 - 1
         j = 1
         do i = 1, 77
            flx(1+i,2,k+1,1) = 0.5*(nf(1+i,2,k+1)+nf(2+i,2,k+1))*ufp1(1+
     1         i,2,k+1,1)
            flx(1+i,2,k+1,2) = 0.5*(nf(1+i,2,k+1)+nf(1+i,3,k+1))*ufp1(1+
     1         i,2,k+1,2)
            flx(1+i,2,k+1,3) = 0.5*(nf(1+i,2,k+1)+nf(1+i,2,k+2))*ufp1(1+
     1         i,2,k+1,3)
         end do
         do j = 1, 8
            do i = 1, 77
               r1 = 0.5*(nf(1+i,j*2+1,k+1)+nf(2+i,j*2+1,k+1))*ufp1(1+i,j
     1            *2+1,k+1,1)
               r2 = 0.5*(nf(1+i,(j+1)*2,k+1)+nf(2+i,(j+1)*2,k+1))*ufp1(1
     1            +i,(j+1)*2,k+1,1)
               flx(1+i,j*2+1,k+1,1) = r1
               flx(1+i,(j+1)*2,k+1,1) = r2
               r3 = 0.5*(nf(1+i,j*2+1,k+1)+nf(1+i,(j+1)*2,k+1))*ufp1(1+i
     1            ,j*2+1,k+1,2)
               r4 = 0.5*(nf(1+i,(j+1)*2,k+1)+nf(1+i,j*2+3,k+1))*ufp1(1+i
     1            ,(j+1)*2,k+1,2)
               flx(1+i,j*2+1,k+1,2) = r3
               flx(1+i,(j+1)*2,k+1,2) = r4
               r5 = 0.5*(nf(1+i,j*2+1,k+1)+nf(1+i,j*2+1,k+2))*ufp1(1+i,j
     1            *2+1,k+1,3)
               r6 = 0.5*(nf(1+i,(j+1)*2,k+1)+nf(1+i,(j+1)*2,k+2))*ufp1(1
     1            +i,(j+1)*2,k+1,3)
               flx(1+i,j*2+1,k+1,3) = r5
               flx(1+i,(j+1)*2,k+1,3) = r6
            end do
         end do
      end do
      return 
      end 
      subroutine vp_trans_nf_laxwend22(vp_mythread, vp_nmthreads, nf1, 
     1   flx, nfp1)
      common /coords/qx, qy, qz, lambda, ri, rj, rk, dz_grid, dz_cell
      common /timestep/ntf, dtsub
      REAL dtsub, ntf
      INTEGER ri, rj, rk
      REAL nf1(79,19,37), flx(79,19,37,3), dz_grid(37), nfp1(79,19,37), 
     1   qx(79), qy(19), qz(37), lambda(37), dz_cell(37)
      INTEGER k, j, i, j9, j10, j11, vp_mythread, vp_nmthreads, j4, j5
 
      j9 = (35 + vp_nmthreads - 1)/vp_nmthreads
      j10 = vp_mythread*j9 + 1
      j11 = min(35 - j10 + 1,j9)
      do k = j10, j11 + j10 - 1
         do j = 1, 17
            do i = 1, 77
               nfp1(1+i,j+1,k+1) = nf1(1+i,j-1+2,k-1+2) + 0.002*(nf1(2+i
     1            ,j-1+2,k-1+2)+nf1(i,j-1+2,k-1+2)+nf1(1+i,j-1+3,k-1+2)+
     2            nf1(1+i,j,k-1+2)+nf1(1+i,j-1+2,k-1+3)+nf1(1+i,j-1+2,k)
     3            -6.0*nf1(1+i,j-1+2,k-1+2)) - (dtsub/10000.0)*(flx(1+i,
     4            j-1+2,k-1+2,1)-flx(i,j-1+2,k-1+2,1)) - (dtsub/10000.0)
     5            *(flx(1+i,j-1+2,k-1+2,2)-flx(1+i,j,k-1+2,2)) - (dtsub/
     6            dz_grid(k-1+2))*(flx(1+i,j-1+2,k-1+2,3)-flx(1+i,j-1+2,
     7            k,3))
            end do
         end do
      end do
      return 
      end 
      subroutine vp_trans_nf_laxwend23(vp_mythread, vp_nmthreads, nfp1)
      REAL nfp1(79,19,37)
      INTEGER j3, j2, j12, j13, j14, vp_mythread, vp_nmthreads, j4, j5
      j12 = (703 + vp_nmthreads - 1)/vp_nmthreads
      j13 = vp_mythread*j12 + 1
      j14 = min(703 - j13 + 1,j12)
      do j3 = j13, j14 + j13 - 1
         nfp1(78,j3,1) = 4E15
         nfp1(79,j3,1) = 4E15
      end do
      return 
      end 
c----------------------------------------------------------------------
 
 
c----------------------------------------------------------------------
      SUBROUTINE trans_pf_LaxWend1(pf,pf1,uf)
c----------------------------------------------------------------------
      PARAMETER (nx = 79, ny = 19, nz = 37)
 
c grid parameters
      PARAMETER (dx = 10000.0, dy = 10000.0)   !units in km
      PARAMETER (delz = 10000.0)          !dz at release coordinates
                                       !this is not the final dz
                                       !d_lamba determines dz
c time stepping parameters
      PARAMETER (tstart = 0.2)   !starting time
      PARAMETER (dt = 1.)     !main time step
      PARAMETER (nt = 5.0)        !number of time steps
c      PARAMETER (dtsub = 0.01) !subcycle time step
      PARAMETER (dtsub_init = 0.1) !subcycle time step
      PARAMETER (ntsub = 10.0)   !number of subcycle time steps
      PARAMETER (nout = 50)      !number of steps to diagnosic output
 
c logical variable for restart
C...Translated by Crescent Bay Software VAST-F    7.0K8 15:17:00   4/ 1/04    
      integer j1, j2, j3
      logical restart
      PARAMETER (restart = .false.)
      PARAMETER (mrestart = 500)      ! use -1 for no save
      PARAMETER (mbegin = 0)      !mstart
 
c inital ambient density at the release point.
      real nf_init,b0_init,nn_coef
      PARAMETER (nf_init = 4e15)   !was 0.01
      PARAMETER (b0_init = 2e-9)    !was 0.2
      PARAMETER (nn_coef = 1e5)
 
c electron ion collision frequency
      real nu_init, eta_init
      PARAMETER (nu_init = 0.001)
      PARAMETER (eta_init = 0.0)
 
c temporary directory for output
      character*9 temp_dir
      PARAMETER (temp_dir = '/tmp/pad/')
 
c total number of neutrals released
      real No
      PARAMETER (No=1.0e28)
 
c neutral cloud expansion characteristics
      PARAMETER(vo = 20.0)
      PARAMETER(vth = 19.3)
      PARAMETER(vsat = 0.0)
      PARAMETER(vsw = 340.0)
 
c total kinetic energy of released neutrals
      PARAMETER(Etot = 2.696e7)       !units of kg m^2 / s^2
 
c density scaling parameter, alpha, and ion particle array dims
 
      real alpha
 
c      PARAMETER (beta =  8.5227271e-18)
      PARAMETER (alpha = 1.9263418e-20) !mH
 
c max number of ion particles to be produced.  This parameter
c is used to dimension the particle arrays.
      integer*4 Ni_max
      PARAMETER (Ni_max = 10000)
 
c earth's magnetic dipole field geometry
c earth radius and equatorial B alt (defining L shell)
      PARAMETER (r_earth=6.38e3)                !earth radius
      PARAMETER (B_max_alt=2.0e3+r_earth)       !ro defining L shell
 
c misc constants
      real mu0,epsilon,pi,rtod,q,mO,mBa,O_to_Ba,km_to_m,kboltz,melec
      real mproton,tempf0
      PARAMETER (pi = 3.14159)
      PARAMETER (rtod = 180.0/pi)  !radians to degreesc
      PARAMETER (mu0 = pi*4.0e-7)  !magnetic permeability of free space
      PARAMETER (epsilon = 8.85e-12) !dielectric constant
      PARAMETER (q = 1.6e-19)        !electron charge (Coulombs)
      PARAMETER (melec = 9.1e-31)  !mass of electron (kg)
      PARAMETER (mproton = 1.67e-27)
      PARAMETER (mO = 1.67e-27)    !mass of H (kg)
c      PARAMETER (mO = 2.3e-25)    !mass of Ba (kg)
      PARAMETER (mBa = 17.0*mO)    !mass of Ba (kg)
      PARAMETER (O_to_Ba = mO/mBa) !convert E and B for particle move
      PARAMETER (km_to_m = 1e3)    !km to meters conversion
      PARAMETER (kboltz = 1.38e-29)   !kg km^2 / s^2 / K
      PARAMETER (tempf0 = 50*11600.)     !K
 
      real Qo, vrad, N_o, Rp, tau_photo
      integer S_radius
      PARAMETER (Qo = 1.0e28)  !neutral source rate
      PARAMETER (vrad = 1.0)   !escape velocity
c      PARAMETER (N_o = 5e34)   !Steady state neutral particle constant
      PARAMETER (N_o = 1e33)   !Steady state neutral particle constant
      PARAMETER (Rp = 20000.0)  !Pluto radius
      PARAMETER (tau_photo = 1.2e6)
c      PARAMETER (dNi = 2500)
      PARAMETER (S_radius = 3) !units of dx
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
c raw grid coordinate data
 
      real qx(nx), qy(ny), qz(nz)
                            !cartisian coordinates of grid points
      real lambda(nz)       !latitude of given grid point
      integer ri, rj, rk    !release coordinates
      real dz_grid(nz)      !dz for grid, variable in z direction
      real dz_cell(nz)      !dz for cell, variable in z direction
 
      common /coords/ qx,qy,qz,lambda,ri,rj,rk,dz_grid,dz_cell
 
 
c Total number of ions produced at a given time during the simulation
      integer*4 Ni_tot, Nn_tot
 
c Location (indices) of particles with the grid
      integer ijkp(Ni_max,3)
      logical ionized(Ni_max)
      common /ion_info/ Ni_tot, Nn_tot, ijkp, ionized
 
      real seed
      common /rndseed/ seed
 
c Indices for boundary cells
      integer ip,im,jp,jm,kp,km
 
c Total input energy (ions)
      real input_E,input_EeP,prev_Etot, bndry_Eflux
      real input_chex, input_bill
      common /energy/ input_E, input_EeP, prev_Etot, bndry_Eflux,
     x                input_chex, input_bill
 
c Momentum conservation
      real surf_tot(3),graduu_tot(3),ugradu_tot(3)
      common /momentum/ surf_tot,graduu_tot,ugradu_tot
 
c Dummy variables for doing vector operations and particle flux
c calculation....etc.
 
      real a(nx,ny,nz,3),c(nx,ny,nz,3),ct(nx,ny,nz,3)
      common /dummy/ a,c,ct
 
c Variable time step for fluid velocity update used in both
c particle update and fluid update
 
      real dtf,ntf,dtsub
      common /timestep/ ntf,dtsub
 
c Weight variables for trilinear interpolation
 
      real wght(Ni_max,8), wquad(Ni_max,3)
      common /weights/ wght, wquad
 
c variable for anomlous resistivity
 
      real etar(nx,ny,nz,3)
      real nuei
      common /resist/ etar,nuei
 
c variable for particle scale
 
      real beta, dNi
      common /scale/ beta, dNi
 
 
c mass array for multi-ion species
      real mrat,m_arr
      common /mass/ mrat(Ni_max),m_arr(Ni_max)
 
 
 
      real pf(nx,ny,nz),
     x     pf1(nx,ny,nz),
     x     uf(nx,ny,nz,3)
 
      real minnf,maxnf
      real t1(nx,ny,nz,3)
      real t1c(nx,ny,nz)
      real t2(nx,ny,nz)
 
      parameter(gamma = 5./3.)
 
      external vp_trans_pf_laxwend11
      external vp_trans_pf_laxwend12
      external vp_trans_pf_laxwend13
      integer j6, j7, j8, j9, j10, j11, j12, j13, j14
      real r1, r2, r3, r4, t1c1
      call periodic_scalar(pf1)
 
      call vp_pcall(vp_trans_pf_laxwend11,2,5,gamma,uf,pf1,t1,t2)
 
      call periodic(t1)
 
      call vp_pcall (vp_trans_pf_laxwend12, 2, 5, t1, t1c, pf1, t2, pf)
 
c      pf = abs(pf)
 
      call periodic_scalar(pf)
      call vp_pcall (vp_trans_pf_laxwend13, 2, 2, tempf0, pf)
 
      return
      end
      subroutine vp_trans_pf_laxwend11(vp_mythread, vp_nmthreads, gamma
     1   , uf, pf1, t1, t2)
      common /coords/qx, qy, qz, lambda, ri, rj, rk, dz_grid, dz_cell
      REAL gamma
      INTEGER ri, rj, rk
      REAL uf(79,19,37,3), pf1(79,19,37), t1(79,19,37,3), dz_grid(37), 
     1   dz_cell(37), t2(79,19,37), qx(79), qy(19), qz(37), lambda(37)
      REAL r1, r2, r3, r4
      INTEGER k, j, i, j6, j7, j8, vp_mythread, vp_nmthreads, j4, j5
 
      j6 = (35 + vp_nmthreads - 1)/vp_nmthreads
      j7 = vp_mythread*j6 + 1
      j8 = min(35 - j7 + 1,j6)
      do k = j7, j8 + j7 - 1
         do j = 1, 17
            r1 = 1./10000.0
            r2 = 1./10000.0
            r3 = 1./dz_grid(k+1)
            r4 = 1./10000.0
            do i = 1, 77
               t1(1+i,j+1,k+1,1) = uf(1+i,j+1,k+1,1)*(pf1(2+i,j+1,k+1)-
     1            pf1(1+i,j+1,k+1))*r1
               t1(1+i,j+1,k+1,2) = uf(1+i,j+1,k+1,2)*(pf1(1+i,j+2,k+1)-
     1            pf1(1+i,j+1,k+1))*r2
               t1(1+i,j+1,k+1,3) = uf(1+i,j+1,k+1,3)*(pf1(1+i,j+1,k+2)-
     1            pf1(1+i,j+1,k+1))*r3
               t2(1+i,j+1,k+1) = gamma*pf1(1+i,j-1+2,k-1+2)*((uf(1+i,j-1
     1            +2,k-1+2,1)-uf(i,j-1+2,k-1+2,1))*r4+(uf(1+i,j-1+2,k-1+
     2            2,2)-uf(1+i,j,k-1+2,2))*r4+(uf(1+i,j-1+2,k-1+2,3)-uf(1
     3            +i,j-1+2,k,3))/dz_cell(k-1+2))
            end do
         end do
      end do
      return 
      end 
      subroutine vp_trans_pf_laxwend12(vp_mythread, vp_nmthreads, t1, 
     1   t1c, pf1, t2, pf)
      common /timestep/ntf, dtsub
      REAL dtsub, ntf
      REAL t1(79,19,37,3), t1c(79,19,37), pf1(79,19,37), t2(79,19,37), 
     1   pf(79,19,37)
      REAL t1c1
      INTEGER k, j, i, j9, j10, j11, vp_mythread, vp_nmthreads, j4, j5
 
      j9 = (35 + vp_nmthreads - 1)/vp_nmthreads
      j10 = vp_mythread*j9 + 1
      j11 = min(35 - j10 + 1,j9)
      do k = j10, j11 + j10 - 1
         do j = 1, 17
            do i = 1, 77
               t1c1 = (1./2.)*(t1(1+i,j-1+2,k-1+2,1)+t1(i,j-1+2,k-1+2,1)
     1            ) + (1./2.)*(t1(1+i,j-1+2,k-1+2,2)+t1(1+i,j,k-1+2,2))
     2             + (1./2.)*(t1(1+i,j-1+2,k-1+2,3)+t1(1+i,j-1+2,k,3))
               pf(1+i,j+1,k+1) = (1.0/12.0)*(pf1(2+i,j-1+2,k-1+2)+pf1(i,
     1            j-1+2,k-1+2)+pf1(1+i,j-1+3,k-1+2)+pf1(1+i,j,k-1+2)+pf1
     2            (1+i,j-1+2,k-1+3)+pf1(1+i,j-1+2,k)+6.0*pf1(1+i,j-1+2,k
     3            -1+2)) - 0.5*dtsub*(t1c1 + t2(1+i,j-1+2,k-1+2))
            end do
 
         end do
      end do
      return 
      end 
      subroutine vp_trans_pf_laxwend13(vp_mythread, vp_nmthreads, tempf0
     1   , pf)
      REAL tempf0
      REAL pf(79,19,37)
      INTEGER j3, j2, j12, j13, j14, vp_mythread, vp_nmthreads, j4, j5
      j12 = (703 + vp_nmthreads - 1)/vp_nmthreads
      j13 = vp_mythread*j12 + 1
      j14 = min(703 - j13 + 1,j12)
      do j3 = j13, j14 + j13 - 1
         pf(78,j3,1) = ((4E15*1.38E-29)*tempf0)
         pf(79,j3,1) = ((4E15*1.38E-29)*tempf0)
      end do
      return 
      end 
c----------------------------------------------------------------------
 
 
c----------------------------------------------------------------------
      SUBROUTINE trans_pf_LaxWend2(pf,pf1,ufp1)
c----------------------------------------------------------------------
      PARAMETER (nx = 79, ny = 19, nz = 37)
 
c grid parameters
      PARAMETER (dx = 10000.0, dy = 10000.0)   !units in km
      PARAMETER (delz = 10000.0)          !dz at release coordinates
                                       !this is not the final dz
                                       !d_lamba determines dz
c time stepping parameters
      PARAMETER (tstart = 0.2)   !starting time
      PARAMETER (dt = 1.)     !main time step
      PARAMETER (nt = 5.0)        !number of time steps
c      PARAMETER (dtsub = 0.01) !subcycle time step
      PARAMETER (dtsub_init = 0.1) !subcycle time step
      PARAMETER (ntsub = 10.0)   !number of subcycle time steps
      PARAMETER (nout = 50)      !number of steps to diagnosic output
 
c logical variable for restart
C...Translated by Crescent Bay Software VAST-F    7.0K8 15:17:00   4/ 1/04    
      integer j1, j2, j3
      logical restart
      PARAMETER (restart = .false.)
      PARAMETER (mrestart = 500)      ! use -1 for no save
      PARAMETER (mbegin = 0)      !mstart
 
c inital ambient density at the release point.
      real nf_init,b0_init,nn_coef
      PARAMETER (nf_init = 4e15)   !was 0.01
      PARAMETER (b0_init = 2e-9)    !was 0.2
      PARAMETER (nn_coef = 1e5)
 
c electron ion collision frequency
      real nu_init, eta_init
      PARAMETER (nu_init = 0.001)
      PARAMETER (eta_init = 0.0)
 
c temporary directory for output
      character*9 temp_dir
      PARAMETER (temp_dir = '/tmp/pad/')
 
c total number of neutrals released
      real No
      PARAMETER (No=1.0e28)
 
c neutral cloud expansion characteristics
      PARAMETER(vo = 20.0)
      PARAMETER(vth = 19.3)
      PARAMETER(vsat = 0.0)
      PARAMETER(vsw = 340.0)
 
c total kinetic energy of released neutrals
      PARAMETER(Etot = 2.696e7)       !units of kg m^2 / s^2
 
c density scaling parameter, alpha, and ion particle array dims
 
      real alpha
 
c      PARAMETER (beta =  8.5227271e-18)
      PARAMETER (alpha = 1.9263418e-20) !mH
 
c max number of ion particles to be produced.  This parameter
c is used to dimension the particle arrays.
      integer*4 Ni_max
      PARAMETER (Ni_max = 10000)
 
c earth's magnetic dipole field geometry
c earth radius and equatorial B alt (defining L shell)
      PARAMETER (r_earth=6.38e3)                !earth radius
      PARAMETER (B_max_alt=2.0e3+r_earth)       !ro defining L shell
 
c misc constants
      real mu0,epsilon,pi,rtod,q,mO,mBa,O_to_Ba,km_to_m,kboltz,melec
      real mproton,tempf0
      PARAMETER (pi = 3.14159)
      PARAMETER (rtod = 180.0/pi)  !radians to degreesc
      PARAMETER (mu0 = pi*4.0e-7)  !magnetic permeability of free space
      PARAMETER (epsilon = 8.85e-12) !dielectric constant
      PARAMETER (q = 1.6e-19)        !electron charge (Coulombs)
      PARAMETER (melec = 9.1e-31)  !mass of electron (kg)
      PARAMETER (mproton = 1.67e-27)
      PARAMETER (mO = 1.67e-27)    !mass of H (kg)
c      PARAMETER (mO = 2.3e-25)    !mass of Ba (kg)
      PARAMETER (mBa = 17.0*mO)    !mass of Ba (kg)
      PARAMETER (O_to_Ba = mO/mBa) !convert E and B for particle move
      PARAMETER (km_to_m = 1e3)    !km to meters conversion
      PARAMETER (kboltz = 1.38e-29)   !kg km^2 / s^2 / K
      PARAMETER (tempf0 = 50*11600.)     !K
 
      real Qo, vrad, N_o, Rp, tau_photo
      integer S_radius
      PARAMETER (Qo = 1.0e28)  !neutral source rate
      PARAMETER (vrad = 1.0)   !escape velocity
c      PARAMETER (N_o = 5e34)   !Steady state neutral particle constant
      PARAMETER (N_o = 1e33)   !Steady state neutral particle constant
      PARAMETER (Rp = 20000.0)  !Pluto radius
      PARAMETER (tau_photo = 1.2e6)
c      PARAMETER (dNi = 2500)
      PARAMETER (S_radius = 3) !units of dx
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
c raw grid coordinate data
 
      real qx(nx), qy(ny), qz(nz)
                            !cartisian coordinates of grid points
      real lambda(nz)       !latitude of given grid point
      integer ri, rj, rk    !release coordinates
      real dz_grid(nz)      !dz for grid, variable in z direction
      real dz_cell(nz)      !dz for cell, variable in z direction
 
      common /coords/ qx,qy,qz,lambda,ri,rj,rk,dz_grid,dz_cell
 
 
c Total number of ions produced at a given time during the simulation
      integer*4 Ni_tot, Nn_tot
 
c Location (indices) of particles with the grid
      integer ijkp(Ni_max,3)
      logical ionized(Ni_max)
      common /ion_info/ Ni_tot, Nn_tot, ijkp, ionized
 
      real seed
      common /rndseed/ seed
 
c Indices for boundary cells
      integer ip,im,jp,jm,kp,km
 
c Total input energy (ions)
      real input_E,input_EeP,prev_Etot, bndry_Eflux
      real input_chex, input_bill
      common /energy/ input_E, input_EeP, prev_Etot, bndry_Eflux,
     x                input_chex, input_bill
 
c Momentum conservation
      real surf_tot(3),graduu_tot(3),ugradu_tot(3)
      common /momentum/ surf_tot,graduu_tot,ugradu_tot
 
c Dummy variables for doing vector operations and particle flux
c calculation....etc.
 
      real a(nx,ny,nz,3),c(nx,ny,nz,3),ct(nx,ny,nz,3)
      common /dummy/ a,c,ct
 
c Variable time step for fluid velocity update used in both
c particle update and fluid update
 
      real dtf,ntf,dtsub
      common /timestep/ ntf,dtsub
 
c Weight variables for trilinear interpolation
 
      real wght(Ni_max,8), wquad(Ni_max,3)
      common /weights/ wght, wquad
 
c variable for anomlous resistivity
 
      real etar(nx,ny,nz,3)
      real nuei
      common /resist/ etar,nuei
 
c variable for particle scale
 
      real beta, dNi
      common /scale/ beta, dNi
 
 
c mass array for multi-ion species
      real mrat,m_arr
      common /mass/ mrat(Ni_max),m_arr(Ni_max)
 
 
 
      real pf(nx,ny,nz),
     x     pf1(nx,ny,nz),
     x     ufp1(nx,ny,nz,3)
 
      real minnf,maxnf
      real t1(nx,ny,nz,3)
      real t1c(nx,ny,nz)
      real t2(nx,ny,nz)
      real pfp1(nx,ny,nz)
 
      parameter(gamma = 5./3.)
 
c      call periodic_scalar(pf1)
 
      external vp_trans_pf_laxwend21
      external vp_trans_pf_laxwend22
      external vp_trans_pf_laxwend23
      integer j6, j7, j8, j9, j10, j11, j12, j13, j14, j15, j16, j17
      real r1, r2, r3, r4, t1c1
      call vp_pcall(vp_trans_pf_laxwend21,2,5,gamma,ufp1,pf,t1,t2)
 
      call periodic(t1)
 
      call vp_pcall(vp_trans_pf_laxwend22,2,5,pfp1,pf1,t1,t1c,t2)
c      pf1 = abs(pf1)
 
      call periodic_scalar(pf1)
      call vp_pcall (vp_trans_pf_laxwend23, 2, 2, tempf0, pf1)
 
      return
      end
      subroutine vp_trans_pf_laxwend21(vp_mythread, vp_nmthreads, gamma
     1   , ufp1, pf, t1, t2)
      common /coords/qx, qy, qz, lambda, ri, rj, rk, dz_grid, dz_cell
      REAL gamma
      INTEGER ri, rj, rk
      REAL ufp1(79,19,37,3), pf(79,19,37), t1(79,19,37,3), dz_grid(37), 
     1   dz_cell(37), t2(79,19,37), qx(79), qy(19), qz(37), lambda(37)
      REAL r1, r2, r3, r4
      INTEGER k, j, i, j6, j7, j8, vp_mythread, vp_nmthreads, j4, j5
      j6 = (35 + vp_nmthreads - 1)/vp_nmthreads
      j7 = vp_mythread*j6 + 1
      j8 = min(35 - j7 + 1,j6)
      do k = j7, j8 + j7 - 1
         do j = 1, 17
            r1 = 1./10000.0
            r2 = 1./10000.0
            r3 = 1./dz_grid(k+1)
            r4 = 1./10000.0
            do i = 1, 77
               t1(1+i,j+1,k+1,1) = ufp1(1+i,j+1,k+1,1)*(pf(2+i,j+1,k+1)-
     1            pf(1+i,j+1,k+1))*r1
               t1(1+i,j+1,k+1,2) = ufp1(1+i,j+1,k+1,2)*(pf(1+i,j+2,k+1)-
     1            pf(1+i,j+1,k+1))*r2
               t1(1+i,j+1,k+1,3) = ufp1(1+i,j+1,k+1,3)*(pf(1+i,j+1,k+2)-
     1            pf(1+i,j+1,k+1))*r3
               t2(1+i,j+1,k+1) = gamma*pf(1+i,j-1+2,k-1+2)*((ufp1(1+i,j-
     1            1+2,k-1+2,1)-ufp1(i,j-1+2,k-1+2,1))*r4+(ufp1(1+i,j-1+2
     2            ,k-1+2,2)-ufp1(1+i,j,k-1+2,2))*r4+(ufp1(1+i,j-1+2,k-1+
     3            2,3)-ufp1(1+i,j-1+2,k,3))/dz_cell(k-1+2))
            end do
         end do
      end do
      return 
      end 
      subroutine vp_trans_pf_laxwend22(vp_mythread, vp_nmthreads, pfp1, 
     1   pf1, t1, t1c, t2)
      common /timestep/ntf, dtsub
      REAL dtsub, ntf
      REAL pfp1(79,19,37), pf1(79,19,37), t1(79,19,37,3), t1c(79,19,37)
     1   , t2(79,19,37)
      REAL t1c1
      INTEGER j3, j2, j1, k, j, i, j9, j10, j11, j12, j13, j14, 
     1   vp_mythread, vp_nmthreads, j4, j5
 
      j9 = (35 + vp_nmthreads - 1)/vp_nmthreads
      j10 = vp_mythread*j9 + 1
      j11 = min(35 - j10 + 1,j9)
      do k = j10, j11 + j10 - 1
         do j = 1, 17
            do i = 1, 77
               t1c1 = (1./2.)*(t1(1+i,j-1+2,k-1+2,1)+t1(i,j-1+2,k-1+2,1)
     1            ) + (1./2.)*(t1(1+i,j-1+2,k-1+2,2)+t1(1+i,j,k-1+2,2))
     2             + (1./2.)*(t1(1+i,j-1+2,k-1+2,3)+t1(1+i,j-1+2,k,3))
               pfp1(1+i,j+1,k+1) = 0.002*(pf1(2+i,j-1+2,k-1+2)+pf1(i,j-1
     1            +2,k-1+2)+pf1(1+i,j-1+3,k-1+2)+pf1(1+i,j,k-1+2)+pf1(1+
     2            i,j-1+2,k-1+3)+pf1(1+i,j-1+2,k)-6.0*pf1(1+i,j-1+2,k-1+
     3            2)) + pf1(1+i,j-1+2,k-1+2) - dtsub*(t1c1 + t2(1+i,j-1+
     4            2,k-1+2))
            end do
         end do
      end do
 
      call vp_barrier (vp_mythread)
      j12 = (55537 + vp_nmthreads - 1)/vp_nmthreads
      j13 = vp_mythread*j12 + 1
      j14 = min(55537 - j13 + 1,j12)
      do j3 = j13, j14 + j13 - 1
         pf1(j3,1,1) = pfp1(j3,1,1)
      end do
      return 
      end 
      subroutine vp_trans_pf_laxwend23(vp_mythread, vp_nmthreads, tempf0
     1   , pf1)
      REAL tempf0
      REAL pf1(79,19,37)
      INTEGER j3, j2, j15, j16, j17, vp_mythread, vp_nmthreads, j4, j5
      j15 = (703 + vp_nmthreads - 1)/vp_nmthreads
      j16 = vp_mythread*j15 + 1
      j17 = min(703 - j16 + 1,j15)
      do j3 = j16, j17 + j16 - 1
         pf1(78,j3,1) = ((4E15*1.38E-29)*tempf0)
         pf1(79,j3,1) = ((4E15*1.38E-29)*tempf0)
      end do
      return 
      end 
c----------------------------------------------------------------------
 
 
 
 
 
 
 
 
 
 
 
 
 
 
