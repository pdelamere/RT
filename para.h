c para.h
c contains simulation parameter list for a barium release

c simulation domain dimensions
      PARAMETER (nx = 209, ny = 3, nz = 459)

c grid parameters
      PARAMETER (dx = 2000.0, dy = 2000.0)   !units in km
      PARAMETER (delz = 2000.0)          !dz at release coordinates
                                       !this is not the final dz
                                       !d_lamba determines dz
c time stepping parameters
      PARAMETER (tstart = 0.2)   !starting timeOrton steadying force, then chaos takes over
      PARAMETER (dt = 0.5)     !main time step
      PARAMETER (nt = 12000)        !number of time steps
c      PARAMETER (dtsub = 0.3) !subcycle time step 
      PARAMETER (dtsub_init = 0.05) !subcycle time step 
      PARAMETER (ntsub = 10.0)   !number of subcycle time steps
      PARAMETER (nout =200)      !number of steps to diagnosic output 

c logical variable for restart
      logical restart
      PARAMETER (restart = .false.)
      PARAMETER (mrestart = 6000)      ! use -1 for no save
      PARAMETER (mbegin = 0)      !mstart


c output directory
      character(50) out_dir
      PARAMETER (out_dir='/Volumes/Scratch/hybrid/RT/run_test/')

c electron ion collision frequency
      real nu_init, eta_init
      PARAMETER (nu_init = 0.002)
      PARAMETER (eta_init = 0.0)

c temporary directory for output
      character*9 temp_dir
      PARAMETER (temp_dir = '/tmp/pad/')

c total number of neutrals released
      real No
      PARAMETER (No=1.0e28)

c neutral cloud expansion characteristics
      real vtop,vbottom,vth,vsw
      PARAMETER(vo = 20.0)
c      PARAMETER(vth = 10.0)
      PARAMETER(vsat = 0.0)
      PARAMETER(vsw = 0.0)
      PARAMETER(vtop = vsw)
      PARAMETER(vbottom = -vsw)

c total kinetic energy of released neutrals 
      PARAMETER(Etot = 2.696e7)       !units of kg m^2 / s^2

c density scaling parameter, alpha, and ion particle array dims
       

c max number of ion particles to be produced.  This parameter
c is used to dimension the particle arrays.
      integer*4 Ni_max
      PARAMETER (Ni_max = 2000000)

c earths magnetic dipole field geometry
c earth radius and equatorial B alt (defining L shell)
      PARAMETER (r_earth=6.38e3)                !earth radius
      PARAMETER (B_max_alt=2.0e3+r_earth)       !ro defining L shell
 
c misc constants
      real mu0,epsilon,pi,rtod,mO,mBa,O_to_Ba,km_to_m,kboltz,melec
      real mproton,tempf0,m_pu
      
      real*8 q

      PARAMETER (pi = 3.14159)
      PARAMETER (rtod = 180.0/pi)  !radians to degreesc
      PARAMETER (mu0 = pi*4.0e-7)  !magnetic permeability of free space
      PARAMETER (epsilon = 8.85e-12) !dielectric constant
      PARAMETER (q = 1.6e-19)        !electron charge (Coulombs)
      PARAMETER (melec = 9.1e-31)  !mass of electron (kg)
      PARAMETER (mproton = 16*1.67e-27)
      PARAMETER (mO = mproton)    !mass of H (kg)
c      PARAMETER (mO = 2.3e-25)    !mass of Ba (kg)
      PARAMETER (m_pu = 1.0)
      PARAMETER (mBa = m_pu*mO)    !mass of Ba (kg)
      PARAMETER (O_to_Ba = mO/mBa) !convert E and B for particle move
      PARAMETER (km_to_m = 1e3)    !km to meters conversion
      PARAMETER (kboltz = 1.38e-29)   !kg km^2 / s^2 / K
      PARAMETER (tempf0 = 50*11600.)     !K

      real*8 alpha

c      PARAMETER (beta =  8.5227271e-18)
c      PARAMETER (alpha = 1.9263418e-20) !mH
      PARAMETER (alpha = pi*4e-10*q*q/mO) !mH


      real Qo, vrad, N_o, Rp, tau_photo
      integer S_radius
      PARAMETER (Qo = 2e27)  !neutral source rate
      PARAMETER (vrad = 0.1)   !escape velocity
c      PARAMETER (N_o = 5e34)   !Steady state neutral particle constant
      PARAMETER (N_o = 1e33)   !Steady state neutral particle constant
      PARAMETER (Rp = 1200.0)  !Pluto radius
      PARAMETER (tau_photo = 1.2e9)
c      PARAMETER (dNi = 2500)
      PARAMETER (S_radius = 30) !units of dx


c inital ambient density at the release point.
      real nf_init,b0_init,nn_coef,np_top,np_bottom
      real b0_top,b0_bottom,Lo,vth_top,vth_bottom,vth_max
      real m_top, m_bottom,m_heavy,np_bottom_proton
      real vth_heavy
      PARAMETER (nf_init = 0.2e15)   !was 0.01
      PARAMETER (m_heavy = 16.0)
      PARAMETER (np_top = nf_init/2.0)
      PARAMETER (np_bottom = nf_init)
      PARAMETER (f_proton_top = 0.8) !fraction relative to top
      PARAMETER (b0_init = 5e-9)    !was 0.2
      PARAMETER (b0_top = 1.0*b0_init)
      PARAMETER (b0_bottom = b0_init)
      PARAMETER (vth_top = 200.00)
      PARAMETER (vth_bottom = 200.00)
      PARAMETER (vth_heavy = 150.0)
      PARAMETER (vth_max = 20*200.0)
      PARAMETER (m_top = mproton)
      PARAMETER (m_bottom = mproton)
      PARAMETER (nn_coef = 1e5)
      PARAMETER (Lo = 5.0*dx) !gradient scale length of boundary


      real grav0
      PARAMETER (grav0 = 0.1*vth_bottom**2/Lo) !km/s^2












