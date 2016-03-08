c----------------------------------------------------------------------
      SUBROUTINE f_update_tlev(uf,uf2,b1,b12,b1p2,bt,b0,bdp)
c loops run 1 to n since values are only being copied
c----------------------------------------------------------------------
      include 'incurv.h'
      
      real uf(nx,ny,nz,3),
     x     uf1(nx,ny,nz,3),
     x     b1(nx,ny,nz,3),
     x     b12(nx,ny,nz,3),
     x     b1p2(nx,ny,nz,3),
     x     bt(nx,ny,nz,3),
     x     b0(nz),
     x     bdp(nx,ny,nz,3)

      do 10 i=1,nx
         do 10 j=1,ny
            do 10 k=1,nz
               bt(i,j,k,1) = b1p2(i,j,k,1) + bdp(i,j,k,1)
               bt(i,j,k,2) = b1p2(i,j,k,2) + bdp(i,j,k,2)
               bt(i,j,k,3) = b1p2(i,j,k,3) + bdp(i,j,k,3) + b0(k)
               do 10 m=1,3
c                  uf2(i,j,k,m) = uf(i,j,k,m)
                  b12(i,j,k,m) = b1(i,j,k,m)
                  b1(i,j,k,m) = b1p2(i,j,k,m)
 10               continue
      
      return
      end
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
      include 'incurv.h'

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

      call periodic(aa)
      call periodic(bbmf)

      do 5 k=1,nz
         zfrc(k) = 0.5*dz_grid(k)/dz_cell(k)
 5       continue

      do 10 i=3,nx      
         do 10 j=3,ny
            do 10 k=3,nz

               im = i-1         !assume daa/dxyz = 0 at boundary
               jm = j-1         !bbmf is given on boundary
               km = k-1

               ax = 0.5*(aa(i,j,k,1) + aa(im,j,k,1))
               bx = 0.5*(bbmf(i,j,k,1) + bbmf(im,j,k,1))

               ay = 0.5*(aa(i,j,k,2) + aa(i,jm,k,2))
               by = 0.5*(bbmf(i,j,k,2) + bbmf(i,jm,k,2))

               az = zfrc(k)*(aa(i,j,k,3) - aa(i,j,km,3)) + aa(i,j,km,3)
               bz = zfrc(k)*(bbmf(i,j,k,3) - bbmf(i,j,km,3))
     x                     + bbmf(i,j,km,3)

               ct(i,j,k,1) = ay*bz - az*by
               ct(i,j,k,2) = az*bx - ax*bz
               ct(i,j,k,3) = ax*by - ay*bx

 10            continue

      call periodic(ct)

c do x side i=2

      i=2
      do 20 j=3,ny
         do 20 k=3,nz

            jm = j-1
            km = k-1

            ax = 2.0*aa(2,j,k,1) - 0.5*(aa(3,j,k,1) + aa(2,j,k,1))
c            bx = 2.0*bbmf(2,j,k,1) - 
c     x           0.5*(bbmf(3,j,k,1) + bbmf(2,j,k,1))
            bx = bbmf(2,j,k,1)

            ay = 0.5*(aa(i,j,k,2) + aa(i,jm,k,2))
            by = 0.5*(bbmf(i,j,k,2) + bbmf(i,jm,k,2))

c            az = 0.5*(aa(i,j,k,3) + aa(i,j,km,3))
c            bz = 0.5*(bbmf(i,j,k,3) + bbmf(i,j,km,3))

            az = zfrc(k)*(aa(i,j,k,3) - aa(i,j,km,3)) + aa(i,j,km,3)
            bz = zfrc(k)*(bbmf(i,j,k,3) - bbmf(i,j,km,3))
     x                     + bbmf(i,j,km,3)

            ct(i,j,k,1) = ay*bz - az*by
            ct(i,j,k,2) = az*bx - ax*bz
            ct(i,j,k,3) = ax*by - ay*bx

 20         continue

c do x edge j,k=2

      j=2
      k=2
      do 25 i=3,nx

         im = i-1

         ax = 0.5*(aa(i,j,k,1) + aa(im,j,k,1))
         bx = 0.5*(bbmf(i,j,k,1) + bbmf(im,j,k,1))

         ay = 2.0*aa(i,2,k,2) - 0.5*(aa(i,3,k,2) + aa(i,2,k,2))
c         by = 2.0*bbmf(i,2,k,2) - 0.5*(bbmf(i,3,k,2) + bbmf(i,2,k,2))
         by = bbmf(i,2,k,2)

         az = 2.0*aa(i,j,2,3) - 
     x        (zfrc(k)*(aa(i,j,3,3) - aa(i,j,2,3)) + aa(i,j,2,3))
c         bz = 2.0*bbmf(i,j,2,3) - 0.5*(bbmf(i,j,3,3) + bbmf(i,j,2,3))
         bz = bbmf(i,j,2,3)

         ct(i,j,k,1) = ay*bz - az*by
         ct(i,j,k,2) = az*bx - ax*bz
         ct(i,j,k,3) = ax*by - ay*bx

 25      continue

c do y side j=2
      
      j=2
      do 30 i=3,nx
         do 30 k=3,nz

            im = i-1
            km = k-1

            ax = 0.5*(aa(i,j,k,1) + aa(im,j,k,1))
            bx = 0.5*(bbmf(i,j,k,1) + bbmf(im,j,k,1))

            ay = 2.0*aa(i,2,k,2) - 0.5*(aa(i,3,k,2) + aa(i,2,k,2))
c            by = 2.0*bbmf(i,2,k,2) - 
c     x           0.5*(bbmf(i,3,k,2) + bbmf(i,2,k,2))
            by = bbmf(i,2,k,2)

c            az = 0.5*(aa(i,j,k,3) + aa(i,j,km,3))
c            bz = 0.5*(bbmf(i,j,k,3) + bbmf(i,j,km,3))

            az = zfrc(k)*(aa(i,j,k,3) - aa(i,j,km,3)) + aa(i,j,km,3)
            bz = zfrc(k)*(bbmf(i,j,k,3) - bbmf(i,j,km,3))
     x                     + bbmf(i,j,km,3)

            ct(i,j,k,1) = ay*bz - az*by
            ct(i,j,k,2) = az*bx - ax*bz
            ct(i,j,k,3) = ax*by - ay*bx

 30         continue

c do y edge i,k=2

      i=2
      k=2
      do 35 j=3,ny

         jm = j-1

         ax = 2.0*aa(2,j,k,1) - 0.5*(aa(3,j,k,1) + aa(2,j,k,1))
c         bx = 2.0*bbmf(2,j,k,1) - 0.5*(bbmf(3,j,k,1) + bbmf(2,j,k,1))
         bx = bbmf(2,j,k,1)

         ay = 0.5*(aa(i,j,k,2) + aa(i,jm,k,2))
         by = 0.5*(bbmf(i,j,k,2) + bbmf(i,jm,k,2))

c         az = 2.0*aa(i,j,2,3) - 0.5*(aa(i,j,3,3) + aa(i,j,2,3))
         az = 2.0*aa(i,j,2,3) - 
     x        zfrc(k)*(aa(i,j,3,3) - aa(i,j,2,3)) + aa(i,j,2,3)
c         bz = 2.0*bbmf(i,j,2,3) - 0.5*(bbmf(i,j,3,3) + bbmf(i,j,2,3))
         bz = bbmf(i,j,2,3)

         ct(i,j,k,1) = ay*bz - az*by
         ct(i,j,k,2) = az*bx - ax*bz
         ct(i,j,k,3) = ax*by - ay*bx

 35      continue


c do z side k=2

      k=2
      do 40 i=3,nx
         do 40 j=3,ny

            im = i-1
            jm = j-1

            ax = 0.5*(aa(i,j,k,1) + aa(im,j,k,1))
            bx = 0.5*(bbmf(i,j,k,1) + bbmf(im,j,k,1))

            ay = 0.5*(aa(i,j,k,2) + aa(i,jm,k,2))
            by = 0.5*(bbmf(i,j,k,2) + bbmf(i,jm,k,2))

c            az = 2.0*aa(i,j,2,3) - 0.5*(aa(i,j,3,3) + aa(i,j,2,3))
            az = 2.0*aa(i,j,2,3) - 
     x           zfrc(k)*(aa(i,j,3,3) - aa(i,j,2,3)) + aa(i,j,2,3)
c            bz = 2.0*bbmf(i,j,2,3) - 
c     x           0.5*(bbmf(i,j,3,3) + bbmf(i,j,2,3))
            bz = bbmf(i,j,2,3)

            ct(i,j,k,1) = ay*bz - az*by
            ct(i,j,k,2) = az*bx - ax*bz
            ct(i,j,k,3) = ax*by - ay*bx

 40         continue

c do z edge i,j=2

      i=2
      j=2
      do 45 k=3,nz

         km=k-1

         ax = 2.0*aa(2,j,k,1) - 0.5*(aa(3,j,k,1) + aa(2,j,k,1))
c         bx = 2.0*bbmf(2,j,k,1) - 0.5*(bbmf(3,j,k,1) + bbmf(2,j,k,1))
         bx = bbmf(2,j,k,1)

         ay = 2.0*aa(i,2,k,2) - 0.5*(aa(i,3,k,2) + aa(i,2,k,2))
c         by = 2.0*bbmf(i,2,k,2) - 0.5*(bbmf(i,3,k,2) + bbmf(i,2,k,2))
         by = bbmf(i,2,k,2)

c         az = 0.5*(aa(i,j,k,3) + aa(i,j,km,3))
c         bz = 0.5*(bbmf(i,j,k,3) + bbmf(i,j,km,3)) 

          az = zfrc(k)*(aa(i,j,k,3) - aa(i,j,km,3)) + aa(i,j,km,3)
          bz = zfrc(k)*(bbmf(i,j,k,3) - bbmf(i,j,km,3))
     x                     + bbmf(i,j,km,3)

         ct(i,j,k,1) = ay*bz - az*by
         ct(i,j,k,2) = az*bx - ax*bz
         ct(i,j,k,3) = ax*by - ay*bx

 45      continue

c do corner i,j,k=2

      i=2
      j=2
      k=2
      ax = 2.0*aa(2,j,k,1) - 0.5*(aa(3,j,k,1) + aa(2,j,k,1))
c      bx = 2.0*bbmf(2,j,k,1) - 0.5*(bbmf(3,j,k,1) + bbmf(2,j,k,1))
      bx = bbmf(2,j,k,1)

      ay = 2.0*aa(i,2,k,2) - 0.5*(aa(i,3,k,2) + aa(i,2,k,2))
c      by = 2.0*bbmf(i,2,k,2) - 0.5*(bbmf(i,3,k,2) + bbmf(i,2,k,2))
      by = bbmf(i,2,k,2)

c      az = 2.0*aa(i,j,2,3) - 0.5*(aa(i,j,3,3) + aa(i,j,2,3))
      az = 2.0*aa(i,j,2,3) - 
     x     zfrc(k)*(aa(i,j,3,3) - aa(i,j,2,3)) + aa(i,j,2,3)
c      bz = 2.0*bbmf(i,j,2,3) - 0.5*(bbmf(i,j,3,3) + bbmf(i,j,2,3))
      bz = bbmf(i,j,2,3)

      ct(i,j,k,1) = ay*bz - az*by
      ct(i,j,k,2) = az*bx - ax*bz
      ct(i,j,k,3) = ax*by - ay*bx


c extrapolate back to main cell contravarient positions.
c ...just average across cells since cell edges are centered
c about the grid points.
      
      do 60 i=2,nx-1
         do 60 j=2,ny-1
            do 60 k=2,nz-1

               ip = i+1
               jp = j+1
               kp = k+1

               cc(i,j,k,1) = 0.5*(ct(i,j,k,1) + ct(ip,j,k,1))
               cc(i,j,k,2) = 0.5*(ct(i,j,k,2) + ct(i,jp,k,2))
               cc(i,j,k,3) = 0.5*(ct(i,j,k,3) + ct(i,j,kp,3))

 60            continue

c do x side, i=nx

      i=nx
      do 70 j=2,ny-1
         do 70 k=2,nz-1
            jp = j+1
            kp = k+1
c            cc(nx-1,j,k,1) = 2.0*ct(nx-1,j,k,1) - 
c     x                     0.5*(ct(nx-1,j,k,1) + ct(nx-2,j,k,1))
            cc(nx,j,k,1) = 2.0*ct(nx,j,k,1) - 
     x                     0.5*(ct(nx,j,k,1) + ct(nx-1,j,k,1))
            cc(i,j,k,2) = 0.5*(ct(i,j,k,2) + ct(i,jp,k,2))
            cc(i,j,k,3) = 0.5*(ct(i,j,k,3) + ct(i,j,kp,3))
 70         continue

c do x edge, j=ny,k=nz

      j=ny
      k=nz
      do 75 i=2,nx-1
         ip=i+1
         cc(i,j,k,1) = 0.5*(ct(i,j,k,1) + ct(ip,j,k,1))
c         cc(i,ny-1,k,2) = 2.0*ct(i,ny-1,k,2) - 
c     x                  0.5*(ct(i,ny-1,k,2) + ct(i,ny-2,k,2))
c         cc(i,j,nz-1,3) = cc(i,j,nz-2,3) +
c     x                  (2.0*dz_cell(nz-1)/dz_grid(nz-1))*
c     x                  (ct(i,j,nz-1,3) - cc(i,j,nz-2,3))
         cc(i,ny,k,2) = 2.0*ct(i,ny,k,2) - 
     x                  0.5*(ct(i,ny,k,2) + ct(i,ny-1,k,2))
         cc(i,j,nz,3) = cc(i,j,nz-1,3) +
     x                  (2.0*dz_cell(nz)/dz_grid(nz))*
     x                  (ct(i,j,nz,3) - cc(i,j,nz-1,3))
 75      continue

c do y side, j=ny

      j=ny
      do 80 i=2,nx-1
         do 80 k=2,nz-1
            ip = i+1
            kp = k+1
            cc(i,j,k,1) = 0.5*(ct(i,j,k,1) + ct(ip,j,k,1))
c            cc(i,ny-1,k,2) = 2.0*ct(i,ny-1,k,2) - 
c     x                     0.5*(ct(i,ny-1,k,2) + ct(i,ny-2,k,2))
            cc(i,ny,k,2) = 2.0*ct(i,ny,k,2) - 
     x                     0.5*(ct(i,ny,k,2) + ct(i,ny-1,k,2))
            cc(i,j,k,3) = 0.5*(ct(i,j,k,3) + ct(i,j,kp,3))
 80         continue
      
c do y egde, i=nx,k=nz

      i=nx
      k=nz
      do 85 j=2,ny-1
         jp=j+1
c         cc(nx-1,j,k,1) = 2.0*ct(nx-1,j,k,1) - 
c     x                  0.5*(ct(nx-1,j,k,1) + ct(nx-2,j,k,1))
         cc(nx,j,k,1) = 2.0*ct(nx,j,k,1) - 
     x                  0.5*(ct(nx,j,k,1) + ct(nx-1,j,k,1))
         cc(i,j,k,2) = 0.5*(ct(i,j,k,2) + ct(i,jp,k,2))
c         cc(i,j,nz-1,3) = cc(i,j,nz-2,3) +
c     x                  (2.0*dz_cell(nz-1)/dz_grid(nz-1))*
c     x                  (ct(i,j,nz-1,3) - cc(i,j,nz-2,3))
         cc(i,j,nz,3) = cc(i,j,nz-1,3) +
     x                  (2.0*dz_cell(nz)/dz_grid(nz))*
     x                  (ct(i,j,nz,3) - cc(i,j,nz-1,3))
 85      continue

c do z side, k=nz

      k=nz
      do 90 i=2,nx-1
         do 90 j=2,ny-1
            ip=i+1
            jp=j+1
            cc(i,j,k,1) = 0.5*(ct(i,j,k,1) + ct(ip,j,k,1))
            cc(i,j,k,2) = 0.5*(ct(i,j,k,2) + ct(i,jp,k,2))
c         cc(i,j,nz-1,3) = cc(i,j,nz-2,3) +
c     x                  (2.0*dz_cell(nz-1)/dz_grid(nz-1))*
c     x                  (ct(i,j,nz-1,3) - cc(i,j,nz-2,3))
         cc(i,j,nz,3) = cc(i,j,nz-1,3) +
     x                  (2.0*dz_cell(nz)/dz_grid(nz))*
     x                  (ct(i,j,nz,3) - cc(i,j,nz-1,3))
 90         continue

c do z edge, i=nx,j=ny

      i=nx
      j=ny
      do 95 k=2,nz-1
         kp=k+1
c         cc(nx-1,j,k,1) = 2.0*ct(nx-1,j,k,1) - 
c     x                  0.5*(ct(nx-1,j,k,1) + ct(nx-2,j,k,1))
c         cc(i,ny-1,k,2) = 2.0*ct(i,ny-1,k,2) - 
c     x                  0.5*(ct(i,ny-1,k,2) + ct(i,ny-2,k,2))
         cc(nx,j,k,1) = 2.0*ct(nx,j,k,1) - 
     x                  0.5*(ct(nx,j,k,1) + ct(nx-1,j,k,1))
         cc(i,ny,k,2) = 2.0*ct(i,ny,k,2) - 
     x                  0.5*(ct(i,ny,k,2) + ct(i,ny-1,k,2))
         cc(i,j,k,3) = 0.5*(ct(i,j,k,3) + ct(i,j,kp,3))
 95      continue

c do corner
      i=nx
      j=ny
      k=nz
c      cc(nx-1,j,k,1) = 2.0*ct(nx-1,j,k,1) - 
c     x               0.5*(ct(nx-1,j,k,1) + ct(nx-2,j,k,1))
c      cc(i,ny-1,k,2) = 2.0*ct(i,ny-1,k,2) - 
c     x               0.5*(ct(i,ny-1,k,2) + ct(i,ny-2,k,2))     
c      cc(i,j,nz,3) = 2.0*ct(i,j,nz,3) - 
c     x               0.5*(ct(i,j,nz,3) + ct(i,j,nz-1,3)) 
c      cc(i,j,nz-1,3) = cc(i,j,nz-2,3) +
c     x               (2.0*dz_cell(nz-1)/dz_grid(nz-1))*
c     x               (ct(i,j,nz-1,3) - cc(i,j,nz-2,3))
      cc(nx,j,k,1) = 2.0*ct(nx,j,k,1) - 
     x               0.5*(ct(nx,j,k,1) + ct(nx-1,j,k,1))
      cc(i,ny,k,2) = 2.0*ct(i,ny,k,2) - 
     x               0.5*(ct(i,ny,k,2) + ct(i,ny-1,k,2))     
c      cc(i,j,nz,3) = 2.0*ct(i,j,nz,3) - 
c     x               0.5*(ct(i,j,nz,3) + ct(i,j,nz-1,3)) 
      cc(i,j,nz,3) = cc(i,j,nz-1,3) +
     x               (2.0*dz_cell(nz)/dz_grid(nz))*
     x               (ct(i,j,nz,3) - cc(i,j,nz-1,3))

c      write(*,*) 'ct1....',ct(nx-1,rj,rk,1),ct(nx-2,rj,rk,1),
c     x                     ct(nx-3,rj,rk,1)
c      write(*,*) 'cc1....',cc(nx-1,rj,rk,1),cc(nx-2,rj,rk,1),
c     x                     cc(nx-3,rj,rk,1)

      call periodic(cc)

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
      include 'incurv.h'

      real bt(nx,ny,nz,3),   !main cell covarient
     x     btmf(nx,ny,nz,3)  !main cell contravarient

      real bx1, bx2, by1, by2, bz1, bz2  !main cell center fields
      real zrat           !ratio for doing linear interpolation
                          !to grid point position.
      real zplus, zminus  !position of main cell edges up and down
      real b_j, b_jm, b_i, b_im !intermediate step in average process

      do 10 i=2,nx-1
         do 10 j=2,ny-1
            do 10 k=2,nz-1

               ip = i+1
               jp = j+1
               kp = k+1
               im = i-1
               jm = j-1
               km = k-1

c The x component of B resides at the k and k-1 edges, so this
c requires the non-uniform grid interpolation

               zplus = (qz(k+1) + qz(k))/2.0
               zminus = (qz(k) + qz(k-1))/2.0
               zrat = (qz(k) - zminus)/(zplus - zminus)
   
               b_j = bt(i,j,km,1) 
     x               + zrat*(bt(i,j,k,1) - bt(i,j,km,1)) 
               b_jm = bt(i,jm,km,1)
     x                + zrat*(bt(i,jm,k,1) - bt(i,jm,km,1))
               bx1 = (b_j + b_jm)/2.0

               b_j = bt(ip,j,km,1) 
     x               + zrat*(bt(ip,j,k,1) - bt(ip,j,km,1)) 
               b_jm = bt(ip,jm,km,1)
     x                + zrat*(bt(ip,jm,k,1) - bt(ip,jm,km,1))
               bx2 = (b_j + b_jm)/2.0

               
               b_i = bt(i,j,km,2) 
     x               + zrat*(bt(i,j,k,2) - bt(i,j,km,2)) 
               b_im = bt(im,j,km,2)
     x                + zrat*(bt(im,j,k,2) - bt(im,j,km,2))           
               by1 = (b_i + b_im)/2.0

               b_i = bt(i,jp,km,2) 
     x               + zrat*(bt(i,jp,k,2) - bt(i,jp,km,2)) 
               b_im = bt(im,jp,km,2)
     x                + zrat*(bt(im,jp,k,2) - bt(im,jp,km,2))
               by2 = (b_i + b_im)/2.0


               bz1 = 0.25*(bt(i,j,k,3) + bt(i,jm,k,3) +
     x                     bt(im,jm,k,3) + bt(im,j,k,3))
               bz2 = 0.25*(bt(i,j,kp,3) + bt(i,jm,kp,3) +
     x                     bt(im,jm,kp,3) + bt(im,j,kp,3))

               btmf(i,j,k,1) = 0.5*(bx1+bx2)
               btmf(i,j,k,2) = 0.5*(by1+by2)
               btmf(i,j,k,3) = 0.5*(bz1+bz2)

 10            continue

c      call boundaries(btmf)
      call periodic(btmf)

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
      include 'incurv.h'

      real b1(nx,ny,nz,3),
     x     nf(nx,ny,nz),
     x     np(nx,ny,nz),
     x     aj(nx,ny,nz,3)

      real curl_B(3)      !dummy for holding curl vector
      real ntot(3)        !total density, np + nf

      call periodic_scalar(np)
      call periodic_scalar(nf)
      call periodic(b1)

      do 10 i=2,nx   
         do 10 j=2,ny
            do 10 k=2,nz

               ip = i+1
               jp = j+1
               kp = k+1

               if (ip .gt. nx) then ip = nx
               if (jp .gt. ny) then jp = ny
               if (kp .gt. nz) then kp = nz

               ntot(1) = 0.5*(nf(i,j,k)+nf(ip,j,k))
     x                 + 0.5*(np(i,j,k)+np(ip,j,k))
               ntot(2) = 0.5*(nf(i,j,k)+nf(i,jp,k))
     x                 + 0.5*(np(i,j,k)+np(i,jp,k))
               ntot(3) = 0.5*(nf(i,j,k)+nf(i,j,kp))
     x                 + 0.5*(np(i,j,k)+np(i,j,kp))

               curl_B(1) = (b1(i,j,k,3)/dy) - (b1(i,j-1,k,3)/dy) 
     x                    + (b1(i,j,k-1,2)/dz_cell(k))  
     x                    - (b1(i,j,k,2)/dz_cell(k))
               curl_B(2) = (b1(i,j,k,1)/dz_cell(k)) 
     x                     - (b1(i,j,k-1,1)/dz_cell(k))
     x                     - (b1(i,j,k,3)/dx) + (b1(i-1,j,k,3)/dx)
               curl_B(3) = (b1(i,j,k,2)/dx) - (b1(i-1,j,k,2)/dx) 
     x                     + (b1(i,j-1,k,1)/dy) - (b1(i,j,k,1)/dy)

               do 10 m=1,3
                  aj(i,j,k,m) = curl_B(m)/(ntot(m)*alpha)
 10            continue

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
      include 'incurv.h'

      real E(nx,ny,nz,3)      !E field, main cell contravarient
      real curl_E(nx,ny,nz,3) !curl of E, main cell covarient
      real lx, ly, lz         !lengths of dual cell edges

      call periodic(E)

      do 10 i=2,nx-1
         do 10 j=2,ny-1
            do 10 k=2,nz-1

               lx = qx(i+1) - qx(i)
               ly = qy(j+1) - qy(j)
               lz = qz(k+1) - qz(k)

               curl_E(i,j,k,1) =  (E(i,j+1,k,3)/ly) - (E(i,j,k,3)/ly)
     x                       + (E(i,j,k,2)/lz) - (E(i,j,k+1,2)/lz)
               curl_E(i,j,k,2) =  (E(i,j,k,3)/lx) - (E(i+1,j,k,3)/lx)
     x                       + (E(i,j,k+1,1)/lz) - (E(i,j,k,1)/lz)
               curl_E(i,j,k,3) =  (E(i,j,k,1)/ly) - (E(i,j+1,k,1)/ly)
     x                       + (E(i+1,j,k,2)/lx) - (E(i,j,k,2)/lx)

 10          continue

      call periodic(curl_E)

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE get_um_dot_BB(u,b,cc)
c uf and btmf are gathered at main cell center and uf.B*B 
c calculated.  Result returned to main cell contravarient
c postion.
c----------------------------------------------------------------------
      include 'incurv.h'

      real u(nx,ny,nz,3),  !main cell contravarient fluid velocity (uf)
     x     b(nx,ny,nz,3),   !main cell contravarient bt (btmf)
     x     cc(nx,ny,nz,3)   !(uf.B)*B

      real ux,uy,uz,bx,by,bz   !temp vars for u and b at cell center
      real temp                !used to vectorize loop
c      real ct(nx,ny,nz,3)      !result are main cell center
      real udotb               !u dot b
      real zfrc(nz)            !0.5*dz_grid(k)/dz_cell(k)

! first gather everything at center

      call periodic(u)
      call periodic(b)

      do 5 k=1,nz
         zfrc(k) = 0.5*dz_grid(k)/dz_cell(k)
 5       continue

      do 10 i=2,nx      
         do 10 j=2,ny
            do 10 k=2,nz

               im = i-1
               jm = j-1     
               km = k-1

               if (i .eq. 2) then 
                  ux = 2.0*u(2,j,k,1) - 
     x                 0.5*(u(3,j,k,1) + u(2,j,k,1))
                  bx = 2.0*b(2,j,k,1) - 
     x                 0.5*(b(3,j,k,1) + b(2,j,k,1))
               else
                  ux = 0.5*(u(i,j,k,1) + u(im,j,k,1))
                  bx = 0.5*(b(i,j,k,1) + b(im,j,k,1))
               endif

               if (j .eq. 2) then 
                  uy = 2.0*u(i,2,k,2) - 
     x                 0.5*(u(i,3,k,2) + u(i,2,k,2))
                  by = 2.0*b(i,2,k,2) - 
     x                 0.5*(b(i,3,k,2) + b(i,2,k,2))
               else
                  uy = 0.5*(u(i,j,k,2) + u(i,jm,k,2))
                  by = 0.5*(b(i,j,k,2) + b(i,jm,k,2))
               endif

               if (k .eq. 2) then
                  uz = 2.0*u(i,j,2,3) - 
     x                 zfrc(k)*(u(i,j,3,3) - u(i,j,2,3)) + u(i,j,2,3) 
c                  uz = 2.0*u(i,j,2,3) - 
c     x                 0.5*(u(i,j,3,3) + u(i,j,2,3))
                  bz = b(i,j,2,3)
c                  bz = 2.0*b(i,j,2,3) - 
c     x                 0.5*(b(i,j,3,3) + b(i,j,2,3))
               else
                  uz = zfrc(k)*(u(i,j,k,3) - u(i,j,km,3)) + u(i,j,km,3)
                  bz = zfrc(k)*(b(i,j,k,3) - b(i,j,km,3)) + b(i,j,km,3)
c                  uz = 0.5*(u(i,j,k,3) + u(i,j,km,3))
c                  bz = 0.5*(b(i,j,k,3) + b(i,j,km,3))            
               endif

               udotb = ux*bx + uy*by + uz*bz

               ct(i,j,k,1) = udotb*bx
               ct(i,j,k,2) = udotb*by
               ct(i,j,k,3) = udotb*bz

 10            continue

      call periodic(ct)

c extrapolate back to main cell contravarient positions.
c ...just average across cells.

      do 60 i=2,nx
         do 60 j=2,ny
            do 60 k=2,nz

               ip = i+1
               jp = j+1
               kp = k+1

               if (i .eq. nx-1) then 
                  cc(nx,j,k,1) = 2.0*ct(nx,j,k,1) - 
     x                           0.5*(ct(nx,j,k,1) + ct(nx-1,j,k,1))
               else
                  cc(i,j,k,1) = 0.5*(ct(i,j,k,1) + ct(ip,j,k,1))
               endif

               if (j .eq. ny-1) then 
                  cc(i,ny,k,2) = 2.0*ct(i,ny,k,2) - 
     x                           0.5*(ct(i,ny,k,2) + ct(i,ny-1,k,2))
               else
                  cc(i,j,k,2) = 0.5*(ct(i,j,k,2) + ct(i,jp,k,2))
               endif
                  
               if (k .eq. nz-1) then
c                  temp = 2.0*ct(i,j,nz,3) - 
c     x                           0.5*(ct(i,j,nz,3) + ct(i,j,nz-1,3))
                   temp = 0.5*(ct(i,j,nz,3) + ct(i,j,nz-1,3)) +
     x                   (2.0*dz_cell(nz)/dz_grid(nz))*(ct(i,j,nz,3) -
     x                    0.5*(ct(i,j,nz,3) + ct(i,j,nz-1,3)))
               else
                  temp = 0.5*(ct(i,j,k,3) + ct(i,j,kp,3))
               endif             
  
               cc(i,j,k,3) = temp

 60            continue

      call periodic(cc)

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE get_ugradu_Lax(uf,ugradu,delta_t)
c----------------------------------------------------------------------
      include 'incurv.h'

      real uf(nx,ny,nz,3),
     x     ugradu(nx,ny,nz,3)

      real ufc(nx,ny,nz,3)
      real ax1,ax2,ay1,ay2,az1,az2       !misc const
      real u1,u2,u3                      !temp vars

      parameter(ad = 0.0)                 !coefficient to add extra
                                          !diffusion
      call periodic(uf)   

      call face_to_center(uf,ufc)

      do 10 i=2,nx-1
         do 10 j=2,ny-1
            do 10 k=2,nz-1
               ip=i+1
               jp=j+1
               kp=k+1
               im=i-1
               jm=j-1
               km=k-1
               if (ip .gt. nx-1) then ip = nx-1
               if (jp .gt. ny-1) then jp = ny-1
               if (kp .gt. nz-1) then kp = nz-1
               if (im .eq. 1) then im = 2
               if (jm .eq. 1) then jm = 2
               if (km .eq. 1) then km = 2

c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

               ax1 = 0.5*ufc(i,j,k,1)/dx
               ax2 = ad*abs(ufc(ip,j,k,1) - ufc(im,j,k,1))
               ay1 = 0.5*ufc(i,j,k,2)/dy
               ay2 = ad*abs(ufc(ip,j,k,2) - ufc(im,j,k,2))
               u1 = ax1*(ufc(ip,j,k,1) - ufc(im,j,k,1)) - 
     x              ax2*(ufc(im,j,k,1) - 2.0*ufc(i,j,k,1) +
     x                      ufc(ip,j,k,1))
               u2 = ay1*(ufc(i,jp,k,1) - ufc(i,jm,k,1)) - 
     x              ay2*(ufc(i,jm,k,1) - 2.0*ufc(i,j,k,1) +
     x                      ufc(i,jp,k,1)) 
               az1 = ufc(i,j,k,3)/(dz_grid(k)+dz_grid(kp))
               az2 = ad*abs(ufc(ip,j,k,3) - ufc(im,j,k,3))
               u3 = az1*(ufc(i,j,kp,1)-ufc(i,j,km,1)) -
     x              az2*(ufc(i,j,km,1) - 2.0*ufc(i,j,k,1) +
     x                   ufc(i,j,kp,1))
               ct(i,j,k,1) = u1 + u2 + u3

c yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy

               ax1 = 0.5*ufc(i,j,k,1)/dx
               ax2 = ad*abs(ufc(i,jp,k,1) - ufc(i,jm,k,1))
               ay1 = 0.5*ufc(i,j,k,2)/dy
               ay2 = ad*abs(ufc(i,jp,k,2) - ufc(i,jm,k,2))
               u1 = ax1*(ufc(ip,j,k,2) - ufc(im,j,k,2)) - 
     x              ax2*(ufc(im,j,k,2) - 2.0*ufc(i,j,k,2) +
     x                      ufc(ip,j,k,2))
               u2 = ay1*(ufc(i,jp,k,2) - ufc(i,jm,k,2)) - 
     x              ay2*(ufc(i,jm,k,2) - 2.0*ufc(i,j,k,2) +
     x                      ufc(i,jp,k,2)) 
               az1 = ufc(i,j,k,3)/(dz_grid(k)+dz_grid(kp))
               az2 = ad*abs(ufc(i,jp,k,3) - ufc(i,jm,k,3))
               u3 = az1*(ufc(i,j,kp,2)-ufc(i,j,km,2)) -
     x              az2*(ufc(i,j,km,2) - 2.0*ufc(i,j,k,2) +
     x                   ufc(i,j,kp,2))
               ct(i,j,k,2) = u1 + u2 + u3

c zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
 
               ax1 = 0.5*ufc(i,j,k,1)/dx
               ax2 = ad*abs(ufc(i,j,kp,1) - ufc(i,j,km,1))
               ay1 = 0.5*ufc(i,j,k,2)/dy
               ay2 = ad*abs(ufc(i,j,kp,2) - ufc(i,j,km,2))
               u1 = ax1*(ufc(ip,j,k,3) - ufc(im,j,k,3)) - 
     x              ax2*(ufc(im,j,k,3) - 2.0*ufc(i,j,k,3) +
     x                      ufc(ip,j,k,3))
               u2 = ay1*(ufc(i,jp,k,3) - ufc(i,jm,k,3)) - 
     x              ay2*(ufc(i,jm,k,3) - 2.0*ufc(i,j,k,3) +
     x                      ufc(i,jp,k,3)) 
               az1 = ufc(i,j,k,3)/(dz_grid(k)+dz_grid(kp))
               az2 = ad*abs(ufc(i,j,kp,3) - ufc(i,j,km,3))
               u3 = az1*(ufc(i,j,kp,3)-ufc(i,j,km,3)) -
     x              az2*(ufc(i,j,km,3) - 2.0*ufc(i,j,k,3) +
     x                   ufc(i,j,kp,3))
               ct(i,j,k,3) = u1 + u2 + u3

 10            continue

      call periodic(ct)

c interpolate back to contravarient positions.

      do 20 i=2,nx-1
         do 20 j=2,ny-1
            do 20 k=2,nz-1
               ip = i+1
               jp = j+1
               kp = k+1
               if (ip .eq. nx) then ip = nx-1
               if (jp .eq. ny) then jp = ny-1
               if (kp .eq. nz) then kp = nz-1
               ugradu(i,j,k,1) = 0.5*(ct(i,j,k,1)+ct(ip,j,k,1))
               ugradu(i,j,k,2) = 0.5*(ct(i,j,k,2)+ct(i,jp,k,2))
               ugradu(i,j,k,3) = 0.5*(ct(i,j,k,3)+ct(i,j,kp,3))
 20         continue

      call periodic(ugradu)


      do 30 i=2,nx
         do 30 j=2,ny
            do 30 m=1,3
               ugradu(i,j,nz-1,m) = ugradu(i,j,nz-2,m)
               ugradu(i,j,nz,m) = ugradu(i,j,nz-2,m)
               ugradu(i,j,2,m) = ugradu(i,j,3,m)
 30         continue

      do 40 i=2,nx
         do 40 k=2,nz
            do 40 m=1,3
               ugradu(i,ny-1,k,m) = ugradu(i,ny-2,k,m)
               ugradu(i,ny,k,m) = ugradu(i,ny-2,k,m)
               ugradu(i,2,k,m) = ugradu(i,3,k,m)
 40         continue

      do 50 j=2,ny
         do 50 k=2,nz
            do 50 m=1,3
               ugradu(nx-1,j,k,m) = ugradu(nx-2,j,k,m)
               ugradu(nx,j,k,m) = ugradu(nx-2,j,k,m)
               ugradu(2,j,k,m) = ugradu(3,j,k,m)
 50         continue

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE get_ugradu_2nd_order(uf,ugradu)
c----------------------------------------------------------------------
      include 'incurv.h'

      real uf(nx,ny,nz,3),
     x     ugradu(nx,ny,nz,3)

      real ax,ay,az,az2       !misc const
      real u1,u2,u3           !temp vars

      do 10 i=2,nx
         do 10 j=2,ny
            do 10 k=2,nz
               ip=i+1
               jp=j+1
               kp=k+1
               im=i-1
               jm=j-1
               km=k-1
               if (ip .gt. nx) then ip = nx
               if (jp .gt. ny) then jp = ny
               if (kp .gt. nz) then kp = nz
               if (im .lt. 2) then im = 2
               if (jm .lt. 2) then jm = 2
               if (km .lt. 2) then km = 2

c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

               ax = 0.5*uf(i,j,k,1)/dx
               ay = 0.5*uf(i,j,k,2)/dy
               u1 = ax*(uf(ip,j,k,1) - uf(im,j,k,1)) - 
     x              0.5*ax*(uf(im,j,k,1) - 2.0*uf(i,j,k,1) +
     x                      uf(ip,j,k,1))
               u2 = ay*(uf(i,jp,k,1) - uf(i,jm,k,1)) - 
     x              0.5*ay*(uf(i,jm,k,1) - 2.0*uf(i,j,k,1) +
     x                      uf(i,jp,k,1)) 
               az = uf(i,j,k,3)/(dz_grid(kp))
               az2 = 0.5*uf(i,j,k,3)/(dz_grid(k) + dz_grid(kp))
               u3 = az2*(uf(i,j,kp,1)-uf(i,j,km,1)) -
     x              az2*(uf(i,j,km,1) - 2.0*uf(i,j,k,1) +
     x                   uf(i,j,kp,1))
               ugradu(i,j,k,1) = u1 + u2 + u3

c yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy
 
               ax = 0.5*uf(i,j,k,1)/dx
               ay = 0.5*uf(i,j,k,2)/dy
               u1 = ax*(uf(ip,j,k,2) - uf(im,j,k,2)) - 
     x              0.5*ax*(uf(im,j,k,2) - 2.0*uf(i,j,k,2) +
     x                      uf(ip,j,k,2))
               u2 = ay*(uf(i,jp,k,2) - uf(i,jm,k,2)) - 
     x              0.5*ay*(uf(i,jm,k,2) - 2.0*uf(i,j,k,2) +
     x                      uf(i,jp,k,2)) 
               az = uf(i,j,k,3)/(dz_grid(kp))
               az2 = 0.5*uf(i,j,k,3)/(dz_grid(k) + dz_grid(kp))
               u3 = az2*(uf(i,j,kp,2)-uf(i,j,km,2)) -
     x              az2*(uf(i,j,km,2) - 2.0*uf(i,j,k,2) +
     x                   uf(i,j,kp,2))
               ugradu(i,j,k,2) = u1 + u2 + u3

c zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
 
               ax = 0.5*uf(i,j,k,1)/dx
               ay = 0.5*uf(i,j,k,2)/dy
               u1 = ax*(uf(ip,j,k,3) - uf(im,j,k,3)) - 
     x              0.5*ax*(uf(im,j,k,3) - 2.0*uf(i,j,k,3) +
     x                      uf(ip,j,k,3))
               u2 = ay*(uf(i,jp,k,3) - uf(i,jm,k,3)) - 
     x              0.5*ay*(uf(i,jm,k,3) - 2.0*uf(i,j,k,3) +
     x                      uf(i,jp,k,3)) 
               az = uf(i,j,k,3)/(dz_grid(kp))
               az2 = 0.5*uf(i,j,k,3)/(dz_grid(k) + dz_grid(kp))
               u3 = az2*(uf(i,j,kp,3)-uf(i,j,km,3)) -
     x              az2*(uf(i,j,km,3) - 2.0*uf(i,j,k,3) +
     x                   uf(i,j,kp,3))
               ugradu(i,j,k,3) = u1 + u2 + u3

 10            continue

      write(*,*) 'ugradu...',ugradu(ri,rj,rk,1),ugradu(ri,rj,rk,2),
     x                       ugradu(ri,rj,rk,3)

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE get_ugradu_leapfrog(uf1,ugradu)
c----------------------------------------------------------------------
      include 'incurv.h'

      real uf1(nx,ny,nz,3),
     x     ugradu(nx,ny,nz,3)

      real ax,ay,az           !misc const
      real u1,u2,u3           !temp vars

      do 10 i=2,nx
         do 10 j=2,ny
            do 10 k=2,nz
               ip=i+1
               jp=j+1
               kp=k+1
               im=i-1
               jm=j-1
               km=k-1
               if (ip .gt. nx) then 
                  ip = nx
                  im = nx-2
                  endif
               if (jp .gt. ny) then 
                  jp = ny
                  im = ny-2
                  endif
               if (kp .gt. nz) then 
                  kp = nz
                  km = nz-2
                  endif
               if (im .lt. 2) then 
                  im = 2
                  ip = 4
                  endif
               if (jm .lt. 2) then 
                  jm = 2
                  jp = 4
                  endif
               if (km .lt. 2) then 
                  km = 2
                  kp = 4
                  endif

c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

               ax = 0.5*uf1(i,j,k,1)/dx
               ay = 0.5*uf1(i,j,k,2)/dy
               u1 = ax*(uf1(ip,j,k,1) - uf1(im,j,k,1)) 
               u2 = ay*(uf1(i,jp,k,1) - uf1(i,jm,k,1)) 
               az = uf1(i,j,k,3)/(dz_grid(k)+dz_grid(kp))
               u3 = az*(uf1(i,j,kp,1)-uf1(i,j,km,1))
               ugradu(i,j,k,1) = u1 + u2 + u3

c yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy
 
               ax = 0.5*uf1(i,j,k,1)/dx
               ay = 0.5*uf1(i,j,k,2)/dy
               u1 = ax*(uf1(ip,j,k,2) - uf1(im,j,k,2)) 
               u2 = ay*(uf1(i,jp,k,2) - uf1(i,jm,k,2)) 
               az = uf1(i,j,k,3)/(dz_grid(k)+dz_grid(kp))
               u3 = az*(uf1(i,j,kp,2)-uf1(i,j,km,2))
               ugradu(i,j,k,2) = u1 + u2 + u3

c zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
 
               ax = 0.5*uf1(i,j,k,1)/dx
               ay = 0.5*uf1(i,j,k,2)/dy
               u1 = ax*(uf1(ip,j,k,3) - uf1(im,j,k,3)) 
               u2 = ay*(uf1(i,jp,k,3) - uf1(i,jm,k,3)) 
               az = uf1(i,j,k,3)/(dz_grid(k)+dz_grid(kp))
               u3 = az*(uf1(i,j,kp,3)-uf1(i,j,km,3))
               ugradu(i,j,k,3) = u1 + u2 + u3

 10            continue

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE get_ugradu_leapfrog_v(uf1,ugradu)
c----------------------------------------------------------------------
      include 'incurv.h'

      real uf1(nx,ny,nz,3),
     x     ugradu(nx,ny,nz,3)

      real ax,ay,az           !misc const
      real u1,u2,u3           !temp vars
      real w1,w2,w3           !ad hoc viscosity terms for stability

      parameter (vc = 50.0)    !viscosity constant

      do 10 i=2,nx
         do 10 j=2,ny
            do 10 k=2,nz
               ip=i+1
               jp=j+1
               kp=k+1
               im=i-1
               jm=j-1
               km=k-1
               if (ip .gt. nx) then 
                  ip = nx
c                  im = nx-2
                  endif
               if (jp .gt. ny) then 
                  jp = ny
c                  im = ny-2
                  endif
               if (kp .gt. nz) then 
                  kp = nz
c                  km = nz-2
                  endif
               if (im .lt. 2) then 
                  im = 2
c                  ip = 4
                  endif
               if (jm .lt. 2) then 
                  jm = 2
c                  jp = 4
                  endif
               if (km .lt. 2) then 
                  km = 2
c                  kp = 4
                  endif

c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

               ax = 0.5*uf1(i,j,k,1)/dx
               ay = 0.5*uf1(i,j,k,2)/dy
               u1 = ax*(uf1(ip,j,k,1) - uf1(im,j,k,1)) 
               u2 = ay*(uf1(i,jp,k,1) - uf1(i,jm,k,1)) 
               az = uf1(i,j,k,3)/(dz_grid(k)+dz_grid(kp))
               u3 = az*(uf1(i,j,kp,1)-uf1(i,j,km,1))
               w1 = 0.5*(uf1(ip,j,k,1) - uf1(im,j,k,1))/dx
               w1 = vc*dx*abs(w1)*w1
               ugradu(i,j,k,1) = u1 + u2 + u3 + w1

c yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy
 
               ax = 0.5*uf1(i,j,k,1)/dx
               ay = 0.5*uf1(i,j,k,2)/dy
               u1 = ax*(uf1(ip,j,k,2) - uf1(im,j,k,2)) 
               u2 = ay*(uf1(i,jp,k,2) - uf1(i,jm,k,2)) 
               az = uf1(i,j,k,3)/(dz_grid(k)+dz_grid(kp))
               u3 = az*(uf1(i,j,kp,2)-uf1(i,j,km,2))
               w2 = 0.5*(uf1(i,jp,k,2) - uf1(i,jm,k,2))/dy
               w2 = vc*dy*abs(w2)*w2
               ugradu(i,j,k,2) = u1 + u2 + u3 + w2

c zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
 
               ax = 0.5*uf1(i,j,k,1)/dx
               ay = 0.5*uf1(i,j,k,2)/dy
               u1 = ax*(uf1(ip,j,k,3) - uf1(im,j,k,3)) 
               u2 = ay*(uf1(i,jp,k,3) - uf1(i,jm,k,3)) 
               az = uf1(i,j,k,3)/(dz_grid(k)+dz_grid(kp))
               u3 = az*(uf1(i,j,kp,3)-uf1(i,j,km,3))
               w3 = (uf1(i,j,kp,3) - uf1(i,j,km,3))/
     x                  (dz_grid(k)+dz_grid(kp))
               w3 = vc*(dz_grid(k)+dz_grid(kp))*abs(w3)*w3
               ugradu(i,j,k,3) = u1 + u2 + u3 + w3

 10            continue


      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE get_ugradu_donor_cell(uf1,ugradu)
c----------------------------------------------------------------------
      include 'incurv.h'

      real uf1(nx,ny,nz,3),
     x     ugradu(nx,ny,nz,3)

      real ax,ay,az1,az2,aaz1,aaz2           !misc const
      real u1,u2,u3                          !temp vars

      do 10 i=2,nx
         do 10 j=2,ny
            do 10 k=2,nz
               ip=i+1
               jp=j+1
               kp=k+1
               im=i-1
               jm=j-1
               km=k-1
               if (ip .gt. nx) then ip = nx
               if (jp .gt. ny) then jp = ny
               if (kp .gt. nz) then kp = nz
               if (im .lt. 2) then im = 2
               if (jm .lt. 2) then jm = 2
               if (km .lt. 2) then km = 2

c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

               ax = 0.5*uf1(i,j,k,1)/dx
               ay = 0.5*uf1(i,j,k,2)/dy
               u1 = ax*(uf1(ip,j,k,1) - uf1(im,j,k,1)) -
     x              abs(ax)*(uf1(ip,j,k,1) - 2.0*uf1(i,j,k,1) +
     x                       uf1(im,j,k,1)) 
               u2 = ay*(uf1(i,jp,k,1) - uf1(i,jm,k,1)) -
     x              abs(ay)*(uf1(i,jp,k,1) - 2.0*uf1(i,j,k,1) +
     x                       uf1(i,jm,k,1))
               az1 = 0.5*uf1(i,j,k,3)/dz_grid(k)
               az2 = 0.5*uf1(i,j,k,3)/dz_grid(kp)
               aaz1 = abs(az1)
               aaz2 = abs(az2)
               u3 = az2*uf1(i,j,kp,1) - az1*uf1(i,j,km,1) -
     x              az2*uf1(i,j,k,1) + az1*uf1(i,j,k,1) -
     x              aaz2*uf1(i,j,kp,1) - aaz1*uf1(i,j,km,1) +
     x              aaz2*uf1(i,j,k,1) + aaz1*uf1(i,j,k,1)
c               if (uf1(i,j,k,3) .gt. 0.0) then
c                  az = uf1(i,j,k,3)/dz_grid(k)
c                  u3 = az*(uf1(i,j,k,1)-uf1(i,j,km,1))
c                  endif
c               if (uf1(i,j,k,3) .le. 0.0) then
c                  az = uf1(i,j,k,3)/dz_grid(kp)
c                  u3 = az*(uf1(i,j,kp,1)-uf1(i,j,k,1))
c                  endif
               ugradu(i,j,k,1) = u1 + u2 + u3

c yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy
 
c               ax = 0.5*uf1(i,j,k,1)/dx
c               ay = 0.5*uf1(i,j,k,2)/dy
               u1 = ax*(uf1(ip,j,k,2) - uf1(im,j,k,2)) -
     x              abs(ax)*(uf1(ip,j,k,2) - 2.0*uf1(i,j,k,2) +
     x                       uf1(im,j,k,2)) 
               u2 = ay*(uf1(i,jp,k,2) - uf1(i,jm,k,2)) -
     x              abs(ay)*(uf1(i,jp,k,2) - 2.0*uf1(i,j,k,2) +
     x                       uf1(i,jm,k,2))
c               az1 = 0.5*uf1(i,j,k,3)/dz_grid(k)
c               az2 = 0.5*uf1(i,j,k,3)/dz_grid(kp)
c               aaz1 = abs(az1)
c               aaz2 = abs(az2)
               u3 = az2*uf1(i,j,kp,2) - az1*uf1(i,j,km,2) -
     x              az2*uf1(i,j,k,2) + az1*uf1(i,j,k,2) -
     x              aaz2*uf1(i,j,kp,2) - aaz1*uf1(i,j,km,2) +
     x              aaz2*uf1(i,j,k,2) + aaz1*uf1(i,j,k,2)
c               if (uf1(i,j,k,3) .gt. 0.0) then
c                  az = uf1(i,j,k,3)/dz_grid(k)
c                  u3 = az*(uf1(i,j,k,2)-uf1(i,j,km,2))
c                  endif
c               if (uf1(i,j,k,3) .le. 0.0) then
c                  az = uf1(i,j,k,3)/dz_grid(kp)
c                  u3 = az*(uf1(i,j,kp,2)-uf1(i,j,k,2))
c                  endif
               ugradu(i,j,k,2) = u1 + u2 + u3

c zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
 
c               ax = 0.5*uf1(i,j,k,1)/dx
c               ay = 0.5*uf1(i,j,k,2)/dy
               u1 = ax*(uf1(ip,j,k,3) - uf1(im,j,k,3)) -
     x              abs(ax)*(uf1(ip,j,k,3) - 2.0*uf1(i,j,k,3) +
     x                       uf1(im,j,k,3)) 
               u2 = ay*(uf1(i,jp,k,3) - uf1(i,jm,k,3)) -
     x              abs(ay)*(uf1(i,jp,k,3) - 2.0*uf1(i,j,k,3) +
     x                       uf1(i,jm,k,3))
c               az1 = 0.5*uf1(i,j,k,3)/dz_grid(k)
c               az2 = 0.5*uf1(i,j,k,3)/dz_grid(kp)
c               aaz1 = abs(az1)
c               aaz2 = abs(az2)
               u3 = az2*uf1(i,j,kp,3) - az1*uf1(i,j,km,3) -
     x              az2*uf1(i,j,k,3) + az1*uf1(i,j,k,3) -
     x              aaz2*uf1(i,j,kp,3) - aaz1*uf1(i,j,km,3) +
     x              aaz2*uf1(i,j,k,3) + aaz1*uf1(i,j,k,3)
c               if (uf1(i,j,k,3) .gt. 0.0) then
c                  az = uf1(i,j,k,3)/dz_grid(k)
c                  u3 = az*(uf1(i,j,k,3)-uf1(i,j,km,3))
c                  endif
c               if (uf1(i,j,k,3) .le. 0.0) then
c                  az = uf1(i,j,k,3)/dz_grid(kp)
c                  u3 = az*(uf1(i,j,kp,3)-uf1(i,j,k,3))
c                  endif
               ugradu(i,j,k,3) = u1 + u2 + u3
c  note the -u3 is to correct a sign error

 10            continue

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE get_Ef(Ef,aj,np,nf,up,uf,btmf,nu,ugradu,delta_t,gradP)
c Need to treat boundaries separately!!
c----------------------------------------------------------------------
      include 'incurv.h'

      real Ef(nx,ny,nz,3),
     x     aj(nx,ny,nz,3),
     x     np(nx,ny,nz),
     x     nf(nx,ny,nz),
     x     up(nx,ny,nz,3),
     x     uf(nx,ny,nz,3),
     x     btmf(nx,ny,nz,3),
     x     nu(nx,ny,nz),
     x     ugradu(nx,ny,nz,3),
     x     gradP(nx,ny,nz,3)
 
      real ntot(3)                !total plasma density
      real fnp(3)                 !fraction, np/n
      real aac(3),bbc(3),ccc(3)
      real cc(nx,ny,nz,3)

      call periodic_scalar(np)
      call periodic_scalar(nf)

      do 10 i=2,nx 
         do 10 j=2,ny
            do 10 k=2,nz

               ip = i+1
               jp = j+1
               kp = k+1

               if (ip .gt. nx) then ip = nx
               if (jp .gt. ny) then jp = ny
               if (kp .gt. nz) then kp = nz

               ntot(1) = 0.5*(nf(i,j,k)+nf(ip,j,k))
     x                 + 0.5*(np(i,j,k)+np(ip,j,k))
               ntot(2) = 0.5*(nf(i,j,k)+nf(i,jp,k))
     x                 + 0.5*(np(i,j,k)+np(i,jp,k))
               ntot(3) = 0.5*(nf(i,j,k)+nf(i,j,kp))
     x                 + 0.5*(np(i,j,k)+np(i,j,kp))

               fnp(1) = 0.5*(np(i,j,k)+np(ip,j,k))/ntot(1)
               fnp(2) = 0.5*(np(i,j,k)+np(i,jp,k))/ntot(2)
               fnp(3) = 0.5*(np(i,j,k)+np(i,j,kp))/ntot(3)
               
               do 10 m=1,3
                  a(i,j,k,m) = aj(i,j,k,m) - fnp(m)*up(i,j,k,m)
 10               continue

      call crossf(a,btmf,c)

c      call get_ugradu_donor_cell(uf,ugradu)
c      call get_ugradu_leapfrog(uf,ugradu)
c      call get_ugradu_leapfrog_v(uf,ugradu)
c      call get_ugradu_2nd_order(uf,ugradu)
      call get_ugradu_Lax(uf,ugradu,delta_t)

      do 20 i=2,nx 
         do 20 j=2,ny
            do 20 k=2,nz

               ip = i+1
               jp = j+1
               kp = k+1

               if (ip .gt. nx) then ip = nx
               if (jp .gt. ny) then jp = ny
               if (kp .gt. nz) then kp = nz

               ntot(1) = 0.5*(nf(i,j,k)+nf(ip,j,k))
     x                 + 0.5*(np(i,j,k)+np(ip,j,k))
               ntot(2) = 0.5*(nf(i,j,k)+nf(i,jp,k))
     x                 + 0.5*(np(i,j,k)+np(i,jp,k))
               ntot(3) = 0.5*(nf(i,j,k)+nf(i,j,kp))
     x                 + 0.5*(np(i,j,k)+np(i,j,kp))

               fnp(1) = 0.5*(np(i,j,k)+np(ip,j,k))/ntot(1)
               fnp(2) = 0.5*(np(i,j,k)+np(i,jp,k))/ntot(2)
               fnp(3) = 0.5*(np(i,j,k)+np(i,j,kp))/ntot(3)

               do 20 m=1,3
                  Ef(i,j,k,m) = c(i,j,k,m) - ugradu(i,j,k,m) 
     x                          + nu(i,j,k)*fnp(m)*up(i,j,k,m)
     x                          + nuei*aj(i,j,k,m) - gradP(i,j,k,m)
c     x                          + etar(i,j,k,m)*aj(i,j,k,m)
c                  Ef(i,j,k,m) = c(i,j,k,m) + 
c     x                          nu(i,j,k)*fnp(m)*up(i,j,k,m)
 20            continue

c      call fix_tangential_E(Ef)
      call periodic(Ef)

      return
      end
c----------------------------------------------------------------------



c----------------------------------------------------------------------
      SUBROUTINE get_uplus_uminus(Ef,btmf,uf,nu,np,nf,uplus,uminus,
     x                            delta_t)
c This is the heart of the fluid velocity update.  It solves eqn. 18
c (Dan's paper) for uf+
c----------------------------------------------------------------------
      include 'incurv.h'

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

      do 10 i=2,nx    
         do 10 j=2,ny
            do 10 k=2,nz
               do 10 m=1,3
                  uminus(i,j,k,m) = uf(i,j,k,m) + 
     x                              0.5*delta_t*Ef(i,j,k,m)
 10            continue

      call crossf(uminus, btmf, um_x_B)
      call get_um_dot_BB(uminus , btmf, um_dot_BB)

      call face_to_center(btmf,btc)

      do 15 i=2,nx 
         do 15 j=2,ny
            do 15 k=2,nz
               bsqrd(i,j,k) =  btc(i,j,k,1)**2 + btc(i,j,k,2)**2 + 
     x               btc(i,j,k,3)**2
 15            continue

      call periodic_scalar(bsqrd)
              
      do 20 i=2,nx
         do 20 j=2,ny
            do 20 k=2,nz

               ip = i+1
               jp = j+1
               kp = k+1

               if (ip .gt. nx) then ip = nx
               if (jp .gt. ny) then jp = ny
               if (kp .gt. nz) then kp = nz

               bsq(1) = 0.5*(bsqrd(i,j,k) + bsqrd(ip,j,k))
               bsq(2) = 0.5*(bsqrd(i,j,k) + bsqrd(i,jp,k))
               bsq(3) = 0.5*(bsqrd(i,j,k) + bsqrd(i,j,kp))

               npave(1) = 0.5*(np(i,j,k)+np(ip,j,k))
               npave(2) = 0.5*(np(i,j,k)+np(i,jp,k))
               npave(3) = 0.5*(np(i,j,k)+np(i,j,kp))

               ntot(1) = 0.5*(nf(i,j,k)+nf(ip,j,k)) + npave(1)
               ntot(2) = 0.5*(nf(i,j,k)+nf(i,jp,k)) + npave(2)
               ntot(3) = 0.5*(nf(i,j,k)+nf(i,j,kp)) + npave(3)

               do 20 m=1,3

                  eta = npave(m)*delta_t/(2.0*ntot(m))
                  QQ = eta/(1.0+nu(i,j,k)*eta)
                  PP = (1.0-nu(i,j,k)*eta)/(1.0+nu(i,j,k)*eta)
                  a1 = 1.0/(1.0 + QQ*QQ*bsq(m))
                  b1 = PP - (QQ*QQ*bsq(m))
                  b2 = (QQ*PP) + QQ
                  b3 = (QQ*QQ*PP) + (QQ*QQ)

                  uplus(i,j,k,m) = a1*(b1*uminus(i,j,k,m) + 
     x                             b2*um_x_B(i,j,k,m) + 
     x                             b3*um_dot_BB(i,j,k,m))

 20            continue

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE predict_uf(Ef,b0,b1,b12,uf,uf2,ufp2,nu,np,nf,uplus, 
     x                      uminus,ugradu,up,gradP,nuin,bdp)
c Calculate the fluid velocity, uf,  at the new time step and replace
c uf1 with the new value, uf, in preparation for the next time step.
c----------------------------------------------------------------------
      include 'incurv.h'

      real Ef(nx,ny,nz,3),
     x     b0(nz),
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
     x     gradP(nx,ny,nz,3),
     x     nuin(nx,ny,nz),
     x     bdp(nx,ny,nz,3)

      real b1h(nx,ny,nz,3)
      real bth(nx,ny,nz,3)
      real btmfh(nx,ny,nz,3)
      real ajh(nx,ny,nz,3)

      real delta_t

      delta_t = 2.0*dtsub

      do 10 i=1,nx
         do 10 j=1,ny
            do 10 k=1,nz
               b1h(i,j,k,1) = 0.5*(b1(i,j,k,1) + b12(i,j,k,1))
               b1h(i,j,k,2) = 0.5*(b1(i,j,k,2) + b12(i,j,k,2))
               b1h(i,j,k,3) = 0.5*(b1(i,j,k,3) + b12(i,j,k,3))
               bth(i,j,k,1) = b1h(i,j,k,1) + bdp(i,j,k,1)
               bth(i,j,k,2) = b1h(i,j,k,2) + bdp(i,j,k,2)
               bth(i,j,k,3) = b1h(i,j,k,3) + b0(k) + bdp(i,j,k,3)
 10            continue

      call cov_to_contra(bth,btmfh)
      call curlB(b1h,nf,np,ajh)

      call get_Ef(Ef,ajh,np,nf,up,uf,btmfh,nu,ugradu,delta_t,gradP)
      call get_uplus_uminus(Ef,btmfh,uf2,nu,np,nf,uplus,uminus,
     x                      delta_t)

      do 20 i=2,nx
         do 20 j=2,ny
            do 20 k=2,nz
               do 20 m=1,3
                  ufp2(i,j,k,m) = uplus(i,j,k,m) + 
     x                            0.5*delta_t*Ef(i,j,k,m) -
     x                        0.5*delta_t*nuin(i,j,k)*uplus(i,j,k,m)
 20            continue

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE correct_uf(Ef,btmf,uf,uf2,ufp2,nu,np,nf,uplus,uminus, 
     x                      ugradu,aj,up,ufp1,gradP,nuin)
c Calculate the fluid velocity, uf,  at the new time step and replace
c uf1 with the new value, uf, in preparation for the next time step.
c----------------------------------------------------------------------
      include 'incurv.h'

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
     x     gradP(nx,ny,nz,3),
     x     nuin(nx,ny,nz)
 
      real delta_t

      delta_t = dtsub
      
      do 10 i=2,nx
         do 10 j=2,ny
            do 10 k=2,nz
               do 10 m=1,3
                  ufp1(i,j,k,m) = 0.5*(uf(i,j,k,m) + ufp2(i,j,k,m))
 10              continue

      call get_Ef(Ef,aj,np,nf,up,ufp1,btmf,nu,ugradu,delta_t,gradP)
      call get_uplus_uminus(Ef,btmf,uf,nu,np,nf,uplus,uminus,
     x                      delta_t)

      do 20 i=2,nx
         do 20 j=2,ny
            do 20 k=2,nz
               do 20 m=1,3
                  uf2(i,j,k,m) = uf(i,j,k,m)
                  uf(i,j,k,m) = uplus(i,j,k,m) + 0.5*dtsub*Ef(i,j,k,m)
     x                          - 0.5*dtsub*nuin(i,j,k)*uplus(i,j,k,m)
 20            continue

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE get_E(E,bt,btmf,aj,up,uf,uf2,np,nf,nu,gradP)
c E must be at time level m. We have uf at levels m-1/2 and m+1/2, so
c the average value is used for uf in the calculation of ui.
c----------------------------------------------------------------------
      include 'incurv.h'

      real E(nx,ny,nz,3),
     x     bt(nx,ny,nz,3),
     x     btmf(nx,ny,nz,3),
     x     aj(nx,ny,nz,3),
     x     up(nx,ny,nz,3),
     x     uf(nx,ny,nz,3),
     x     uf2(nx,ny,nz,3),
     x     np(nx,ny,nz),
     x     nf(nx,ny,nz),
     x     nu(nx,ny,nz),
     x     gradP(nx,ny,nz,3)

      real ntot(3)         !total density np + nf
      real fnp(3),fnf(3)   !fraction np and nf of n
      real npave(3)

c      real a(nx,ny,nz,3), 
c     x     c(nx,ny,nz,3)  !dummy vars for doing cross product

      call periodic_scalar(np)
      call periodic_scalar(nf)

      do 10 i=2,nx    
         do 10 j=2,ny
            do 10 k=2,nz

               ip = i+1
               jp = j+1
               kp = k+1

               if (ip .gt. nx) then ip = nx
               if (jp .gt. ny) then jp = ny
               if (kp .gt. nz) then kp = nz

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

               do 10 m=1,3
                  a(i,j,k,m) = aj(i,j,k,m) - fnp(m)*up(i,j,k,m) - 
     x                         fnf(m)*0.5*(uf2(i,j,k,m)+uf(i,j,k,m))
c                  a(i,j,k,m) = - fnp(m)*up(i,j,k,m) - 
c     x                         fnf(m)*0.5*(uf2(i,j,k,m)+uf(i,j,k,m))
 10               continue

      call crossf(a,btmf,c)
               
      do 20 i=2,nx      
         do 20 j=2,ny   
            do 20 k=2,nz
               do 20 m=1,3 
                  E(i,j,k,m) = c(i,j,k,m) + nu(i,j,k)*aj(i,j,k,m)
     x                         + nuei*aj(i,j,k,m) - gradP(i,j,k,m)
c     x                         + etar(i,j,k,m)*aj(i,j,k,m)
 20               continue

c      call fix_tangential_E(E)
      call periodic(E)

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE predict_B(b1,b12,b1p2,bt,btmf,E,aj,up,uf,uf2,np,nf,nu,
     x                     gradP)
c Predictor step in magnetic field update.
c----------------------------------------------------------------------
      include 'incurv.h'

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
     x     nu(nx,ny,nz),
     x     gradP(nx,ny,nz,3)

      real curl_E(nx,ny,nz,3)   !curl of E

      call get_E(E,bt,btmf,aj,up,uf,uf2,np,nf,nu,gradP)  !E at time level m 

      call curlE(E,curl_E)
c      call fix_tangential_E(E)
      call periodic(E)

      do 10 i=2,nx-1
         do 10 j=2,ny-1
            do 10 k=2,nz-1
               do 10 m=1,3
                  b1p2(i,j,k,m) = b12(i,j,k,m) - 
     x                            2.0*dtsub*curl_E(i,j,k,m)
 10               continue

c      call boundaries(b1p2)
c      call damp(b1p2)
      call periodic(b1p2)

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE get_Ep1(E,b0,b1,b1p2,aj,up,uf,np,nf,nu,gradP,bdp)
c The main feature here is that E must be calculated at time level
c m + 1/2.  That means that we need B at m + 1/2.  So b1p1 is
c calculated as 0.5*(b1 + b1p2).  uf and np are already at time level
c m + 1/2, so they are used as is. 
c----------------------------------------------------------------------
      include 'incurv.h'

      real E(nx,ny,nz,3),
     x     b0(nz),
     x     b1(nx,ny,nz,3),
     x     b1p2(nx,ny,nz,3),
     x     aj(nx,ny,nz,3),
     x     up(nx,ny,nz,3),
     x     uf(nx,ny,nz,3),
     x     np(nx,ny,nz),
     x     nf(nx,ny,nz),
     x     nu(nx,ny,nz),
     x     gradP(nx,ny,nz,3),
     x     bdp(nx,ny,nz,3)

      real b1p1(nx,ny,nz,3)   !b1 at time level m + 1/2
      real btp1(nx,ny,nz,3)   !bt at time level m + 1/2
      real btp1mf(nx,ny,nz,3) !btp1 at contravarient position
      real ntot(3)            !total density np + nf
      real fnp(3),fnf(3)      !fraction np and nf of n
      real npave(3)

c      real a(nx,ny,nz,3),
c     x     c(nx,ny,nz,3)    !dummy vars for doing cross product


      do 5 i=1,nx
         do 5 j=1,ny
            do 5 k=1,nz
               btp1(i,j,k,1) = 0.5*(b1p2(i,j,k,1) + b1(i,j,k,1)) +
     x              bdp(i,j,k,1)
               b1p1(i,j,k,1) = 0.5*(b1p2(i,j,k,1) + b1(i,j,k,1))
               btp1(i,j,k,2) = 0.5*(b1p2(i,j,k,2) + b1(i,j,k,2)) +
     x              bdp(i,j,k,2)
               b1p1(i,j,k,2) = 0.5*(b1p2(i,j,k,2) + b1(i,j,k,2))
               btp1(i,j,k,3) = b0(k) + bdp(i,j,k,3) +
     x              0.5*(b1p2(i,j,k,3)+b1(i,j,k,3))
               b1p1(i,j,k,3) = 0.5*(b1p2(i,j,k,3)+b1(i,j,k,3))
 5             continue

      call curlB(b1p1,nf,np,aj)

      call periodic_scalar(np)
      call periodic_scalar(nf)

      do 10 i=2,nx       
         do 10 j=2,ny
            do 10 k=2,nz

               ip = i+1
               jp = j+1
               kp = k+1

               if (ip .gt. nx) then ip = nx
               if (jp .gt. ny) then jp = ny
               if (kp .gt. nz) then kp = nz

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

               do 10 m=1,3
                  a(i,j,k,m) = aj(i,j,k,m) - fnp(m)*up(i,j,k,m) - 
     x                                       fnf(m)*uf(i,j,k,m)
c                  a(i,j,k,m) = - fnp(m)*up(i,j,k,m) - 
c     x                           fnf(m)*uf(i,j,k,m)
 10               continue

      call cov_to_contra(btp1,btp1mf)
      call crossf(a,btp1mf,c)
               
      do 20 i=2,nx       
         do 20 j=2,ny     
            do 20 k=2,nz  
               do 20 m=1,3 
                  E(i,j,k,m) = c(i,j,k,m) + nu(i,j,k)*aj(i,j,k,m)
     x                         + nuei*aj(i,j,k,m) - gradP(i,j,k,m)
c     x                         + etar(i,j,k,m)*aj(i,j,k,m)
 20               continue

c      call fix_tangential_E(E)
      call periodic(E)

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE correct_B(b0,b1,b1p2,E,aj,up,uf,np,nf,nu,gradP,bdp)
c Corrector step in magnetic field update.
c----------------------------------------------------------------------
      include 'incurv.h'

      real b0(nz),
     x     b1(nx,ny,nz,3),
     x     b1p2(nx,ny,nz,3),
     x     E(nx,ny,nz,3),
     x     aj(nx,ny,nz,3),
     x     up(nx,ny,nz,3),
     x     uf(nx,ny,nz,3),
     x     np(nx,ny,nz),
     x     nf(nx,ny,nz),
     x     nu(nx,ny,nz),
     x     gradP(nx,ny,nz,3),
     x     bdp(nx,ny,nz,3)

      real curl_E(nx,ny,nz,3)            !curl of E

      call get_Ep1(E,b0,b1,b1p2,aj,up,uf,np,nf,nu,gradP,bdp)  
                                                   !E at time level m 

      call curlE(E,curl_E)
c      call fix_tangential_E(E)
      call periodic(E)

      do 10 i=2,nx-1
         do 10 j=2,ny-1
            do 10 k=2,nz-1
               do 10 m=1,3
                  b1p2(i,j,k,m) = b1(i,j,k,m) - 
     x                            dtsub*curl_E(i,j,k,m)
 10               continue

c      call boundaries(b1p2)
c      call damp(b1p2)
      call periodic(b1p2)

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE predict_nf(nf,nf1,nf3,nfp1,uf,divu,b1)
c----------------------------------------------------------------------
      include 'incurv.h'

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
 
      minnf = 10.0e20
      maxnf = 10.0e20

      call periodic_scalar(nf1)

      do 5 i=2,nx-1
         do 5 j=2,ny-1
            do 5 k=2,nz-1
               flx(i,j,k,1) = 0.5*(nf1(i,j,k)+nf1(i+1,j,k))*uf(i,j,k,1)
               flx(i,j,k,2) = 0.5*(nf1(i,j,k)+nf1(i,j+1,k))*uf(i,j,k,2)
               flx(i,j,k,3) = 0.5*(nf1(i,j,k)+nf1(i,j,k+1))*uf(i,j,k,3)
 5             continue

      do 20 j=1,ny
         do 20 k=1,nz
            do 20 m=1,3
               flx(1,j,k,m) = flx(2,j,k,m)
               flx(nx,j,k,m) = flx(nx-1,j,k,m)
 20         continue

      do 30 i=1,nx
         do 30 k=1,nz
            do 30 m=1,3
               flx(i,1,k,m) = flx(i,2,k,m)
               flx(i,ny,k,m) = flx(i,ny-1,k,m)
 30         continue

      do 40 i=1,nx
         do 40 j=1,ny
            do 40 m=1,3
               flx(i,j,1,m) = flx(i,j,2,m)
               flx(i,j,nz,m) = flx(i,j,nz-1,m)
 40         continue
     
      call periodic(flx)

      do 10 i=2,nx-1
         do 10 j=2,ny-1
            do 10 k=2,nz-1
               nfp1(i,j,k) = nf3(i,j,k) 
     x         - (2.0*dtsub/dx)*(flx(i,j,k,1)-flx(i-1,j,k,1))
     x         - (2.0*dtsub/dy)*(flx(i,j,k,2)-flx(i,j-1,k,2))
     x         - (2.0*dtsub/(dz_cell(k)))*
     x           (flx(i,j,k,3)-flx(i,j,k-1,3))
               divu(i,j,k) = (uf(i,j,k,1) - uf(i-1,j,k,1))/dx +
     x                       (uf(i,j,k,2) - uf(i,j-1,k,2))/dy +
     x                       (uf(i,j,k,3) - uf(i,j,k-1,3))/dz_cell(k)
 10            continue

      do 50 i=2,nx-1
         do 50 j=2,ny-1
            do 50 k=2,nz-1
               nf(i,j,k) = 0.5*(nfp1(i,j,k) + nf1(i,j,k))
               nf3(i,j,k) = nf1(i,j,k)
 50            continue

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE correct_nf(nf,nf1,ufp1)
c----------------------------------------------------------------------
      include 'incurv.h'

      real nf(nx,ny,nz),
     x     nf1(nx,ny,nz),
     x     ufp1(nx,ny,nz,3) 

      real flx(nx,ny,nz,3)
      real minnf,maxnf
      real ufp1c(nx,ny,nz,3)

      minnf = 10.0000000000000000e20
      maxnf = 10.0000000000000000e20

c      call face_to_center(ufp1,ufp1c)

      do 5 i=2,nx-1
         do 5 j=2,ny-1
            do 5 k=2,nz-1
             flx(i,j,k,1) = 0.5*(nf(i,j,k)+nf(i+1,j,k))*ufp1(i,j,k,1)
             flx(i,j,k,2) = 0.5*(nf(i,j,k)+nf(i,j+1,k))*ufp1(i,j,k,2)
             flx(i,j,k,3) = 0.5*(nf(i,j,k)+nf(i,j,k+1))*ufp1(i,j,k,3)
 5           continue

      do 20 j=1,ny
         do 20 k=1,nz
            do 20 m=1,3
               flx(1,j,k,m) = flx(2,j,k,m)
               flx(nx,j,k,m) = flx(nx-1,j,k,m)
 20         continue

      do 30 i=1,nx
         do 30 k=1,nz
            do 30 m=1,3
               flx(i,1,k,m) = flx(i,2,k,m)
               flx(i,ny,k,m) = flx(i,ny-1,k,m)
 30         continue

      do 40 i=1,nx
         do 40 j=1,ny
            do 40 m=1,3
               flx(i,j,1,m) = flx(i,j,2,m)
               flx(i,j,nz,m) = flx(i,j,nz-1,m)
 40         continue

      do 10 i=2,nx-1
         do 10 j=2,ny-1
            do 10 k=2,nz-1
               nf1(i,j,k) = nf1(i,j,k) 
     x            - (dtsub/dx)*(flx(i,j,k,1)-flx(i-1,j,k,1))
     x            - (dtsub/dy)*(flx(i,j,k,2)-flx(i,j-1,k,2))
     x            - (dtsub/(dz_cell(k)))*
     x              (flx(i,j,k,3)-flx(i,j,k-1,3))
               if (nf(i,j,k) .lt. minnf) then minnf = nf(i,j,k)
               if (nf(i,j,k) .ge. maxnf) then maxnf = nf(i,j,k)
 10            continue


c      write(*,*) 'Min nf.....',minnf
c      write(*,*) 'Max nf.....',maxnf

      call periodic_scalar(nf1)

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE trans_nf_Lax(nf,nf1,nfp1,uf)
c----------------------------------------------------------------------
      include 'incurv.h'

      real nf(nx,ny,nz),
     x     nf1(nx,ny,nz),
     x     nfp1(nx,ny,nz),
     x     uf(nx,ny,nz,3)
 
      real minnf,maxnf
      real flx(nx,ny,nz,3)
      real ufc(nx,ny,nz,3)

c      call face_to_center(uf,ufc)
 
      minnf = 10.0e20
      maxnf = 10.0e20

      do 3 i=2,nx-1
         do 3 j=2,ny-1
            do 3 k=2,nz-1
               nf1(i,j,k) = nfp1(i,j,k)
 3             continue      

      call periodic_scalar(nf1)

      do 5 i=2,nx-1
         do 5 j=2,ny-1
            do 5 k=2,nz-1
               flx(i,j,k,1) = 0.5*(nf1(i,j,k)+nf1(i+1,j,k))*uf(i,j,k,1)
               flx(i,j,k,2) = 0.5*(nf1(i,j,k)+nf1(i,j+1,k))*uf(i,j,k,2)
               flx(i,j,k,3) = 0.5*(nf1(i,j,k)+nf1(i,j,k+1))*uf(i,j,k,3)
 5             continue

      do 20 j=1,ny
         do 20 k=1,nz
            do 20 m=1,3
               flx(1,j,k,m) = flx(2,j,k,m)
               flx(nx,j,k,m) = flx(nx-1,j,k,m)
 20         continue

      do 30 i=1,nx
         do 30 k=1,nz
            do 30 m=1,3
               flx(i,1,k,m) = flx(i,2,k,m)
               flx(i,ny,k,m) = flx(i,ny-1,k,m)
 30         continue

      do 40 i=1,nx
         do 40 j=1,ny
            do 40 m=1,3
               flx(i,j,1,m) = flx(i,j,2,m)
               flx(i,j,nz,m) = flx(i,j,nz-1,m)
 40         continue

      call periodic(flx)

      do 10 i=2,nx-1
         do 10 j=2,ny-1
            do 10 k=2,nz-1
               nfp1(i,j,k) = (1.0/6.0)*(nf1(i+1,j,k)+nf1(i-1,j,k)+
     x                               nf1(i,j+1,k)+nf1(i,j-1,k)+
     x                               nf1(i,j,k+1)+nf1(i,j,k-1))  
     x         - 0.5*(dtsub/dx)*(flx(i+1,j,k,1)-flx(i-1,j,k,1))
     x         - 0.5*(dtsub/dy)*(flx(i,j+1,k,2)-flx(i,j-1,k,2))
     x         - (dtsub/(dz_grid(k)+dz_grid(k+1)))*
     x           (flx(i,j,k+1,3)-flx(i,j,k-1,3))
 10            continue

      do 50 i=2,nx-1
         do 50 j=2,ny-1
            do 50 k=2,nz-1
               nf(i,j,k) = 0.5*(nf1(i,j,k) + nfp1(i,j,k))
 50            continue

      call periodic_scalar(nf)

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE trans_nf_LaxWend1(nf,nf1,nfp1,uf)
c----------------------------------------------------------------------
      include 'incurv.h'

      real nf(nx,ny,nz),
     x     nf1(nx,ny,nz),
     x     nfp1(nx,ny,nz),
     x     uf(nx,ny,nz,3)
 
      real minnf,maxnf
      real flx(nx,ny,nz,3)

      do 3 i=2,nx-1
         do 3 j=2,ny-1
            do 3 k=2,nz-1
               nf1(i,j,k) = nfp1(i,j,k)
 3             continue   
   
      call periodic_scalar(nf1)

      do 5 i=2,nx-1
         do 5 j=2,ny-1
            do 5 k=2,nz-1
               flx(i,j,k,1) = 0.5*(nf1(i,j,k)+nf1(i+1,j,k))*uf(i,j,k,1)
               flx(i,j,k,2) = 0.5*(nf1(i,j,k)+nf1(i,j+1,k))*uf(i,j,k,2)
               flx(i,j,k,3) = 0.5*(nf1(i,j,k)+nf1(i,j,k+1))*uf(i,j,k,3)
 5             continue

c      do 20 j=1,ny
c         do 20 k=1,nz
c            do 20 m=1,3
c               flx(1,j,k,m) = flx(2,j,k,m)
c               flx(nx,j,k,m) = flx(nx-1,j,k,m)
c 20         continue

c      do 30 i=1,nx
c         do 30 k=1,nz
c            do 30 m=1,3
c               flx(i,1,k,m) = flx(i,2,k,m)
c               flx(i,ny,k,m) = flx(i,ny-1,k,m)
c 30         continue

c      do 40 i=1,nx
c         do 40 j=1,ny
c            do 40 m=1,3
c               flx(i,j,1,m) = flx(i,j,2,m)
c               flx(i,j,nz,m) = flx(i,j,nz-1,m)
c 40         continue

      call periodic(flx)
      call periodic_scalar(nf1)

      do 10 i=2,nx-1
         do 10 j=2,ny-1
            do 10 k=2,nz-1
               nf(i,j,k) = (1.0/6.0)*(nf1(i+1,j,k)+nf1(i-1,j,k)+
     x                               nf1(i,j+1,k)+nf1(i,j-1,k)+
     x                               nf1(i,j,k+1)+nf1(i,j,k-1))  
     x         - 0.5*(0.5*dtsub/dx)*(flx(i+1,j,k,1)-flx(i-1,j,k,1))
     x         - 0.5*(0.5*dtsub/dy)*(flx(i,j+1,k,2)-flx(i,j-1,k,2))
     x         - (0.5*dtsub/(dz_grid(k)+dz_grid(k+1)))*
     x           (flx(i,j,k+1,3)-flx(i,j,k-1,3))
 10            continue

      call periodic_scalar(nf)

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE trans_nf_LaxWend2(nf,nf1,nfp1,ufp1)
c----------------------------------------------------------------------
      include 'incurv.h'

      real nf(nx,ny,nz),
     x     nf1(nx,ny,nz),
     x     nfp1(nx,ny,nz),
     x     ufp1(nx,ny,nz,3)
 
      real minnf,maxnf
      real flx(nx,ny,nz,3)

      call periodic_scalar(nf)

      do 5 i=2,nx-1
         do 5 j=2,ny-1
            do 5 k=2,nz-1
               flx(i,j,k,1) = 0.5*(nf(i,j,k)+nf(i+1,j,k))*ufp1(i,j,k,1)
               flx(i,j,k,2) = 0.5*(nf(i,j,k)+nf(i,j+1,k))*ufp1(i,j,k,2)
               flx(i,j,k,3) = 0.5*(nf(i,j,k)+nf(i,j,k+1))*ufp1(i,j,k,3)
 5             continue

c      do 20 j=1,ny
c         do 20 k=1,nz
c            do 20 m=1,3
c               flx(1,j,k,m) = flx(2,j,k,m)
c               flx(nx,j,k,m) = flx(nx-1,j,k,m)
c 20         continue

c      do 30 i=1,nx
c         do 30 k=1,nz
c            do 30 m=1,3
c               flx(i,1,k,m) = flx(i,2,k,m)
c               flx(i,ny,k,m) = flx(i,ny-1,k,m)
c 30         continue

c      do 40 i=1,nx
c         do 40 j=1,ny
c            do 40 m=1,3
c               flx(i,j,1,m) = flx(i,j,2,m)
c               flx(i,j,nz,m) = flx(i,j,nz-1,m)
c 40         continue

      call periodic(flx)
      call periodic_scalar(nf1)

      do 10 i=2,nx-1
         do 10 j=2,ny-1
            do 10 k=2,nz-1
               nfp1(i,j,k) = (1.0/6.0)*(nf1(i+1,j,k)+nf1(i-1,j,k)+
     x                               nf1(i,j+1,k)+nf1(i,j-1,k)+
     x                               nf1(i,j,k+1)+nf1(i,j,k-1))  
     x         - 0.5*(dtsub/dx)*(flx(i+1,j,k,1)-flx(i-1,j,k,1))
     x         - 0.5*(dtsub/dy)*(flx(i,j+1,k,2)-flx(i,j-1,k,2))
     x         - (dtsub/(dz_grid(k)+dz_grid(k+1)))*
     x           (flx(i,j,k+1,3)-flx(i,j,k-1,3))
 10            continue

      call periodic_scalar(nfp1)


      return
      end
c----------------------------------------------------------------------













