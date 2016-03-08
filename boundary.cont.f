c---------------------------------------------------------------------
      SUBROUTINE periodic(aa)
c---------------------------------------------------------------------
      include 'incurv.h'

      real aa(nx,ny,nz,3)

c x direction

      do 10 j=2,ny-1
         do 10 k=2,nz-1
            do 10 m=1,3
               aa(1,j,k,m) = aa(nx-1,j,k,m)
               aa(nx,j,k,m) = aa(2,j,k,m)
 10            continue

c y direction

      do 20 i=2,nx-1
         do 20 k=2,nz-1
            do 20 m=1,3
               aa(i,1,k,m) = aa(i,ny-1,k,m)
               aa(i,ny,k,m) = aa(i,2,k,m)
 20            continue

c z direction

      do 30 i=2,nx-1
         do 30 j=2,ny-1
            do 30 m=1,3
               aa(i,j,1,m) = aa(i,j,nz-1,m)
               aa(i,j,nz,m) = aa(i,j,2,m)
 30            continue

      return
      end
c---------------------------------------------------------------------


c---------------------------------------------------------------------
      SUBROUTINE damp(aa)
c---------------------------------------------------------------------
      include 'incurv.h'

      real aa(nx,ny,nz,3)

      do 10 j=1,ny
         do 10 k=1,nz
            do 10 m=1,3

               aa(2,j,k,m) = aa(2,j,k,m)*0.8 
               aa(3,j,k,m) = aa(3,j,k,m)*0.85
               aa(4,j,k,m) = aa(4,j,k,m)*0.9
               aa(5,j,k,m) = aa(5,j,k,m)*0.95
               aa(nx,j,k,m) = aa(nx,j,k,m)*0.8
               aa(nx-1,j,k,m) = aa(nx-1,j,k,m)*0.85 
               aa(nx-2,j,k,m) = aa(nx-2,j,k,m)*0.9 
               aa(nx-3,j,k,m) = aa(nx-3,j,k,m)*0.95 

 10            continue


      do 20 i=2,nx-1
         do 20 k=1,nz
            do 20 m=1,3

               aa(i,2,k,m) = aa(i,2,k,m)*0.8 
               aa(i,3,k,m) = aa(i,3,k,m)*0.85
               aa(i,4,k,m) = aa(i,4,k,m)*0.9
               aa(i,5,k,m) = aa(i,5,k,m)*0.95
               aa(i,ny,k,m) = aa(i,ny,k,m)*0.8
               aa(i,ny-1,k,m) = aa(i,ny-1,k,m)*0.85 
               aa(i,ny-2,k,m) = aa(i,ny-2,k,m)*0.9 
               aa(i,ny-3,k,m) = aa(i,ny-3,k,m)*0.95 

 20            continue


      do 30 i=2,nx-1
         do 30 j=2,ny-1
            do 30 m=1,3

               aa(i,j,2,m) = aa(i,j,2,m)*0.8 
               aa(i,j,3,m) = aa(i,j,3,m)*0.85
               aa(i,j,4,m) = aa(i,j,4,m)*0.9
               aa(i,j,5,m) = aa(i,j,5,m)*0.95
               aa(i,j,nz-1,m) = aa(i,j,nz-1,m)*0.8
               aa(i,j,nz-2,m) = aa(i,j,nz-2,m)*0.85 
               aa(i,j,nz-3,m) = aa(i,j,nz-3,m)*0.9 
               aa(i,j,nz-4,m) = aa(i,j,nz-4,m)*0.95 

 30            continue

      return
      end
c---------------------------------------------------------------------


c---------------------------------------------------------------------
      SUBROUTINE extrapolated(aa)
c---------------------------------------------------------------------
      include 'incurv.h'

      real aa(nx,ny,nz,3)

      do 40 j=2,ny-1 
         do 40 k=2,nz-1
            do 40 m=1,3
               aa(1,j,k,m) = 3.0*aa(2,j,k,m) - 3.0*aa(3,j,k,m)  
     x                          + aa(4,j,k,m)
               aa(nx,j,k,m) = 3.0*aa(nx-1,j,k,m) - 3.0*aa(nx-2,j,k,m)  
     x                           + aa(nx-3,j,k,m)
 40            continue


      do 50 i=2,nx-1
         do 50 k=2,nz-1
            do 50 m=1,3
               aa(i,1,k,m) = -(3.0*aa(i,2,k,m) - 3.0*aa(i,3,k,m) 
     x                          + aa(i,4,k,m))
               aa(i,ny,k,m) = 3.0*aa(i,ny-1,k,m) - 3.0*aa(i,ny-2,k,m)
     x                           + aa(i,ny-3,k,m)
 50            continue    

      do 60 i=2,nx-1
         do 60 j=2,ny-1
            do 60 m=1,3

c               aa(i,j,1,m) = -(3.0*aa(i,j,2,m) - 3.0*aa(i,j,3,m)
c     x                          + aa(i,j,4,m))
c               aa(i,j,nz,m) = 3.0*aa(i,j,nz-1,m) - 3.0*aa(i,j,nz-2,m)
c     x                           + aa(i,j,nz-3,m)
 60            continue

      return
      end
c---------------------------------------------------------------------


c---------------------------------------------------------------------
      SUBROUTINE edges_corners_1(aa)
c---------------------------------------------------------------------
      include 'incurv.h'

      real aa(nx,ny,nz,3)
      real f2,f3

      f2=(1.0/2.0)    !damp edges slightly...should be 1/2
      f3=(1.0/3.0)    !damp corners slightly...should be 1/3

c x edge 

      do 10 i=2,nx-1
         do 10 m=1,3
            aa(i,1,1,m)   = f2*(aa(i,2,1,m)     + aa(i,1,2,m))
            aa(i,1,nz,m)  = f2*(aa(i,2,nz,m)    + aa(i,1,nz-1,m))
            aa(i,ny,1,m)  = f2*(aa(i,ny,2,m)    + aa(i,ny-1,1,m))
            aa(i,ny,nz,m) = f2*(aa(i,ny,nz-1,m) + aa(i,ny-1,nz,m))
 10      continue

c y edge

      do 20 j=2,ny-1
         do 20 m=1,3
            aa(1,j,1,m)   = f2*(aa(2,j,1,m)     + aa(1,j,2,m))
            aa(1,j,nz,m)  = f2*(aa(2,j,nz,m)    + aa(1,j,nz-1,m))
            aa(nx,j,1,m)  = f2*(aa(nx,j,2,m)    + aa(nx-1,j,1,m))
            aa(nx,j,nz,m) = f2*(aa(nx,j,nz-1,m) + aa(nx-1,j,nz,m))
 20      continue

c z edge

      do 30 k=2,nz-1
         do 30 m=1,3
            aa(1,1,k,m)   = f2*(aa(2,1,k,m)     + aa(1,2,k,m))
            aa(1,ny,k,m)  = f2*(aa(2,ny,k,m)    + aa(1,ny-1,k,m))
            aa(nx,1,k,m)  = f2*(aa(nx,2,k,m)    + aa(nx-1,1,k,m))
            aa(nx,ny,k,m) = f2*(aa(nx,ny-1,k,m) + aa(nx-1,ny,k,m))
 30      continue

c corners

      do 40 m=1,3
         aa(1,1,1,m) = f3*(aa(1,1,2,m) + aa(1,2,1,m) + aa(2,1,1,m))
         aa(1,1,nz,m) = f3*(aa(1,1,nz-1,m) + aa(1,2,nz,m) 
     x                      + aa(2,1,nz,m))
         aa(1,ny,1,m) = f3*(aa(1,ny,2,m) + aa(1,ny-1,1,m)
     x                      + aa(2,ny,1,m))
         aa(nx,1,1,m) = f3*(aa(nx,1,2,m) + aa(nx,2,1,m) 
     x                      + aa(nx-1,1,1,m))
         aa(1,ny,nz,m) = f3*(aa(1,ny,nz-1,m) + aa(1,ny-1,nz,m)
     x                       + aa(2,ny,nz,m))
         aa(nx,1,nz,m) = f3*(aa(nx,1,nz-1,m) + aa(nx,2,nz,m)
     x                       + aa(nx-1,1,nz,m))
         aa(nx,ny,1,m) = f3*(aa(nx,ny,2,m) + aa(nx,ny-1,1,m)
     x                       + aa(nx-1,ny,1,m))
         aa(nx,ny,nz,m) = f3*(aa(nx,ny,nz-1,m) + aa(nx,ny-1,nz,m)
     x                        + aa(nx-1,ny,nz,m))

 40      continue

      return
      end
c---------------------------------------------------------------------


c---------------------------------------------------------------------
      SUBROUTINE edges_corners_2(aa)
c---------------------------------------------------------------------
      include 'incurv.h'

      real aa(nx,ny,nz,3)
      real f3

      f3=(1.0/3.0)

c x edge 

      do 10 i=3,nx-1
         do 10 m=1,3
            aa(i,2,2,m)   = 0.5*(aa(i,3,2,m)     + aa(i,2,3,m))
            aa(i,2,nz,m)  = 0.5*(aa(i,3,nz,m)    + aa(i,2,nz-1,m))
            aa(i,ny,2,m)  = 0.5*(aa(i,ny,3,m)    + aa(i,ny-1,2,m))
            aa(i,ny,nz,m) = 0.5*(aa(i,ny,nz-1,m) + aa(i,ny-1,nz,m))
 10      continue

c y edge

      do 20 j=3,ny-1
         do 20 m=1,3
            aa(2,j,2,m)   = 0.5*(aa(3,j,2,m)     + aa(2,j,3,m))
            aa(2,j,nz,m)  = 0.5*(aa(3,j,nz,m)    + aa(2,j,nz-1,m))
            aa(nx,j,2,m)  = 0.5*(aa(nx,j,3,m)    + aa(nx-1,j,2,m))
            aa(nx,j,nz,m) = 0.5*(aa(nx,j,nz-1,m) + aa(nx-1,j,nz,m))
 20      continue

c z edge

      do 30 k=3,nz-1
         do 30 m=1,3
            aa(2,2,k,m)   = 0.5*(aa(3,2,k,m)     + aa(2,3,k,m))
            aa(2,ny,k,m)  = 0.5*(aa(3,ny,k,m)    + aa(2,ny-1,k,m))
            aa(nx,2,k,m)  = 0.5*(aa(nx,3,k,m)    + aa(nx-1,2,k,m))
            aa(nx,ny,k,m) = 0.5*(aa(nx,ny-1,k,m) + aa(nx-1,ny,k,m))
 30      continue

c corners

      do 40 m=1,3
         aa(2,2,2,m) = f3*(aa(2,2,3,m) + aa(2,3,2,m) + aa(3,2,2,m))
         aa(2,2,nz,m) = f3*(aa(2,2,nz-1,m) + aa(2,3,nz,m) 
     x                      + aa(3,2,nz,m))
         aa(2,ny,2,m) = f3*(aa(2,ny,3,m) + aa(2,ny-1,2,m)
     x                      + aa(3,ny,2,m))
         aa(nx,2,2,m) = f3*(aa(nx,2,3,m) + aa(nx,3,2,m) 
     x                      + aa(nx-1,2,2,m))
         aa(2,ny,nz,m) = f3*(aa(2,ny,nz-1,m) + aa(2,ny-1,nz,m)
     x                       + aa(3,ny,nz,m))
         aa(nx,2,nz,m) = f3*(aa(nx,2,nz-1,m) + aa(nx,3,nz,m)
     x                       + aa(nx-1,2,nz,m))
         aa(nx,ny,2,m) = f3*(aa(nx,ny,3,m) + aa(nx,ny-1,2,m)
     x                       + aa(nx-1,ny,2,m))
         aa(nx,ny,nz,m) = f3*(aa(nx,ny,nz-1,m) + aa(nx,ny-1,nz,m)
     x                        + aa(nx-1,ny,nz,m))

 40      continue

      return
      end
c---------------------------------------------------------------------


c---------------------------------------------------------------------
      SUBROUTINE normal_B(b)  !divergence B
c---------------------------------------------------------------------
      include 'incurv.h'

      real b(nx,ny,nz,3)

c copy normal components to boundaries first
c x surfaces
      do 4 j=2,ny-1
         do 4 k=2,nz-1

            b(2,j,k,1) = b(3,j,k,1)       !normal
c            b(1,j,k,1) = b(3,j,k,1)
            b(nx,j,k,1) = b(nx-1,j,k,1)

 4       continue

c y surfaces
       do 6 i=2,nx-1
          do 6 k=2,nz-1

             b(i,2,k,2) = b(i,3,k,2)      !normal
c             b(i,1,k,2) = b(i,3,k,2)
             b(i,ny,k,2) = b(i,ny-1,k,2)

 6        continue

c z surfaces
       do 8 i=2,nx-1
          do 8 j=2,ny-1

             b(i,j,2,3) = b(i,j,3,3)      !normal
c             b(i,j,1,3) = b(i,j,3,3)
             b(i,j,nz,3) = b(i,j,nz-1,3)

 8        continue

      call edges_corners_2(b)

c normal x components
      do 10 j=2,ny-1
         do 10 k=2,nz-1
            b(2,j,k,1) = dx*(b(2,j+1,k,2) - b(2,j,k,2))/dy +
     x                    dx*(b(2,j,k+1,3) - b(2,j,k,3))/dz_grid(k) +
     x                    b(3,j,k,1)

            b(nx,j,k,1) = b(nx-1,j,k,1) -
     x                  dx*(b(nx-1,j+1,k,2) - b(nx-1,j,k,2))/dy -
     x                  dx*(b(nx-1,j,k+1,3) - b(nx-1,j,k,3))/dz_grid(k)
 10         continue

c normal y components
      do 20 i=2,nx-1
         do 20 k=2,nz-1
            b(i,2,k,2) = dy*(b(i+1,2,k,1) - b(i,2,k,1))/dx + 
     x                    dy*(b(i,2,k+1,3) - b(i,2,k,3))/dz_grid(k) + 
     x                    b(i,3,k,2)

            b(i,ny,k,2) = b(i,ny-1,k,2) -
     x                  dy*(b(i+1,ny-1,k,1) - b(i,ny-1,k,1))/dx -
     x                  dy*(b(i,ny-1,k+1,3) - b(i,ny-1,k,3))/dz_grid(k)
 20         continue

c normal z components
      do 30 i=2,nx-1
         do 30 j=2,ny-1
            b(i,j,2,3) = dz_grid(2)*(b(i+1,j,2,1) - b(i,j,2,1))/dx + 
     x                   dz_grid(2)*(b(i,j+1,2,2) - b(i,j,2,2))/dy +
     x                   b(i,j,3,3)

            b(i,j,nz,3) = b(i,j,nz-1,3) -
     x               dz_grid(nz)*(b(i+1,j,nz-1,1) - b(i,j,nz-1,1))/dx -
     x               dz_grid(nz)*(b(i,j+1,nz-1,2) - b(i,j,nz-1,2))/dy

 30         continue

      return 
      end
c---------------------------------------------------------------------


c---------------------------------------------------------------------
      SUBROUTINE tangent_B_zero(b) !normal derivative = 0
c The normal derivative of the tangential components is used to 
c determine the tangential boundary values.  ALSO, the normal
c derivative of the normal components temporarily sets the boundary
c values of the normal components.  These values are undated later
c by the requirement that divB=0.  This helps with the values on
c corners and edges.
c---------------------------------------------------------------------
      include 'incurv.h'

      real b(nx,ny,nz,3)

c x surfaces
      do 10 j=2,ny-1
         do 10 k=2,nz-1

            b(1,j,k,1) = b(2,j,k,1)       !normal
c            b(1,j,k,1) = b(3,j,k,1)
c            b(nx,j,k,1) = b(nx-1,j,k,1)

            b(1,j,k,2) = b(2,j,k,2)       !tangential
            b(1,j,k,3) = b(2,j,k,3)
            b(nx,j,k,2) = b(nx-1,j,k,2)
            b(nx,j,k,3) = b(nx-1,j,k,3)


 10         continue

c y surfaces
       do 20 i=2,nx-1
          do 20 k=2,nz-1

             b(i,1,k,2) = b(i,2,k,2)      !normal
c             b(i,1,k,2) = b(i,3,k,2)
c             b(i,ny,k,2) = b(i,ny-1,k,2)

             b(i,1,k,1) = b(i,2,k,1)      !tangential
             b(i,1,k,3) = b(i,2,k,3)
             b(i,ny,k,1) = b(i,ny-1,k,1)
             b(i,ny,k,3) = b(i,ny-1,k,3)

 20          continue

c z surfaces
       do 30 i=2,nx-1
          do 30 j=2,ny-1

             b(i,j,1,3) = b(i,j,2,3)      !normal
c             b(i,j,1,3) = b(i,j,3,3)
c             b(i,j,nz,3) = b(i,j,nz-1,3)

             b(i,j,1,1) = b(i,j,2,1)      !tangential
             b(i,j,1,2) = b(i,j,2,2)
             b(i,j,nz,1) = b(i,j,nz-1,1)
             b(i,j,nz,2) = b(i,j,nz-1,2)
             
 30          continue


      return
      end
c---------------------------------------------------------------------


c---------------------------------------------------------------------
      SUBROUTINE copy_to_boundary(b)
c---------------------------------------------------------------------
      include 'incurv.h'

      real b(nx,ny,nz,3)

c x surfaces
      do 10 j=1,ny
         do 10 k=1,nz
            do 10 m=1,3
               b(1,j,k,m) = b(2,j,k,m)       !tangential
               b(nx-1,j,k,m) = b(nx-2,j,k,m)
               b(nx,j,k,m) = b(nx-1,j,k,m)
 10         continue

c y surfaces
       do 20 i=1,nx
          do 20 k=1,nz
             do 20 m=1,3
                b(i,1,k,m) = b(i,2,k,m)      !tangential
                b(i,ny-1,k,m) = b(i,ny-2,k,m)
                b(i,ny,k,m) = b(i,ny-1,k,m)
 20          continue

c z surfaces
       do 30 i=1,nx
          do 30 j=1,ny
             do 30 m=1,3
                b(i,j,1,m) = b(i,j,2,m)      !tangential
                b(i,j,nz-1,m) = b(i,j,nz-2,m)
                b(i,j,nz,m) = b(i,j,nz-1,m)
 30          continue

      return
      end
c---------------------------------------------------------------------


c---------------------------------------------------------------------
      SUBROUTINE fix_normal_b(b)
c---------------------------------------------------------------------
      include 'incurv.h'

      real b(nx,ny,nz,3)

c normal x components
      do 10 j=2,ny-1
         do 10 k=2,nz-1
            b(2,j,k,1) = dx*(b(2,j+1,k,2) - b(2,j,k,2))/dy +
     x                    dx*(b(2,j,k+1,3) - b(2,j,k,3))/dz_grid(k) +
     x                    b(3,j,k,1)

            b(nx,j,k,1) = b(nx-1,j,k,1) -
     x                  dx*(b(nx-1,j+1,k,2) - b(nx-1,j,k,2))/dy -
     x                  dx*(b(nx-1,j,k+1,3) - b(nx-1,j,k,3))/dz_grid(k)
 10         continue

c normal y components
      do 20 i=2,nx-1
         do 20 k=2,nz-1
            b(i,2,k,2) = dy*(b(i+1,2,k,1) - b(i,2,k,1))/dx + 
     x                    dy*(b(i,2,k+1,3) - b(i,2,k,3))/dz_grid(k) + 
     x                    b(i,3,k,2)

            b(i,ny,k,2) = b(i,ny-1,k,2) -
     x                  dy*(b(i+1,ny-1,k,1) - b(i,ny-1,k,1))/dx -
     x                  dy*(b(i,ny-1,k+1,3) - b(i,ny-1,k,3))/dz_grid(k)
 20         continue

c normal z components
      do 30 i=2,nx-1
         do 30 j=2,ny-1
            b(i,j,2,3) = dz_grid(2)*(b(i+1,j,2,1) - b(i,j,2,1))/dx + 
     x                   dz_grid(2)*(b(i,j+1,2,2) - b(i,j,2,2))/dy +
     x                   b(i,j,3,3)

            b(i,j,nz,3) = b(i,j,nz-1,3) -
     x               dz_grid(nz)*(b(i+1,j,nz-1,1) - b(i,j,nz-1,1))/dx -
     x               dz_grid(nz)*(b(i,j+1,nz-1,2) - b(i,j,nz-1,2))/dy

 30         continue


      return
      end
c---------------------------------------------------------------------

c---------------------------------------------------------------------
      SUBROUTINE smooth_boundary(b)
c---------------------------------------------------------------------
      include 'incurv.h'

      real b(nx,ny,nz,3)

c x surfaces
      do 10 j=2,ny-1
         do 10 k=2,nz-1
            do 10 m=1,3
               b(1,j,k,m) = b(2,j,k,m)       !tangential
               b(nx,j,k,m) = b(nx-1,j,k,m)
 10         continue

c y surfaces
       do 20 i=1,nx
          do 20 k=1,nz
             do 20 m=1,3
                b(i,1,k,m) = b(i,2,k,m)      !tangential
                b(i,ny,k,m) = b(i,ny-1,k,m)
 20          continue

c z surfaces
       do 30 i=1,nx
          do 30 j=1,ny
             do 30 m=1,3
                b(i,j,1,m) = b(i,j,2,m)      !tangential
                b(i,j,nz,m) = b(i,j,nz-1,m)
 30          continue

      return
      end
c---------------------------------------------------------------------


c---------------------------------------------------------------------
      SUBROUTINE fix_tangential_E(E)
c---------------------------------------------------------------------
      include 'incurv.h'

      real E(nx,ny,nz,3)

c     i = 2 & i = nx
      do 10 j=2,ny
         do 10 k=2,nz
            E(2,j,k,2) = E(3,j,k,2)
            E(2,j,k,3) = E(3,j,k,3)
c            E(nx-1,j,k,2) = E(nx-2,j,k,2)
c            E(nx-1,j,k,3) = E(nx-2,j,k,3)
            E(nx,j,k,2) = E(nx-1,j,k,2)
            E(nx,j,k,3) = E(nx-1,j,k,3)
 10         continue

c     j = 2 & j = ny
      do 20 i=2,nx
         do 20 k=2,nz
            E(i,2,k,1) = E(i,3,k,1)
            E(i,2,k,3) = E(i,3,k,3)
c            E(i,ny-1,k,1) = E(i,ny-2,k,1)
c            E(i,ny-1,k,3) = E(i,ny-2,k,3)
            E(i,ny,k,1) = E(i,ny-1,k,1)
            E(i,ny,k,3) = E(i,ny-1,k,3)
 20         continue

c     k = 2 & k = nz
      do 30 i=2,nx
         do 30 j=2,ny
c            E(i,j,2,1) = E(i,j,3,1)
c            E(i,j,2,2) = E(i,j,3,2)
c            E(i,j,nz,1) = E(i,j,nz-1,1)
c            E(i,j,nz,2) = E(i,j,nz-1,2)
            E(i,j,2,1) = 0.0
            E(i,j,2,2) = 0.0
c            E(i,j,nz-1,1) = 0.0
c            E(i,j,nz-1,2) = 0.0
            E(i,j,nz,1) = 0.0
            E(i,j,nz,2) = 0.0
 30         continue


cc i = 2 & i = nx
c      do 40 j=2,ny 
c         do 40 k=2,nz
c               E(2,j,k,2) = 3.0*E(3,j,k,2) - 3.0*E(4,j,k,2)  
c     x                          + E(5,j,k,2)
c               E(nx,j,k,2) = 3.0*E(nx-1,j,k,2) - 3.0*E(nx-2,j,k,2)  
c     x                           + E(nx-3,j,k,2)
c               E(2,j,k,3) = 3.0*E(3,j,k,3) - 3.0*E(4,j,k,3)  
c     x                          + E(5,j,k,3)
c               E(nx,j,k,3) = 3.0*E(nx-1,j,k,3) - 3.0*E(nx-2,j,k,3)  
c     x                           + E(nx-3,j,k,3)
c 40            continue


cc j = 2 & j = ny
c      do 50 i=2,nx
c         do 50 k=2,nz
c               E(i,2,k,1) = 3.0*E(i,3,k,1) - 3.0*E(i,4,k,1) 
c     x                          + E(i,5,k,1)
c               E(i,ny,k,1) = 3.0*E(i,ny-1,k,1) - 3.0*E(i,ny-2,k,1)
c     x                           + E(i,ny-3,k,1)
c               E(i,2,k,3) = 3.0*E(i,3,k,3) - 3.0*E(i,4,k,3) 
c     x                          + E(i,5,k,3)
c               E(i,ny,k,3) = 3.0*E(i,ny-1,k,3) - 3.0*E(i,ny-2,k,3)
c     x                           + E(i,ny-3,k,3)
c 50            continue    

cc k = 2 & k = nz
c      do 60 i=2,nx
c         do 60 j=2,ny
c               E(i,j,2,1) = 3.0*E(i,j,3,1) - 3.0*E(i,j,4,1)
c     x                          + E(i,j,5,1)
c               E(i,j,nz,1) = 3.0*E(i,j,nz-1,1) - 3.0*E(i,j,nz-2,1)
c     x                           + E(i,j,nz-3,1)
c               E(i,j,2,2) = 3.0*E(i,j,3,2) - 3.0*E(i,j,4,2)
c     x                          + E(i,j,5,2)
c               E(i,j,nz,2) = 3.0*E(i,j,nz-1,2) - 3.0*E(i,j,nz-2,2)
c     x                           + E(i,j,nz-3,2)
c 60            continue

      return
      end
c---------------------------------------------------------------------


c---------------------------------------------------------------------
      SUBROUTINE boundaries(b)
c---------------------------------------------------------------------
      include 'incurv.h'

      real b(nx,ny,nz,3)

c      call periodic(aa)
c      call extrapolated(aa)
c      call edges_corners(aa)

c      call normal_B(b)
c      call tangent_B_zero(b)
cc      call edges_corners_2(b)
c      call edges_corners_1(b)

      call copy_to_boundary(b)
      call fix_normal_b(b)
c      call smooth_boundary(b)

      return
      end
c---------------------------------------------------------------------

















