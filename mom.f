c----------------------------------------------------------------------
      SUBROUTINE get_graduu_intgl(uf,nf)
c----------------------------------------------------------------------
      include 'incurv.h'

      real uf(nx,ny,nz,3)         !uf at cell face
      real nf(nx,ny,nz)
      real ufc(nx,ny,nz,3)        !uf at cell center

      real ax1,ay1,az1                   !misc const
      real u1,u2,u3                      !temp vars

      call face_to_center(uf,ufc)

      graduu_tot(1) = 0.0
      graduu_tot(2) = 0.0
      graduu_tot(3) = 0.0

      do 10 i=3,nx-1
         do 10 j=3,ny-1
            do 10 k=3,nz-1
c               ip=i+1
c               jp=j+1
c               kp=k+1
               im=i-1
               jm=j-1
               km=k-1
c               if (ip .gt. nx-1) then ip = nx-1
c               if (jp .gt. ny-1) then jp = ny-1
c               if (kp .gt. nz-1) then kp = nz-1
c               if (im .eq. 1) then im = 2
c               if (jm .eq. 1) then jm = 2
c               if (km .eq. 1) then km = 2
               vol = dx*dy*dz_cell(k)

c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

               ax1 = ufc(i,j,k,1)/dx
               ay1 = ufc(i,j,k,1)/dy
               u1 = ax1*(uf(i,j,k,1) - uf(im,j,k,1)) 
               u2 = ay1*(uf(i,j,k,2) - uf(i,jm,k,2))  
               az1 = ufc(i,j,k,1)/dz_cell(k)
               u3 = az1*(uf(i,j,k,3)-uf(i,j,km,3)) 
               
               graduu_tot(1) = graduu_tot(1) + 
     x                         (mO*nf(i,j,k)*vol)*(u1 + u2 + u3)

c yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy

               ax1 = ufc(i,j,k,2)/dx
               ay1 = ufc(i,j,k,2)/dy
               u1 = ax1*(uf(i,j,k,1) - uf(im,j,k,1)) 
               u2 = ay1*(uf(i,j,k,2) - uf(i,jm,k,2))  
               az1 = ufc(i,j,k,2)/dz_cell(k)
               u3 = az1*(uf(i,j,k,3)-uf(i,j,km,3)) 

               graduu_tot(2) = graduu_tot(2) + 
     x                         (mO*nf(i,j,k)*vol)*(u1 + u2 + u3)


c zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz

               ax1 = ufc(i,j,k,3)/dx
               ay1 = ufc(i,j,k,3)/dy
               u1 = ax1*(uf(i,j,k,1) - uf(im,j,k,1)) 
               u2 = ay1*(uf(i,j,k,2) - uf(i,jm,k,2))  
               az1 = ufc(i,j,k,3)/dz_cell(k)
               u3 = az1*(uf(i,j,k,3)-uf(i,j,km,3)) 

               graduu_tot(3) = graduu_tot(3) + 
     x                         (mO*nf(i,j,k)*vol)*(u1 + u2 + u3)


 10            continue

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE get_surface_intgl(b1,uf,nf)
c----------------------------------------------------------------------
      include 'incurv.h'

      real b1(nx,ny,nz,3)
      real uf(nx,ny,nz,3)
      real nf(nx,ny,nz)

      real b1c(nx,ny,nz,3)
      real ufc(nx,ny,nz,3)

      real Tm(3,3)
      real B2
      real ax,ay,az

      call face_to_center(b1,b1c)
      call face_to_center(uf,ufc)

      surf_tot(1) = 0.0
      surf_tot(2) = 0.0
      surf_tot(3) = 0.0

c x sides  

      do 10 j=2,ny-1
         do 10 k=2,nz-1

            i = 2

            ax = dy*dz_cell(k)
            B2 = b1(i,j,k,1)**2 + b1(i,j,k,2)**2 + b1(i,j,k,3)**2

            Tm(1,1) = (mO/alpha)*(b1(i,j,k,1)*b1(i,j,k,1) - 0.5*B2)
     x               - mO*nf(i,j,k)*uf(i,j,k,1)*uf(i,j,k,1)
            Tm(1,2) = (mO/alpha)*(b1(i,j,k,1)*b1(i,j,k,2))
     x               - mO*nf(i,j,k)*uf(i,j,k,1)*uf(i,j,k,2)
            Tm(1,3) = (mO/alpha)*(b1(i,j,k,1)*b1(i,j,k,3))
     x               - mO*nf(i,j,k)*uf(i,j,k,1)*uf(i,j,k,3)

            surf_tot(1) = surf_tot(1) - ax*Tm(1,1)                    
            surf_tot(2) = surf_tot(2) - ax*Tm(1,2)                    
            surf_tot(3) = surf_tot(3) - ax*Tm(1,3)                    

            i = nx-1

c            ax = dy*dz_cell(k)
c            B2 = b1(i,j,k,1)**2 + b1(i,j,k,2)**2 + b1(i,j,k,3)**2

            Tm(1,1) = (mO/alpha)*(b1(i,j,k,1)*b1(i,j,k,1) - 0.5*B2)
     x               - mO*nf(i,j,k)*uf(i,j,k,1)*uf(i,j,k,1)
            Tm(1,2) = (mO/alpha)*(b1(i,j,k,1)*b1(i,j,k,2))
     x               - mO*nf(i,j,k)*uf(i,j,k,1)*uf(i,j,k,2)
            Tm(1,3) = (mO/alpha)*(b1(i,j,k,1)*b1(i,j,k,3))
     x               - mO*nf(i,j,k)*uf(i,j,k,1)*uf(i,j,k,3)

            surf_tot(1) = surf_tot(1) + ax*Tm(1,1)                    
            surf_tot(2) = surf_tot(2) + ax*Tm(1,2)                    
            surf_tot(3) = surf_tot(3) + ax*Tm(1,3)                    


 10         continue

c y sides  

      do 20 i=2,nx-1
         do 20 k=2,nz-1

            j = 2

            ay = dx*dz_cell(k)
            B2 = b1(i,j,k,1)**2 + b1(i,j,k,2)**2 + b1(i,j,k,3)**2

            Tm(2,1) = (mO/alpha)*(b1(i,j,k,2)*b1(i,j,k,1))
     x               - mO*nf(i,j,k)*uf(i,j,k,2)*uf(i,j,k,1)
            Tm(2,2) = (mO/alpha)*(b1(i,j,k,2)*b1(i,j,k,2) - 0.5*B2)
     x               - mO*nf(i,j,k)*uf(i,j,k,2)*uf(i,j,k,2)
            Tm(2,3) = (mO/alpha)*(b1(i,j,k,2)*b1(i,j,k,3))
     x               - mO*nf(i,j,k)*uf(i,j,k,2)*uf(i,j,k,3)

            surf_tot(1) = surf_tot(1) - ay*Tm(2,1)                    
            surf_tot(2) = surf_tot(2) - ay*Tm(2,2)                    
            surf_tot(3) = surf_tot(3) - ay*Tm(2,3)                    

            j = ny-1

c            ay = dx*dz_cell(k)
c            B2 = b1(i,j,k,1)**2 + b1(i,j,k,2)**2 + b1(i,j,k,3)**2

            Tm(2,1) = (mO/alpha)*(b1(i,j,k,2)*b1(i,j,k,1))
     x               - mO*nf(i,j,k)*uf(i,j,k,2)*uf(i,j,k,1)
            Tm(2,2) = (mO/alpha)*(b1(i,j,k,2)*b1(i,j,k,2) - 0.5*B2)
     x               - mO*nf(i,j,k)*uf(i,j,k,2)*uf(i,j,k,2)
            Tm(2,3) = (mO/alpha)*(b1(i,j,k,2)*b1(i,j,k,3))
     x               - mO*nf(i,j,k)*uf(i,j,k,2)*uf(i,j,k,3)

            surf_tot(1) = surf_tot(1) + ay*Tm(2,1)                    
            surf_tot(2) = surf_tot(2) + ay*Tm(2,2)                    
            surf_tot(3) = surf_tot(3) + ay*Tm(2,3)                    

 20         continue

c z sides  

      do 30 i=2,nx-1
         do 30 j=2,ny-1

            k = 2

            az = dx*dy
            B2 = b1(i,j,k,1)**2 + b1(i,j,k,2)**2 + b1(i,j,k,3)**2

            Tm(3,1) = (mO/alpha)*(b1(i,j,k,3)*b1(i,j,k,1))
     x               - mO*nf(i,j,k)*uf(i,j,k,3)*uf(i,j,k,1)
            Tm(3,2) = (mO/alpha)*(b1(i,j,k,3)*b1(i,j,k,2))
     x               - mO*nf(i,j,k)*uf(i,j,k,3)*uf(i,j,k,2)
            Tm(3,3) = (mO/alpha)*(b1(i,j,k,3)*b1(i,j,k,3) - 0.5*B2)
     x               - mO*nf(i,j,k)*uf(i,j,k,3)*uf(i,j,k,3)

            surf_tot(1) = surf_tot(1) - az*Tm(3,1)                    
            surf_tot(2) = surf_tot(2) - az*Tm(3,2)                    
            surf_tot(3) = surf_tot(3) - az*Tm(3,3)                    

            k = nz-1

c            az = dx*dy
c            B2 = b1(i,j,k,1)**2 + b1(i,j,k,2)**2 + b1(i,j,k,3)**2

            Tm(3,1) = (mO/alpha)*(b1(i,j,k,3)*b1(i,j,k,1))
     x               - mO*nf(i,j,k)*uf(i,j,k,3)*uf(i,j,k,1)
            Tm(3,2) = (mO/alpha)*(b1(i,j,k,3)*b1(i,j,k,2))
     x               - mO*nf(i,j,k)*uf(i,j,k,3)*uf(i,j,k,2)
            Tm(3,3) = (mO/alpha)*(b1(i,j,k,3)*b1(i,j,k,3) - 0.5*B2)
     x               - mO*nf(i,j,k)*uf(i,j,k,3)*uf(i,j,k,3)

            surf_tot(1) = surf_tot(1) + az*Tm(3,1)                    
            surf_tot(2) = surf_tot(2) + az*Tm(3,2)                    
            surf_tot(3) = surf_tot(3) + az*Tm(3,3)                    

 30         continue

      return
      end
c----------------------------------------------------------------------



c----------------------------------------------------------------------
      SUBROUTINE check_momentum(uf,nf,b1,ugradu)
c----------------------------------------------------------------------
      include 'incurv.h'

      real uf(nx,ny,nz,3),
     x     nf(nx,ny,nz),
     x     b1(nx,ny,nz,3),
     x     ugradu(nx,ny,nz,3)

      ugradu_tot(1) = 0.0
      ugradu_tot(2) = 0.0
      ugradu_tot(3) = 0.0

      do 10 i=2,nx
         do 10 j=2,ny
            do 10 k=2,nz
               vol = dx*dy*dz_cell(k)
               do 10 m=1,3 
                  ugradu_tot(m) = ugradu_tot(m) + 
     x                            mO*nf(i,j,k)*vol*ugradu(i,j,k,m)
 10               continue

      call get_graduu_intgl(uf,nf)
      call get_surface_intgl(b1,uf,nf)

c      write(*,*) 'Momentum conservation....'
c      write(*,*) '   Surface intgl T.da....',surf_tot(1),surf_tot(2),
c     x                                      surf_tot(3)
c      write(*,*) '   Volume intgl graduu...',graduu_tot(1),
c     x                  graduu_tot(2),graduu_tot(3)
c      write(*,*) '   Total.................',surf_tot(1)+graduu_tot(1),
c     x            surf_tot(2)+graduu_tot(2),surf_tot(3)+graduu_tot(3)
c      write(*,*) '   Volume intgl ugradu...',ugradu_tot(1),
c     x                 ugradu_tot(2),ugradu_tot(3)

      return
      end
     
c----------------------------------------------------------------------
