
c----------------------------------------------------------------------
      SUBROUTINE grd7()
c----------------------------------------------------------------------
      include 'incurv.h'

      parameter(nrgrd = 5)

      rk=nz/2
      rj=ny/2
      ri=nx-15

      do 10 i=1,nx
         qx(i) = i*dx
 10            continue

      do 20 j=1,ny
         qy(j) = j*dy
 20            continue

c up from release
      do 32 k = rk,rk+nrgrd
         dz_grid(k) = delz
 32   continue
c CDIR@ NEXTSCALAR
      do 34 k = rk+nrgrd+1,nz
         dz_grid(k) = delz +
     x     0.0*(k-(rk+nrgrd+1))/(nz-(rk+nrgrd+1)) 
c     x     2.0*sin((k-(rk+nrgrd+1))*0.5*pi/(nz-(rk+nrgrd+1)))**2 
c                                !dz_grid(k-1) + 0.01*delz 
 34   continue

c down from release
      do 36 k = rk-nrgrd,rk-1
         dz_grid(k) = delz
 36      continue
c CDIR@ NEXTSCALAR
      do 37 k = 1,rk-nrgrd-1
         ind = rk-nrgrd-k
         dz_grid(ind) = delz + 
     x     0.0*(rk-nrgrd-1-ind)/(rk-nrgrd-1)
c     x     2.0*sin((rk-nrgrd-1-ind)*(-0.5*pi)/(rk-nrgrd-1))**2 
c                                !dz_grid(ind+1) + 0.01*delz
 37   continue

      qz(1) = 0.0
      do 39 k=2,nz
c         write(*,*) 'dz_grid...',k,dz_grid(k)
         qz(k) = qz(k-1)+dz_grid(k)
 39   continue

      dz_cell(1) = dz_grid(1)
      dz_cell(nz) = dz_grid(nz)
      do 40 k=2,nz-1
         dz_cell(k) = ((qz(k+1) + qz(k))/2.0) -
     x                ((qz(k) + qz(k-1))/2.0)
 40            continue

c      call assign('assign -F system -N ultrix f:' //'c.coord.dat')
      open(40,file='c.coord.dat',status='unknown',
     x         form='unformatted')

c      open(40,file='coord.dat',status='unknown',form='unformatted')
      write(40) nx
      write(40) ny
      write(40) nz
      write(40) qx
      write(40) qy
      write(40) qz
      write(40) dz_grid
      write(40) dz_cell
      close(40)

      return
      end
c----------------------------------------------------------------------


      PROGRAM write_coords

      include 'incurv.h'

      call grd7()

      stop 
      end
