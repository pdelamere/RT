 
;mixing normalization

dx = 300.
nx = 279.
Am = 4.0*nx
device,decompose=0
c1 = fsc_color('red')
c2 = fsc_color('green')
c3 = fsc_color('blue')
c4 = fsc_color('cyan')
c5 = fsc_color('yellow')
c6 = fsc_color('purple')
c7 = fsc_color('orange')
c8 = fsc_color('Maroon')

t1 = 75
t2 = 275


files = ['27','28','29','30','31','32','33','34','35','36','37','38','39','40']
;files = ['41','42','43','44','45','46','47','48','49','50','51','52','53','54','55','56','57']
;files = ['61','test1','63','64','65','46','47','48','49','50','51','52','53','54','55','56','57']
mass = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]
mass = [1,2,3,4,5,6,7,8,9,10,11,12,13,14]


@gops
!p.multi=[0,1,1]
device,/inches,xsize=6.0,ysize=6.0,xoffset=1

device,filename='case1_mix_profiles.ps'
;device,/color


;mass 1
close,1
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_'+files(0)  
f_read_coord,rundir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz

nfrm = 20

lnvz = fltarr(nfrm)
mix_arr = fltarr(nfrm)

nout=50
dt = 0.5
tm = dt*nout+findgen(nfrm)*nout*dt

for i = 1,nfrm do begin 
   f_read_3d_m_32,rundir+'/mixed_1',i,mixed
   f_read_3d_vec_m_32,rundir+'/upall_1',i,up

   wh = where((mixed le 0.75) and (mixed ge 0.25))
   mix_arr(i-1) = n_elements(wh)/Am

   lnvz(i-1) = alog(max(abs(smooth(up(*,1,*,2),2))))
endfor

plot,tm,mix_arr,linestyle=0,xtitle='t (s)',ytitle='Mix!dtot!n/A!do!n',xrange=[50,500],$
   /xsty

;mass 2
close,1
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_'+files(1)  
f_read_coord,rundir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz

nfrm = 20

lnvz = fltarr(nfrm)
mix_arr = fltarr(nfrm)

nout=50
dt = 0.5
tm = dt*nout+findgen(nfrm)*nout*dt

for i = 1,nfrm do begin 
   f_read_3d_m_32,rundir+'/mixed_1',i,mixed
   f_read_3d_vec_m_32,rundir+'/upall_1',i,up

   wh = where((mixed le 0.75) and (mixed ge 0.25))
   mix_arr(i-1) = n_elements(wh)/Am

   lnvz(i-1) = alog(max(abs(smooth(up(*,1,*,2),2))))
endfor

;oplot,tm,mix_arr,linestyle=1;,color=c1

;mass 4
close,1
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_'+files(3)  
f_read_coord,rundir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz

nfrm = 20

lnvz = fltarr(nfrm)
mix_arr = fltarr(nfrm)

nout=50
dt = 0.5
tm = dt*nout+findgen(nfrm)*nout*dt

for i = 1,nfrm do begin 
   f_read_3d_m_32,rundir+'/mixed_1',i,mixed
   f_read_3d_vec_m_32,rundir+'/upall_1',i,up

   wh = where((mixed le 0.75) and (mixed ge 0.25))
   mix_arr(i-1) = n_elements(wh)/Am

   lnvz(i-1) = alog(max(abs(smooth(up(*,1,*,2),2))))
endfor

oplot,tm,mix_arr,linestyle=4;,color=c2


;mass 6
close,1
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_'+files(5)  
f_read_coord,rundir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz

nfrm = 20

lnvz = fltarr(nfrm)
mix_arr = fltarr(nfrm)

nout=50
dt = 0.5
tm = dt*nout+findgen(nfrm)*nout*dt

for i = 1,nfrm do begin 
   f_read_3d_m_32,rundir+'/mixed_1',i,mixed
   f_read_3d_vec_m_32,rundir+'/upall_1',i,up

   wh = where((mixed le 0.75) and (mixed ge 0.25))
   mix_arr(i-1) = n_elements(wh)/Am

   lnvz(i-1) = alog(max(abs(smooth(up(*,1,*,2),2))))
endfor

;oplot,tm,mix_arr,linestyle=2;,color=c3

;mass 8
close,1
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_'+files(7)  
f_read_coord,rundir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz

nfrm = 20

lnvz = fltarr(nfrm)
mix_arr = fltarr(nfrm)

nout=50
dt = 0.5
tm = dt*nout+findgen(nfrm)*nout*dt

for i = 1,nfrm do begin 
   f_read_3d_m_32,rundir+'/mixed_1',i,mixed
   f_read_3d_vec_m_32,rundir+'/upall_1',i,up

   wh = where((mixed le 0.75) and (mixed ge 0.25))
   mix_arr(i-1) = n_elements(wh)/Am

   lnvz(i-1) = alog(max(abs(smooth(up(*,1,*,2),2))))
endfor

oplot,tm,mix_arr,linestyle=2;,color=c4

;mass 10
close,1
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_'+files(9)  
f_read_coord,rundir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz

nfrm = 20

lnvz = fltarr(nfrm)
mix_arr = fltarr(nfrm)

nout=50
dt = 0.5
tm = dt*nout+findgen(nfrm)*nout*dt

for i = 1,nfrm do begin 
   f_read_3d_m_32,rundir+'/mixed_1',i,mixed
   f_read_3d_vec_m_32,rundir+'/upall_1',i,up

   wh = where((mixed le 0.75) and (mixed ge 0.25))
   mix_arr(i-1) = n_elements(wh)/Am

   lnvz(i-1) = alog(max(abs(smooth(up(*,1,*,2),2))))
endfor

;oplot,tm,mix_arr,linestyle=4;,color=c5


;mass 12
close,1
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_'+files(11)  
f_read_coord,rundir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz

nfrm = 20

lnvz = fltarr(nfrm)
mix_arr = fltarr(nfrm)

nout=50
dt = 0.5
tm = dt*nout+findgen(nfrm)*nout*dt

for i = 1,nfrm do begin 
   f_read_3d_m_32,rundir+'/mixed_1',i,mixed
   f_read_3d_vec_m_32,rundir+'/upall_1',i,up

   wh = where((mixed le 0.75) and (mixed ge 0.25))
   mix_arr(i-1) = n_elements(wh)/Am

   lnvz(i-1) = alog(max(abs(smooth(up(*,1,*,2),2))))
endfor

oplot,tm,mix_arr,linestyle=3;,color=c5



legend,['mass ratio 1','4','8','12'],linestyle=[0,4,2,3],$
  /bottom,/right


device,/close
set_plot,'x

mproton=1.67e-27
m1 = 1.*mproton
m2 = 1.*mproton
n1 = 0.1e6
n2 = 0.1e6

rho1 = m1*n1
rho2 = m2*n2

a1 = rho1/(rho1+rho2)
a2 = rho2/(rho1+rho2)

vrel=2*250e3
b = 5e-9
muo = !pi*4e-7

va1 = b/sqrt(muo*rho1)
va2 = b/sqrt(muo*rho2)

kref = 2*!pi*3/(279*300e3) ;m = 3 mode
phi = 2*!dtor

q = kref*sqrt(a1*a2*vrel^2*cos(phi)^2 - (a1*va1^2 + a2*va2^2)*sin(phi)^2)

print,q,q*2*4*300e3/(2*250e3)



end
