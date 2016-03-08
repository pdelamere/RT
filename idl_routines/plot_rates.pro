 
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


;mass 1
close,1
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_18'  
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

!x.title='time (s)'
!p.multi=[0,1,1]
;plot,tm,mix_arr,ytitle='Mixing/(Lo*nx)';,/ylog,yrange=[100,1e5]
;p1 = !P & x1 = !X & y1 = !Y
plot,tm,lnvz,yrange=[3.0,6],xrange=[25,500],/ysty,/xsty,$
    ytitle='ln(v!dmax!n)'
p2 = !P & x2 = !X & y2 = !Y

wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)

f = fit(1)
s = sigma(1)

oplot,tm,fit(0)+fit(1)*tm

;mass 2
close,1
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_17'
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

;!p.multi=[0,1,2]
;!P = p1 & !X = x1 & !Y = y1
;oplot,tm,mix_arr,linestyle=1
!P = p2 & !X = x2 & !Y = y2
oplot,tm,lnvz,linestyle=1,color=c1


wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)

oplot,tm,fit(0)+fit(1)*tm,color=c1

f = [f,fit(1)]
s = [s,sigma(1)]


;mass 3
close,1
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_23'
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

;!p.multi=[0,1,2]
;!P = p1 & !X = x1 & !Y = y1
;oplot,tm,mix_arr,linestyle=1
!P = p2 & !X = x2 & !Y = y2
oplot,tm,lnvz,linestyle=1,color=c1


wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)

oplot,tm,fit(0)+fit(1)*tm,color=c1

f = [f,fit(1)]
s = [s,sigma(1)]


;mass 4

close,1
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_16'
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



;!p.multi=[0,1,2]
;!P = p1 & !X = x1 & !Y = y1
;oplot,tm,mix_arr,linestyle=2
!P = p2 & !X = x2 & !Y = y2
oplot,tm,lnvz,linestyle=2,color=c2


wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)

oplot,tm,fit(0)+fit(1)*tm,color=c2

f = [f,fit(1)]
s = [s,sigma(1)]

;mass 5

close,1
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_24'
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



;!p.multi=[0,1,2]
;!P = p1 & !X = x1 & !Y = y1
;oplot,tm,mix_arr,linestyle=2
!P = p2 & !X = x2 & !Y = y2
oplot,tm,lnvz,linestyle=2,color=c2


wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)

oplot,tm,fit(0)+fit(1)*tm,color=c2

f = [f,fit(1)]
s = [s,sigma(1)]


;mass 6
close,1
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_19'
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

;!p.multi=[0,1,2]
;!P = p1 & !X = x1 & !Y = y1
;oplot,tm,mix_arr,linestyle=3
!P = p2 & !X = x2 & !Y = y2
oplot,tm,lnvz,linestyle=3,color=c3

wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)

oplot,tm,fit(0)+fit(1)*tm,color=c3

f = [f,fit(1)]
s = [s,sigma(1)]


;mass 7

close,1
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_20'
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



;!p.multi=[0,1,2]
;!P = p1 & !X = x1 & !Y = y1
;oplot,tm,mix_arr,linestyle=4
!P = p2 & !X = x2 & !Y = y2
oplot,tm,lnvz,linestyle=4,color=c4

wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)

oplot,tm,fit(0)+fit(1)*tm,color=c4

f = [f,fit(1)]
s = [s,sigma(1)]


;mass 8

close,1
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_15'
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

wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)

;!p.multi=[0,1,2]
;!P = p1 & !X = x1 & !Y = y1
;oplot,tm,mix_arr,linestyle=5
!P = p2 & !X = x2 & !Y = y2
oplot,tm,lnvz,linestyle=5,color=c5

oplot,tm,fit(0)+fit(1)*tm,color=c5
;oplot,tm,fit(0)+(fit(1)+sigma(1))*tm,color=c5
;oplot,tm,fit(0)+(fit(1)-sigma(1))*tm,color=c5

f = [f,fit(1)]
s = [s,sigma(1)]

;mass 9

close,1
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_22'
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

wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)

;!p.multi=[0,1,2]
;!P = p1 & !X = x1 & !Y = y1
;oplot,tm,mix_arr,linestyle=5
!P = p2 & !X = x2 & !Y = y2
oplot,tm,lnvz,linestyle=5,color=c8

oplot,tm,fit(0)+fit(1)*tm,color=c8
;oplot,tm,fit(0)+(fit(1)+sigma(1))*tm,color=c5
;oplot,tm,fit(0)+(fit(1)-sigma(1))*tm,color=c5

f = [f,fit(1)]
s = [s,sigma(1)]


;mass 10

close,1
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_21'
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

   lnvz(i-1) = alog(max(abs(smooth(reform(up(*,1,*,2)),2))))

endfor



;!p.multi=[0,1,2]
;!P = p1 & !X = x1 & !Y = y1
;oplot,tm,mix_arr,linestyle=6
!P = p2 & !X = x2 & !Y = y2
oplot,tm,lnvz,linestyle=7,color=c7

wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)
;fit = ladfit(tm(wh),lnvz(wh),absdev=sigma)
oplot,tm,fit(0)+fit(1)*tm,color=c7

f = [f,fit(1)]
s = [s,sigma(1)]


;mass 11

close,1
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_25'
f_read_coord,rundir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz

nfrm = 10

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

   lnvz(i-1) = alog(max(abs(smooth(reform(up(*,1,*,2)),2))))

endfor



;!p.multi=[0,1,2]
;!P = p1 & !X = x1 & !Y = y1
;oplot,tm,mix_arr,linestyle=6
!P = p2 & !X = x2 & !Y = y2
oplot,tm,lnvz,linestyle=7,color=c7

wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)
;fit = ladfit(tm(wh),lnvz(wh),absdev=sigma)
oplot,tm,fit(0)+fit(1)*tm,color=c7

f = [f,fit(1)]
s = [s,sigma(1)]

;mass 12

close,1
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_14'
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

   lnvz(i-1) = alog(max(abs(smooth(reform(up(*,1,*,2)),2))))

endfor



;!p.multi=[0,1,2]
;!P = p1 & !X = x1 & !Y = y1
;oplot,tm,mix_arr,linestyle=6
!P = p2 & !X = x2 & !Y = y2
oplot,tm,lnvz,linestyle=6,color=c6

wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)
;fit = ladfit(tm(wh),lnvz(wh),absdev=sigma)
oplot,tm,fit(0)+fit(1)*tm,color=c6

f = [f,fit(1)]
s = [s,sigma(1)]

legend,['1','2','4','6','7','8','12'],linestyle=[0,1,2,3,4,5,6],/bottom,/right

;stop

ff = (2*4*300./500)
mass = [1,2,3,4,5,6,7,8,9,10,11,12]
plot,mass,f*ff,psym=1,xrange=[0,13],/xsty,yrange=ff*[0.006,0.015],/ysty
errplot,mass,(f+s)*ff,(f-s)*ff

fit = poly_fit(mass,f,1)

oplot,mass,fit(0)*ff+fit(1)*ff*mass

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
