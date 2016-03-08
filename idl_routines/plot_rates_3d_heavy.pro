 
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


;files = ['27','28','29','30','31','32','33','34','35','36','37','38','39','40']
;files = ['41','42','43','44','45','46','47','48','49','50','51','52','53','54','55','56','57']
files = ['61','63','64','65','66','67','68','69','70','71']
mass = [1,2,3,4,5,6,7,8,12,16]

;mass 1
close,1
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_'+files(0)  
f_read_coord,rundir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz

nfrm = 30

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

!x.title='mass ratio'
!p.multi=[0,4,4]
!y.range=[3.0,6.0]
!x.range=[25,1200]
;plot,tm,mix_arr,ytitle='Mixing/(Lo*nx)';,/ylog,yrange=[100,1e5]
;p1 = !P & x1 = !X & y1 = !Y
plot,tm,lnvz,yrange=[3.0,6],/ysty,/xsty,$
    ytitle='ln(v!dmax!n)',linestyle=2,title='1'
;p2 = !P & x2 = !X & y2 = !Y

wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)

f = fit(1)
s = sigma(1)

oplot,tm,fit(0)+fit(1)*tm

;mass 2
t1=150
t2=250
close,1
;rundir = '/Volumes/MacD72-2/hybrid/KHI/run_28'
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_'+files(1)  
f_read_coord,rundir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz

nfrm = 30

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
;!P = p2 & !X = x2 & !Y = y2
plot,tm,lnvz,linestyle=2,title='2';,color=c1


wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)

oplot,tm,fit(0)+fit(1)*tm,color=c1

f = [f,fit(1)]
s = [s,sigma(1)]


;mass 3
close,1
;rundir = '/Volumes/MacD72-2/hybrid/KHI/run_29'
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_'+files(2)  
f_read_coord,rundir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz

nfrm = 30

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
;!P = p2 & !X = x2 & !Y = y2
plot,tm,lnvz,linestyle=2,title='3';,color=c3


wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)

oplot,tm,fit(0)+fit(1)*tm,color=c3

f = [f,fit(1)]
s = [s,sigma(1)]


;mass 4
t1 = 200
t2 = 425
close,1
;rundir = '/Volumes/MacD72-2/hybrid/KHI/run_30'
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_'+files(3)  
f_read_coord,rundir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz

nfrm = 30

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
;!P = p2 & !X = x2 & !Y = y2
plot,tm,lnvz,linestyle=2,title='4';,color=c4


wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)

oplot,tm,fit(0)+fit(1)*tm,color=c4

f = [f,fit(1)]
s = [s,sigma(1)]

;mass 5
t1 = 200
t2 = 425
close,1
;rundir = '/Volumes/MacD72-2/hybrid/KHI/run_31'
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_'+files(4)  
f_read_coord,rundir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz

nfrm = 30

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
;!P = p2 & !X = x2 & !Y = y2
plot,tm,lnvz,linestyle=2,title='5';,color=c5


wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)

oplot,tm,fit(0)+fit(1)*tm,color=c5

f = [f,fit(1)]
s = [s,sigma(1)]


;mass 6
t1 = 200
t2 = 425
close,1
;rundir = '/Volumes/MacD72-2/hybrid/KHI/run_31'
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_'+files(5)  
f_read_coord,rundir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz

nfrm = 30

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
;!P = p2 & !X = x2 & !Y = y2
plot,tm,lnvz,linestyle=2,title='5';,color=c5


wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)

oplot,tm,fit(0)+fit(1)*tm,color=c5

f = [f,fit(1)]
s = [s,sigma(1)]


;mass 7
t1 = 200
t2 = 425
close,1
;rundir = '/Volumes/MacD72-2/hybrid/KHI/run_31'
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_'+files(6)  
f_read_coord,rundir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz

nfrm = 30

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
;!P = p2 & !X = x2 & !Y = y2
plot,tm,lnvz,linestyle=2,title='5';,color=c5


wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)

oplot,tm,fit(0)+fit(1)*tm,color=c5

f = [f,fit(1)]
s = [s,sigma(1)]


;mass 8
t1 = 200
t2 = 475
close,1
;rundir = '/Volumes/MacD72-2/hybrid/KHI/run_31'
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_'+files(7)  
f_read_coord,rundir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz

nfrm = 40

lnvz = fltarr(nfrm)
mix_arr = fltarr(nfrm)

nout=50
dt = 0.4
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
;!P = p2 & !X = x2 & !Y = y2
plot,tm,lnvz,linestyle=2,title='5';,color=c5


wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)

oplot,tm,fit(0)+fit(1)*tm,color=c5

f = [f,fit(1)]
s = [s,sigma(1)]

;mass 12
t1 = 250
t2 = 700
close,1
;rundir = '/Volumes/MacD72-2/hybrid/KHI/run_31'
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_'+files(8)  
f_read_coord,rundir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz

nfrm = 80

lnvz = fltarr(nfrm)
mix_arr = fltarr(nfrm)

nout=50
dt = 0.2
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
;!P = p2 & !X = x2 & !Y = y2
plot,tm,lnvz,linestyle=2,title='5';,color=c5


wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)

oplot,tm,fit(0)+fit(1)*tm,color=c5

f = [f,fit(1)]
s = [s,sigma(1)]



;mass 16
t1 = 275
t2 = 700
close,1
;rundir = '/Volumes/MacD72-2/hybrid/KHI/run_31'
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_'+files(9)  
f_read_coord,rundir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz

nfrm = 117

lnvz = fltarr(nfrm)
mix_arr = fltarr(nfrm)

nout=50
dt = 0.2
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
;!P = p2 & !X = x2 & !Y = y2
plot,tm,lnvz,linestyle=2,title='5';,color=c5


wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)

oplot,tm,fit(0)+fit(1)*tm,color=c5

f = [f,fit(1)]
s = [s,sigma(1)]


;;mass 32
;t1 = 500
;t2 = 1000
;close,1
;;rundir = '/Volumes/MacD72-2/hybrid/KHI/run_31'
;rundir = '/Volumes/MacD72-2/hybrid/KHI/run_'+files(10)  
;f_read_coord,rundir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz

;nfrm = 112

;lnvz = fltarr(nfrm)
;mix_arr = fltarr(nfrm)

;nout=100
;dt = 0.1
;tm = dt*nout+findgen(nfrm)*nout*dt

;for i = 1,nfrm do begin 
;   f_read_3d_m_32,rundir+'/mixed_1',i,mixed
;   f_read_3d_vec_m_32,rundir+'/upall_1',i,up

;   wh = where((mixed le 0.75) and (mixed ge 0.25))
;   mix_arr(i-1) = n_elements(wh)/Am

;   lnvz(i-1) = alog(max(abs(smooth(up(*,1,*,2),2))))

;endfor

;;!p.multi=[0,1,2]
;;!P = p1 & !X = x1 & !Y = y1
;;oplot,tm,mix_arr,linestyle=1
;;!P = p2 & !X = x2 & !Y = y2
;plot,tm,lnvz,linestyle=2,title='5';,color=c5


;wh = where((tm ge t1) and (tm le t2))
;fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)

;oplot,tm,fit(0)+fit(1)*tm,color=c5

;f = [f,fit(1)]
;s = [s,sigma(1)]


;legend,['1','2','4','6','7','8','12'],linestyle=[0,,2,3,4,5,6],/bottom,/right

;stop

ff = (2*4*300./500)

plot,mass,f*ff,psym=1,xrange=[0,17],/xsty,yrange=ff*[0.002,0.012],/ysty
errplot,mass,(f+s)*ff,(f-s)*ff

fit = poly_fit(mass,f,1,sigma=sigma)

oplot,mass,fit(0)*ff+fit(1)*ff*mass
oplot,mass,fit(0)*ff-sigma(0)*ff+(fit(1)*ff+sigma(1)*ff)*mass,linestyle=1
oplot,mass,fit(0)*ff+sigma(0)*ff+(fit(1)*ff-sigma(1)*ff)*mass,linestyle=1
print,fit*ff,sigma*ff

@gops
!p.multi=[0,1,1]
device,filename='case3_growth_rates.ps'
;ff = (2*4*300./500)
device,/inches,xsize=6.0,ysize=6.0,xoffset=1
f = f*1e3
s = s*1e3

save,filename='heavy.sav',mass,f

plot,mass,f,psym=1,xrange=[0,17],/xsty,yrange=[0.002,0.012]*1e3,/ysty,$
  ytitle='q (10!u-3!n s!u-1!n)',xtitle='mass (amu)'
errplot,mass,(f+s),(f-s)

fit = poly_fit(mass,f,1,sigma=sigma)

oplot,mass,fit(0)+fit(1)*mass
oplot,mass,fit(0)-sigma(0)+(fit(1)+sigma(1))*mass,linestyle=1
oplot,mass,fit(0)+sigma(0)+(fit(1)-sigma(1))*mass,linestyle=1

;fit = poly_fit(mass,alog(f),1,sigma=sigma)
;fit = fit
;sigma = sigma
;;mass = [mass,17]

;m = findgen(64)/2.

;oplot,m,exp(fit(0))*exp(fit(1)*m)
;oplot,m,exp(fit(0)-sigma(0))*exp((fit(1)+sigma(1))*m),linestyle=1
;oplot,m,exp(fit(0)+sigma(0))*exp((fit(1)-sigma(1))*m),linestyle=1
;;oplot,mass,fit(0)+sigma(0)+(fit(1)-sigma(1))*mass,linestyle=1
device,/close


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
