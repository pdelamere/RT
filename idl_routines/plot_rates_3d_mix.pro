 
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
files = ['41','42','43','44','45','46','47','48','49','50','51','52','53','54','55','56','57']
;files = ['61','test1','63','64','65','46','47','48','49','50','51','52','53','54','55','56','57']
mass = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]
;mass = [1,2,3,4,5,6,7,8,9,10,11,12,13,14]

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

!x.title='mass ratio'
!p.multi=[0,4,4]
!y.range=[3.0,6.0]
!x.range=[25,500]
;plot,tm,mix_arr,ytitle='Mixing/(Lo*nx)';,/ylog,yrange=[100,1e5]
;p1 = !P & x1 = !X & y1 = !Y
plot,tm,lnvz,yrange=[3.0,6],xrange=[25,500],/ysty,/xsty,$
    ytitle='ln(v!dmax!n)',linestyle=2,title='1'
;p2 = !P & x2 = !X & y2 = !Y

wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)

f = fit(1)
s = sigma(1)

oplot,tm,fit(0)+fit(1)*tm

;mass 2
close,1
;rundir = '/Volumes/MacD72-2/hybrid/KHI/run_28'
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
;!P = p2 & !X = x2 & !Y = y2
plot,tm,lnvz,linestyle=2,title='3';,color=c3


wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)

oplot,tm,fit(0)+fit(1)*tm,color=c3

f = [f,fit(1)]
s = [s,sigma(1)]


;mass 4
close,1
;rundir = '/Volumes/MacD72-2/hybrid/KHI/run_30'
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
close,1
;rundir = '/Volumes/MacD72-2/hybrid/KHI/run_31'
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_'+files(4)  
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
;!P = p2 & !X = x2 & !Y = y2
plot,tm,lnvz,linestyle=2,title='5';,color=c5


wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)

oplot,tm,fit(0)+fit(1)*tm,color=c5

f = [f,fit(1)]
s = [s,sigma(1)]

;mass 6
close,1
;rundir = '/Volumes/MacD72-2/hybrid/KHI/run_32'
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
;   lnvz(i-1) = alog(max(abs(up(*,1,*,2))))

endfor

;!p.multi=[0,1,2]
;!P = p1 & !X = x1 & !Y = y1
;oplot,tm,mix_arr,linestyle=1
;!P = p2 & !X = x2 & !Y = y2
plot,tm,lnvz,linestyle=2,title='6';,color=c6


wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)

oplot,tm,fit(0)+fit(1)*tm,color=c6

f = [f,fit(1)]
s = [s,sigma(1)]


;mass 7
close,1
;rundir = '/Volumes/MacD72-2/hybrid/KHI/run_33'
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_'+files(6)  
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
;   lnvz(i-1) = alog(max(abs(up(*,1,*,2))))

endfor

;!p.multi=[0,1,2]
;!P = p1 & !X = x1 & !Y = y1
;oplot,tm,mix_arr,linestyle=1
;!P = p2 & !X = x2 & !Y = y2
plot,tm,lnvz,linestyle=2,title='7';,color=c7


wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)

oplot,tm,fit(0)+fit(1)*tm,color=c7

f = [f,fit(1)]
s = [s,sigma(1)]


;mass 8
close,1
;rundir = '/Volumes/MacD72-2/hybrid/KHI/run_34'
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
;   lnvz(i-1) = alog(max(abs(up(*,1,*,2))))

endfor

;!p.multi=[0,1,2]
;!P = p1 & !X = x1 & !Y = y1
;oplot,tm,mix_arr,linestyle=1
;!P = p2 & !X = x2 & !Y = y2
plot,tm,lnvz,linestyle=2,title='8';,color=c7


wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)

oplot,tm,fit(0)+fit(1)*tm,color=c7

f = [f,fit(1)]
s = [s,sigma(1)]


;mass 9
close,1
;rundir = '/Volumes/MacD72-2/hybrid/KHI/run_35'
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_'+files(8)  
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
;   lnvz(i-1) = alog(max(abs(up(*,1,*,2))))

endfor

;!p.multi=[0,1,2]
;!P = p1 & !X = x1 & !Y = y1
;oplot,tm,mix_arr,linestyle=1
;!P = p2 & !X = x2 & !Y = y2
plot,tm,lnvz,linestyle=2,title='9';,color=c7


wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)

oplot,tm,fit(0)+fit(1)*tm,color=c7

f = [f,fit(1)]
s = [s,sigma(1)]


;mass 10
close,1
;rundir = '/Volumes/MacD72-2/hybrid/KHI/run_36'
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
;   lnvz(i-1) = alog(max(abs(up(*,1,*,2))))

endfor

;!p.multi=[0,1,2]
;!P = p1 & !X = x1 & !Y = y1
;oplot,tm,mix_arr,linestyle=1
;!P = p2 & !X = x2 & !Y = y2
plot,tm,lnvz,linestyle=2,title='10';,color=c7


wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)

oplot,tm,fit(0)+fit(1)*tm,color=c7

f = [f,fit(1)]
s = [s,sigma(1)]


;mass 11
close,1
;rundir = '/Volumes/MacD72-2/hybrid/KHI/run_37'
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_'+files(10)  
f_read_coord,rundir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz

nfrm = 11

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
;   lnvz(i-1) = alog(max(abs(up(*,1,*,2))))

endfor

;!p.multi=[0,1,2]
;!P = p1 & !X = x1 & !Y = y1
;oplot,tm,mix_arr,linestyle=1
;!P = p2 & !X = x2 & !Y = y2
plot,tm,lnvz,linestyle=2,title='10';,color=c7


wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)

oplot,tm,fit(0)+fit(1)*tm,color=c7

f = [f,fit(1)]
s = [s,sigma(1)]


;mass 12
close,1
;rundir = '/Volumes/MacD72-2/hybrid/KHI/run_38'
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
;   lnvz(i-1) = alog(max(abs(up(*,1,*,2))))

endfor

;!p.multi=[0,1,2]
;!P = p1 & !X = x1 & !Y = y1
;oplot,tm,mix_arr,linestyle=1
;!P = p2 & !X = x2 & !Y = y2
plot,tm,lnvz,linestyle=2,title='10';,color=c7


wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)

oplot,tm,fit(0)+fit(1)*tm,color=c7

f = [f,fit(1)]
s = [s,sigma(1)]

;mass 13
close,1
;rundir = '/Volumes/MacD72-2/hybrid/KHI/run_39'
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_'+files(12)  
f_read_coord,rundir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz

nfrm = 11

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
;   lnvz(i-1) = alog(max(abs(up(*,1,*,2))))

endfor

;!p.multi=[0,1,2]
;!P = p1 & !X = x1 & !Y = y1
;oplot,tm,mix_arr,linestyle=1
;!P = p2 & !X = x2 & !Y = y2
plot,tm,lnvz,linestyle=2,title='10';,color=c7


wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)

oplot,tm,fit(0)+fit(1)*tm,color=c7

f = [f,fit(1)]
s = [s,sigma(1)]



;mass 14
;t1 = 100
;t2 = 250
close,1
;rundir = '/Volumes/MacD72-2/hybrid/KHI/run_40'
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_'+files(13)  
f_read_coord,rundir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz

nfrm = 11

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
;   lnvz(i-1) = alog(max(abs(up(*,1,*,2))))

endfor

;!p.multi=[0,1,2]
;!P = p1 & !X = x1 & !Y = y1
;oplot,tm,mix_arr,linestyle=1
;!P = p2 & !X = x2 & !Y = y2
plot,tm,lnvz,linestyle=2,title='10';,color=c7


wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)

oplot,tm,fit(0)+fit(1)*tm,color=c7

f = [f,fit(1)]
s = [s,sigma(1)]

;mass 15
;t1 = 100
;t2 = 250
close,1
;rundir = '/Volumes/MacD72-2/hybrid/KHI/run_40'
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_'+files(14)  
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
;   lnvz(i-1) = alog(max(abs(up(*,1,*,2))))

endfor

;!p.multi=[0,1,2]
;!P = p1 & !X = x1 & !Y = y1
;oplot,tm,mix_arr,linestyle=1
;!P = p2 & !X = x2 & !Y = y2
plot,tm,lnvz,linestyle=2,title='10';,color=c7


wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)

oplot,tm,fit(0)+fit(1)*tm,color=c7

f = [f,fit(1)]
s = [s,sigma(1)]


;mass 16
;t1 = 100
;t2 = 250
close,1
;rundir = '/Volumes/MacD72-2/hybrid/KHI/run_40'
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_'+files(15)  
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
;   lnvz(i-1) = alog(max(abs(up(*,1,*,2))))

endfor

;!p.multi=[0,1,2]
;!P = p1 & !X = x1 & !Y = y1
;oplot,tm,mix_arr,linestyle=1
;!P = p2 & !X = x2 & !Y = y2
plot,tm,lnvz,linestyle=2,title='10';,color=c7


wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)

oplot,tm,fit(0)+fit(1)*tm,color=c7

f = [f,fit(1)]
s = [s,sigma(1)]



;mass 17
;t1 = 100
;t2 = 250
close,1
;rundir = '/Volumes/MacD72-2/hybrid/KHI/run_40'
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_'+files(16)  
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
;   lnvz(i-1) = alog(max(abs(up(*,1,*,2))))

endfor

;!p.multi=[0,1,2]
;!P = p1 & !X = x1 & !Y = y1
;oplot,tm,mix_arr,linestyle=1
;!P = p2 & !X = x2 & !Y = y2
plot,tm,lnvz,linestyle=2,title='10';,color=c7


wh = where((tm ge t1) and (tm le t2))
fit = poly_fit(tm(wh),lnvz(wh),1,sigma=sigma)

oplot,tm,fit(0)+fit(1)*tm,color=c7

f = [f,fit(1)]
s = [s,sigma(1)]



;legend,['1','2','4','6','7','8','12'],linestyle=[0,1,2,3,4,5,6],/bottom,/right

;stop

ff = (2*4*300./500)

plot,mass,f*ff,psym=1,xrange=[0,18],/xsty,yrange=ff*[0.004,0.015],/ysty
errplot,mass,(f+s)*ff,(f-s)*ff

fit = poly_fit(mass,f,1,sigma=sigma)

oplot,mass,fit(0)*ff+fit(1)*ff*mass
oplot,mass,fit(0)*ff-sigma(0)*ff+(fit(1)*ff+sigma(1)*ff)*mass,linestyle=1
oplot,mass,fit(0)*ff+sigma(0)*ff+(fit(1)*ff-sigma(1)*ff)*mass,linestyle=1
print,fit*ff,sigma*ff
@gops
!p.multi=[0,1,1]
device,filename='case2_growth_rates.ps'
;ff = (2*4*300./500)

f = f*1e3
s = s*1e3

device,/inches,xsize=6.0,ysize=6.0,xoffset=1

plot,mass,f,psym=1,xrange=[0.5,17.5],/xsty,yrange=[0.004,0.012]*1e3,/ysty,$
  ytitle='q (10!u-3!n s!u-1!n)',xtitle='mass ratio'
errplot,mass,(f+s),(f-s)

fit = poly_fit(mass,f,1,sigma=sigma)
fit = fit
sigma = sigma

oplot,mass,fit(0)+fit(1)*mass
oplot,mass,fit(0)-sigma(0)+(fit(1)+sigma(1))*mass,linestyle=1
oplot,mass,fit(0)+sigma(0)+(fit(1)-sigma(1))*mass,linestyle=1
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
