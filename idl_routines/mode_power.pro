pro img_cont, a, x, y,dt, WINDOW_SCALE = window_scale, $
                ASPECT = aspect, INTERP = interp,POSTSCRIPT=postscript
;+
; NAME:
;	IMAGE_CONT
; PURPOSE:
;	Overlay an image and a contour plot.
; CATEGORY:
;	General graphics.
; CALLING SEQUENCE:
;	IMAGE_CONT, A
; INPUTS:
;	A = 2 dimensional array to display.
; KEYWORD PARAMETERS:
;	/WINDOW_SCALE = set to scale the window size to the image size,
;		otherwise the image size is scaled to the window size.
;		Ignored when outputting to devices with scalable pixels.
;	/ASPECT = set to retain image's aspect ratio.  Assumes square
;		pixels.  If /WINDOW_SCALE is set, the aspect ratio is
;		retained.
;	/INTERP = set to bi-linear interpolate if image is resampled.
; OUTPUTS:
;	No explicit outputs.
; COMMON BLOCKS:
;	none.
; SIDE EFFECTS:
;	The currently selected display is affected.
; RESTRICTIONS:
;	None that are obvious.
; PROCEDURE:
;	If the device has scalable pixels then the image is written over
;	the plot window.
; MODIFICATION HISTORY:
;	DMS, May, 1988.
;-

nclr = !d.n_colors

if keyword_set(postscript) then begin
set_plot,'ps
device,filename='img.eps
device,/encapsulated
!p.font=0
device,/palatino
device,bits=8
device,/color
endif

;f_read_coord,'coord.dat',x,y,z,dzg,dzc,nx,ny,nz
;ax = y
;az = z

;a = bytscl(a)
;a = a<255
maxa=max(a)
mina=min(a)
;help,a,maxa,mina
;print,!d.n_colors,nclr
;a = (not(bytscl(a)))*nclr/256
;print,max(a)
;stop
sz=size(a)
;ab = bytscl(a)
;a=rebin(a,4*sz(1),4*sz(2))
;sz=size(a)
 

;!x.title = 'y (10!u3!n km)' 
;!y.title = 'z (10!u3!n km)'
!x.title = 'time (s)' 
;!y.title = 'mode'

!p.charsize=1.4
!y.margin = [4,4]
!x.margin = [10,10]
;clrb = findgen(sz(2))*nclr/max(findgen(sz(2)))
;;clrb = findgen(sz(2))*255/max(findgen(sz(2)))
;for i=sz(1)-15,sz(1)-1 do begin
;   a(i,*) = max(clrb)-clrb
;;   a(i,*) = clrb
;endfor

on_error,2                      ;Return to caller if an error occurs
sz = size(a)			;Size of image
if sz(0) lt 2 then message, 'Parameter not 2D'

	;set window used by contour

contour,[[0,0],[1,1]],/nodata, xstyle=1, ystyle = 1,ytitle='mode',$
            xrange=[min(x),max(x)+dt],yrange=[min(y),max(y)+1]
p1 = !P & x1 = !X & y1= !Y

px = !x.window * !d.x_vsize	;Get size of window in device units
py = !y.window * !d.y_vsize
swx = px(1)-px(0)		;Size in x in device units
swy = py(1)-py(0)		;Size in Y
six = float(sz(1))		;Image sizes
siy = float(sz(2))
aspi = six / siy		;Image aspect ratio
aspw = swx / swy		;Window aspect ratio
f = aspi / aspw			;Ratio of aspect ratios

if (!d.flags and 1) ne 0 then begin	;Scalable pixels?
  if keyword_set(aspect) then begin	;Retain aspect ratio?
				;Adjust window size
	if f ge 1.0 then swy = swy / f else swx = swx * f
	endif
  print,max(a)
  tv,a,px(0),py(0),xsize = swx, ysize = swy, /device


endif else begin	;Not scalable pixels	
   if keyword_set(window_scale) then begin ;Scale window to image?
	tv,a,px(0),py(0)	;Output image
	swx = six		;Set window size from image
	swy = siy
    endif else begin		;Scale window
	if keyword_set(aspect) then begin
		if f ge 1.0 then swy = swy / f else swx = swx * f
		endif		;aspect
        print,max(a)
	tv,poly_2d((a),$	;Have to resample image
                   [[0,0],[six/swx,0]], [[0,siy/swy],[0,0]],$
                   keyword_set(interp),swx,swy), $
           px(0),py(0)

	endelse			;window_scale
  endelse			;scalable pixels

     axis,xaxis=2,xstyle=1
     axis,yaxis=2,ystyle=1
;     axis,yaxis=1,ystyle=1,$
;;     yticks=2,$
;;     ytickv=[min(y),(min(y)+max(y))/2,max(y)], $
;;     ytickname=[string(min(a)),string((min(a)+max(a))/2),string(max(a))], $
;;     ytitle='E!dz!n (mV/m)',$
;;     yrange = [mina*(2.3e-25/1.6e-19)*1e6,maxa*(2.3e-25/1.6e-19*1e6)]
;     ytitle='(u!de!n)!dz!n (km/s)',$
;     yrange = [mina,maxa]


mx = !d.n_colors-1		;Brightest color
colors = [mx,mx,mx,0,0,0]	;color vectors
if !d.name eq 'PS' then colors = mx - colors ;invert line colors for pstscrp
m = max(a)
;contour,[[0,0],[1,1]],/nodata, xstyle=1, ystyle = 1
;contour,a,/noerase,/nodata,/yst,$	;Do the contour
;	  pos = [px(0),py(0), px(0)+swx,py(0)+swy],/dev,$
;          c_color =  colors, charsize=1.5, $
;          xtitle='x (km)', ytitle='y (km)', $
;          levels=[0.02*m,0.04*m,0.06*m,0.2*m,0.4*m,0.5*m,0.7*m,0.9*m]
;;	  xrange=[min(ax),max(ax)], xstyle=1
;dx = !x.crange(1) - !x.crange(0)
;dy = !y.crange(1) - !y.crange(0)
;xyouts,0.05*dx+!x.crange(0),0.85*dy+!y.crange(0),tit,/data
;!p.multi=[0,1,2]

!P = p1 & !X = x1 & !Y = y1

;contour,[[0,0],[1,1]],/nodata, xstyle=1, ystyle = 1,/noerase

if keyword_set(postscript) then begin
device,/close
set_plot,'x'
!p.font=-1
endif

return
end
;---------------------------------------------------------------------------;
;mass 1
nav = 1
mmax = 40
mmin = 1
nplot=6
!p.multi=[0,3,2]

close,1
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_18'
f_read_coord,rundir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz

nfrm = 20

loadct,0
s_max = 105
s_min = 0
nnx = nx
nny = nz

zm = 1

img = fltarr(nnx,nny,nfrm)
 
for i = 1,nfrm do begin 
;   f_read_3d_vec_m_32,rundir+'/upall_1',i,np
   f_read_3d_vec_m_32,rundir+'/upall_1',i,np
   img(*,*,i-1)=bytscl(reform(np(*,1,*,2)))
endfor


dx = 300.0

Ni = nx+1
Ti = dx
k21 = Ni/2 + 1
k = indgen(Ni)
k(k21) = k21 - Ni + findgen(k21 - 2)
k = 2*!pi*k/(Ni*Ti)

m = (shift(k(1:*),-k21))*(Ni*Ti)/(2*!pi)
print,m

wp = fltarr(nfrm,nx-1)
im = fltarr(nx)   

for j = 0,nfrm-1 do begin

   for i = 0,nx-1 do begin
;      im(i) = total(reform(img(i,nz/2-0:nz/2+5,j)))
      im(i) = total(reform((img(i,nz/2-0:nz/2+nav,j))))/(nav +1)
   endfor

   ftarr = fft(im)
   ftarr = (abs(ftarr)^2)

   wp(j,*) = shift(abs(ftarr(1:*)),-k21)

endfor



wh = where(m ge mmin)
dt = 0.5
nout = 50
tm = indgen(nfrm)*nout*dt


device,decompose=0
loadct,39

;wp = wp/max(wp)
;print,max(wp^2)
;stop


lwp = alog(wp(*,wh(0:mmax)))

sz = size(lwp)

wpt = fltarr(sz(1),nplot)

for i = 0,sz(1)-1 do begin
   wpt(i,0) = total(wp(i,*))
endfor

lwp_arr = fltarr(sz(1),sz(2),nplot)

lwp_arr(*,*,0) = lwp

;img_cont,alog(wp(*,wh(0:mmax))),tm+dt*nout/2,m(wh(0:mmax))-0.5,dt*nout

;-------------------------------------
;mass 2
;nav = 5
;mmax = 40
;!p.multi=[0,2,2]

close,1
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_17'
f_read_coord,rundir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz

nfrm = 20

loadct,0
s_max = 105
s_min = 0
nnx = nx
nny = nz

zm = 1

img = fltarr(nnx,nny,nfrm)
 
for i = 1,nfrm do begin 
;   f_read_3d_vec_m_32,rundir+'/upall_1',i,np
   f_read_3d_vec_m_32,rundir+'/upall_1',i,np
   img(*,*,i-1)=bytscl(reform(np(*,1,*,2)))
endfor

dx = 300.0

Ni = nx+1
Ti = dx
k21 = Ni/2 + 1
k = indgen(Ni)
k(k21) = k21 - Ni + findgen(k21 - 2)
k = 2*!pi*k/(Ni*Ti)

m = (shift(k(1:*),-k21))*(Ni*Ti)/(2*!pi)
print,m

wp = fltarr(nfrm,nx-1)
im = fltarr(nx)   

for j = 0,nfrm-1 do begin

   for i = 0,nx-1 do begin
;      im(i) = total(reform(img(i,nz/2-0:nz/2+5,j)))
      im(i) = total(reform((img(i,nz/2-0:nz/2+nav,j))))/(nav +1)
   endfor

   ftarr = fft(im)
   ftarr = (abs(ftarr)^2)

   wp(j,*) = shift(abs(ftarr(1:*)),-k21)

endfor



wh = where(m ge mmin)
dt = 0.5
nout = 50
tm = indgen(nfrm)*nout*dt


device,decompose=0
loadct,39

;wp = wp/max(wp)
;print,max(wp^2)
;stop

lwp = alog(wp(*,wh(0:mmax)))

for i = 0,sz(1)-1 do begin
   wpt(i,1) = total(wp(i,*))
endfor


lwp_arr(*,*,1) = lwp

;img_cont,alog(wp(*,wh(0:mmax))),tm+dt*nout/2,m(wh(0:mmax))-0.5,dt*nout


;------------------------------------
; mass 4

close,1
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_16'
f_read_coord,rundir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz

nfrm = 20

loadct,0
s_max = 105
s_min = 0
nnx = nx
nny = nz

zm = 1

img = fltarr(nnx,nny,nfrm)
 
for i = 1,nfrm do begin 
   f_read_3d_vec_m_32,rundir+'/upall_1',i,np
   img(*,*,i-1)=bytscl(reform(np(*,1,*,2)))
endfor

dx = 300.0

Ni = nx+1
Ti = dx
k21 = Ni/2 + 1
k = indgen(Ni)
k(k21) = k21 - Ni + findgen(k21 - 2)
k = 2*!pi*k/(Ni*Ti)

m = (shift(k(1:*),-k21))*(Ni*Ti)/(2*!pi)
print,m

wp = fltarr(nfrm,nx-1)
im = fltarr(nx)   

for j = 0,nfrm-1 do begin

   for i = 0,nx-1 do begin
;      im(i) = total(reform(img(i,nz/2-0:nz/2+5,j)))
      im(i) = total(reform((img(i,nz/2-0:nz/2+nav,j))))/(nav +1)
   endfor

   ftarr = fft(im)
   ftarr = (abs(ftarr)^2)

   wp(j,*) = shift(abs(ftarr(1:*)),-k21)

endfor

wh = where(m ge mmin)
dt = 0.5
nout = 50
tm = indgen(nfrm)*nout*dt

device,decompose=0
loadct,39

lwp = alog(wp(*,wh(0:mmax)))

for i = 0,sz(1)-1 do begin
   wpt(i,2) = total(wp(i,*))
endfor


lwp_arr(*,*,2) = lwp

;img_cont,alog(wp(*,wh(0:mmax))),tm+dt*nout/2,m(wh(0:mmax))-0.5,dt*nout


;------------------------------------
; mass 6

close,1
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_19'
f_read_coord,rundir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz

nfrm = 20

loadct,0
s_max = 105
s_min = 0
nnx = nx
nny = nz

zm = 1

img = fltarr(nnx,nny,nfrm)
 
for i = 1,nfrm do begin 
   f_read_3d_vec_m_32,rundir+'/upall_1',i,np
   img(*,*,i-1)=bytscl(reform(np(*,1,*,2)))
endfor

dx = 300.0

Ni = nx+1
Ti = dx
k21 = Ni/2 + 1
k = indgen(Ni)
k(k21) = k21 - Ni + findgen(k21 - 2)
k = 2*!pi*k/(Ni*Ti)

m = (shift(k(1:*),-k21))*(Ni*Ti)/(2*!pi)
print,m

wp = fltarr(nfrm,nx-1)
im = fltarr(nx)   

for j = 0,nfrm-1 do begin

   for i = 0,nx-1 do begin
;      im(i) = total(reform(img(i,nz/2-0:nz/2+5,j)))
      im(i) = total(reform((img(i,nz/2-0:nz/2+nav,j))))/(nav +1)
   endfor

   ftarr = fft(im)
   ftarr = (abs(ftarr)^2)

   wp(j,*) = shift(abs(ftarr(1:*)),-k21)

endfor

wh = where(m ge mmin)
dt = 0.5
nout = 50
tm = indgen(nfrm)*nout*dt

lwp = alog(wp(*,wh(0:mmax)))

for i = 0,sz(1)-1 do begin
   wpt(i,3) = total(wp(i,*))
endfor


lwp_arr(*,*,3) = lwp

;img_cont,alog(wp(*,wh(0:mmax))),tm+dt*nout/2,m(wh(0:mmax))-0.5,dt*nout

;----------------------------
; mass 8

close,1
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_15'
f_read_coord,rundir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz

nfrm = 20

loadct,0
s_max = 105
s_min = 0
nnx = nx
nny = nz

zm = 1

img = fltarr(nnx,nny,nfrm)
 
for i = 1,nfrm do begin 
   f_read_3d_vec_m_32,rundir+'/upall_1',i,np
   img(*,*,i-1)=bytscl(reform(np(*,1,*,2)))
endfor

dx = 300.0

Ni = nx+1
Ti = dx
k21 = Ni/2 + 1
k = indgen(Ni)
k(k21) = k21 - Ni + findgen(k21 - 2)
k = 2*!pi*k/(Ni*Ti)

m = (shift(k(1:*),-k21))*(Ni*Ti)/(2*!pi)
print,m

wp = fltarr(nfrm,nx-1)
im = fltarr(nx)   

for j = 0,nfrm-1 do begin

   for i = 0,nx-1 do begin
      im(i) = total(reform((img(i,nz/2-0:nz/2+nav,j))))/(nav +1)
   endfor

   ftarr = fft(im)
   ftarr = (abs(ftarr)^2)
   
   wp(j,*) = shift(abs(ftarr(1:*)),-k21)
   
endfor

wh = where(m ge mmin)
dt = 0.5
nout = 50
tm = indgen(nfrm)*nout*dt

device,decompose=0
loadct,39

lwp = alog(wp(*,wh(0:mmax)))

for i = 0,sz(1)-1 do begin
   wpt(i,4) = total(wp(i,*))
endfor


lwp_arr(*,*,4) = lwp


;---------------------------------
; mass 12
close,1
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_14'
f_read_coord,rundir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz

nfrm = 20

loadct,0
s_max = 105
s_min = 0
nnx = nx
nny = nz

zm = 1

img = fltarr(nnx,nny,nfrm)
 
for i = 1,nfrm do begin 
   f_read_3d_vec_m_32,rundir+'/upall_1',i,np
   img(*,*,i-1)=bytscl(reform(np(*,1,*,2)))
endfor


dx = 300.0

Ni = nx+1
Ti = dx
k21 = Ni/2 + 1
k = indgen(Ni)
k(k21) = k21 - Ni + findgen(k21 - 2)
k = 2*!pi*k/(Ni*Ti)

m = (shift(k(1:*),-k21))*(Ni*Ti)/(2*!pi)
print,m

wp = fltarr(nfrm,nx-1)
im = fltarr(nx)   

for j = 0,nfrm-1 do begin

   for i = 0,nx-1 do begin
      im(i) = total(reform((img(i,nz/2-0:nz/2+nav,j))))/(nav +1)
   endfor
   
   ftarr = fft(im)
   ftarr = (abs(ftarr)^2)

   wp(j,*) = shift(abs(ftarr(1:*)),-k21)

endfor

wh = where(m ge mmin)
dt = 0.5
nout = 50
tm = indgen(nfrm)*nout*dt

lwp = alog(wp(*,wh(0:mmax)))

for i = 0,sz(1)-1 do begin
   wpt(i,5) = total(wp(i,*))
endfor


lwp_arr(*,*,5) = lwp

print,max(lwp_arr(*,*,0)),min(lwp_arr(*,*,0))

lwp_arr = bytscl(lwp_arr)
loadct,39
for i = 0,nplot-1 do begin
   img_cont,lwp_arr(*,*,i)<254,tm+dt*nout/2,m(wh(0:mmax))-0.5,dt*nout
   colorbar,position=[!y.window(0), !x.window(1)+0.01, !y.window(1), $
                      !x.window(1)+0.03],ncolors=255,/vertical,/right,$
            maxrange=[1],divisions=4,title='log(wave power)',$
            charsize=1.4
endfor
stop
!p.multi=[0,1,1]

wpt=(wpt/max(wpt))
device,decompose=0
plot,tm+dt*nout/2,wpt(*,0),linestyle=0,/ylog,yrange=[0.03,1.0],/ysty,$
   ytitle='normalized log(wave power)',/nodata
oplot,tm+dt*nout/2,wpt(*,0),linestyle=0,color=fsc_color('red')
oplot,tm+dt*nout/2,wpt(*,1),linestyle=1,color=fsc_color('blue')
oplot,tm+dt*nout/2,wpt(*,2),linestyle=2,color=fsc_color('green')
oplot,tm+dt*nout/2,wpt(*,3),linestyle=3,color=fsc_color('yellow')
oplot,tm+dt*nout/2,wpt(*,4),linestyle=4,color=fsc_color('cyan')
oplot,tm+dt*nout/2,wpt(*,5),linestyle=5,color=fsc_color('purple')

legend,['1','2','4','6','8','12'],linestyle=[0,1,2,3,4,5],/bottom,/right

help,lwp_arr



;plot,tm+dt*nout/2,wp(0,*)


end
