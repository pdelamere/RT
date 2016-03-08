pro img_cont, a, x, y,nfrm, WINDOW_SCALE = window_scale, $
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
;!x.title = 'time (s)' 
;!y.title = 'mode'

!p.charsize=1.2
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

contour,[[0,0],[1,1]],/nodata, xstyle=1, ystyle = 1,$
            ytitle='z (10!u3!n km)',xtitle='x (10!u3!n km)',$
            xrange=[min(x),max(x)],yrange=[min(y),max(y)+1],$
;            title='t = '+strmid(strtrim(string(0.5*nfrm*100.),2),0,6) + ' (s)',$
            /isotropic
;            xrange=[min(x),800],yrange=[min(y),max(y)+1]
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

pro img_cont_ylog, a, x, y,nfrm, WINDOW_SCALE = window_scale, $
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

maxa=max(a)
mina=min(a)

sz=size(a)

!p.charsize=1.2
!y.margin = [4,4]
!x.margin = [10,10]

on_error,2                      ;Return to caller if an error occurs
sz = size(a)			;Size of image
if sz(0) lt 2 then message, 'Parameter not 2D'

	;set window used by contour

contour,[[0,0],[1,1]],/nodata, xstyle=1, ystyle = 1,$
            ytitle='Energy (eV)',xtitle='z (10!u3!n km)',$
            xrange=[min(x),max(x)],yrange=[min(y),max(y)],/ylog
;            xrange=[min(x),max(x)],yrange=[10,5000],/ylog
;            title='t = '+strmid(strtrim(string(0.5*nfrm*100.),2),0,6) + ' (s)',$
;            /isotropic
;            xrange=[min(x),800],yrange=[min(y),max(y)+1]
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
  ;print,max(a)
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
        ;print,max(a)
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



;----------------------------------------------------------------
PRO read_part,file,nfrm,xp
;----------------------------------------------------------------

Ni_max=long(0)
nt=0l
ntout=0l
frm=0l

file = file+'.dat'
print,' reading...',file
openr,1,file,/f77_unformatted
readu,1,nt
readu,1,ntout
readu,1,Ni_max
print,nt,ntout,Ni_max

xp=fltarr(Ni_max,3,/nozero)

readu,1,frm
print,'  image #.....',frm
readu,1,xp
frmcnt = 1

while (frmcnt lt nfrm) do begin

   readu,1,frm
   print,'  image #.....',frm
   readu,1,xp
   frmcnt = frmcnt + 1

endwhile

close,1

return
end

;----------------------------------------------------------------


;----------------------------------------------------------------
PRO read_part_scalar,file,nfrm,xp
;----------------------------------------------------------------

Ni_max=long(0)
nt=0
ntout=0
frm=0

file = file+'.dat'
print,' reading...',file
openr,1,file,/f77_unformatted
readu,1,nt
readu,1,ntout
readu,1,Ni_max
print,nt,ntout,Ni_max

xp=fltarr(Ni_max,/nozero)

readu,1,frm
print,'  image #.....',frm
readu,1,xp
frmcnt = 1

while (frmcnt lt nfrm) do begin

   readu,1,frm
   print,'  image #.....',frm
   readu,1,xp
   frmcnt = frmcnt + 1

endwhile

close,1

return
end
;----------------------------------------------------------------


;----------------------------------------------------------------
pro get_dist,xyz_traj,lxyz,nlxyz,t
;----------------------------------------------------------------
close,1

mp=1.67e-27
nfile = '1'
nfrm = 15
;define look direction of instrument
;theta = 180.*!dtor
;phi = !pi/4
;dtheta= 20.*!dtor
;dphi = 20.*!dtor

rundir = '/Volumes/MacD97-2/hybrid/KHI/run_101/'

part_dir = rundir+'part_1/'
fluid_fields_dir = rundir

print,part_dir

f_read_coord,fluid_fields_dir+'coord.dat',x,y,z,dzc,dzg,nx,ny,nz
f_read_3d_m,fluid_fields_dir+'npall_'+strtrim(string(nfile),2),nfrm,npt
f_read_3d_m,fluid_fields_dir+'np_b_'+strtrim(string(nfile),2),nfrm,npb
f_read_3d_vec_m_32,fluid_fields_dir+'b1all_'+strtrim(string(nfile),2),nfrm,b1
f_read_3d_vec_m_32,fluid_fields_dir+'upall_'+strtrim(string(nfile),2),nfrm,up
f_read_3d_m,fluid_fields_dir+'temp_p_'+strtrim(string(nfile),2),nfrm,ti

up2d = fltarr(nx,nz,3)
up2d(*,*,*) = reform(up(*,1,*,*))

np = npt+npb
;np = ti

xx = x
z = z - z(nz/2)

z1=nz/2
z2 = nz/2

x1 = 0
x2 = nx-1

zz = z(nz/2-z1:nz/2+z2)

xx  = x(x1:x2)

np_img = reform(np(x1:x2,1,nz/2-z1:nz/2+z2))
b1_img = reform(sqrt(b1(*,1,*,0)^2 + b1(*,1,*,1)^2 + b1(*,1,*,2)^2))

window,0
!p.multi=[0,1,1]
device,decompose=0
loadct,39
xx = indgen(nx)
zz = indgen(nz)
;xx = x
;zz = z
window,0
x1 = 0
x2 = nx-1
z1 = nz/2-100
z2 = nz/2+100
img_cont,bytscl(b1_img(x1:x2,z1:z2)),xx(x1:x2),zz(z1:z2),nfrm

consx = 8
consz = consx
sz = size(up2d(x1:x2,z1:z2))

xxx = congrid(xx,sz(1)/consz)
zzz = congrid(zz(z1:z2),sz(2)/consz)

vx = congrid(reform(up2d(x1:x2,z1:z2,0)),sz(1)/consz,sz(2)/consz)
vz = congrid(reform(up2d(x1:x2,z1:z2,2)),sz(1)/consz,sz(2)/consz)
   
   
velovect,vx(1:*,1:*),vz(1:*,1:*),xxx(1:*),zzz(1:*),/overplot,color=fsc_color("white"),thick=1.5,length=2.0
   
cursor,xcur,zcur,/down,/data

window,1
window,2
;window,3
window,4


;write to postscrip
;@gops
;device,filename='caps_synth_den_vec.ps'
;device,/color
;img_cont,bytscl(np_img(x1:x2,z1:z2)),xx(x1:x2),zz(z1:z2),nfrm
;velovect,vx(1:*,1:*),vz(1:*,1:*),xxx(1:*),zzz(1:*),/overplot,color=fsc_color("w;hite"),thick=1.5,length=2.0
;device,/close
;set_plot,'x'


xcur = fix(xcur)
zcur = fix(zcur)
ycur=1


print,xcur,ycur,zcur

ztemp = fltarr(35)

dv = 10.0
nn = 2000.
alpha=1.9263418e-20
d3v = dv^3
vxyp = fltarr(nn/dv,nn/dv)
vxzp = fltarr(nn/dv,nn/dv)
vyzp = fltarr(nn/dv,nn/dv)
;fvp = fltarr(nn/dv,nn/dv,nn/dv)
;vxyi = fltarr(nn/dv,nn/dv)
;vxzi = fltarr(nn/dv,nn/dv)
;vyzi = fltarr(nn/dv,nn/dv)
;fvi = fltarr(nn/dv,nn/dv,nn/dv)
vx = -(nn/(2*dv))*dv + findgen(nn/dv)*dv
vy = -(nn/(2*dv))*dv + findgen(nn/dv)*dv
vz = -(nn/(2*dv))*dv + findgen(nn/dv)*dv

;xcur = fix(nx/2
;ycur = ny/2 
;zcur = nz/2-30

;xcur = xyz_traj(0)
;ycur = xyz_traj(1)
;zcur = xyz_traj(2)

dx = 250.0
ndx = 2.0
;Rp = 1200.
;svol = (2.0*ndx*dx*1e5)^3
b = b1(xcur,ycur,zcur,*)
bx = b1(xcur,ycur,zcur,0)
by = b1(xcur,ycur,zcur,1)
bz = b1(xcur,ycur,zcur,2)

;bmag = sqrt(bx^2 + by^2 + bz^2)
;bxhat = bx/bmag
;byhat = by/bmag
;bzhat = bz/bmag

bxy = b(1)/b(0)
bzx = b(2)/b(0)
bzy = b(2)/b(1)

;;look direction xyz
;lxhat = lxyz(0)
;lyhat = lxyz(1)
;lzhat = lxyz(2)

;lmag = sqrt(lxhat^2 + lyhat^2 + lzhat^2)
;lxhat = lxhat/lmag
;lyhat = lyhat/lmag
;lzhat = lzhat/lmag

;nlxhat = nlxyz(0)  ;axis perpedicular to the look direction
;nlyhat = nlxyz(1)
;nlzhat = nlxyz(2)

;look direction perp b
;nlxhat = bxhat
;nlyhat = byhat
;nlzhat = bzhat

;lxhat = 1.0
;lyhat = -lxhat*bxhat/byhat
;lzhat = 0.0

;lmag = sqrt(lxhat^2 + lyhat^2 + lzhat^2)
;lxhat = lxhat/lmag
;lyhat = lyhat/lmag
;lzhat = lzhat/lmag

!p.multi=[0,1,1]
f_read_coord,fluid_fields_dir+'coord.dat',x,y,z,dzc,dzg,nx,ny,nz
for nfil = 0,5 do begin
    xfile = part_dir+'xp_'+strtrim(string(nfil),2)
    vfile = part_dir+'vp_'+strtrim(string(nfil),2)
    mratfile = part_dir+'mrat_'+strtrim(string(nfil),2)

    read_part,xfile,nfrm,xp
    read_part,vfile,nfrm,vp
    read_part_scalar,mratfile,nfrm,mrat
 
endfor

save,filename='part_arr.sav',xp,vp,mrat
;restore,'part_arr.sav'

i = xcur
j = ycur
k = zcur

dE = 20.
hmin = 1.0
hmax = 10000.

nh = fix((hmax-hmin)/dE)+1

zcnt=100
nnz = (nz/2+zcnt) - (nz/2 - zcnt)
evst = fltarr(nnz+1,nh)

wset,0
plots,[xx(i),xx(i)],[zz(nz/2-zcnt),zz(nz/2+zcnt)],/data,linestyle=0,thick=4


cnt = 0
for zind = nz/2-zcnt,nz/2+zcnt do begin


   i = xcur
   j = 1
   k = zind
   print,zind,np_img(i,k)



   wset,0
   plots,xx(i),zz(k),/data,psym=6

   wh = where((xp(*,0) ge x(i)-ndx*dx) and (xp(*,0) le x(i)+ndx*dx) and $
              (xp(*,1) ge y(j)-ndx*dx) and (xp(*,1) le y(j)+ndx*dx) and $ 
              (xp(*,2) ge z(k)-ndx*dx) and (xp(*,2) le z(k)+ndx*dx))
   print,n_elements(wh),z(k)

   if (wh(0) gt -1) then begin
      e_arr = 0
      cnt_arr = 0
      for l = 0ll,n_elements(wh)-1 do begin
;         ii = fix(vp(wh(l),0)/dv) + (nn/(2*dv))
;         jj = fix(vp(wh(l),1)/dv) + (nn/(2*dv))
;         kk = fix(vp(wh(l),2)/dv) + (nn/(2*dv))
;         vxyp(ii,jj) = 1.0 + vxyp(ii,jj)
;         vxzp(ii,kk) = 1.0 + vxzp(ii,kk)
;         vyzp(jj,kk) = 1.0 + vyzp(jj,kk)
         vpp = reform(vp(wh(l),*))
;      vpp2 = sqrt(vpp(0)^2 + vpp(1)^2 + vpp(2)^2)
;      vpp = vpp/vpp2
         vdotl = transpose(vpp(*))#[0,0,-1]

         if (vdotl gt cos(80*!dtor)) then begin
            e_arr = [e_arr,(vpp(0)^2 + vpp(1)^2 + vpp(2)^2)/mrat(wh(l))]
            cnt_arr = [cnt_arr,vdotl]
         endif
         
      endfor
   endif
   
   help,e_arr
   e_arr = 0.5*1.67e-27*e_arr*1e6/1.6e-19


   h = histogram(e_arr,binsize=dE,min = hmin, max = hmax,reverse_indices=ri)
   h1 = fltarr(n_elements(h))
   for i = 0,n_elements(h)-2 do begin
      if (ri[i] ne ri[i+1]) then begin
         h1(i) = total(cnt_arr(ri(ri(i):ri(i+1)-1)))
      endif
   endfor

   h = h1

   xE = dE*(findgen(n_elements(h)) + hmin)

;rebin to log scale
   s = {energy_bin, e_min: 0.0, e_max: 0.0}
   close,2
   openr,2,'caps_e_bin.dat'

   levst = 0
   lxE = 10
   while not(eof(2)) do begin
      readf,2,s
      emin = s.e_min
      emax = s.e_max      
      if ((emin ge 10.0) and (emax le 10000)) then begin
         wh = where((xE gt emin) and (xE le emax))
                                ; if (wh(0) ge 0) then levst = [levst,total(h(wh))]
         if (wh(0) ge 0) then levst = [total(h(wh)),levst]
         if (wh(0) eq -1) then levst = [0,levst]
         lxE = [lxE,(emax+emin)/2]
      endif
   endwhile

   emin = 10
   emax = 10

;   while (emax lt 10000.0) do begin
;      emax = emin*1.2
;      ;print,emax
;      wh = where((xE gt emin) and (xE le emax))
;      if (wh(0) ge 0) then levst = [levst,total(h(wh))]
;      if (wh(0) eq -1) then levst = [levst,0]
;      lxE = [lxE,(emax+emin)/2]
;      emin=emax
;   endwhile

   wset,1
   lxE = lxE(0:n_elements(lxE)-1)
   plot,lxE,levst,psym=10,/xlog,xrange=[10,10000],/xsty
;   plot,xE,h,psym=10,xrange=[10,10000],/xsty
;   stop
   if (cnt eq 0) then evst = fltarr(nnz+1,n_elements(lxE))
;   whzero = where(levst lt 0.01)
;   levst(whzero) = 0.1

   evst(cnt,*) = levst

   wset,2
   img_cont_ylog,bytscl(evst)<254,findgen(nnz+1),lxE
;   contour,evst,findgen(nnz+1),xE,/ylog,yrange=[10,10000],/ysty,nlev=255,/fill
;   tvscl,evst
   colorbar,position=[!y.window(0), !x.window(1)+0.01, !y.window(1), $
                      !x.window(1)+0.03],ncolors=255,/vertical,/right,$
            maxrange=[max(evst)],divisions=4,title='counts',$
            charsize=1.4
   cnt = cnt+1
endfor

wset,4
plot,b1(xcur,1,*,1),/nodata,yrange=[-0.2,0.6]
oplot,b1(xcur,1,*,1),color=fsc_color('red')
oplot,b1(xcur,1,*,0),color=fsc_color('green')
oplot,b1(xcur,1,*,2),color=fsc_color('blue')


;write postscript
;@gops
;device,/color
;device,filename='caps_synth_espec.ps'
;   img_cont_ylog,bytscl(evst),findgen(nnz+1),lxE
;;   contour,evst,findgen(nnz+1),xE,/ylog,yrange=[10,10000],/ysty,nlev=255,/fill
;;   tvscl,evst
;   colorbar,position=[!y.window(0), !x.window(1)+0.01, !y.window(1), $
;                      !x.window(1)+0.03],ncolors=255,/vertical,/right,$
;            maxrange=[max(evst)],divisions=4,title='counts',$
;            charsize=1.4
;device,/close
;set_plot,'x'



xx = reverse((x - x(55)))
xx = ((x - x(nx-55)))
yy = (y - y(ny/2))

;window,0,xsize=800,ysize=800
;contlevs = [0.02,0.05,0.1,0.2,0.5,1.0]*1e15
;contour,reform(np(*,*,zcur)),levels=contlevs,xx/1e3,yy/1e3,$
;    /isotropic,xrange=[-140,40],yrange=[-70,70],xtitle='x 10!u3!n km',$
;    ytitle='y 10!u3!n km',/xsty,/ysty
;plots,[xx(xcur)-ndx*dx,xx(xcur)+ndx*dx]/1e3,([y(ycur)-ndx*dx,y(ycur)-ndx*dx]-y(ny/2))/1e3,thick=3
;plots,[xx(xcur)-ndx*dx,xx(xcur)+ndx*dx]/1e3,([y(ycur)+ndx*dx,y(ycur)+ndx*dx]-y(ny/2))/1e3,thick=3
;plots,[xx(xcur)-ndx*dx,xx(xcur)-ndx*dx]/1e3,([y(ycur)-ndx*dx,y(ycur)+ndx*dx]-y(ny/2))/1e3,thick=3
;plots,[xx(xcur)+ndx*dx,xx(xcur)+ndx*dx]/1e3,([y(ycur)-ndx*dx,y(ycur)+ndx*dx]-y(ny/2))/1e3,thick=3

;t1 = 25e3/cos(15*!dtor) - atan(15*!dtor)*x
t2 = -13.5e3/cos(15*!dtor) - atan(15*!dtor)*xx
;t3 = 7.5e3/cos(15*!dtor) - atan(15*!dtor)*x
 
;oplot,x/1e3,t1/1e3,linestyle=1,color=255,thick=4
oplot,xx/1e3,t2/1e3,linestyle=2,thick=1
;oplot,x/1e3,t3/1e3,linestyle=1,color=255,thick=4
;oplot,t.x/1e3,t.y/1e3,linestyle=3,thick=1

xrng = 600

device,decompose=0
loadct,39
;szf = size(fvp)
;fvp = smooth(fvp,2)
sz = size(vxyp)
lev=[1/exp(4),1/exp(3),1/exp(2),1/exp(1),1.0]*1.0*max(vxyp)

contour,smooth(vxyp,2),vx,vy,/isotropic,levels=lev,$
  xtitle='vx (km/s)',ytitle='vy (km/s)',$
  xrange=[-xrng,xrng],yrange=[-xrng,xrng],/xsty,/ysty

contour,smooth(vxzp,2),vx,vz,/isotropic,levels=lev,$
 xtitle='vx (km/s)',ytitle='vz (km/s)',$
  xrange=[-xrng,xrng],yrange=[-xrng,xrng],/xsty,/ysty

;xrng=400
vxzp_sm = vxzp 
vxzp_sm = smooth(vxzp,2)
vyzp_sm = median(vyzp,2)

;lev=[1/exp(4),1/exp(3),1/exp(2),1/exp(1),1.0]*1.1*max(vyzp_sm)
contour,smooth(vyzp,2),vy,vz,/isotropic,levels=lev,$
  xtitle='vy (km/s)',ytitle='vz (km/s)',$
  xrange=[-xrng,xrng],yrange=[-xrng,xrng],/xsty,/ysty

;vxy = reform(fvp(*,*,szf(3)/2))
;vxy_sm = smooth(vxy,10)

vxyp_sm = smooth(vxyp,20)
vxzp_sm = smooth(vxzp,20)

plot,vx,vxyp_sm(*,sz(2)/2),xrange=[-600,600],/xsty,xtitle='vx'
oplot,vx,1.0*max(vxyp_sm(*,sz(2)/2))*exp(-(vx+250)^2/139.^2),linestyle=1
plot,vz,vxzp_sm(sz(1)/2,*),xrange=[-600,600],/xsty,xtitle='vz'
oplot,vz,1.0*max(vxzp_sm(sz(1)/2,*))*exp(-vz^2/139.^2),linestyle=1
;contour,smooth(vyzi,2),vy,vz,/isotropic,nlev=10,xtitle='vy (km/s)',ytitle='vz (km/s)',$
;  xrange=[-xrng,xrng],yrange=[-xrng,xrng],/xsty,/ysty,/overplot,$
;  c_color=fsc_color('green')

;wh = where((abs(bzy*vy) le xrng) and (abs(vy) lt xrng))
;plots,vy(wh),bzy(wh)*vy(wh),linestyle=1
;arrow,0,0,lyhat*xrng/2,lzhat*xrng/2,color=fsc_color('red'),/data
;arrow,0,0,nlyhat*xrng/2,nlzhat*xrng/2,color=fsc_color('blue'),/data

;sz = size(fvp)
;scale3,xrange=[0,sz(1)],yrange=[0,sz(2)],zrange=[0,sz(3)]

;shade_volume, fvp,0.1*max(fvp),v,p  
;tv,polyshade(v,p,/t3d)           

;protons

wh = where(fvp gt 0.0)
ijkarr = array_indices(fvp,wh)
iarr = reform(ijkarr(0,*))
jarr = reform(ijkarr(1,*))
karr = reform(ijkarr(2,*))
vxarr = vx(iarr)
vyarr = vy(jarr)
vzarr = vz(karr)
varr = sqrt(vxarr^2 + vyarr^2 + vzarr^2)
vxhat = vxarr/varr
vyhat = vyarr/varr
vzhat = vzarr/varr

vhat = [vxhat,vyhat,vzhat]
lhat = [lxhat,lyhat,lzhat]
nlhat = [nlxhat,nlyhat,nlzhat]

lhat = [lxhat,lyhat,lzhat]
Earr = 0.0
cntarr = 0.0
flg = 'false'
for i = 0,n_elements(wh)-1 do begin
   vhat = [vxhat(i),vyhat(i),vzhat(i)]
   vdotl=transpose(vhat)#lhat
   vdotnl = transpose(vhat)#nlhat
;   print,vdotl
   if ((vdotl lt cos(45*!dtor)) and (abs(vdotnl) lt cos(80*!dtor))) then begin
      Earr = [Earr,0.5*mp*(vxarr(i)^2 + vyarr(i)^2 + vzarr(i)^2)*1e6/1.6e-19] 
;      print,Earr
      flg = 'true'
   endif
endfor

if (flg eq 'false') then begin
   print, 'No proton counts...'
   plot,[10,10000],[0,1],/nodata,title='Protons',/xlog
   goto, SKIP
endif
Earr = Earr(1:*)

dE = 10.
hmin = 1.0
hmax = 10000.

h = histogram(Earr,binsize=dE,min = hmin, max = hmax)
xE = dE*(findgen(n_elements(h)) + hmin)


;!p.multi=[0,1,1]
;window,1,xsize=400,ysize=400
plot,xE,h,psym=10,xtitle='Energy (eV)',ytitle='Counts',/xlog,$
   xrange=[10,10000],yrange=[0,max(h)*1.1],/ysty,/xsty,title='Protons'


SKIP:
;heavy ions

wh = where(fvi gt 0.0)
ijkarr = array_indices(fvi,wh)
iarr = reform(ijkarr(0,*))
jarr = reform(ijkarr(1,*))
karr = reform(ijkarr(2,*))
vxarr = vx(iarr)
vyarr = vy(jarr)
vzarr = vz(karr)
varr = sqrt(vxarr^2 + vyarr^2 + vzarr^2)
vxhat = vxarr/varr
vyhat = vyarr/varr
vzhat = vzarr/varr

vhat = [vxhat,vyhat,vzhat]
lhat = [lxhat,lyhat,lzhat]
nlhat = [nlxhat,nlyhat,nlzhat]

lhat = [lxhat,lyhat,lzhat]
Earr = 0.0
cntarr = 0.0
flg = 'false'
for i = 0,n_elements(wh)-1 do begin
   vhat = [vxhat(i),vyhat(i),vzhat(i)]
   vdotl=transpose(vhat)#lhat
   vdotnl = transpose(vhat)#nlhat
;   print,vdotl
   if ((vdotl lt 0.0) and (abs(vdotnl) lt cos(80*!dtor))) then begin
      Earr = [Earr,0.5*28.0*mp*(vxarr(i)^2 + vyarr(i)^2 + vzarr(i)^2)*1e6/1.6e-19] 
;      print,Earr
      flg = 'true'
   endif
endfor

if (flg eq 'false') then begin
   print, 'No pickup counts...'
   plot,[10,10000],[0,1],/nodata,title='Pickup ions',/xlog
   goto,SKIP_PICKUP
endif
Earr = Earr(1:*)

dE = 10
hmin = 1.0
hmax = 7500.

h = histogram(Earr,binsize=dE,min = hmin, max = hmax)

xE = dE*(findgen(n_elements(h))+hmin)

;!p.multi=[0,1,1]
;window,2,xsize=400,ysize=400
;wh1 = where(xE1 le 100)
;wh2 = where(xE2 gt 100)
plot,xE,h,psym=10,xtitle='Energy (eV)',ytitle='Counts',/xlog,$
   xrange=[10,10000],yrange=[0,max(h)*1.1],/ysty,/xsty,title='Pickup ions'


SKIP_PICKUP:

;device,/close

return
end
;----------------------------------------------------------------------


;----------------------------------------------------------------------
pro read_traj,t
;----------------------------------------------------------------------

t = {trajectory,Julian: 0.0, x: 0.0, y: 0.0, z: 0.0, vx: 0.0, vy: 0.0, $
     vz: 0.0,sx: 0.0, sy: 0.0, sz: 0.0, px: 0.0, py: 0.0, pz: 0.0} 

d = {trajectory,Julian: 0.0, x: 0.0, y: 0.0, z: 0.0, vx: 0.0, vy: 0.0, $
     vz: 0.0,sx: 0.0, sy: 0.0, sz: 0.0, px: 0.0, py: 0.0, pz: 0.0} 

close,1
openr, 1, 'trajectory1.csv'

junk = ' ' 
readf, 1, junk
print,junk


while not(eof(1)) do begin
   readf,1,d
   t = [t,d]
endwhile

t = t(1:*)

close,1

return
end
;----------------------------------------------------------------------


;----------------------------------------------------------------------
;main program
;----------------------------------------------------------------------

;set_plot,'ps'
;;device,/landscape
;!p.font=0
;device,filename='vdist.ps'
;!p.thick=2.0
;!x.thick=2.0
;!y.thick=2.0
!p.multi=[0,2,3]
;!p.charsize=1.5
;@x6x9
;device,/color


;read_traj,t
rundir = '/Volumes/MacD97-2/hybrid/KHI/run_102/'
f_read_coord,rundir+'coord.dat',x,y,z,dzc,dzg,nx,ny,nz

xyz_traj=[nx/2,ny/2,nz/2]
lxyz = [1.0,0.0,0.0]
nlxyz = [0.0,1.0,0.0]

get_dist,xyz_traj,lxyz,nlxyz

;stop

;t2 = -13.5e3/cos(15*!dtor) - atan(15*!dtor)*xx
;dx = 2200.

;xx = x - x(55)
;yy = y - y(ny/2)

;nstep = 10
;for i = 0,nstep-1 do begin 
;   isc = (nstep-i)*10
;   jsc = ny/2 - round((13.5e3/dx)/cos(15*!dtor)) - atan(15*!dtor)*(isc-(nx-55))
;   print,isc,jsc
;   xyz_traj=[isc,jsc,nz/2-35]
;   get_dist,xyz_traj,lxyz,nlxyz,t
;endfor

;evice,/close

stop
end
;----------------------------------------------------------------------
