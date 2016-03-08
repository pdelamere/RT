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

close,1
rundir = '/Volumes/Macd97-2/hybrid/KHI/run_test'
f_read_coord,rundir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz
;@gops
;device,/color
;device,/landscape
;@x6x9

;set_plot,'ps'
;!p.font=0
;;!p.thick=3
;!x.thick=3
;!y.thick=3
;device,filename='z_slice.eps'
;device,/color
;;device,decompose=0
!p.charsize=1.5
!x.charsize=1.5
!y.charsize=1.5

nframe= 2

ninit = 0.4e6
mb = 1
So = ninit*250.0*1.6e-19/(1.67e-27*ninit)^(5./3)
print,So
;stop

;im = tvrd()
;sz = size(im)
nnx = 1200
nny = 700


b12d = fltarr(nx,nz,nframe)
s2d = fltarr(nx,nz,nframe)
p2d = fltarr(nx,nz,nframe)
npt2d = fltarr(nx,nz,nframe)
npb2d = fltarr(nx,nz,nframe)
ti2d = fltarr(nx,nz,nframe)
up2d = fltarr(nx,nz,3,nframe)
rho2d = fltarr(nx,nz,nframe)

for i = 0,nframe-1 do begin 
                                ;f_read_3d_m_32,rundir+'/npall_1',i,np
   f_read_3d_m_32,rundir+'/mixed_1',i+1,mix
   f_read_3d_m_32,rundir+'/npall_1',i+1,np
   f_read_3d_m_32,rundir+'/mnp_1',i+1,mnp
   f_read_3d_m_32,rundir+'/np_b_1',i+1,npb
   f_read_3d_m_32,rundir+'/temp_p_1',i+1,ti
                                ;f_read_3d_m_32,rundir+'/np_b_1',i,np
   f_read_3d_vec_m_32,rundir+'/upall_1',i+1,up
   f_read_3d_vec_m_32,rundir+'/b1all_1',i+1,b1
   b1 = b1*1.67e-27/1.6e-19
;   ti=reform(ti(*,1,*))
   b12d(*,*,i)=reform((1.0*(b1(*,1,*,2))^2 + 1.0*(b1(*,1,*,0))^2) + 0.0*(b1(*,1,*,1))^2)
   mp = 1.67e-27
   gamma= 5./3.
   np = np/1e9                  ;convert km to m
   npb = npb/1e9
   ti = ti*1.6e-19              ;convert eV to J
   ti(*,1,nz-2) = ti(*,1,nz-3)
   ti(*,1,nz-1) = ti(*,1,nz-3)

   mnp = mnp/1e9

   p = (np+npb)*ti/1e-11        ;convert to nPa'
   s =  double((np+npb)*ti)/double(mnp)^gamma

      
   b12d(*,*,i) = b12d(*,*,i)/(2*!pi*4e-7)/1e-11
   s2d(*,*,i) = smooth((reform(s(*,1,*))),2)/So   ;normalize to initial entropy
   p2d(*,*,i) = reform(p(*,1,*))
   npt2d(*,*,i) = reform(np(*,1,*))/1e6
   npb2d(*,*,i) = reform(npb(*,1,*))/1e6
   ti2d(*,*,i) = reform(ti(*,1,*))
   up2d(*,*,*,i) = reform(up(*,1,*,*))
   rho2d(*,*,i) = reform(npt2d(*,*,i) + mb*npb2d(*,*,i))

endfor



bimg = bytscl(b12d)
simg = bytscl(s2d)
pimg = bytscl(p2d)
nptimg = bytscl(npt2d)
npbimg = bytscl(npb2d)
tiimg =  bytscl(ti2d)
rhoimg = bytscl(rho2d)
   
xx = x
z = z - z(nz/2)

z1=120/2
z2 = 120/2

x1 = 0
x2 = nx-1

zz = z(nz/2-z1:nz/2+z2)

xx  = x(x1:x2)

bimg_arr=bimg(x1:x2,nz/2-z1:nz/2+z2,*)
simg_arr=simg(x1:x2,nz/2-z1:nz/2+z2,*)
pimg_arr=pimg(x1:x2,nz/2-z1:nz/2+z2,*)
tiimg_arr = tiimg(x1:x2,nz/2-z1:nz/2+z2,*)
npt_arr = nptimg(x1:x2,nz/2-z1:nz/2+z2,*)
npb_arr = npbimg(x1:x2,nz/2-z1:nz/2+z2,*)
rho_arr = rhoimg(x1:x2,nz/2-z1:nz/2+z2,*)


sz = size(bimg(x1:x2,nz/2-z1:nz/2+z2,0))
consx =5
consz =consx

xxx = congrid(xx,sz(1)/consz)
zzz = congrid(zz,sz(2)/consz)



zm =1
nstartup = 0
nf_tot = nframe-nstartup
XINTERANIMATE, SET=[zm*nnx,zm*nny, nf_tot], /SHOWLOAD 


for nf = nstartup,nframe-1 do begin
   nfrm = nf+1

   loadct,39
;s_max = 155
;s_min = 0
   nnx = nx
   nny = nz
;stretch,s_min,s_max
   device,decompose=0
   
   zm = 1
   xslc=60
   zslc = 295
;@x6x9
   !p.multi=[0,2,3]
   
   keden = 0.5*1.67e-27*0.4e6*250e3^2
   
   loadct,39
   img_cont,bimg_arr(*,*,nf),xx/1e3,zz/1e3,nfrm
   
   mxb12d = max(b12d(x1:x2,nz/2 - z1: nz/2 + z2,*))
;   print,mxb12d
;stop
   colorbar,minrange=0,maxrange=mxb12d,format='(f5.2)',$
            position=[!y.window(0), !x.window(1)+0.01, !y.window(1), $
                      !x.window(1)+0.03],ncolors=255,/vertical,/right,$
            title='magnetic energy density (10!u-11!n Pa)',charsize=1.2,divisions=5,minor = 0,yaxis = 1

vx = congrid(reform(up2d(x1:x2,nz/2-z1:nz/2+z2,0,nf)),sz(1)/consz,sz(2)/consz)
vz = congrid(reform(up2d(x1:x2,nz/2-z1:nz/2+z2,2,nf)),sz(1)/consz,sz(2)/consz)
   
   
   velovect,vx(1:*,1:*),vz(1:*,1:*),xxx(1:*)/1e3,zzz(1:*)/1e3,/overplot,color=fsc_color("white"),thick=1.5,length=2.0
   
   
   
   img_cont,pimg_arr(*,*,nf),xx/1e3,zz/1e3,nfrm
   
   mxp = max(p2d(x1:x2,nz/2 - z1: nz/2 + z2,*))
   colorbar,minrange=0,maxrange=mxp,format='(f5.2)',$
            position=[!y.window(0), !x.window(1)+0.01, !y.window(1), $
                      !x.window(1)+0.03],ncolors=255,/vertical,/right,$
            title='Pressure (10!u-11!n Pa)',charsize=1.2,divisions=5,minor = 0,yaxis = 1
   
   velovect,vx(1:*,1:*),vz(1:*,1:*),xxx(1:*)/1e3,zzz(1:*)/1e3,/overplot,color=fsc_color("white"),thick=1.5,length=2.0
   
   
;   img_cont,tiimg_arr(*,*,nf),xx/1e3,zz/1e3,nfrm
   img_cont,rho_arr(*,*,nf),xx/1e3,zz/1e3,nfrm
   
;   mxti = max(ti2d(x1:x2,nz/2 - z1: nz/2 + z2,*)/1e-19/1e2)
  mxti = max(rho2d(x1:x2,nz/2 - z1: nz/2 + z2,*))
   colorbar,minrange=0,maxrange=mxti,format='(f5.2)',$
            position=[!y.window(0), !x.window(1)+0.01, !y.window(1), $
                      !x.window(1)+0.03],ncolors=255,/vertical,/right,$
            title='mass density (AMU/cm!u-3!d)',charsize=1.2,divisions=5,minor = 0,yaxis = 1
   
   velovect,vx(1:*,1:*),vz(1:*,1:*),xxx(1:*)/1e3,zzz(1:*)/1e3,/overplot,color=fsc_color("white"),thick=1.5,length=2.0


;   img_cont,tiimg_arr(*,*,nf),xx/1e3,zz/1e3,nfrm

;   mxti = max(ti2d(x1:x2,nz/2 - z1: nz/2 + z2,*)/1e-19/1e2)
;   colorbar,minrange=0,maxrange=mxti,format='(f5.2)',$
;            position=[!y.window(0), !x.window(1)+0.01, !y.window(1), $
;                      !x.window(1)+0.03],ncolors=255,/vertical,/right,$
;            title='Temperature (100 eV)',charsize=1.2,divisions=5,minor = 0,yaxis = 1
   
;   velovect,vx(1:*,1:*),vz(1:*,1:*),xxx(1:*)/1e3,zzz(1:*)/1e3,/overplot,color=fsc_color("white"),thick=1.5,length=2.0
   
      
   img_cont,npt_arr(*,*,nf),xx/1e3,zz/1e3,nfrm
   mxnpt = max(npt2d(x1:x2,nz/2 - z1: nz/2 + z2,*))
   colorbar,minrange=0,maxrange=mxnpt,format='(f5.2)',$
            position=[!y.window(0), !x.window(1)+0.01, !y.window(1), $
                      !x.window(1)+0.03],ncolors=255,/vertical,/right,$
            title='n!di!n top (cm!u-3!n)',charsize=1.2,divisions=5,minor = 0,yaxis = 1
   
   velovect,vx(1:*,1:*),vz(1:*,1:*),xxx(1:*)/1e3,zzz(1:*)/1e3,/overplot,color=fsc_color("white"),thick=1.5,length=2.0
   
   
   img_cont,npb_arr(*,*,nf),xx/1e3,zz/1e3,nfrm
   mxnpb = max(npb2d(x1:x2,nz/2 - z1: nz/2 + z2,*))
   colorbar,minrange=0,maxrange=mxnpb,format='(f5.2)',$
            position=[!y.window(0), !x.window(1)+0.01, !y.window(1), $
                      !x.window(1)+0.03],ncolors=255,/vertical,/right,$
            title='n!di!n bottom (cm!u-3!n)',charsize=1.2,divisions=5,minor = 0,yaxis = 1
   
   velovect,vx(1:*,1:*),vz(1:*,1:*),xxx(1:*)/1e3,zzz(1:*)/1e3,/overplot,color=fsc_color("white"),thick=1.5,length=2.0
   
   s2d = s2d
   img_cont,simg_arr(*,*,nf),xx/1e3,zz/1e3,nfrm
   
   mxs2d = max(s2d(x1:x2,nz/2 - z1: nz/2 + z2,*))
   mins2d = min(s2d(x1:x2,nz/2 - z1: nz/2 + z2,*))
   colorbar,minrange=mins2d,maxrange=mxs2d,format='(f5.2)',$
            position=[!y.window(0), !x.window(1)+0.01, !y.window(1), $
                      !x.window(1)+0.03],ncolors=255,/vertical,/right,$
            title='normalized (P/rho^gamma)',charsize=1.2,divisions=5,minor = 0,yaxis = 1
   
   velovect,vx(1:*,1:*),vz(1:*,1:*),xxx(1:*)/1e3,zzz(1:*)/1e3,/overplot,color=fsc_color("white"),thick=1.5,length=2.0
   
   im = tvrd(true=1)
   tv,im

   xinteranimate, frame = nf-nstartup, image = im
;stop   
   
endfor
;device,/close

xinteranimate,/keep_pixmaps







end
