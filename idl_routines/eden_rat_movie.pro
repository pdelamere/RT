close,1

f_read_coord,'run1/coord.dat',x,y,z,dzc,dzg,nx,ny,nz

nfrm = 134

loadct,33
s_max = 255
s_min = 0
nnx = nx*4
nny = 141*4
stretch,s_min,s_max
device,decompose=0

eden = 0.5*11e6*1.67e-27*200e3^2 ;+ (1./3.)*1.67e-27*11e6*100e3^2

bflx = 16e-9*sin(10*!dtor)
L = x(nx-1)

print,(bflx^2/(2*!pi*4e-7))*(L/(70.)),eden
stop

img = fltarr(2*nnx,nny,nfrm)
 
; Initialize XINTERANIMATE: 
XINTERANIMATE, SET=[2*nnx,nny, nfrm], /SHOWLOAD 
 
lnvy = fltarr(nfrm)

for i = 1,nfrm do begin 
;   f_read_3d_m_32,'run1/npall_1',i,np

   f_read_3d_vec_m_32,'run1/b1all_1',i,b1

   b1 = b1*1.67e-27/1.6e-19

;   lnvy(i-1) = alog(max(abs(up(*,1,30:120,2))))
   
;   h=bytscl(rebin(reform(np(*,1,30:130)),nnx,nny))
   h=rebin(reform((0*b1(*,1,10:150,0)^2+b1(*,1,10:150,2)^2)/(2*!pi*4e-7)),nnx,nny)
   img(0:nnx-2,*,i-1) = h(0:nnx-2,*)
   img(nnx-1:2*nnx-3,*,i-1) = h(1:*,*)
;   xinteranimate, frame = i-1, image = img<s_max
endfor


img = img/eden
print,max(img)
;stop

img = bytscl(img)

for i = 0,nfrm-1 do begin
   xinteranimate, frame = i, image = img(*,*,i)
;   contour,img(*,*,i),levels=[0.1,0.2,0.5,1.0,1.2,1.5]
endfor


;plot,lnvy

xinteranimate,/keep_pixmaps



end
