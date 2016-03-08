close,1
rundir = '/Volumes/Scratch/hybrid/KHI/run_test'

f_read_coord,rundir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz

nfrm = 17

loadct,33
s_max = 255
s_min = 20
nnx = nx
nny = nz
stretch,s_min,s_max
device,decompose=0

zm = 2

img = fltarr(nnx,nny,nfrm)
 
; Initialize XINTERANIMATE: 
XINTERANIMATE, SET=[zm*nnx,zm*nny, nfrm], /SHOWLOAD 
 
lnvy = fltarr(nfrm)

for i = 1,nfrm do begin 
;   f_read_3d_m_32,'run1/npall_1',i,np

   f_read_3d_vec_m_32,rundir+'/b1all_1',i,b1

;   lnvy(i-1) = alog(max(abs(up(*,1,30:120,2))))
   
;   h=bytscl(rebin(reform(np(*,1,30:130)),nnx,nny))
   h=reform(sqrt(b1(*,1,*,0)^2+0.*b1(*,1,*,1)^2)+b1(*,1,*,2)^2)
   print,i

   img(0:nx-1,*,i-1) = h(*,*)
;   img(nx-1:2*nx-2,*,i-1) = h(*,*)
;   xinteranimate, frame = i-1, image = img<s_max
endfor

img = bytscl(img)
img(*,nz/2-96/2,*) = 255
img(*,nz/2+96/2,*) = 255

for i = 0,nfrm-1 do begin
   xinteranimate, frame = i, image = rebin(img(*,*,i),nnx*zm,nny*zm)
endfor

;plot,lnvy

xinteranimate,/keep_pixmaps



end
