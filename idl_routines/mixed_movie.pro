close,1
rundir = '/Volumes/Scratch/hybrid/KHI/'+'run_RT'
;rundir = './'
f_read_coord,rundir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz

nfrm = 10

loadct,33
s_max = 205
s_min = 0
nnx = nx*2
nny = nz
stretch,s_min,s_max
device,decompose=0

zm = 2.0

img = fltarr(nnx,nny,nfrm)
 
; Initialize XINTERANIMATE: 
XINTERANIMATE, SET=[zm*nnx,zm*nny, nfrm], /SHOWLOAD 
 
lnvy = fltarr(nfrm)

for i = 1,nfrm do begin 
   f_read_3d_m_32,rundir+'/mixed_1',i,np

;   f_read_3d_vec_m_32,'run2/b1all_2',i,b1

;   lnvy(i-1) = alog(max(abs(up(*,1,30:120,2))))
   h=bytscl(reform(np(*,1,*)<13e15))

;   h=rebin(reform(sqrt(b1(*,1,*,0)^2+b1(*,1,*,2)^2)),nnx,nny)
   img(0:nx-1,*,i-1) = h(*,*)
   img(nx-1:2*nx-3,*,i-1) = h(1:*,*)
;   xinteranimate, frame = i-1, image = img<s_max
endfor

;img = bytscl(img)

img = bytscl(img)
;img(*,nz/2-96/2,*) = 255
;img(*,nz/2+96/2,*) = 255


for i = 0,nfrm-1 do begin
   xinteranimate, frame = i, image = rebin(img(*,*,i)<255,zm*nnx,zm*nny)
endfor

;plot,lnvy

xinteranimate,/keep_pixmaps



end
