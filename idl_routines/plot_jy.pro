
nfrm = 90
device,decompose=0
LOADCT, 39, /SILENT 

f_read_coord,'run1/coord.dat',x,y,z,dzg,dzc,nx,ny,nz
zsz = 60
xx = x(*)
zz = z(nz/2-zsz:nz/2+zsz)



f_read_3d_vec_m_32,'run1/ajall_5',0,up
jy = reform(up(*,1,nz/2-zsz:nz/2+zsz,1))

sz = size(jy)
jyarr = fltarr(sz(1),sz(2),nfrm)
jyarr(*,*,0) =jy

for i = 1,nfrm-1 do begin
   f_read_3d_vec_m_32,'run1/ajall_5',i,up

   jy = reform(up(*,1,nz/2-zsz:nz/2+zsz,1))

   jyarr(*,*,i) = jy(*,*) 

endfor

   xx = x(*)
   zz = z(nz/2-zsz:nz/2+zsz)

jyarr = bytscl(abs(jyarr))


for i = 1,nfrm-1 do begin

   contour,smooth(jyarr(*,*,i),2),xx,zz,/isotropic,/fill,nlev=255,/xsty,/ysty

endfor

end
