nfile = '1'
nfrm = 1

for l = 1,nfile do begin

nfile = string(l)

part_dir = 'part_'+strtrim(string(nfile),2)
fluid_fields_dir = '/Volumes/MacD72-2/hybrid/KHI/run_test'

f_read_coord,fluid_fields_dir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz
f_read_3d_vec_m_32,fluid_fields_dir+'/up_b_'+strtrim(string(nfile),2),nfrm,up

curlv = fltarr(nx,nz)

up2d = reform(up(*,ny/2,*,*))
dx = 1800.


for i = 1,nx-2 do begin
   for k = 1,nz-2 do begin
      km = k-1
      im = i-1
      curlv(i,k) = (up2d(i,k,2)-up2d(im,k,2))/dx - (up2d(i,k,0)-up2d(i,km,0))/dx
   endfor
endfor

loadct,14
device,decompose=0
contour,curlv,nlev=255,/fill,/isotropic

endfor

end
