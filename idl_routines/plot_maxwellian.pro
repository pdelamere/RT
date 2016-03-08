nx = 100

dx = 300
f = fltarr(nx,nx,nx)
rl = fltarr(nx,nx,nx)
b = 5e-9
m = 2*1.67e-27


vsw = 250.0
vo = 1000.0
vth = 200.0

vx = findgen(nx)*vo/nx - vo/2
vy = findgen(nx)*vo/nx - vo/2
vz = findgen(nx)*vo/nx - vo/2


for i= 0,nx-1 do begin
   for j = 0,nx-1 do begin
      for k = 0,nx-1 do begin
         f(i,j,k) = exp(-((vx(i)-vsw)^2 + vy(j)^2 + vz(k)^2)/vth^2)
         v = sqrt((vx(i))^2 + vy(j)^2 + vz(k)^2)
         rL(i,j,k)=m*v/(1.6e-19*b)
      endfor
   endfor
endfor

contour,f(*,*,nx/2),vx,vy,levels=[exp(-2),exp(-1),0.9],/isotropic,$
   xrange=[-200,600],yrange=[-200,200]


rd = fsc_color('red')
contour,rL(*,*,nx/2),vx,vy,levels=[1,2,3,4,5]*dx,/isotropic,$
   c_label=[1,1,1,1,1],$
   xrange=[-200,600],yrange=[-200,200],/noerase,c_color=rd

wh = where(rL le 1.0*dx)

print,'fraction of under-resolved gyromotion (dx)....',total(f(wh))/total(f)

wh = where(rL le 2.0*dx)

print,'fraction of under-resolved gyromotion (2dx)....',total(f(wh))/total(f)

end
