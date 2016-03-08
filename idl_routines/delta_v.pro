nx = 100
ny = 100

x = findgen(nx)-nx/2
y = findgen(ny)-ny/2

v1 = 1.0

v = fltarr(nx,ny,2)
Lx = nx
Lo = 2.0

plot,y,tanh(y/Lo)*cosh(y/Lo)^(-2)
plot,y,cosh(y/Lo)^(-2)

for i = 0,nx-1 do begin
   for j = 0,ny-1 do begin
      v(i,j,1) = v1*sin(!pi*x(i)/(nx/2))*(cosh(y(j)/Lo))^(-2)
      v(i,j,0) = -2.0*v1*cos(!pi*x(i)/(nx/2))*tanh(y(j)/Lo)*(cosh(y(j)/Lo))^(-2)
   endfor
endfor

velovect,v(*,*,0),v(*,*,1)




end

