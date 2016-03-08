b = 16e-9

theta = 18.0

bx = b*sin(theta*!dtor)

nsw = 7.0e6
vsw = 200e3
vrel = vsw*2
msw = 1.67e-27



muo = !pi*4e-7

rj = 7.14e4

dx = 70.0
nx =305.0

L = dx*nx

bEden = bx^2/(2*muo)
kEden = 0.5*nsw*msw*vrel^2

print,'max vortex size (nz)...',kEden/bEden
print,'max wave length (nx)...',kEden/((bx/sqrt(2))^2/(2*muo))

print,((bx/sqrt(2))^2/(2*muo))*L/dx
print,0.5*nsw*msw*vrel^2


print,'vortex size....',(L/2)/rj


end
