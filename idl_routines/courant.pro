;pro courant

no=double(0.4e15/17)   ;km^-3  pluto 0.01
Bo=double(5e-9)  ;pluto 0.2

q=double(1.6e-19)
mo = double(17*1.67e-27)
mel = double(9.1e-31)
muo = !pi*4e-7/1e3
c = 3e8
ep = 8.85e-12
wpi = sqrt((no/1e9)*q*q/(ep*mo))
wpe = sqrt((no/1e9)*q*q/(ep*mel))

print,'wpi...',wpi

alpha = muo*q*q/mo
print,'alpha...',alpha

Bo = Bo*q/mo     

dx = 250.0     ;km
 
k = !pi/(dx)

a1 = k^2*Bo/(alpha*no)
a2 = (k*Bo)^2/(alpha*no)

w = 0.5*(a1 + sqrt(a1^2 + 4*a2))

phi = w/k

dt = dx/phi

vsw = 175.0
rL = mo*vsw/(q*Bo*mo/q)

print,'phi,dt....',phi,dt
print,'Whister mode....',phi
print,'Alfven mode.....',(Bo*mo/q)/sqrt(muo*no*mo)
print,'Magnetic field..',Bo,(Bo*mo/q)
print,'proton qyroradius (km)...',rL,rL/dx
print,'pickup qyroradius (km)...',28*rL,28*rL/dx
print,'ion inertial length (km)...',(c/wpi)/1e3,(c/wpi)/1e3/dx
print,'electron inertial length (km)...',(c/wpe)/1e3


; source density
nsrc = 4000e15  ;km^3
bsrc = 300e-9  ;T
vsrc = 20.0 ;km/s
va_src = (bsrc)/sqrt(muo*nsrc*mo)

print,'source, va/vth...',va_src/vsrc,va_src,vsrc

vol = dx*16.*dx*16.*dx*55.
beta = (150000/vol)/no

print,'beta....',beta


end


