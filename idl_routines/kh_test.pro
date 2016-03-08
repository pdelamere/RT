bsh = 16e-9*0
bms = 16e-9*0

nsh = 11.0e6
nms = 11.0e6


mp = 1.67e-27
msh = mp
mms = mp

muo = !pi*4e-7

va_sh = bsh/sqrt(muo*nsh*msh)

va_ms = bms/sqrt(muo*nms*mms)

rho_sh = msh*nsh
rho_ms = mms*nms

vrel = sqrt(((rho_sh + rho_ms)/(rho_ms*rho_sh))*(rho_sh*va_sh^2 + rho_ms*va_ms^2))

print,'vrel...',vrel/1e3

c = 3e8
q = double(1.6e-19)
wpi = sqrt(nsh*q*q/(8.85e-12*mp))
print,c/wpi

lambda = 3*70e3  
lambda = findgen(100)*c/wpi
vr = 400e3
gamma = (2*!pi/lambda)*sqrt(rho_ms*rho_sh*vr^2/(rho_ms+rho_sh)^2 + $
                            (rho_ms*va_ms^2 + rho_sh*va_sh^2)/(rho_ms+rho_sh))

plot,lambda/1e3,1/gamma

print,'t1...',rho_ms*rho_sh*vr^2/(rho_ms+rho_sh)^2
print,'t2...',(rho_ms*va_ms^2 + rho_sh*va_sh^2)/(rho_ms+rho_sh)


end
