
vo = 250e3

dt = 0.5
nx = 279.0
dx = 300e3  ;m

m1 = double(1.67e-27)
m2 = double(1.67e-27)

n1 = 0.1e6 ;m^-3
n2 = 0.1e6 ;m^-3


rho1 = m1*n1
rho2 = m2*n2

v1 = 250e3
v2 = -250e3

vrel = v2-v1
   
phi = 2.0*!pi/180.
b1 = (5e-9)*sin(phi)
b2 = (5e-9)*sin(phi)

muo = !pi*4e-7

m = [1,2,3,4,5,6,7,8]
lambda = nx*dx/m
  
va1 = b1/sqrt(muo*rho1)
va2 = b2/sqrt(muo*rho2)

print,'alfven velocity...',va1/1e3,va2/1e3

vrel2 = ((rho1 + rho2)/(muo*rho1*rho2))*(b1^2 + b2^2)
   
print,sqrt(vrel2)/1e3,vrel
stop

t1 = (rho1*rho2*vrel^2/(rho1+rho2)^2)
t2 = rho1*va1^2/(rho1+rho2)
t3 = rho2*va2^2/(rho1+rho2)
   
;print,t1,t2,t3
;stop

gamma = (2*!pi/lambda)*sqrt(t1 + t2 + t2)
!p.multi=[0,1,2]
plot,m,gamma*2*4*300./500.,psym=6,xtitle='m',ytitle='growth rate (1/s)'
plot,m,(1/gamma),psym=6,xtitle='m',ytitle='growth time (s)'
;stop
print,m,gamma,(1/gamma)

  
print,'growth rate, tau...',gamma,1/gamma
stop
   
;   print,'gyroperiod...',1/(1.6e-19*0.2e-9/1.67e-27)
   
   
;   print,'vrel^2...',vrel^2,(((rho1+rho2)/(rho1*rho2))*(rho1*va1^2 +
;   rho2*va2^2))

rhs = ((rho1+rho2)/(rho1*rho2))*(rho1*va1^2 + rho2*va2^2)

print,'vel^2 = ',vrel^2,rhs
print,'rhs....',vrel^2/rhs
   
end
