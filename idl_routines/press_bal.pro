muo = !pi*4e-7
mi = 1.67e-27

n1 = 0.4e6
n2 = 0.07e6

T1 = 200.0*1.6e-19 ;input 1-d thermal velocity
T2 = 1000.0*1.6e-19

b1 = 5e-9
b2 = 5e-9


eb1 = b1^2/(2*muo)
eb2 = b2^2/(2*muo)

ep1 = (3./2.)*n1*T1
ep2 = (3./2.)*n2*T2

etot1 = eb1 + ep1
etot2 = eb2 + ep2


print,'etot1...',etot1,eb1,ep1
print,'etot2...',etot2

;ep2 = c - eb2
;print,'ep2...',ep2,ep1
T2 = n1*T1/n2
print,'T1...',T1/1.6e-19
print,'T2...',T2/1.6e-19

print,'beta 1...',ep1/eb1,ep1+eb1
print,'beta 2...',ep2/eb2,ep2+eb2

;print,T1/1.6e-19,T2/1.6e-19

vth1 = sqrt(3*T1/mi)
vth2 = sqrt(3*T2/(8*mi))


print,'thermal velocity...',vth1/1e3,vth2/1e3


z = findgen(100)

n = 0.5*(n1 + n2) + $
      0.5*(n1 - n2)*tanh((z-50)/10)

vth = (0.5*(vth1 + vth2) + $
      0.5*(vth1 - vth2)*tanh((z-50)/10))

T = (0.5*(T1 + T2) + $
      0.5*(T1 - T2)*tanh((z-50)/10))

!p.multi=[0,2,2]
;plot,z,vth


;T = (1./3.)*mi*vth^2


;plot,z,T
;plot,z,n
;plot,z,n*t
;oplot,z,n^2, linestyle=2


;solve numerically for pressure balance

P = n1*T1
n0 = P/T

erase
plot,z,T
plot,z,n0
oplot,z,n,linestyle=1
plot,z,n0*T


end
