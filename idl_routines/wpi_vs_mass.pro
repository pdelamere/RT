@gops
device,filename='mode_suppress.ps'
!p.charsize=1.0
!x.charsize=1.0
!y.charsize=1.0
c = 3e8

ep = 8.85e-12

mp = 1.67e-27

nx = 279.0
dx = 300.0

m = 1+indgen(40)

rho = 0.1e6*mp

mass = (1+indgen(40))*mp

ni = rho/mass
q = 1.6e-19
b = 5e-9
vth = 200e3

wpi = sqrt(ni*q*q/(ep*mass))
lpi = c/wpi

;set_plot,'x'
;plot,mass/mp,lpi/1e6
;stop

!p.multi=[0,2,1]
!x.title='mass (amu)'

lpi_arr = fltarr(n_elements(mass),n_elements(lpi))

for j=0,n_elements(lpi)-1 do begin
   lpi_arr(*,j) = nx*dx*1e3/lpi(j)
endfor

lev=[2,4,6,8,10,20,50,100]

contour,lpi_arr,mass/mp,lpi/1e6,levels=lev,$,
   c_label=[1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,1,0,1,0,1],$
   c_linestyle=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],$
   xrange=[1,40],/xsty,ytitle='ion inertial length (10!u3!n km)'
;stop

oplot,mass/mp,lpi/1e6

;for i = 1,20 do begin
;   oplot,[!x.crange(0),!x.crange(1)],[dx*nx/i,dx*nx/i],linestyle=1
;endfor

contour,lpi_arr,mass/mp,lpi/1e6,levels=lev,$,
   c_label=[1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,1,0,1,0,1],$
   c_linestyle=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],$
   xrange=[1,40],/xsty,ytitle='ion thermal gyroradius (10!u3!n km)'

oplot,mass/mp,mass*vth/(q*b)/1e6

;for i = 1,20 do begin
;   oplot,[!x.crange(0),!x.crange(1)],[dx*nx/i,dx*nx/i],linestyle=1
;endfor

device,/close


end
