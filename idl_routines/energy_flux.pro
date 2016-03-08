close,1
rundir = '/Volumes/MacD72-2/hybrid/KHI/'+'run_80'
f_read_coord,rundir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz

nfrm = 6

zslc = nz/2 +20

qoverm = 1.6e-19/1.67e-27
muo = !pi*4e-7
mp = 1.67e-27

vsw = 2*250e3
vth = 200e3
n = 0.1e6
m = mp
bo = 5e-9

Efo = 0.5*n*m*vsw*vsw*vsw + 1.5*n*m*vth^2*vsw + vsw*bo^2/muo

print,Efo 

eflx = 0.0

for i = 1,nfrm do begin 
   f_read_3d_m_32,rundir+'/temp_p_1',i,t
   f_read_3d_m_32,rundir+'/npall_1',i,np
   np = np/1e9
   f_read_3d_vec_m_32,rundir+'/upall_1',i,up
   up = up*1e3
   f_read_3d_vec_m_32,rundir+'/b1all_1',i,b
   b = b/qoverm
   f_read_3d_vec_m_32,rundir+'/Eall_1',i,E
   E = E*1e3/qoverm

   S = (E(*,1,*,0)*b(*,1,*,1) - E(*,1,*,1)*b(*,1,*,0))/muo
   Ek = 0.5*np(*,1,*)*mp*up(*,1,*,2)^2*up(*,1,*,2)
   Ei = (1.5)*np(*,1,*)*t(*,1,*)*1.6e-19*up(*,1,*,2)

   S = reform(S)
   Ek = reform(Ek)
   Ei = reform(Ei)


   dx = 300e3
   Efo_tot = Efo*(nz/2)*dx*dx
   A = dx*dx
   Ef_tot = (total(S(*,zslc))+ total(Ek(*,zslc)) + total(Ei(*,zslc)))*A

   eflx = [eflx,Ef_tot]

endfor

S = reform(S)
Ek = reform(Ek)
Ei = reform(Ei)



!p.multi=[0,2,2]

plot,S(*,zslc)/Efo
plot,Ek(*,zslc)/Efo
plot,Ei(*,zslc)/Efo

dx = 300e3
Efo_tot = Efo*(nz/2)*dx*dx
A = dx*dx
Abox = 279*dx*dx
Asat = 4*6e7*dx*45
Asat = 4*40*6e7^2
Ef_tot = (total(S(*,zslc))+ total(Ek(*,zslc)) + total(Ei(*,zslc)))*A

Esw = 0.5*0.05e6*1.67e-27*450e3^3*(4*40*6e7^2)

print,'total energy flux in...',Efo_tot
print,'total S................',total(S(*,zslc))*A
print,'total Ek...............',total(Ek(*,zslc))*A
print,'total Ei...............',total(Ei(*,zslc))*A
print,'total Ef...............',Ef_tot*Asat/Abox,100*Ef_tot/Efo_tot,$
             100*Ef_tot*Asat/Abox/Esw


plot,eflx/Efo_tot

end

