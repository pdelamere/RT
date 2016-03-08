close,1
rundir = '/Volumes/Scratch/hybrid/KHI/'+'run_test'
f_read_coord,rundir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz


nfrm=2

f_read_3d_m_32,rundir+'/temp_p_1',nfrm,tp
f_read_3d_m_32,rundir+'/npall_1',nfrm,np_t
f_read_3d_m_32,rundir+'/np_b_1',nfrm,np_b

np = np_t+np_b

nt = np*tp/1e15 

!p.multi=[0,1,2]
surface,reform(nt(*,1,*)),charsize=2.0
im = image(reform(nt(*,1,*)), rgb_table=33)
plot,smooth(nt(10,1,*),2)

end
