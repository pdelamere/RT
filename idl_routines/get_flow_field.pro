
rundir = '/Volumes/MacD72-2/hybrid/KHI/run_6/'
f_read_coord,rundir+'coord.dat',x,y,z,dzg,dzc,nx,ny,nz
f_read_3d_vec_m_32,rundir+'upall_1',20,up

;vx = reform(up(*,1,nz/2-20:nz/2+20,0))
;vz = reform(up(*,1,nz/2-20:nz/2+20,2))

zsz = 40
xsz = 75
vx = reform(up(*,1,nz/2-zsz:nz/2+zsz,0))
vz = reform(up(*,1,nz/2-zsz:nz/2+zsz,2))

vx = smooth(vx,2)
vz = smooth(vz,2)

device,decompose=0
LOADCT, 39, /SILENT 
TVLCT, r, g, b, /GET 
rgbTable = [[r],[g],[b]] 


ivector,vx,vz,scale_isotropic=1,streamline_stepsize=1.0,$
        /streamlines,$
        streamline_nsteps=10,x_streamparticles=50,y_streamparticles=50,$
        view_number=1,auto_color=1,rgb_table=rgbtable,subsample_method=1,$
        view_zoom=1.2,head_size=1.0,dimensions= [1200,600]
;ivector,vx,vz,scale_isotropic=1,identifier = 1.0,data_location=1,$
;    view_zoom=1.0,/no_saveprompt,dimensions=[1200,600]

;stop

f_read_3d_vec_m_32,rundir+'upall_1',1,up



;ivector,vx,vz,scale_isotropic=1,streamline_stepsize=0.05,$
;      /overplot,/streamlines,$
;     streamline_nsteps=100,x_streamparticles=25,y_streamparticles=25


xx = x(*)
zz = z(nz/2-zsz:nz/2+zsz)


;vx(1,1,0) = 2.0*up(nx/2,ny/2,nz/2,1)

ivector,vx,vz,xx,zz,scale_isotropic=1,streamline_stepsize=10,$
        /streamlines,$
        streamline_nsteps=10,x_streamparticles=50,y_streamparticles=50,$
        view_number=1,auto_color=1,rgb_table=rgbtable,subsample_method=1,$
        view_zoom=1.2,head_size=0.75

;stop

tool = ITGETCURRENT(TOOL=oTool) 


for i = 0,19 do begin
   f_read_3d_vec_m_32,rundir+'b1all_1',i,up



;   vx = smooth(reform(up(*,1,nz/2-zsz:nz/2+zsz,0)),2)
;   vz = smooth(reform(up(*,1,nz/2-zsz:nz/2+zsz,2)),2)

   vx = reform(up(*,1,nz/2-zsz:nz/2+zsz,0))
   vz = reform(up(*,1,nz/2-zsz:nz/2+zsz,2))

   vx = smooth(vx,2)
   vz = smooth(vz,2)


;   vx (0,0,0) = up(0,0,0,1)

   xx = x(*)
   zz = z(nz/2-zsz:nz/2+zsz)

ivector,vx,vz,scale_isotropic=1,streamline_stepsize=5,$
        /streamlines,$
        streamline_nsteps=100,x_streamparticles=50,y_streamparticles=50,$
        view_number=1,auto_color=1,rgb_table=rgbtable,subsample_method=1,$
        view_zoom=1.2,head_size=1.0,dimensions= [1200,600]

;   ivector,vx,vz,xx,zz,scale_isotropic=1,streamline_stepsize=1.0,$
;           /streamlines,$
;           streamline_nsteps=20,x_streamparticles=50,y_streamparticles=50,$
;           view_number=1,auto_color=1,rgb_table=rgbtable,subsample_method=1,$
;           view_zoom=1.2,head_size=1.0

   idTool=itgetcurrent(tool=oTool)
                                ; Export the window to a JPEG file.
                                ; First disable the Export wizard dialog.
   void = oTool->DoSetProperty('Operations/File/Export', 'SHOW_EXECUTION_UI', 0)
                                ; Export the entire window.
   void = oTool->DoSetProperty('Operations/File/Export', 'SOURCE', 1)
   void = oTool->DoSetProperty('Operations/File/Export', 'FILENAME', $
                               'iplot'+string(i,format='(I3)')+'.jpg')
   void = oTool->DoAction('Operations/File/Export')
endfor



end
 
