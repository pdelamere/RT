;f_read_coord,rundir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz

openr,1,'dc.dat',/f77_unformatted
;readu,1,m
readu,1,dc

dcarr = dc

while not(eof(1)) do begin
;   readu,1,mll
   readu,1,dc
   print,dc
   dcarr = [dcarr,dc]
endwhile

plot,dcarr*1e6

end
