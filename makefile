#F77 = mpif90 -i4 -r4 -O2 -byteswapio
F77 = mpif90 -i4 -real-size 32 -O4

FILES = maind.f gutsf.f gutsp.f misc.f boundary.f fluid_init.f part_init.f 
INCLUDE = incurv.h para.h
OBJECTS = maind.o gutsf.o gutsp.o misc.o boundary.o part_init.o initial.o
OBJECTS_F = maind_f.o gutsf.o misc.o boundary.o fluid_init.o part_init.o initial.o gutsp.o mom.o
OBJECTS_P = maind_p.o gutsp.o misc.o part_init.o gutsf.o boundary.o initial.o

hybrid:	$(OBJECTS) 
	$(F77) -o hybrid_2dkh $(OBJECTS) 

fluid:	$(OBJECTS_F)
	$(F77) -o fluid $(OBJECTS_F)

part:	$(OBJECTS_P)
	$(F77) -o part $(OBJECTS_P)

clean:
	rm *.o hybrid fluid part

maind.o:maind.f $(INCLUDE);$(F77) -c maind.f
maind_f.o:maind_f.f $(INCLUDE);$(F77) -c maind_f.f
maind_p.o:maind_p.f $(INCLUDE);$(F77) -c maind_p.f
gutsf.o:gutsf.f $(INCLUDE);$(F77) -c gutsf.f
gutsp.o:gutsp.f $(INCLUDE);$(F77) -c gutsp.f
misc.o:misc.f $(INCLUDE);$(F77) -c misc.f
boundary.o:boundary.f $(INCLUDE);$(F77) -c boundary.f
fluid_init.o:fluid_init.f $(INCLUDE);$(F77) -c fluid_init.f
part_init.o:part_init.f $(INCLUDE);$(F77) -c part_init.f
initial.o:initial.f $(INCLUDE);$(F77) -c initial.f
mom.o:mom.f $(INCLUDE);$(F77) -c mom.f


