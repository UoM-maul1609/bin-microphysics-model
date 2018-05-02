DEBUG = -fbounds-check -g
MPI    =#-DMPI1
OPT    =-O3

# these three lines should be edited for your system. On systems 
# that do not have separate fortran and c libraries, set NETCDF_FOR and NETCDF_C
# to the same, and set NETCDF_LIB to -lnetcdf (i.e. without the extra f)
#NETCDF_FOR=/Users/mccikpc2/Dropbox/programming/netcdf-4.4.4-mac/
#NETCDF_C=/Users/mccikpc2/Dropbox/programming/netcdf-4.4.1.1-mac/
NETCDF_LIB=-lnetcdff 

NETCDFLIB=-L ${NETCDF_FOR}/lib/  \
          -L ${NETCDF_C}/lib/
NETCDFMOD= ${NETCDF_FOR}/include/


FOR = gfortran -c  
FOR2 = gfortran  

AR = ar 
RANLIB = ranlib 
OBJ = o
FFLAGS = $(OPT)  $(DEBUG)  -o 
FFLAGS2 =  $(DEBUG) -O3 -o 



main.exe	:  b_micro_lib.a  main.$(OBJ) bin_microphysics_module.$(OBJ) 
	$(FOR2) $(FFLAGS2)main.exe main.$(OBJ) bin_microphysics_module.$(OBJ) -lm b_micro_lib.a \
		${NETCDFLIB} -I ${NETCDFMOD} ${NETCDF_LIB} $(DEBUG)
b_micro_lib.a	:   nrtype.$(OBJ) nr.$(OBJ) nrutil.$(OBJ) zbrent.$(OBJ) trapzd.$(OBJ) \
                    qsimp.$(OBJ) polint.$(OBJ) locate.$(OBJ) qromb.$(OBJ) brent.$(OBJ) \
                    midpnt.$(OBJ) dfridr.$(OBJ) fmin.$(OBJ) broydn.$(OBJ) fdjac.$(OBJ) \
                    qrdcmp.$(OBJ) qrupdt.$(OBJ) rsolv.$(OBJ) lnsrch.$(OBJ) rotate.$(OBJ) \
                    pythag.$(OBJ) d1mach.$(OBJ) vode.$(OBJ) dlinpk.$(OBJ) acdc.$(OBJ) \
                    colmod.$(OBJ) opkdmain.$(OBJ) opkda1.$(OBJ) opkda2.$(OBJ) \
                    rkqs.$(OBJ) rkck.$(OBJ) odeint.$(OBJ)
	$(AR) rc b_micro_lib.a nrtype.$(OBJ) nr.$(OBJ) nrutil.$(OBJ) zbrent.$(OBJ) \
	                trapzd.$(OBJ) qsimp.$(OBJ) polint.$(OBJ) locate.$(OBJ) qromb.$(OBJ) \
	                brent.$(OBJ) midpnt.$(OBJ) dfridr.$(OBJ) fmin.$(OBJ) broydn.$(OBJ) \
	                fdjac.$(OBJ) qrdcmp.$(OBJ) qrupdt.$(OBJ) rsolv.$(OBJ) lnsrch.$(OBJ) \
	                rotate.$(OBJ) pythag.$(OBJ) d1mach.$(OBJ) vode.$(OBJ) dlinpk.$(OBJ) \
	                acdc.$(OBJ) colmod.$(OBJ) opkdmain.$(OBJ) opkda1.$(OBJ) opkda2.$(OBJ) \
                    rkqs.$(OBJ) rkck.$(OBJ) odeint.$(OBJ)
rkqs.$(OBJ)	: nr/rkqs.f90
	$(FOR) nr/rkqs.f90 $(FFLAGS)rkqs.$(OBJ)	
rkck.$(OBJ)	: nr/rkck.f90
	$(FOR) nr/rkck.f90 $(FFLAGS)rkck.$(OBJ)	
odeint.$(OBJ)	: nr/odeint.f90
	$(FOR) nr/odeint.f90 $(FFLAGS)odeint.$(OBJ)	
zbrent.$(OBJ)	: nr/zbrent.f90
	$(FOR) nr/zbrent.f90 $(FFLAGS)zbrent.$(OBJ)
trapzd.$(OBJ)	: nr/trapzd.f90
	$(FOR) nr/trapzd.f90 $(FFLAGS)trapzd.$(OBJ)
qsimp.$(OBJ)	: nr/qsimp.f90
	$(FOR) nr/qsimp.f90 $(FFLAGS)qsimp.$(OBJ)
polint.$(OBJ)	: nr/polint.f90
	$(FOR) nr/polint.f90 $(FFLAGS)polint.$(OBJ)
locate.$(OBJ)	: nr/locate.f90
	$(FOR) nr/locate.f90 $(FFLAGS)locate.$(OBJ)
qromb.$(OBJ)	: nr/qromb.f90
	$(FOR) nr/qromb.f90 $(FFLAGS)qromb.$(OBJ)
brent.$(OBJ)	: nr/brent.f90
	$(FOR) nr/brent.f90 $(FFLAGS)brent.$(OBJ)
midpnt.$(OBJ)	: nr/midpnt.f90
	$(FOR) nr/midpnt.f90 $(FFLAGS)midpnt.$(OBJ)
fmin.$(OBJ)	: nr/fmin.f90
	$(FOR) nr/fmin.f90 $(FFLAGS)fmin.$(OBJ)
broydn.$(OBJ)	: nr/broydn.f90
	$(FOR) nr/broydn.f90 $(FFLAGS)broydn.$(OBJ)
fdjac.$(OBJ)	: nr/fdjac.f90
	$(FOR) nr/fdjac.f90 $(FFLAGS)fdjac.$(OBJ)
qrdcmp.$(OBJ)	: nr/qrdcmp.f90
	$(FOR) nr/qrdcmp.f90 $(FFLAGS)qrdcmp.$(OBJ)
qrupdt.$(OBJ)	: nr/qrupdt.f90
	$(FOR) nr/qrupdt.f90 $(FFLAGS)qrupdt.$(OBJ)
rsolv.$(OBJ)	: nr/rsolv.f90
	$(FOR) nr/rsolv.f90 $(FFLAGS)rsolv.$(OBJ)
lnsrch.$(OBJ)	: nr/lnsrch.f90
	$(FOR) nr/lnsrch.f90 $(FFLAGS)lnsrch.$(OBJ)
rotate.$(OBJ)	: nr/rotate.f90
	$(FOR) nr/rotate.f90 $(FFLAGS)rotate.$(OBJ)
pythag.$(OBJ)	: nr/pythag.f90
	$(FOR) nr/pythag.f90 $(FFLAGS)pythag.$(OBJ)
dfridr.$(OBJ)	: nr/dfridr.f90
	$(FOR) nr/dfridr.f90 $(FFLAGS)dfridr.$(OBJ)
nrtype.$(OBJ)	: nr/nrtype.f90
	$(FOR) nr/nrtype.f90 $(FFLAGS)nrtype.$(OBJ)
nr.$(OBJ)	: nr/nr.f90 
	$(FOR) nr/nr.f90 $(FFLAGS)nr.$(OBJ)
nrutil.$(OBJ)	: nr/nrutil.f90
	$(FOR) nr/nrutil.f90 $(FFLAGS)nrutil.$(OBJ)
d1mach.$(OBJ) 	: netlib/d1mach.f 
	$(FOR) netlib/d1mach.f $(FFLAGS)d1mach.$(OBJ)  
opkdmain.$(OBJ) : netlib/opkdmain.f 
	$(FOR) netlib/opkdmain.f $(FFLAGS)opkdmain.$(OBJ)
opkda1.$(OBJ) : netlib/opkda1.f 
	$(FOR) netlib/opkda1.f -w $(FFLAGS)opkda1.$(OBJ)
opkda2.$(OBJ) : netlib/opkda2.f 
	$(FOR) netlib/opkda2.f -w $(FFLAGS)opkda2.$(OBJ)
vode.$(OBJ) : netlib/vode.f 
	$(FOR) 	netlib/vode.f $(FFLAGS)vode.$(OBJ)
acdc.$(OBJ) : netlib/acdc.f 
	$(FOR) 	netlib/acdc.f $(FFLAGS)acdc.$(OBJ)
dlinpk.$(OBJ) : netlib/dlinpk.f 
	$(FOR) 	netlib/dlinpk.f -w $(FFLAGS)dlinpk.$(OBJ)
colmod.$(OBJ) : netlib/colmod.f 
	$(FOR) 	netlib/colmod.f -Wno-all $(FFLAGS)colmod.$(OBJ)
hygfx.$(OBJ) : fun/hygfx.for 
	$(FOR) fun/hygfx.for $(FFLAGS)hygfx.$(OBJ) 
bin_microphysics_module.$(OBJ)	: bin_microphysics_module.f90
	$(FOR) bin_microphysics_module.f90 -I ${NETCDFMOD} \
	     $(FFLAGS)bin_microphysics_module.$(OBJ)
main.$(OBJ)   : main.f90 bin_microphysics_module.$(OBJ) 
	$(FOR)  main.f90 -I ${NETCDFMOD} $(FFLAGS)main.$(OBJ) 


clean :
	rm *.exe  *.o *.mod *~ \
	b_micro_lib.a

