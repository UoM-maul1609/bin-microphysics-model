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



main.exe	:  osnf_lib.a  main.$(OBJ) 
	$(FOR2) $(FFLAGS2)main.exe main.$(OBJ) -lm osnf_lib.a \
		${NETCDFLIB} -I ${NETCDFMOD} ${NETCDF_LIB} $(DEBUG)
osnf_lib.a	:  numerics.$(OBJ) zeroin.$(OBJ) sfmin.$(OBJ) \
            fmin.$(OBJ) r1mach.$(OBJ) d1mach.$(OBJ) \
             dfsid1.$(OBJ) \
            poly_int.$(OBJ) find_pos.$(OBJ) numerics_type.$(OBJ) \
            svode.$(OBJ) slinpk.$(OBJ) vode.$(OBJ) dlinpk.$(OBJ) \
            vode_integrate.$(OBJ) erfinv.$(OBJ) tridiagonal.$(OBJ) \
            hygfx.$(OBJ) random.$(OBJ)
	$(AR) rc osnf_lib.a numerics.$(OBJ) zeroin.$(OBJ) sfmin.$(OBJ) \
	    fmin.$(OBJ) r1mach.$(OBJ) \
	    d1mach.$(OBJ) dfsid1.$(OBJ) poly_int.$(OBJ) find_pos.$(OBJ) numerics_type.$(OBJ) \
            svode.$(OBJ) slinpk.$(OBJ) vode.$(OBJ) dlinpk.$(OBJ) vode_integrate.$(OBJ) \
            erfinv.$(OBJ) tridiagonal.$(OBJ) hygfx.$(OBJ) random.$(OBJ)
numerics_type.$(OBJ)	: numerics_type.f90 
	$(FOR) numerics_type.f90 $(FFLAGS)numerics_type.$(OBJ)
numerics.$(OBJ)	: numerics.f90 numerics_type.$(OBJ)
	$(FOR) numerics.f90 $(FFLAGS)numerics.$(OBJ)
zeroin.$(OBJ)	: zeroin.f
	$(FOR) zeroin.f $(FFLAGS)zeroin.$(OBJ)
sfmin.$(OBJ)	: sfmin.f
	$(FOR) sfmin.f $(FFLAGS)sfmin.$(OBJ)
fmin.$(OBJ)	: fmin.f
	$(FOR) fmin.f $(FFLAGS)fmin.$(OBJ)
r1mach.$(OBJ) 	: r1mach.f 
	$(FOR) r1mach.f $(FFLAGS)r1mach.$(OBJ)  
d1mach.$(OBJ) 	: d1mach.f 
	$(FOR) d1mach.f $(FFLAGS)d1mach.$(OBJ)  
svode.$(OBJ) 	: svode.f
	$(FOR) svode.f -cpp -DDUAL $(FFLAGS)svode.$(OBJ)  
slinpk.$(OBJ) 	: slinpk.f 
	$(FOR) slinpk.f -w $(FFLAGS)slinpk.$(OBJ)  
vode.$(OBJ) 	: vode.f 
	$(FOR) vode.f $(FFLAGS)vode.$(OBJ)  
dlinpk.$(OBJ) 	: dlinpk.f 
	$(FOR) dlinpk.f -w $(FFLAGS)dlinpk.$(OBJ)  
vode_integrate.$(OBJ) 	: vode_integrate.f90
	$(FOR) vode_integrate.f90 $(FFLAGS)vode_integrate.$(OBJ)  
dfsid1.$(OBJ) 	: dfsid1.f 
	$(FOR) dfsid1.f $(FFLAGS)dfsid1.$(OBJ)  
find_pos.$(OBJ)	: find_pos.f90 numerics_type.$(OBJ)
	$(FOR) find_pos.f90 $(FFLAGS)find_pos.$(OBJ)
poly_int.$(OBJ)	: poly_int.f90 numerics_type.$(OBJ)
	$(FOR) poly_int.f90 $(FFLAGS)poly_int.$(OBJ)
tridiagonal.$(OBJ)	: tridiagonal.f90 numerics_type.$(OBJ)
	$(FOR) tridiagonal.f90 $(FFLAGS)tridiagonal.$(OBJ)
erfinv.$(OBJ)	: erfinv.f90 numerics_type.$(OBJ)
	$(FOR) erfinv.f90 $(FFLAGS)erfinv.$(OBJ)
random.$(OBJ)	: random.f90 numerics_type.$(OBJ)
	$(FOR) random.f90 $(FFLAGS)random.$(OBJ)
hygfx.$(OBJ) : hygfx.for 
	$(FOR) hygfx.for $(FFLAGS)hygfx.$(OBJ) 
main.$(OBJ)   : main.f90 
	$(FOR)  main.f90 -I ${NETCDFMOD} $(FFLAGS)main.$(OBJ) 
	
clean :
	rm *.exe *.o *.mod *~ \
	osnf_lib.a;rm -R *.dSYM

