OSNF_DIR = osnf

.PHONY: osnf_code cleanall
CLEANDIRS = $(OSNF_DIR) ./


DEBUG = -g
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
VAR_TYPE = 1 # 0 single, 1 double



main.exe	:  sce_micro_lib.a  osnf_code main.$(OBJ) sce_module.$(OBJ) 
	$(FOR2) $(FFLAGS2)main.exe main.$(OBJ) sce_module.$(OBJ) \
	        -lm sce_micro_lib.a $(OSNF_DIR)/osnf_lib.a \
		${NETCDFLIB} -I ${NETCDFMOD} ${NETCDF_LIB} $(DEBUG)
sce_micro_lib.a	:  osnf_code 
	cp $(OSNF_DIR)/osnf_lib.a sce_micro_lib.a 
sce_module.$(OBJ)	: sce_module.f90
	$(FOR) sce_module.f90 -I ${NETCDFMOD} -I${OSNF_DIR}\
	     $(FFLAGS)sce_module.$(OBJ)
main.$(OBJ)   : main.f90 sce_module.$(OBJ) 
	$(FOR)  main.f90 -I ${NETCDFMOD}  -I${OSNF_DIR} $(FFLAGS)main.$(OBJ) 
	
osnf_code:
	$(MAKE) -C $(OSNF_DIR)
clean :
	rm *.exe *.o *.mod *~ \
	*.a;rm -R *.dSYM

cleanall:
	for i in $(CLEANDIRS); do \
		$(MAKE) -C $$i clean; \
	done


