OSNF_DIR = osnf
SCE_DIR = sce
SCE_OSNF_DIR = sce/osnf

.PHONY: osnf_code sce_code cleanall
CLEANDIRS = $(OSNF_DIR) $(SCE_DIR)  $(SCE_OSNF_DIR) ./


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



main.exe	:  b_micro_lib.a  sce_code \
        osnf_code main.$(OBJ) bin_microphysics_module.$(OBJ) 
	$(FOR2) $(FFLAGS2)main.exe main.$(OBJ) bin_microphysics_module.$(OBJ) \
	        -lm b_micro_lib.a $(OSNF_DIR)/osnf_lib.a $(SCE_DIR)/sce_micro_lib.a \
	        $(SCE_DIR)/sce_module.$(OBJ) \
		${NETCDFLIB} -I ${NETCDFMOD} ${NETCDF_LIB} $(DEBUG)
b_micro_lib.a	:  osnf_code 
	cp $(OSNF_DIR)/osnf_lib.a b_micro_lib.a 
bin_microphysics_module.$(OBJ)	: bin_microphysics_module.f90
	$(FOR) bin_microphysics_module.f90 -I ${NETCDFMOD} -I${OSNF_DIR} -I${SCE_DIR}\
	     $(FFLAGS)bin_microphysics_module.$(OBJ)
main.$(OBJ)   : main.f90 bin_microphysics_module.$(OBJ) $(SCE_DIR)/sce_module.$(OBJ) 
	$(FOR)  main.f90 -I ${NETCDFMOD}  -I${OSNF_DIR} -I${SCE_DIR} $(FFLAGS)main.$(OBJ) 
	
osnf_code:
	$(MAKE) -C $(OSNF_DIR)
sce_code:
	$(MAKE) -C $(SCE_DIR)
clean :
	rm *.exe *.o *.mod *~ \
	*.a;rm -R *.dSYM

cleanall:
	for i in $(CLEANDIRS); do \
		$(MAKE) -C $$i clean; \
	done


