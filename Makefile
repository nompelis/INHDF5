HDF_INCLUDE = -I $(HDF5_HOME)/include
HDF_LIB = -L $(HDF5_HOME)/lib

MPICC = mpicc
 COPTS = -Wall -O0
#COPTS += -D  _DEBUG_

MPIF90 = mpif90
 FOPTS = -Wall -ffixed-line-length-132 -O0 -fbounds-check -fcheck=all
#FOPTS += -D  _DEBUG_


all:
	$(MPIF90) $(FOPTS) $(HDF_INCLUDE) -c inUtils_HDF5_fortran.F
	$(MPIF90) $(COPTS) $(HDF_INCLUDE) -c inUtils_HDF5.c

exam1:
	$(MPIF90) $(FOPTS) $(HDF_INCLUDE) main.f  inUtils_HDF5_fortran.o \
           $(HDF_LIB) -lhdf5 -lhdf5_fortran

exam2:
	$(MPIF90) $(COPTS) $(HDF_INCLUDE) main.c inUtils_HDF5.o \
           $(HDF_LIB) -lhdf5

clean:
	rm -f *.o a.out

