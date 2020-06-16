HDF_INCLUDE = -I $(HDF5_HOME)/include
HDF_LIB = -L $(HDF5_HOME)/lib

MPICC = mpicc
 COPTS = -Wall -O0

MPIF90 = mpif90
 FOPTS = -Wall -ffixed-line-length-132 -O0 -fbounds-check -fcheck=all
#FOPTS += -D  _DEBUG_


all:
	$(MPIF90) $(FOPTS) $(HDF_INCLUDE) -c inUtils_HDF5.F

exam:
	$(MPIF90) $(FOPTS) $(HDF_INCLUDE) main.f  inUtils_HDF5.o \
           $(HDF_LIB) -lhdf5 -lhdf5_fortran

clean:
	rm -f *.o a.out

