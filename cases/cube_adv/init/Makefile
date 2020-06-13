EXEC       = VOF
FC         = mpif90
DEBUG      = true
VERBOSE    = false
TEST_NUM   = 2
DIMENSION  = TWO_DIMENSIONAL
PRECISION  = DOUBLE_PRECISION
HDF5_DIR   = $(HOME)/hdf5/gnu
# HDF5_DIR = $(HOME)/HPC/gnu/hdf5

SRC_DIR = ../../../src
SRCS = $(SRC_DIR)/basics.f90 $(SRC_DIR)/tools.f90 init.f90 main.f90
CPP = -cpp -D$(DIMENSION) -D$(PRECISION) -DTEST$(TEST_NUM)
ifeq ($(DEBUG), true)
        FC_OPT = -O0 -g -Wall -Wextra  -Wrealloc-lhs-all -fstack-check  -fbounds-check -fbacktrace -Wargument-mismatch
#        FC_OPT = -O3
        EXEC:= TEST-$(TEST_NUM)
else
        EXEC:= $(EXEC)_$(DIMENSION)
        FC_OPT = -O3
endif

ifeq ($(VERBOSE), true)
	CPP := $(CPP) -DVERBOSE
endif

LIBS = -L$(HDF5_DIR)/lib/libhdf5hl_fortran.a $(HDF5_DIR)/lib/libhdf5_hl.a $(HDF5_DIR)/lib/libhdf5_fortran.a $(HDF5_DIR)/lib/libhdf5.a  -lz -ldl -lm 
INCS = -I$(HDF5_DIR)/include -I$(SRC_DIR) 

# rule for *.f90 to *.o
.SUFFIXES: .o .f90
OBJS = $(SRCS:.f90=.o)
OBJS_TESTS = $(SRCS_TESTS:.f90=.o)
.f90.o:
	$(FC) -c $(CPP) $(FC_OPT) $(FFLAGS) $(INCS) $*.f90

# after all .o is gathered, link and generate the executable
$(EXEC): $(OBJS)
	$(FC)  $(FC_OPT) $(INCS) -o $(EXEC)  *.o $(LIBS)

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *mod.f90
	rm -f $(EXEC)
