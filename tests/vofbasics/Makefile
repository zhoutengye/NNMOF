EXEC = VOF
FC = mpif90
DEBUG = true
TEST_NUM = 2
DIMENSION = TWO_DIMENSIONAL
PRECISION = DOUBLE_PRECISION
HDF5_DIR = $(HOME)/hdf5/gnu
# HDF5_DIR = $(HOME)/HPC/gnu/hdf5

# Choose VOF METHOD
# NORM = NORM_PY
# NORM = NORM_CS
NORM = NORM_MYCS

SRC_DIR = ../../src
SRCS = $(SRC_DIR)/basics.f90 $(SRC_DIR)/tools.f90 $(SRC_DIR)/vof_norm.f90 $(SRC_DIR)/vof_flood.f90 $(SRC_DIR)/vof_adv.f90  test.f90
CPP = -cpp -D$(DIMENSION) -D$(PRECISION) -DTEST$(TEST_NUM) -D$(NORM)
ifeq ($(DEBUG), true)
        FC_OPT = -O0 -g -Wall -fstack-check
#        FC_OPT = -O0 -g -traceback -check all -fp-stack-check -warn
        EXEC:= TEST-$(TEST_NUM)
else
        EXEC:= $(EXEC)_$(DIMENSION)
        FC_OPT = -O2
endif

LIBS = -L$(HDF5_DIR)/lib/libhdf5hl_fortran.a $(HDF5_DIR)/lib/libhdf5_hl.a $(HDF5_DIR)/lib/libhdf5_fortran.a $(HDF5_DIR)/lib/libhdf5.a -lz -ldl -lm
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
	rm -f *.mod.f90
	rm -f $(EXEC)
