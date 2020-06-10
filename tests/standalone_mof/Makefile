EXEC = VOF
FC = mpif90
DEBUG = true
TEST_NUM = 1
HDF5_DIR = $(HOME)/hdf5/gnu
# HDF5_DIR = $(HOME)/HPC/gnu/hdf5


SRCS = standalone_mof.f90 main.f90
CPP = -cpp -DTEST$(TEST_NUM)
ifeq ($(DEBUG), true)
        FC_OPT = -O0 -g -Wall -Wextra  -Wrealloc-lhs-all -fstack-check  -fbounds-check -fbacktrace -Wargument-mismatch
        FC_OPT = -O3
        EXEC:= TEST-$(TEST_NUM)
else
        EXEC:= $(EXEC)_$(DIMENSION)
        FC_OPT = -O3
endif


# rule for *.f90 to *.o
.SUFFIXES: .o .f90
OBJS = $(SRCS:.f90=.o)
OBJS_TESTS = $(SRCS_TESTS:.f90=.o)
.f90.o:
	$(FC) -c $(CPP) $(FC_OPT) $(FFLAGS)  $*.f90

# after all .o is gathered, link and generate the executable
$(EXEC): $(OBJS)
	$(FC)  $(FC_OPT) -o $(EXEC)  *.o

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *mod.f90
	rm -f $(EXEC)
