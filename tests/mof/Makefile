compiler   = gnu
EXEC       = MOF
FC         = mpif90
DEBUG      = false
VERBOSE    = false
TEST_NUM   = 1
PRECISION  = DOUBLE_PRECISION
HDF5_DIR   = $(HOME)/hdf5/$(compiler)
# HDF5_DIR = $(HOME)/HPC/$(compiler)/hdf5

SRC_DIR = ../../src
LEMOINE_DIR = ../../src/MOFLemoine/$(compiler)
SUSSMAN_DIR = ../../src/MOFSussman/$(compiler)
SRCS = $(SRC_DIR)/basics.f90 $(SRC_DIR)/tools.f90  $(SRC_DIR)/vof_func.f90 $(SRC_DIR)/vof_adv.f90 test.f90
CPP = -cpp -D$(PRECISION) -DTEST$(TEST_NUM)
ifeq ($(DEBUG), true)
        FC_OPT = -O0 -g -Wall -Wextra  -Wrealloc-lhs-all -fstack-check  -fbounds-check -fbacktrace -Wargument-mismatch
#        FC_OPT = -O3
	CPP = -cpp -D$(PRECISION) -DTEST$(TEST_NUM) -DDEBUG
        EXEC:= TEST-$(TEST_NUM)
else
        FC_OPT = -O3
endif

ifeq ($(VERBOSE), true)
	CPP := $(CPP) -DVERBOSE
endif

MOFLIB = $(LEMOINE_DIR)/lemoinemof.a $(SUSSMAN_DIR)/libmof.a

LIBS =  -L$(HDF5_DIR)/lib/libhdf5hl_fortran.a $(HDF5_DIR)/lib/libhdf5_hl.a $(HDF5_DIR)/lib/libhdf5_fortran.a $(HDF5_DIR)/lib/libhdf5.a  -lz -ldl -lm  
INCS = -I$(HDF5_DIR)/include -I$(SRC_DIR)  -I$(LEMOINE_DIR) -I$(SUSSMAN_DIR)
MOFS = -J$(SRC_DIR)


# rule for *.f90 to *.o
.SUFFIXES: .o .f90
OBJS = $(SRCS:.f90=.o)
OBJS_TESTS = $(SRCS_TESTS:.f90=.o)
.f90.o:
	$(FC) -c $(CPP) $(FC_OPT) $(FFLAGS) $(MODS) $(INCS) $*.f90

# after all .o is gathered, link and generate the executable
$(EXEC): $(OBJS)
	$(FC)  $(FC_OPT) $(MODS) -o $(EXEC) $(INCS) *.o $(MOFLIB) $(LIBS) 

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *mod.f90
	rm -f $(EXEC)
