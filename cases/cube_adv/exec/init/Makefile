MOF_DIR    = ../../../../
compiler   = intel
EXEC       = Zalesak-init
FC         = mpif90
DEBUG      = true
PRECISION  = DOUBLE_PRECISION

# Uncomment to customize the path of HDF5 and hypre path
HDF5_DIR   = $(HOME)/HPC/hdf5/$(compiler)
HYPRE_DIR  = $(HOME)/HPC/hypre/$(compiler)

CASE_SRCS = init.f90 main.f90

#-----------Include essential files---------------
include $(MOF_DIR)/Make-Essential


#-----------Make the executable---------------
$(EXEC): $(OBJS)
	$(FC)  $(FC_OPT) $(INCS) -o $(EXEC)  $(OBJS) $(LIBS)

#-----------Clean redundent files---------------
clean:
	rm -f *.o
	rm -f $(SRC_DIR)/*.o
	rm -f *.mod
	rm -f *mod.f90

clobber:
	rm -f *.o
	rm -f $(SRC_DIR)/*.o
	rm -f *.mod
	rm -f *mod.f90
	rm -f $(EXEC)


#--------------Print the value of Make variable-------
print-%  : ; @echo $* = $($*)

#--------------Check the compiler and mpi version-----
check-env:
	mpif90 --version
	mpirun --version

