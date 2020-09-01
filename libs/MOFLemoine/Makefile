INSTALL_DIR = intel
MOFLIB = lemoinemof.a
FC         = ifort
PRECISION  = DOUBLE_PRECISION

SRCS_GEOMETRY = transformations.f90 points.f90 polyhedron.f90 complete_polyhedron_structure.f90  flood_polyhedron.f90
SRCS_MOF = variables.f90 mof3d_analytic_centroid.f90 mof3d_gradient.f90 mof3d_bfgs.f90
SRCS = $(SRCS_GEOMETRY) $(SRCS_MOF)
FC_OPT = -O3

# rule for *.f90 to *.o
.SUFFIXES: .o .f90
OBJS = $(SRCS:.f90=.o)
.f90.o:
	$(FC) -c $(FC_OPT)  $(INCS) $*.f90

# after all .o is gathered, link and generate the executable
$(MOFLIB): $(OBJS)
	mkdir -p $(INSTALL_DIR)
	ar -crv $(MOFLIB) $(OBJS)
	mv *.a $(INSTALL_DIR)
	mv *.mod $(INSTALL_DIR)
	make clean

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *mod.f90
	rm -f $(MOFLIB)
