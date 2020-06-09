FC = gfortran
DIM = 3

OPTFLAG = -g -O3 -fimplicit-none
CPPFLAG = -DBL_LANG_FORT -DBL_SPACEDIM=$(DIM) -DBL_FORT_USE_UNDERSCORE -DBL_LINUX


.SUFFIXES: .o .F90 
SRCS = decision_tree.F90 PROBCOMMON.F90 tetrahedron_keast_rule.F90 triangle_fekete_rule.F90 GLOBALUTIL.F90 MOF.F90 MOFSussman.F90 
OBJS = $(SRCS:.F90=.o)

# after all .o is gathered, link and generate the executable
libmof: $(OBJS)
	ar rcs libmof.a $(OBJS)

.F90.o:
	$(FC) -c -cpp $(CPPFLAG) $(OPTFLAG) -I. $*.F90

clean:
	rm -f *.o
	rm -f *.mod
	rm -f $(EXEC)
