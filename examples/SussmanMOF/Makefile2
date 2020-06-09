EXEC = MOF.ex
FC = gfortran
FCOPT = -O3
DIM = 3

OPTFLAG = -g -O3 -fimplicit-none
CPPFLAG = -DBL_LANG_FORT -DBL_SPACEDIM=$(DIM) -DBL_FORT_USE_UNDERSCORE -DBL_LINUX


.SUFFIXES: .o .F90 
SRCS = decision_tree.F90 PROBCOMMON.F90 tetrahedron_keast_rule.F90 triangle_fekete_rule.F90 GLOBALUTIL.F90 MOF.F90 MOF_ML.F90 
OBJS = $(SRCS:.F90=.o)

# after all .o is gathered, link and generate the executable
$(EXEC): libmof
	$(FC) -c -cpp $(OPTFLAG) main_MOF_standalone.F90
	$(FC) $(OPTFLAG) -o $(EXEC)  main_MOF_standalone.o -L. -lmof 

libmof: $(OBJS)
	ar rcs libmof.a $(OBJS)

$dt: $(OBJS)

.F90.o:
	$(FC) -c -cpp $(CPPFLAG) $(OPTFLAG) -I. $*.F90

clean:
	rm -f *.o
	rm -f *.mod
	rm -f $(EXEC)
