#--------- paris-ParisSimulator main Makefile --------------------------

# Change the following line to force openmpi to use another compiler
OMPI_FC = gfortran-mp-8

FC = mpifort
CC = mpicxx  # could be mpic++
CUDAC = nvcc

LDFLAGS = -M

FLAGS = -O3 -cpp -fimplicit-none -ffpe-summary=invalid,zero,overflow # -fbounds-check # -DOLD_BDRY_COND -fimplicit-none -fbounds-check # -g 
#FLAGS = -O0 -cpp  -fimplicit-none -fimplicit-none -fbounds-check -g 
FLAGSDEBUG = -cpp -g -gstabs -O2  -fimplicit-none # -ffpe-summary=invalid,zero,overflow
HARDFLAGS = -fbounds-check -ffpe-trap=invalid,zero,overflow -Wall # -ffpe-trap=invalid,zero,overflow,underflow,precision,denormal 

CFLAGS = -O -Wall # -g -gstabs
BINDIR = $(HOME)/bin

# select option for hypre
# babbage
# HYPRE_DIR = /share/apps/hypre
# Local
SILO_DIR = $(HOME)/cfd/libs/silo-4.10.2

HYPRE_DIR = $(HOME)/cfd/libs/hypre-2.10.0b/src/lib
HYPRE_LIBS =  -L$(HYPRE_DIR) -lHYPRE 
INC_SILO = -I$(SILO_DIR)/include
#DIR_MPI = /usr/local/include
#INC_MPI = -I$(DIR_MPI)
INC = $(INC_SILO) # $(INC_MPI)


ifdef HAVE_VOFI
FLAGSA = -DHAVE_VOFI $(FLAGS) 
VOFI_DIR =  $(HOME)/lib
VOFI_LIBS = -L$(VOFI_DIR) -lvofi
else
FLAGSA = $(FLAGS)
endif

ifdef PHASE_CHANGE
FLAGSB = -DPHASE_CHANGE $(FLAGSA)
else
FLAGSB = $(FLAGSA)
endif

ifdef HAVE_SILO
FLAGSC = -DHAVE_SILO $(FLAGSB) $(INC) 
SILO_LIB =  -L$(SILO_DIR)/lib -lsilo -lm -lstdc++
else
FLAGSC = $(FLAGSB)
endif

ifdef HAVE_CUDA
FLAGSD = -DHAVE_CUDA $(FLAGSC) $(INC) -L/usr/local/cuda/lib64 -I/usr/local/cuda/include  -lcudart -lcuda -lstdc++
CFLAGS = -O -Wall # -g -gstabs
CUDAFLAGS = -ccbin $(CC) -L/usr/local/cuda/lib -I/usr/local/cuda/include  -lcudart -lcuda  -lstdc++
CUDAOBJ = module_CUDA.o poissonCUDA.o cudaFun.o
INSTLIBS = -Xlinker  -L/usr/local/cuda/lib64 -I/usr/local/cuda/include  -lcudart -lcuda -lstdc++
else
FLAGSD = $(FLAGSC)
endif

FFLAGS = $(FLAGSD)
LIBS = $(HYPRE_LIBS) $(VOFI_LIBS) $(SILO_LIB)

#------------------------No changes needed beyond this line----------------------------------------------

OBJ = paris.o solids.o modules.o vofmodules.o front.o surface_tension.o lppmodules.o st_testing.o newsolver.o MGsolver.o $(CUDAOBJ) freesurface.o boiling.o vofnonmodule.o vof_functions.o rudman_modules.o rudman_procedures.o rudman_advection.o rudman_bc.o 

SRC = $(wildcard  *.f90) 

install: $(OBJ)
	@echo compiling with FLAGS = $(FFLAGS)  linking with FOPTS = $(FOPTS)
	@echo compiler is FC = $(FC), mpi override is OMPI_FC = $(OMPI_FC)
	$(FC) -o paris $(FOPTS) $(OBJ) $(LIBS) $(INSTLIBS)
	@if [ ! -d $(BINDIR) ] ; then echo "directory bin does not exist creating it" ; mkdir $(BINDIR) ; fi 
	mv paris $(BINDIR)/paris
	@find .  -name "*.sh" -exec chmod +x  {} \; 

all: tags install pariscompare parisdeconv
	@echo making tags executable utility programs.

clean:
	@echo Cleaning all output files, .o files, tmp files, report.
	@rm -fR *.o *.mod paris stats *~ track out* errftc tmp* *.tmp fort.* *.visit core.* *stats*
	@cd Tests; sh ./clean.sh; cd ..
	@cd Documentation; sh ./cleandoc.sh; cd ..
	@find . -type l -exec /bin/rm {} \;

distclean: clean
	@echo Cleaning xml files, visit session files.
	@rm -fR  session* *.xml TAGS tags input

test:  install pariscompare parisdeconv pariscompare3D
	@echo "The test suite takes less than 4 minutes on a 4-core intel i7 MacBook Pro"
	@cd Tests; chmod +x ./runtests.sh; ./runtests.sh

hardtest:  pariscompare parisdeconv pariscompare3D
	make clean install FLAGS:="$(FLAGSDEBUG) $(HARDFLAGS)"
	@echo "The test suite takes less than 4 minutes on an 4-core MacBook Pro mid-2015"
	@cd Tests; chmod +x ./runtests.sh; ./runtests.sh

longtest:  install pariscompare parisdeconv
	@echo "This test takes 28 minutes on an 4-core MacBook Pro mid-2015"
	@cd Tests; chmod +x ./runlongtests.sh; ./runlongtests.sh

# single processor test
minitest: install
	@cd Tests/Mini; chmod +x ./run.sh; ./run.sh

tags:	$(SRC)
# @SZ Create a tags file named TAGS for use by emacs
	@etags $(SRC)
# @SZ Create a tags file named tags for use by vi or textwrangler
# @SZ On MacOS tags and TAGS are identical ! 
# @SZ	ctags paris.f90 

paris.o:  paris.f90 solids.o modules.o vofmodules.o front.o surface_tension.o lppmodules.o st_testing.o newsolver.o MGsolver.o $(CUDAOBJ) freesurface.o  boiling.o rudman_modules.o 
	$(FC) -c $(FFLAGS) $<

vofmodules.o: vofmodules.f90 modules.o rudman_modules.o
	$(FC) -c $(FFLAGS) $<

lppmodules.o: lppmodules.f90 vofmodules.o modules.o
	$(FC) -c $(FFLAGS) $<

surface_tension.o: surface_tension.f90 vofmodules.o modules.o newsolver.o MGsolver.o $(CUDAOBJ)
	$(FC) -c $(FFLAGS) $<

st_testing.o: st_testing.f90 vofmodules.o modules.o surface_tension.o boiling.o
	$(FC) -c $(FFLAGS) $<

solids.o:  solids.f90 modules.o
	$(FC) -c $(FFLAGS) $<

front.o:  front.f90 modules.o
	$(FC) -c $(FFLAGS) $<

newsolver.o:  newsolver.f90 $(CUDAOBJ) modules.o MGsolver.o
	$(FC) -c $(FFLAGS) $<

MGsolver.o:  MGsolver.f90 modules.o
	$(FC) -c $(FFLAGS) $<

freesurface.o: freesurface.f90 modules.o
	$(FC) -c $(FFLAGS) $<

boiling.o: boiling.f90 modules.o
	$(FC) -c $(FFLAGS) $<

averages.o: averages.f90 modules.o solids.o vofmodules.o
	$(FC) -c $(FFLAGS) $<

rudman_modules.o:  rudman_modules.f90 modules.o
	$(FC) -c $(FFLAGS) $<

%.o : %.f90
	$(FC) -c $(FFLAGS) $<

ifdef HAVE_CUDA
poissonCUDA.o: poissonCUDA.cu cudaFun.o 
	$(CUDAC) -c $(CUDAFLAGS) $<
	
module_CUDA.o: module_CUDA.f90  modules.o poissonCUDA.o
	$(FC) -c $(FFLAGS) $<
	
cudaFun.o: cudaFun.f90  modules.o 
	$(FC) -c -ffree-form $(FFLAGS) $<
endif
	
pariscompare: compare.o
	$(CC) -o pariscompare compare.o -lm
	mv pariscompare $(BINDIR)


pariscompare3D: compare_4cols.o
	$(CC) -o pariscompare3D compare_4cols.o -lm
	mv pariscompare3D $(BINDIR)

parisdeconv: parisdeconv.o
	$(CC) -o parisdeconv parisdeconv.o -lm
	mv parisdeconv $(BINDIR)

.c.o:   $< 
	$(CC) -c $(CFLAGS) $< 
