program = openmc

src = main.f90 types.f90 global.f90 fileio.f90 output.f90 \
	string.f90 geometry.f90 mcnp_random.f90 source.f90
objects = $(src:.f90=.o)

#--------------------------------------------------------------------
# Compiler Options
F90 = gfortran
F90FLAGS = -g

#--------------------------------------------------------------------
# Targets
all:	$(program)
$(program): $(objects)
	$(F90) $(F90FLAGS) $(objects) -o $@
clean:
	@rm -f *.o *.mod $(program)
neat:
	@rm -f *.o *.mod

#--------------------------------------------------------------------
# Rules & misc
.SUFFIXES: .f90 .o
.PHONY: all clean neat

%.o: %.f90
	$(F90) -c $(F90FLAGS) $(MPI_COMPILE_FLAGS) $<

#--------------------------------------------------------------------
# Dependencies
types.o:
global.o: types.o
string.o: global.o output.o
output.o: global.o
source.o: global.o mcnp_random.o
geometry.o: types.o global.o output.o string.o
fileio.o: types.o global.o string.o output.o
main.o: global.o fileio.o output.o geometry.o mcnp_random.o \
		source.o
