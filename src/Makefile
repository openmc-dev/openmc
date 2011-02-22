program = openmc

main = main.f90
test = unittest.f90
modules = 	ace.f90 \
		cross_section.f90 \
		data_structures.f90 \
		endf.f90 \
		energy_grid.f90 \
		fileio.f90 \
		geometry.f90 \
		global.f90 \
		mcnp_random.f90 \
		output.f90 \
		physics.f90 \
		search.f90 \
		source.f90 \
		string.f90 \
		types.f90

main_objects = $(modules:.f90=.o) $(main:.f90=.o)
test_objects = $(modules:.f90=.o) $(test:.f90=.o)

#--------------------------------------------------------------------
# Compiler Options
F90 = gfortran
F90FLAGS = -g

#--------------------------------------------------------------------
# Targets
all:	$(program)
$(program): $(main_objects)
	$(F90) $(F90FLAGS) $(main_objects) -o $@
clean:
	@rm -f *.o *.mod $(program)
neat:
	@rm -f *.o *.mod
unittest: $(test_objects)
	$(F90) $(F90FLAGS) $(test_objects) -o $@

#--------------------------------------------------------------------
# Rules & misc
.SUFFIXES: .f90 .o
.PHONY: all clean neat

%.o: %.f90
	$(F90) -c $(F90FLAGS) $(MPI_COMPILE_FLAGS) $<

#--------------------------------------------------------------------
# Dependencies
ace.o: global.o output.o string.o fileio.o string.o endf.o
cross_section.o: global.o string.o data_structures.o output.o
data_structures.o: global.o
endf.o: global.o
energy_grid.o: global.o output.o data_structures.o
fileio.o: types.o global.o string.o output.o data_structures.o
geometry.o: types.o global.o output.o string.o data_structures.o
global.o: types.o
main.o: global.o fileio.o output.o geometry.o mcnp_random.o \
		source.o physics.o cross_section.o data_structures.o \
		ace.o energy_grid.o
output.o: global.o
physics.o: types.o global.o mcnp_random.o geometry.o output.o \
		search.o endf.o
search.o: output.o
source.o: global.o mcnp_random.o
string.o: global.o output.o
types.o:
unittest.o: global.o energy_grid.o
