source = $(wildcard *.f90)
test_objects = $(source:.f90=.o)
executables = $(basename $(source))

objects = ../search.o ../constants.o ../error.o ../global.o
objects += ../mesh.o ../mesh_header.o

#===============================================================================
# Compiler Options
#===============================================================================

F90 = ifort
F90FLAGS = -I..
LDFLAGS = 

#===============================================================================
# Targets
#===============================================================================

all: $(executables)
$(executables): $(test_objects) $(objects)
	$(F90) $(LDFLAGS) $@.o $(objects) -o $@
clean:
	@rm -f *.o *.mod $(tests)
neat:
	@rm -f *.o *.mod

#===============================================================================
# Rules
#===============================================================================

.SUFFIXES: .f90 .o
.PHONY: all clean neat

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<
