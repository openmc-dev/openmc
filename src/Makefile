program = openmc

#===============================================================================
# Object Files
#===============================================================================

include OBJECTS

#===============================================================================
# User Options
#===============================================================================

DEBUG       = no
PROFILE     = no
OPTIMIZE    = no
USE_MPI     = no
USE_OPENMP  = no
USE_COARRAY = no

#===============================================================================
# Compiler Options
#===============================================================================

F90      = ifort
F90FLAGS = -cpp -warn
LDFLAGS  =

ifeq ($(DEBUG),yes)
  F90FLAGS += -g -traceback -ftrapuv -fp-stack-check -check all
  LDFLAGS  += -g
endif

ifeq ($(PROFILE),yes)
  F90FLAGS += -pg
  LDFLAGS  += -pg
endif

ifeq ($(OPTIMIZE),yes)
  F90FLAGS += -ipo -O3
endif

ifeq ($(USE_MPI),yes)
  MPI = /opt/mpich2/1.4-intel
  F90 = $(MPI)/bin/mpif90
  F90FLAGS += -DMPI
endif

ifeq ($(USE_OPENMP),yes)
  F90FLAGS += -openmp
  LDFLAGS += -openmp
endif

ifeq ($(USE_COARRAY),yes)
  F90FLAGS += -coarray
endif

#===============================================================================
# Targets
#===============================================================================

all:	$(program)
$(program): $(objects)
	$(F90) $(objects) -o $@ $(LDFLAGS)
clean:
	@rm -f *.o *.mod $(program)
neat:
	@rm -f *.o *.mod

#===============================================================================
# Rules
#===============================================================================

.SUFFIXES: .f90 .o
.PHONY: all clean neat

%.o: %.f90
	$(F90) -c $(F90FLAGS) $<

#===============================================================================
# Dependencies
#===============================================================================

include DEPENDENCIES
