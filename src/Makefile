program = openmc

#===============================================================================
# Object Files
#===============================================================================

include OBJECTS

#===============================================================================
# User Options
#===============================================================================

COMPILER    = intel
DEBUG       = no
PROFILE     = no
OPTIMIZE    = no
USE_MPI     = no
USE_OPENMP  = no
USE_COARRAY = no

#===============================================================================
# Compiler Options
#===============================================================================

# Intel Fortran compiler options

ifeq ($(COMPILER),intel)
  F90 = ifort
  F90FLAGS = -I../lib/intel -fpp -warn
  LDFLAGS = -L../lib/intel -lxmlparse
endif

# GNU Fortran compiler options

ifeq ($(COMPILER),gfortran)
  F90 = gfortran
  F90FLAGS = -I../lib/gfortran -cpp -Wall
  LDFLAGS = -L../lib/gfortran -lxmlparse
endif

# Set compiler flags for debugging

ifeq ($(DEBUG),yes)
  F90FLAGS += -g 
  LDFLAGS  += -g
  ifeq ($(COMPILER),intel)
    F90FLAGS += -traceback -ftrapuv -fp-stack-check -check all
  endif
  ifeq ($(COMPILER),gfortran)
    F90FLAGS += -pedantic -std=f2008 -fbacktrace -fbounds-check \
                -ffpe-trap=invalid,zero,overflow,underflow
  endif
endif

# Set compiler flags for profiling

ifeq ($(PROFILE),yes)
  F90FLAGS += -pg
  LDFLAGS  += -pg
endif

# Set compiler flags for high optimization

ifeq ($(OPTIMIZE),yes)
  F90FLAGS += -O3
  ifeq ($(COMPILER),intel)
    F90FLAGS += -ipo
  endif
endif

# Use MPI for parallelism

ifeq ($(USE_MPI),yes)
  MPI = /opt/mpich2/1.4.1-intel
  F90 = $(MPI)/bin/mpif90
  F90FLAGS += -DMPI
endif

# Use OpenMP for shared-memory parallelism

ifeq ($(USE_OPENMP),yes)
  F90FLAGS += -openmp
  LDFLAGS += -openmp
endif

# Use Fortran 2008 Coarrays for parallelism

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
	$(F90) $(F90FLAGS) -c $<

#===============================================================================
# Dependencies
#===============================================================================

include DEPENDENCIES
