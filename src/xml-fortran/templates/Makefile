templates = $(wildcard *.xml)
objects = $(templates:.xml=.o)

#===============================================================================
# User Options
#===============================================================================

COMPILER = intel
DEBUG    = no

#===============================================================================
# Compiler Options
#===============================================================================

# Set compiler defaults

F90 = gfortran
F90FLAGS = -I..
LDFLAGS =

# Intel Fortran compiler options

ifeq ($(COMPILER),intel)
  F90 = ifort
  F90FLAGS += -fpp
endif

# GNU Fortran compiler options

ifeq ($(COMPILER),gfortran)
  F90 = gfortran
  F90FLAGS += -cpp -Wall
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

#===============================================================================
# Targets
#===============================================================================

all: $(objects)
clean:
	@rm -f *.o *.mod *.out *.f90

#===============================================================================
# Rules
#===============================================================================

.SUFFIXES: .f90 .o .xml
.PHONY: all clean
.PRECIOUS: %.f90

%.f90: %.xml
	../xmlreader $(basename $@)

%.o: %.f90
	$(F90) $(F90FLAGS) -c $(basename $@).f90
