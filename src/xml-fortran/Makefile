reader = xmlreader
xml_fortran_source  = $(wildcard *.f90)
xml_fortran_objects = $(xml_fortran_source:.f90=.o)

#===============================================================================
# User Options
#===============================================================================

COMPILER = intel
DEBUG    = no

#===============================================================================
# Compiler Options
#===============================================================================

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

all: $(reader)
$(reader): $(xml_fortran_objects)
	$(F90) $(xml_fortran_objects) -o $@ $(LDFLAGS)
clean:
	@rm -f *.o *.mod $(reader)
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

read_xml_primitives.o: xmlparse.o
write_xml_primitives.o: xmlparse.o
xmlreader.o: xmlparse.o
