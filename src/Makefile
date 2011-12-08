program = openmc

templates = $(wildcard xml-fortran/templates/*.o)
xml_fort = xml-fortran/xmlparse.o \
           xml-fortran/read_xml_primitives.o \
           xml-fortran/write_xml_primitives.o

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

#===============================================================================
# Intel Fortran compiler options
#===============================================================================

ifeq ($(COMPILER),intel)
  F90 = ifort
  F90FLAGS := -fpp -warn -assume byterecl
  LDFLAGS =

  # Debugging options
  ifeq ($(DEBUG),yes)
    F90FLAGS += -g -traceback -ftrapuv -fp-stack-check -check all
    LDFLAGS  += -g
  endif

  # Profiling options
  ifeq ($(PROFILE),yes)
    F90FLAGS += -pg
    LDFLAGS  += -pg
  endif

  # Optimization
  ifeq ($(OPTIMIZE),yes)
    F90FLAGS += -O3 # -ipo
  endif
endif

#===============================================================================
# GNU Fortran compiler options
#===============================================================================

ifeq ($(COMPILER),gfortran)
  F90 = gfortran
  F90FLAGS := -cpp -Wall
  LDFLAGS =

  # Debugging options
  ifeq ($(DEBUG),yes)
    F90FLAGS += -g -pedantic -std=f2008 -fbacktrace -fbounds-check \
                -ffpe-trap=invalid,zero,overflow,underflow
    LDFLAGS  += -g
  endif

  # Profiling options
  ifeq ($(PROFILE),yes)
    F90FLAGS += -pg
    LDFLAGS  += -pg
  endif

  ifeq ($(OPTIMIZE),yes)
    F90FLAGS += -O3
  endif
endif

#===============================================================================
# PGI compiler options
#===============================================================================

ifeq ($(COMPILER),pgi)
  F90 = pgf90
  F90FLAGS := -Mpreprocess -Minform=inform
  LDFLAGS =

  # Since ERFC is currently not supported in the PGI compiler and this intrinsic
  # is used in the doppler module (which is not actually used yet), it is
  # removed from the source objects
  objects := $(subst doppler.o, ,$(objects))

  # Debugging options
  ifeq ($(DEBUG),yes)
    F90FLAGS += -g -Mbounds -Mchkptr -Mchkstk -traceback
    LDFLAGS  += -g
  endif

  # Profiling options
  ifeq ($(PROFILE),yes)
    F90FLAGS += -pg
    LDFLAGS  += -pg
  endif

  ifeq ($(OPTIMIZE),yes)
    F90FLAGS += -fast -Mipa
  endif
endif

#===============================================================================
# Cray compiler options
#===============================================================================

ifeq ($(COMPILER),cray)
  F90 = ftn
  F90FLAGS := -e Z -m 0

  # Debugging options
  ifeq ($(DEBUG),yes)
    F90FLAGS += -g -R abcnsp -O0
    LDFLAGS  += -g
  endif
endif

#===============================================================================
# Miscellaneous compiler options
#===============================================================================

# Use MPI for parallelism

ifeq ($(USE_MPI),yes)
  MPI = /opt/mpich2/1.4.1-intel
  F90 = $(MPI)/bin/mpif90
  F90FLAGS += -DMPI
endif

#===============================================================================
# Special options for ORNL Jaguar supercomputer
#===============================================================================

HOSTNAME = $(shell hostname)
ifneq (,$(findstring jaguar,$(HOSTNAME)))
  F90 = ftn
  F90FLAGS += -DMPI
endif

#===============================================================================
# Targets
#===============================================================================

all: xml-fortran $(program)
xml-fortran:
	cd xml-fortran; make F90=$(F90) F90FLAGS="$(F90FLAGS)"
	cd xml-fortran/templates; make F90=$(F90) F90FLAGS="$(F90FLAGS)"
$(program): $(objects)
	$(F90) $(objects) $(templates) $(xml_fort) -o $@ $(LDFLAGS)
distclean: clean
	cd xml-fortran; make clean
	cd xml-fortran/templates; make clean
clean:
	@rm -f *.o *.mod $(program)
neat:
	@rm -f *.o *.mod

#===============================================================================
# Rules
#===============================================================================

.SUFFIXES: .F90 .o
.PHONY: all xml-fortran clean neat distclean 

%.o: %.F90
	$(F90) -Ixml-fortran -Ixml-fortran/templates $(F90FLAGS) -c $<

#===============================================================================
# Dependencies
#===============================================================================

include DEPENDENCIES
