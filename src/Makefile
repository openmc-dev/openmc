program = openmc
prefix = /usr/local

templates = $(wildcard templates/*.o)
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

COMPILER = gnu
DEBUG    = no
PROFILE  = no
OPTIMIZE = no
MPI      = no
HDF5     = no
PETSC    = no

#===============================================================================
# External Library Paths
#===============================================================================

MPI_DIR   = /opt/mpich/3.0.4-$(COMPILER)
HDF5_DIR  = /opt/hdf5/1.8.11-$(COMPILER)
PHDF5_DIR = /opt/phdf5/1.8.11-$(COMPILER)
PETSC_DIR = /opt/petsc/3.3-p6-$(COMPILER)

#===============================================================================
# Add git SHA-1 hash
#===============================================================================

GIT_SHA1 = $(shell git log -1 2>/dev/null | head -n 1 | awk '{print $$2}')

#===============================================================================
# GNU Fortran compiler options
#===============================================================================

ifeq ($(COMPILER),gnu)
  FC = gfortran
  FCFLAGS := -cpp -fbacktrace
  LDFLAGS =

  # Debugging
  ifeq ($(DEBUG),yes)
    FCFLAGS += -g -Wall -pedantic -std=f2008 -fbounds-check \
                -ffpe-trap=invalid,overflow,underflow
    LDFLAGS  += -g
  endif

  # Profiling
  ifeq ($(PROFILE),yes)
    FCFLAGS += -pg
    LDFLAGS  += -pg
  endif

  # Optimization
  ifeq ($(OPTIMIZE),yes)
    FCFLAGS += -O3
  endif
endif

#===============================================================================
# Intel Fortran compiler options
#===============================================================================

ifeq ($(COMPILER),intel)
  FC = ifort
  FCFLAGS := -fpp -warn -assume byterecl -traceback
  LDFLAGS =

  # Debugging
  ifeq ($(DEBUG),yes)
    FCFLAGS += -g -ftrapuv -fp-stack-check -check all -fpe0
    LDFLAGS  += -g
  endif

  # Profiling
  ifeq ($(PROFILE),yes)
    FCFLAGS += -pg
    LDFLAGS  += -pg
  endif

  # Optimization
  ifeq ($(OPTIMIZE),yes)
    FCFLAGS += -O3
  endif
endif

#===============================================================================
# PGI compiler options
#===============================================================================

ifeq ($(COMPILER),pgi)
  FC = pgf90
  FCFLAGS := -Mpreprocess -DNO_F2008 -Minform=inform -traceback
  LDFLAGS =

  # Debugging
  ifeq ($(DEBUG),yes)
    FCFLAGS += -g -Mbounds -Mchkptr -Mchkstk
    LDFLAGS  += -g
  endif

  # Profiling
  ifeq ($(PROFILE),yes)
    FCFLAGS += -pg
    LDFLAGS  += -pg
  endif

  # Optimization
  ifeq ($(OPTIMIZE),yes)
    FCFLAGS += -fast -Mipa
  endif
endif

#===============================================================================
# IBM XL compiler options
#===============================================================================

ifeq ($(COMPILER),ibm)
  FC = xlf2003
  FCFLAGS := -WF,-DNO_F2008 -O2

  # Debugging
  ifeq ($(DEBUG),yes)
    FCFLAGS += -g -C -qflag=i:i -u
    LDFLAGS  += -g
  endif

  # Profiling
  ifeq ($(PROFILE),yes)
    FCFLAGS += -p
    LDFLAGS  += -p
  endif

  # Optimization
  ifeq ($(OPTIMIZE),yes)
    FCFLAGS += -O3
  endif
endif

#===============================================================================
# Cray compiler options
#===============================================================================

ifeq ($(COMPILER),cray)
  FC = ftn
  FCFLAGS := -e Z -m 0

  # Debugging
  ifeq ($(DEBUG),yes)
    FCFLAGS += -g -R abcnsp -O0
    LDFLAGS  += -g
  endif
endif

#===============================================================================
# Setup External Libraries
#===============================================================================

# MPI for distributed-memory parallelism and HDF5 for I/O

ifeq ($(MPI),yes)
  ifeq ($(HDF5),yes)
    FC = $(PHDF5_DIR)/bin/h5pfc
    FCFLAGS += -DHDF5
  else
    FC = $(MPI_DIR)/bin/mpif90
  endif
  FCFLAGS += -DMPI
else
  ifeq ($(HDF5),yes)
    FC = $(HDF5_DIR)/bin/h5fc
    FCFLAGS += -DHDF5
  endif
endif

# PETSC for CMFD functionality

ifeq ($(PETSC),yes)
  # Check to make sure MPI is set
  ifneq ($(MPI),yes)
    $(error MPI must be enabled to compile with PETSC!)
  endif

  # Set up PETSc environment
  include $(PETSC_DIR)/conf/petscvariables
  FCFLAGS += -I$(PETSC_DIR)/include -DPETSC
  LDFLAGS += $(PETSC_LIB)
endif

#===============================================================================
# Machine-specific setup
#===============================================================================

# IBM Blue Gene/P ANL supercomputer

ifeq ($(MACHINE),bluegene)
  FC = /bgsys/drivers/ppcfloor/comm/xl/bin/mpixlf2003
  FCFLAGS = -WF,-DNO_F2008,-DMPI -O3
  LDFLAGS = -lmpich.cnkf90
endif

# Cray XK6 ORNL Titan supercomputer

ifeq ($(MACHINE),crayxk6)
  FC = ftn
  FCFLAGS += -DMPI
endif

# IBM Blue Gene/Q ANL supercomputer

ifeq ($(MACHINE),bluegeneq)
  FC = mpixlf2003
  FCFLAGS = -WF,-DNO_F2008,-DMPI -O5
endif

#===============================================================================
# Targets
#===============================================================================

all: xml-fortran $(program)
xml-fortran:
	cd xml-fortran; make MACHINE=$(MACHINE) FC=$(FC) FCFLAGS="$(FCFLAGS)"
	cd templates; make FC=$(FC) FCFLAGS="$(FCFLAGS)"
$(program): $(objects)
	$(FC) $(objects) $(templates) $(xml_fort) $(LDFLAGS) -o $@
install:
	@install -D $(program) $(DESTDIR)$(prefix)/bin/$(program)
	@install -D utils/statepoint_cmp.py $(DESTDIR)$(prefix)/bin/statepoint_cmp
	@install -D utils/statepoint_histogram.py $(DESTDIR)$(prefix)/bin/statepoint_histogram
	@install -D utils/statepoint_meshplot.py $(DESTDIR)$(prefix)/bin/statepoint_meshplot
	@install -D ../man/man1/openmc.1 $(DESTDIR)$(prefix)/share/man/man1/openmc.1
	@install -D ../LICENSE $(DESTDIR)$(prefix)/share/doc/$(program)/copyright
uninstall:
	@rm $(DESTDIR)$(prefix)/bin/$(program)
	@rm $(DESTDIR)$(prefix)/bin/statepoint_cmp
	@rm $(DESTDIR)$(prefix)/bin/statepoint_histogram
	@rm $(DESTDIR)$(prefix)/bin/statepoint_meshplot
	@rm $(DESTDIR)$(prefix)/share/man/man1/openmc.1
	@rm $(DESTDIR)$(prefix)/share/doc/$(program)/copyright
distclean: clean
	cd xml-fortran; make clean
	cd templates; make clean
clean:
	@rm -f *.o *.mod $(program)
neat:
	@rm -f *.o *.mod

#===============================================================================
# Rules
#===============================================================================

.SUFFIXES: .F90 .o
.PHONY: all xml-fortran install uninstall clean neat distclean 

%.o: %.F90
	$(FC) $(FCFLAGS) -DGIT_SHA1="\"$(GIT_SHA1)\"" -Ixml-fortran -Itemplates -c $<

#===============================================================================
# Dependencies
#===============================================================================

include DEPENDENCIES
