templates = $(wildcard *.xml)
objects = $(templates:.xml=.o)

#===============================================================================
# Compiler Options
#===============================================================================

# Add include for subdirectory

override FCFLAGS += -I../xml-fortran

# Ignore unused variables

ifeq ($(FC),ifort)
  override FCFLAGS += -warn nounused
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
	../xml-fortran/xmlreader $(basename $@)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $(basename $@).f90
