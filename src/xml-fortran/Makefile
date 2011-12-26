reader = xmlreader
source  = $(wildcard *.f90)
objects = $(source:.f90=.o)

#===============================================================================
# Compiler Options
#===============================================================================

# Ignore unused variables

ifeq ($(MACHINE),bluegene)
  override F90 = xlf2003
endif

ifeq ($(F90),ifort)
  override F90FLAGS += -warn nounused
endif

#===============================================================================
# Targets
#===============================================================================

all: $(reader)
$(reader): $(objects)
	$(F90) $(objects) -o $@
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
