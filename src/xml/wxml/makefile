source  = $(wildcard *.F90)
objects = $(source:.F90=.o)

#===============================================================================
# Compiler Options
#===============================================================================

# Ignore unusd variables

ifeq ($(MACHINE),bluegene)
  override F90 = xlf2003
endif

ifeq ($(F90),ifort)
  override F90FLAGS += -warn nounused
endif

#===============================================================================
# Targets
#===============================================================================

all: $(objects)

clean:
	@rm -f *.o *.mod
neat:
	@rm -f *.o *.mod

#===============================================================================
# Rules
#===============================================================================

.SUFFIXES: .F90 .o

.PHONY: clean neat

%.o: %.F90
	$(F90) $(F90FLAGS) -c -I../fsys -I../common -I../utils $<

#===============================================================================
# Dependencies
#===============================================================================

FoX_wxml.o: m_wxml_core.o m_wxml_overloads.o
m_wxml_core.o: m_wxml_escape.o
m_wxml_overloads.o: m_wxml_core.o
