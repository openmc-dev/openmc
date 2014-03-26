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
	mv *.mod ../include
	mv *.o ../lib
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
	$(F90) $(F90FLAGS) -c -I../include $<

#===============================================================================
# Dependencies
#===============================================================================

FoX_dom.o: m_dom_dom.o m_dom_error.o m_dom_extras.o m_dom_parse.o m_dom_utils.o
m_dom_dom.o: m_dom_error.o
m_dom_extras.o: m_dom_dom.o m_dom_error.o
m_dom_parse.o: m_dom_dom.o m_dom_error.o
m_dom_utils.o: m_dom_dom.o m_dom_error.o
