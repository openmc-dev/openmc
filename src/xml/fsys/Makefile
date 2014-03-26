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
	$(F90) $(F90FLAGS) -I../include -c $<

#===============================================================================
# Dependencies
#===============================================================================

fox_m_fsys_format.o: fox_m_fsys_realtypes.o fox_m_fsys_abort_flush.o
fox_m_fsys_parse_input.o: fox_m_fsys_realtypes.o
fox_m_fsys_count_parse_input.o: fox_m_fsys_realtypes.o
fox_m_fsys_string_list.o: fox_m_fsys_array_str.o
