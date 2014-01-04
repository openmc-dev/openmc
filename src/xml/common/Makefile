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

FoX_common.o: m_common_attrs.o
m_common_attrs.o: m_common_element.o m_common_error.o
m_common_buffer.o: m_common_charset.o m_common_error.o
m_common_element.o: m_common_charset.o m_common_content_model.o m_common_error.o m_common_namecheck.o
m_common_elstack.o: m_common_content_model.o m_common_error.o
m_common_entities.o: m_common_charset.o m_common_error.o
m_common_entity_expand.o: m_common_entities.o m_common_error.o
m_common_entity_expand.o: m_common_namecheck.o m_common_struct.o
m_common_io.o: m_common_error.o
m_common_namecheck.o: m_common_charset.o
m_common_namespaces.o: m_common_attrs.o m_common_charset.o m_common_error.o
m_common_namespaces.o: m_common_namecheck.o m_common_struct.o
m_common_notations.o: m_common_error.o
m_common_struct.o: m_common_charset.o m_common_element.o m_common_entities.o
m_common_struct.o: m_common_notations.o
