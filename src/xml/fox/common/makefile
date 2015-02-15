#
default: install
#
OBJFILES=m_common_attrs.$(OBJEXT) m_common_buffer.$(OBJEXT) \
         m_common_charset.$(OBJEXT) m_common_namespaces.$(OBJEXT) \
         m_common_error.$(OBJEXT) m_common_elstack.$(OBJEXT) \
         m_common_io.$(OBJEXT) FoX_common.$(OBJEXT) \
         m_common_namecheck.$(OBJEXT) m_common_entities.$(OBJEXT) \
         m_common_notations.$(OBJEXT) m_common_element.$(OBJEXT) \
         m_common_struct.$(OBJEXT) m_common_entity_expand.$(OBJEXT) \
         m_common_content_model.$(OBJEXT)
#
LIBRARY=libFoX_common.$(LIBEXT)
INCFLAGS=$(MOD_PREFIX)../objs/finclude
#
#
include ../arch.make

.m4.F90:
	m4 -I ../m4 $< > $@

# Uncomment one of the following lines for debugging/profiling support
#
#FFLAGS=$(FFLAGS_DEBUG)
#FFLAGS=$(FFLAGS_PROFILE)
#FFLAGS=$(FFLAGS_CHECK)
#
CP=cp -p
install: $(OBJFILES)
	@echo "  ==> Creating $(LIBRARY) with $(OBJFILES)"
	@$(AR) $(ARFLAGS_EXTRA) cru $(LIBRARY) $(OBJFILES)
	@$(RANLIB) $(LIBRARY)
	$(CP) $(LIBRARY) $(LIB_DIR) 
	for i in *.$(MOD_EXT); do \
      	    $(CP) $$i $(MOD_DIR); done
#
check:
	(cd test; make) | tee -a ../common_lib_check.out
#
clean:
	rm -f *.$(OBJEXT)  *.$(MOD_EXT) $(LIBRARY)

# DO NOT DELETE THIS LINE - used by make depend
FoX_common.$(OBJEXT): m_common_attrs.$(OBJEXT)
m_common_attrs.$(OBJEXT): m_common_element.$(OBJEXT) m_common_error.$(OBJEXT)
m_common_buffer.$(OBJEXT): m_common_charset.$(OBJEXT) m_common_error.$(OBJEXT)
m_common_element.$(OBJEXT): m_common_charset.$(OBJEXT) m_common_content_model.$(OBJEXT) m_common_error.$(OBJEXT) m_common_namecheck.$(OBJEXT)
m_common_elstack.$(OBJEXT): m_common_content_model.$(OBJEXT) m_common_error.$(OBJEXT)
m_common_entities.$(OBJEXT): m_common_charset.$(OBJEXT) m_common_error.$(OBJEXT)
m_common_entity_expand.$(OBJEXT): m_common_entities.$(OBJEXT) m_common_error.$(OBJEXT)
m_common_entity_expand.$(OBJEXT): m_common_namecheck.$(OBJEXT) m_common_struct.$(OBJEXT)
m_common_io.$(OBJEXT): m_common_error.$(OBJEXT)
m_common_namecheck.$(OBJEXT): m_common_charset.$(OBJEXT)
m_common_namespaces.$(OBJEXT): m_common_attrs.$(OBJEXT) m_common_charset.$(OBJEXT) m_common_error.$(OBJEXT)
m_common_namespaces.$(OBJEXT): m_common_namecheck.$(OBJEXT) m_common_struct.$(OBJEXT)
m_common_notations.$(OBJEXT): m_common_error.$(OBJEXT)
m_common_struct.$(OBJEXT): m_common_charset.$(OBJEXT) m_common_element.$(OBJEXT) m_common_entities.$(OBJEXT)
m_common_struct.$(OBJEXT): m_common_notations.$(OBJEXT)

