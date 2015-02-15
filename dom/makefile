#
default: install
#
M4FILES = m_dom_types.m4 m_dom_node.m4 m_dom_nodelist.m4 \
	  m_dom_namednodemap.m4 m_dom_implementation.m4 \
  	  m_dom_document.m4 m_dom_document_type.m4 \
   	  m_dom_element.m4 m_dom_attribute.m4 \
       	  m_dom_character_data.m4 m_dom_entity.m4 \
          m_dom_processing_instruction.m4 m_dom_text.m4 \
 	  m_dom_common.m4 m_dom_exception.m4 m_dom_dom.m4 \
          m_dom_namespaces.m4 m_dom_configuration.m4 \
          m_dom_treewalk.m4 m_dom_exception.m4 m_dom_object.m4
OBJFILES= m_dom_error.$(OBJEXT) m_dom_parse.$(OBJEXT) m_dom_utils.$(OBJEXT) \
	  m_dom_extras.$(OBJEXT) m_dom_dom.$(OBJEXT) FoX_dom.$(OBJEXT)
#
LIBRARY=libFoX_dom.$(LIBEXT)
INCFLAGS=$(MOD_PREFIX)../objs/finclude
#
include ../arch.make
#
m_dom_extras.F90: m_dom_extras.m4
	m4 -I ../m4 $< > $@
m_dom_parse.f90: m_dom_parse.m4
	m4 -I ../m4 $< > $@
m_dom_utils.f90: m_dom_utils.m4
	m4 -I ../m4 $< > $@
#
m_dom_dom.F90: $(M4FILES)
	m4 -I ../m4 m_dom_dom.m4 > $@ 
#
CP=cp -p
install: $(OBJFILES)
	@echo "  ==> Creating $(LIBRARY) with $(OBJFILES)"
	$(AR) $(ARFLAGS_EXTRA) cru $(LIBRARY) $(OBJFILES)
	$(RANLIB) $(LIBRARY)
	$(CP) $(LIBRARY) $(LIB_DIR)
	for i in *.$(MOD_EXT) ; do  \
           $(CP) $$i $(MOD_DIR); done
#
check:
	(cd test; make) | tee -a ../dom_lib_check.out
	
clean:
	rm -f *.$(OBJEXT)  *.$(MOD_EXT) $(LIBRARY)
# DO NOT DELETE THIS LINE - used by make depend
FoX_dom.$(OBJEXT): m_dom_dom.$(OBJEXT) m_dom_error.$(OBJEXT) m_dom_extras.$(OBJEXT) m_dom_parse.$(OBJEXT) m_dom_utils.$(OBJEXT)
m_dom_dom.$(OBJEXT): m_dom_error.$(OBJEXT)
m_dom_extras.$(OBJEXT): m_dom_dom.$(OBJEXT) m_dom_error.$(OBJEXT)
m_dom_parse.$(OBJEXT): m_dom_dom.$(OBJEXT) m_dom_error.$(OBJEXT)
m_dom_utils.$(OBJEXT): m_dom_dom.$(OBJEXT) m_dom_error.$(OBJEXT)
