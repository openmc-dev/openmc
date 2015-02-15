#
default: install
#
OBJFILES= m_sax_types.$(OBJEXT) m_sax_tokenizer.$(OBJEXT) \
          m_sax_reader.$(OBJEXT) m_sax_parser.$(OBJEXT) \
          m_sax_operate.$(OBJEXT) m_sax_xml_source.$(OBJEXT) \
          FoX_sax.$(OBJEXT)
#
LIBRARY=libFoX_sax.$(LIBEXT)
INCFLAGS=$(MOD_PREFIX)../objs/finclude
#
include ../arch.make
#
CP=cp -p
install: $(OBJFILES)
	@echo "  ==> Creating $(LIBRARY) with $(OBJFILES)"
	@$(AR) $(ARFLAGS_EXTRA) cru $(LIBRARY) $(OBJFILES)
	@$(RANLIB) $(LIBRARY)
	$(CP) $(LIBRARY) $(LIB_DIR) 
	for i in *.$(MOD_EXT) ; do  \
           $(CP) $$i $(MOD_DIR);done
#
check:
	(cd test; make; make check) | tee -a ../sax_lib_check.out
#
clean:
	rm -f *.$(OBJEXT)  *.$(MOD_EXT) $(LIBRARY)
# DO NOT DELETE THIS LINE - used by make depend
FoX_sax.$(OBJEXT): m_sax_operate.$(OBJEXT)
m_sax_operate.$(OBJEXT): m_sax_parser.$(OBJEXT) m_sax_reader.$(OBJEXT) m_sax_types.$(OBJEXT)
m_sax_parser.$(OBJEXT): m_sax_reader.$(OBJEXT) m_sax_tokenizer.$(OBJEXT) m_sax_types.$(OBJEXT)
m_sax_reader.$(OBJEXT): m_sax_xml_source.$(OBJEXT)
m_sax_tokenizer.$(OBJEXT): m_sax_reader.$(OBJEXT) m_sax_types.$(OBJEXT)
m_sax_types.$(OBJEXT): m_sax_reader.$(OBJEXT)
