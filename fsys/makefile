#
default: install
#
OBJFILES=fox_m_fsys_abort_flush.$(OBJEXT) \
         fox_m_fsys_array_str.$(OBJEXT) \
         fox_m_fsys_format.$(OBJEXT) \
         fox_m_fsys_parse_input.$(OBJEXT) \
         fox_m_fsys_count_parse_input.$(OBJEXT) \
         fox_m_fsys_string.$(OBJEXT) \
         fox_m_fsys_string_list.$(OBJEXT) \
         fox_m_fsys_realtypes.$(OBJEXT) \
         fox_m_fsys_varstr.$(OBJEXT)
#
LIBRARY=libFoX_fsys.$(LIBEXT)
INCFLAGS=

.m4.F90:
	m4 -I ../m4 $< > $@
#
include ../arch.make
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
clean:
	rm -f *.$(OBJEXT) *.$(MOD_EXT) $(LIBRARY)

# DO NOT DELETE THIS LINE - used by make depend
fox_m_fsys_format.$(OBJEXT): fox_m_fsys_realtypes.$(OBJEXT) fox_m_fsys_abort_flush.$(OBJEXT)
fox_m_fsys_parse_input.$(OBJEXT): fox_m_fsys_realtypes.$(OBJEXT)
fox_m_fsys_count_parse_input.$(OBJEXT): fox_m_fsys_realtypes.$(OBJEXT)
fox_m_fsys_string_list.$(OBJEXT): fox_m_fsys_array_str.$(OBJEXT)
