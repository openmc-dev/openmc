#
default: install
#
OBJFILES= m_wxml_escape.$(OBJEXT) m_wxml_core.$(OBJEXT) \
          m_wxml_overloads.$(OBJEXT) FoX_wxml.$(OBJEXT)
#
LIBRARY=libFoX_wxml.$(LIBEXT)
INCFLAGS=$(MOD_PREFIX)../objs/finclude
#
.m4.F90:
	m4 -I ../m4 $< > $@
#
include ../arch.make
#
CP=cp -p
install: $(OBJFILES)
	@echo "==> Creating $(LIBRARY) with $(OBJFILES)"
	$(AR) $(ARFLAGS_EXTRA) cru $(LIBRARY) $(OBJFILES)
	$(RANLIB) $(LIBRARY)
	$(CP) $(LIBRARY) $(LIB_DIR)
	for i in *.$(MOD_EXT); do \
	  $(CP) $$i $(MOD_DIR); done
#
check:
	(cd test; make) | tee -a ../wxml_lib_check.out
#
clean:
	rm -f *.$(OBJEXT) *.$(MOD_EXT) $(LIBRARY)

m_wxml_overloads.F90: m_wxml_overloads.m4
# DO NOT DELETE THIS LINE - used by make depend
FoX_wxml.$(OBJEXT): m_wxml_core.$(OBJEXT) m_wxml_overloads.$(OBJEXT)
m_wxml_core.$(OBJEXT): m_wxml_escape.$(OBJEXT)
m_wxml_overloads.$(OBJEXT): m_wxml_core.$(OBJEXT)
