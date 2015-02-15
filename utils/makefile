#
default: install
#
OBJFILES=FoX_utils.$(OBJEXT) fox_m_utils_mtprng.$(OBJEXT) \
         fox_m_utils_uuid.$(OBJEXT) fox_m_utils_uri.$(OBJEXT)
#
LIBRARY=libFoX_utils.$(LIBEXT)
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
	for i in *.$(MOD_EXT); do \
      	    $(CP) $$i $(MOD_DIR); done
#
clean:
	rm -f *.$(OBJEXT)  *.$(MOD_EXT) $(LIBRARY)
#
check:
	(cd test; make) | tee -a ../utils_lib_check.out
# DO NOT DELETE THIS LINE - used by make depend
FoX_utils.$(OBJEXT): fox_m_utils_uuid.$(OBJEXT) fox_m_utils_uri.$(OBJEXT)
fox_m_utils_uuid.$(OBJEXT): fox_m_utils_mtprng.$(OBJEXT)
