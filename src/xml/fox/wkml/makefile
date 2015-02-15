#
default: install
#
OBJFILES=FoX_wkml.$(OBJEXT) m_wkml_lowlevel.$(OBJEXT) m_wkml_color.$(OBJEXT) m_wkml_styling.$(OBJEXT) m_wkml_features.$(OBJEXT) m_wkml_coverage.$(OBJEXT) m_wkml_core.$(OBJEXT) m_wkml_contours.$(OBJEXT) m_contours.$(OBJEXT) m_wkml_color_def.$(OBJEXT) m_wkml_chart.$(OBJEXT)
#
LIBRARY=libFoX_wkml.$(LIBEXT)
INCFLAGS=$(MOD_PREFIX)../objs/finclude
#
.m4.F90:
	m4 -I ../m4 $< > $@
#
include ../arch.make

# Uncomment one of the following lines for debugging/profiling support
#
#FFLAGS=$(FFLAGS_DEBUG)
#FFLAGS=$(FFLAGS_PROFILE)
#FFLAGS=$(FFLAGS_CHECK)
#
CP=cp -p
install: $(OBJFILES)
	@echo "  ==> Updating $(LIBRARY) with $(OBJFILES)"
	@$(AR) $(ARFLAGS_EXTRA) cru $(LIBRARY) $(OBJFILES)
	@$(RANLIB) $(LIBRARY)
	$(CP) $(LIBRARY) $(LIB_DIR) ;
	for i in *.$(MOD_EXT) ; do  \
           $(CP) $$i $(MOD_DIR); done
#
check:
	(cd test; make) | tee -a ../wkml_lib_check.out

#
clean:
	rm -f *.$(OBJEXT)  *.$(MOD_EXT) $(LIBRARY)
# DO NOT DELETE THIS LINE - used by make depend
FoX_wkml.$(OBJEXT): m_wkml_color.$(OBJEXT) m_wkml_contours.$(OBJEXT) m_wkml_core.$(OBJEXT) m_wkml_coverage.$(OBJEXT)
FoX_wkml.$(OBJEXT): m_wkml_features.$(OBJEXT) m_wkml_lowlevel.$(OBJEXT) m_wkml_styling.$(OBJEXT) m_wkml_chart.$(OBJEXT)
m_wkml_color.$(OBJEXT): m_wkml_color_def.$(OBJEXT)
m_wkml_contours.$(OBJEXT): m_contours.$(OBJEXT) m_wkml_color.$(OBJEXT) m_wkml_features.$(OBJEXT)
m_wkml_contours.$(OBJEXT): m_wkml_lowlevel.$(OBJEXT) m_wkml_styling.$(OBJEXT)
m_wkml_core.$(OBJEXT): m_wkml_color.$(OBJEXT) m_wkml_lowlevel.$(OBJEXT) m_wkml_styling.$(OBJEXT)
m_wkml_coverage.$(OBJEXT): m_wkml_color.$(OBJEXT) m_wkml_core.$(OBJEXT) m_wkml_lowlevel.$(OBJEXT)
m_wkml_coverage.$(OBJEXT): m_wkml_styling.$(OBJEXT) m_wkml_chart.$(OBJEXT) m_wkml_features.$(OBJEXT)
m_wkml_features.$(OBJEXT): m_wkml_color.$(OBJEXT) m_wkml_core.$(OBJEXT) m_wkml_lowlevel.$(OBJEXT)
m_wkml_features.$(OBJEXT): m_wkml_styling.$(OBJEXT) m_wkml_chart.$(OBJEXT)
m_wkml_lowlevel.$(OBJEXT): m_wkml_color.$(OBJEXT)
m_wkml_styling.$(OBJEXT): m_wkml_color.$(OBJEXT) m_wkml_lowlevel.$(OBJEXT)
m_wkml_chart.$(OBJEXT): m_wkml_lowlevel.$(OBJEXT)
