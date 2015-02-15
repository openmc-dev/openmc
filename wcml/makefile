#
default: install
#
OBJFILES=FoX_wcml.$(OBJEXT) m_wcml_coma.$(OBJEXT) m_wcml_core.$(OBJEXT) m_wcml_stml.$(OBJEXT) \
         m_wcml_parameter.$(OBJEXT) m_wcml_property.$(OBJEXT) m_wcml_metadata.$(OBJEXT) \
         m_wcml_lattice.$(OBJEXT) m_wcml_geometry.$(OBJEXT) m_wcml_molecule.$(OBJEXT) \
         m_wcml_lists.$(OBJEXT) m_wcml_inputdec.$(OBJEXT)
#
LIBRARY=libFoX_wcml.$(LIBEXT)
INCFLAGS=$(MOD_PREFIX)../objs/finclude
#
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
	$(CP) $(LIBRARY) $(LIB_DIR) ;
	for i in *.$(MOD_EXT) ; do  \
           $(CP) $$i $(MOD_DIR); done
#
check:
	(cd test; make) | tee -a ../wcml_lib_check.out
#
clean:
	rm -f *.$(OBJEXT)  *.$(MOD_EXT) $(LIBRARY)
#
m_wcml_coma.F90: m_wcml_coma.m4
m_wcml_geometry.F90: m_wcml_geometry.m4
m_wcml_lattice.F90: m_wcml_lattice.m4
m_wcml_lists.F90: m_wcml_lists.m4
m_wcml_molecule.F90: m_wcml_molecule.m4
m_wcml_parameter.F90: m_wcml_parameter.m4 ../m4/quantity.m4 ../m4/foreach.m4 ../m4/datatypes.m4 ../m4/common.m4
m_wcml_property.F90: m_wcml_property.m4 ../m4/quantity.m4 ../m4/foreach.m4 ../m4/datatypes.m4 ../m4/common.m4 
m_wcml_stml.F90: m_wcml_stml.m4 ../m4/foreach.m4 ../m4/common.m4 ../m4/datatypes.m4
# DO NOT DELETE THIS LINE - used by make depend
FoX_wcml.$(OBJEXT): m_wcml_coma.$(OBJEXT) m_wcml_core.$(OBJEXT) m_wcml_geometry.$(OBJEXT) m_wcml_lattice.$(OBJEXT)
FoX_wcml.$(OBJEXT): m_wcml_lists.$(OBJEXT) m_wcml_metadata.$(OBJEXT) m_wcml_molecule.$(OBJEXT)
FoX_wcml.$(OBJEXT): m_wcml_parameter.$(OBJEXT) m_wcml_property.$(OBJEXT) m_wcml_stml.$(OBJEXT)
FoX_wcml.$(OBJEXT): m_wcml_inputdec.$(OBJEXT)
m_wcml_coma.$(OBJEXT): m_wcml_stml.$(OBJEXT)
m_wcml_core.$(OBJEXT): m_wcml_metadata.$(OBJEXT)
m_wcml_parameter.$(OBJEXT): m_wcml_stml.$(OBJEXT)
m_wcml_property.$(OBJEXT): m_wcml_stml.$(OBJEXT)
m_wcml_inputdec.$(OBJEXT): m_wcml_metadata.$(OBJEXT) m_wcml_lists.$(OBJEXT)
