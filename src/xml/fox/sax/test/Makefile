include ../../arch.make

INCFLAGS=`../../FoX-config --fcflags`

testdrivers: sax_well_formed.ns.yes.exe sax_well_formed.ns.no.exe sax_valid.ns.yes.exe sax_valid.ns.no.exe

sax_well_formed.ns.yes.exe: sax_well_formed.ns.yes.f90
	$(FC) $(FFLAGS) $(INCFLAGS) $(LDFLAGS) $(FCFLAGS_free_f90) $(LINK_O_FLAG) $@ $< `../../FoX-config --libs --sax`
sax_well_formed.ns.no.exe: sax_well_formed.ns.no.f90
	$(FC) $(FFLAGS) $(INCFLAGS) $(LDFLAGS) $(FCFLAGS_free_f90) $(LINK_O_FLAG) $@ $< `../../FoX-config --libs --sax`
sax_valid.ns.yes.exe: sax_valid.ns.yes.f90
	$(FC) $(FFLAGS) $(INCFLAGS) $(LDFLAGS) $(FCFLAGS_free_f90) $(LINK_O_FLAG) $@ $< `../../FoX-config --libs --sax`
sax_valid.ns.no.exe: sax_valid.ns.no.f90
	$(FC) $(FFLAGS) $(INCFLAGS) $(LDFLAGS) $(FCFLAGS_free_f90) $(LINK_O_FLAG) $@ $< `../../FoX-config --libs --sax`
.f90.exe:
	$(FC) $(FFLAGS) $(INCFLAGS) $(LDFLAGS) $(FCFLAGS_free_f90) $(LINK_O_FLAG) $@ $< m_handlers.$(OBJEXT) `../../FoX-config --libs --sax`

check: clean m_handlers.$(OBJEXT)
	./run_tests.sh
clean:
	rm -f *.$(OBJEXT) *.exe
