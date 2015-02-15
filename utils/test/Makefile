include ../../arch.make

INCFLAGS=`../../FoX-config --fcflags`

.f90.exe: ../libFoX_utils.$(LIBEXT)
	$(FC) $(FFLAGS) $(INCFLAGS) $(LDFLAGS) $(FCFLAGS_free_f90) $(LINK_O_FLAG) $@ $< `../../FoX-config --libs`

check: clean
	OBJEXT=$(OBJEXT) ./run_tests.sh

clean:
	rm -f *.$(OBJEXT) *.exe
