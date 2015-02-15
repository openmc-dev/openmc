include ../../arch.make

INCFLAGS=`../../FoX-config --fcflags`

.f90.exe: ../libFoX_common.$(LIBEXT)
	$(FC) $(FFLAGS) $(INCFLAGS) $(LDFLAGS) $(FCFLAGS_free_f90) $(LINK_O_FLAG) $@ $< `../../FoX-config --libs --common`

check:
	OBJEXT=$(OBJEXT) ./run_tests.sh

clean:
	rm -f *.$(OBJEXT) *.exe
