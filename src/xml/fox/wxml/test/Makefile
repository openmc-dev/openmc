include ../../arch.make

INCFLAGS=`../../FoX-config --fcflags`

.f90.exe:
	$(FC) $(FFLAGS) $(INCFLAGS) $(LDFLAGS) $(FCFLAGS_free_f90) $(LINK_O_FLAG) $@ $< `../../FoX-config --libs --wxml`

check: clean
	./run_tests.sh
	
clean:
	rm -f *.$(OBJEXT) *.exe
