include ../../arch.make

INCFLAGS=`../../FoX-config --fcflags`

.f90.exe:
	$(FC) $(FFLAGS) $(INCFLAGS) $(LDFLAGS) $(FCFLAGS_free_f90) $(LINK_O_FLAG) $@ $< `../../FoX-config --libs --wkml`

check: clean
	./run_tests.sh

validate: clean
	wget http://schemas.opengis.net/kml/2.2.0/ogckml22.xsd
	wget http://schemas.opengis.net/kml/2.2.0/atom-author-link.xsd
	xmllint --noout --schema ./ogckml22.xsd *.xml
	rm -f *.xsd

clean:
	rm -f *.$(MOD_EXT) *.$(OBJEXT) *.exe 
	rm -rf test_kmlContours_*_real.*
	rm -f test.xml test.out
