#!/usr/bin/make

include ../arch.make

default: $(BUILD_TARGETS)

INCFLAGS=`../FoX-config --fcflags`

EXAMPLES=dom_example_2$(EXEEXT) dom_example_3$(EXEEXT) sax_example$(EXEEXT) sax_example_2$(EXEEXT) wcml_example$(EXEEXT) wxml_example$(EXEEXT) wkml_example$(EXEEXT) wkml_example_2$(EXEEXT)

dom_lib: sax_lib wxml_lib dom_example_2$(EXEEXT) dom_example_3$(EXEEXT) 
sax_lib: sax_example$(EXEEXT) sax_example_2$(EXEEXT)
wcml_lib: wxml_lib wcml_example$(EXEEXT)
wkml_lib: wxml_lib wkml_example$(EXEEXT) wkml_example_2$(EXEEXT)
wxml_lib: wxml_example$(EXEEXT)

times: time_test.xml dom_canonicalize.ns.yes.exe
	@echo "Timing test is to cannonicalize file: " `du -h test.xml`
	./dom_canonicalize.ns.yes.exe

time_test.xml: test.xml

test.xml:
	@echo "Creating data file test.xml ..."
	@echo "<?xml version='1.0' encoding='utf-8'?>" > test.xml 
	@echo "<data xmlns='http://www.xml-cml.org/schema'>" >> test.xml 
	@echo "<line>" >> test.xml ;
	@for i in {1..10000}; \
	do \
	  for i in {1..10} ;  \
	    do \
	        echo "$j 1.206365800459885      0.1074916670129148     0.02144660316977548    0.005740756862772227\c " >> test.xml ;\
	  done; \
	done
	@echo "</line>" >> test.xml ;
	@echo "</data>" >> test.xml
	@echo "                            ... done"

time_clean:
	rm -f test.xml

# Note that we cannot use $< below since IBM make doesn't understand it.

dom_canonicalize.ns.no.exe: dom_canonicalize.ns.no.$(OBJEXT) ../objs/lib/libFoX_dom.$(LIBEXT) ../objs/lib/libFoX_sax.$(LIBEXT) ../objs/lib/libFoX_wxml.$(LIBEXT) ../objs/lib/libFoX_common.$(LIBEXT) ../objs/lib/libFoX_fsys.$(LIBEXT)
	$(FC) $(FFLAGS) $(LDFLAGS) $(LINK_O_FLAG) $@ dom_canonicalize.ns.no.$(OBJEXT) $$(../FoX-config --libs --dom)
dom_canonicalize.ns.yes.exe: dom_canonicalize.ns.yes.$(OBJEXT) ../objs/lib/libFoX_dom.$(LIBEXT) ../objs/lib/libFoX_sax.$(LIBEXT) ../objs/lib/libFoX_wxml.$(LIBEXT) ../objs/lib/libFoX_common.$(LIBEXT) ../objs/lib/libFoX_fsys.$(LIBEXT)
	$(FC) $(FFLAGS) $(LDFLAGS) $(LINK_O_FLAG) $@ dom_canonicalize.ns.yes.$(OBJEXT) $$(../FoX-config --libs --dom)
dom_configuration.exe: dom_configuration.$(OBJEXT) ../objs/lib/libFoX_dom.$(LIBEXT) ../objs/lib/libFoX_sax.$(LIBEXT) ../objs/lib/libFoX_wxml.$(LIBEXT) ../objs/lib/libFoX_common.$(LIBEXT) ../objs/lib/libFoX_fsys.$(LIBEXT) 
	$(FC) $(FFLAGS) $(LDFLAGS) $(LINK_O_FLAG) $@ dom_configuration.$(OBJEXT) $$(../FoX-config --libs --dom)
dom_example_2$(EXEEXT): dom_example_2.$(OBJEXT) ../objs/lib/libFoX_dom.$(LIBEXT) ../objs/lib/libFoX_sax.$(LIBEXT) ../objs/lib/libFoX_wxml.$(LIBEXT) ../objs/lib/libFoX_common.$(LIBEXT) ../objs/lib/libFoX_fsys.$(LIBEXT)
	$(FC) $(FFLAGS) $(LDFLAGS) $(LINK_O_FLAG) $@ dom_example_2.$(OBJEXT) $$(../FoX-config --libs --dom)
dom_example_3$(EXEEXT): dom_example_3.$(OBJEXT) ../objs/lib/libFoX_dom.$(LIBEXT) ../objs/lib/libFoX_sax.$(LIBEXT) ../objs/lib/libFoX_wxml.$(LIBEXT) ../objs/lib/libFoX_common.$(LIBEXT) ../objs/lib/libFoX_fsys.$(LIBEXT)
	$(FC) $(FFLAGS) $(LDFLAGS) $(LINK_O_FLAG) $@ dom_example_3.$(OBJEXT) $$(../FoX-config --libs --dom)
#
sax_example$(EXEEXT): sax_example.$(OBJEXT) ../objs/lib/libFoX_sax.$(LIBEXT) ../objs/lib/libFoX_common.$(LIBEXT) ../objs/lib/libFoX_fsys.$(LIBEXT)
	$(FC) $(FFLAGS) $(LDFLAGS) $(LINK_O_FLAG) $@ sax_example.$(OBJEXT) $$(../FoX-config --libs --sax)
sax_example_2$(EXEEXT): sax_example_2.$(OBJEXT) ../objs/lib/libFoX_sax.$(LIBEXT) ../objs/lib/libFoX_common.$(LIBEXT) ../objs/lib/libFoX_fsys.$(LIBEXT)
	$(FC) $(FFLAGS) $(LDFLAGS) $(LINK_O_FLAG) $@ sax_example_2.$(OBJEXT) $$(../FoX-config --libs --sax)
#
wxml_example$(EXEEXT): wxml_example.$(OBJEXT) ../objs/lib/libFoX_wxml.$(LIBEXT) ../objs/lib/libFoX_common.$(LIBEXT) ../objs/lib/libFoX_fsys.$(LIBEXT)
	$(FC) $(FFLAGS) $(LDFLAGS) $(LINK_O_FLAG) $@ wxml_example.$(OBJEXT) $$(../FoX-config --libs --wxml)
#
wcml_example$(EXEEXT): wcml_example.$(OBJEXT)  ../objs/lib/libFoX_wcml.$(LIBEXT) ../objs/lib/libFoX_wxml.$(LIBEXT) ../objs/lib/libFoX_common.$(LIBEXT) ../objs/lib/libFoX_fsys.$(LIBEXT)
	$(FC) $(FFLAGS) $(LDFLAGS) $(LINK_O_FLAG) $@ wcml_example.$(OBJEXT) $$(../FoX-config --libs --wcml)
#

wkml_example$(EXEEXT): wkml_example.$(OBJEXT)  ../objs/lib/libFoX_wkml.$(LIBEXT) ../objs/lib/libFoX_wxml.$(LIBEXT) ../objs/lib/libFoX_common.$(LIBEXT) ../objs/lib/libFoX_fsys.$(LIBEXT)
	$(FC) $(FFLAGS) $(LDFLAGS) $(LINK_O_FLAG) $@ wkml_example.$(OBJEXT) $$(../FoX-config --libs --wkml)

wkml_example_2$(EXEEXT): wkml_example_2.$(OBJEXT)  ../objs/lib/libFoX_wkml.$(LIBEXT) ../objs/lib/libFoX_wxml.$(LIBEXT) ../objs/lib/libFoX_common.$(LIBEXT) ../objs/lib/libFoX_fsys.$(LIBEXT)
	$(FC) $(FFLAGS) $(LDFLAGS) $(LINK_O_FLAG) $@ wkml_example_2.$(OBJEXT) $$(../FoX-config --libs --wkml)

clean: time_clean
	rm -f *.$(OBJEXT) *.$(MOD_EXT) $(EXAMPLES) wkml_example.kml wkml_example_2.kml dom_example_3.xml
