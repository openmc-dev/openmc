#!/bin/make
#
include arch.make
#
default: objsdir $(BUILD_TARGETS) examples_build
	@if [ -z "$(BUILD_TARGETS)" ]; then echo "FoX is not configured!"; else touch .FoX; fi
#
objsdir:
	mkdir -p objs/lib objs/finclude
#
# Note the hackery to fix the prefix in FoX-config when installed without breaking
# use from the local directory (which would just need a one-line fix in 
# FoX-config.in). We restore the origional FoX-config so local version still works.
install: objsdir $(BUILD_TARGETS)
	$(MKDIR_P) $(install_prefix)/lib $(install_prefix)/finclude $(install_prefix)/bin
	$(INSTALL) objs/lib/* $(install_prefix)/lib
	$(INSTALL) -m 644 objs/finclude/* $(install_prefix)/finclude
	sed -e s#comp_prefix=.*#comp_prefix=$(install_prefix)# FoX-config > FoX-config.tmp
	mv FoX-config FoX-config.old ; mv FoX-config.tmp FoX-config
	$(INSTALL) FoX-config $(install_prefix)/bin
	mv FoX-config.old FoX-config
#
examples_build: $(BUILD_TARGETS)
	if test -d examples; then (cd examples; $(MAKE) VPATH=$(VPATH)/examples) fi
#
#---------------------------
#
# Recursive make for each module
#
dom_lib: objsdir sax_lib wxml_lib
	(cd dom; $(MAKE) VPATH=$(VPATH)/dom)
dom_lib_clean:
	if test -d dom; then (cd dom; $(MAKE) VPATH=$(VPATH)/dom clean) fi
dom_lib_check: $(BUILD_TARGETS) common_lib utils_lib
	if test -d examples && ! grep DUMMYLIB arch.make > /dev/null ; \
	then \
	    rm -f dom_lib_check ; \
	    rm -f dom_lib_check.out ; \
	    (cd dom; $(MAKE) VPATH=$(VPATH)/dom check) ; \
	    touch dom_lib_check ; \
	    touch dom_lib_check.out ; \
	fi
#
sax_lib: objsdir common_lib utils_lib fsys_lib
	(cd sax; $(MAKE) VPATH=$(VPATH)/sax)
sax_lib_clean:
	if test -d sax; then (cd sax; $(MAKE) VPATH=$(VPATH)/sax clean) fi
sax_lib_check: $(BUILD_TARGETS) common_lib utils_lib
	if test -d examples && ! grep DUMMYLIB arch.make > /dev/null ; \
	then \
	    rm -f sax_lib_check ; \
	    rm -f sax_lib_check.out ; \
	    (cd sax; $(MAKE) VPATH=$(VPATH)/sax check) ; \
	    touch sax_lib_check ; \
	    touch sax_lib_check.out ; \
	fi
#
wxml_lib: objsdir common_lib fsys_lib 
	(cd wxml; $(MAKE) VPATH=$(VPATH)/wxml)
wxml_lib_clean:
	if test -d wxml; then (cd wxml; $(MAKE) VPATH=$(VPATH)/wxml clean) fi
wxml_lib_check: $(BUILD_TARGETS) common_lib utils_lib
	if test -d examples && ! grep DUMMYLIB arch.make > /dev/null ; \
	then \
	    rm -f wxml_lib_check ; \
	    rm -f wxml_lib_check.out ; \
	    (cd wxml; $(MAKE) VPATH=$(VPATH)/wxml check) ; \
	    touch wxml_lib_check ; \
	    touch wxml_lib_check.out ; \
	fi
#
wcml_lib: objsdir utils_lib wxml_lib
	(cd wcml; $(MAKE) VPATH=$(VPATH)/wcml)
wcml_lib_clean: 
	if test -d wcml; then (cd wcml; $(MAKE) VPATH=$(VPATH)/wcml clean) fi
wcml_lib_check: $(BUILD_TARGETS) common_lib utils_lib
	if test -d examples && ! grep DUMMYLIB arch.make > /dev/null ; \
	then \
	    rm -f wcml_lib_check ; \
	    rm -f wcml_lib_check.out ; \
	    (cd wcml; $(MAKE) VPATH=$(VPATH)/wcml check) ; \
	    touch wcml_lib_check ; \
	    touch wcml_lib_check.out ; \
	fi
#
wkml_lib: objsdir utils_lib wxml_lib
	(cd wkml; $(MAKE) VPATH=$(VPATH)/wkml)
wkml_lib_clean:
	if test -d wkml; then (cd wkml; $(MAKE) VPATH=$(VPATH)/wkml clean) fi
wkml_lib_check: wkml_lib
	if test -d examples && ! grep DUMMYLIB arch.make > /dev/null ; \
	then \
	    rm -f wkml_lib_check ; \
	    rm -f wkml_lib_check.out ; \
	    (cd wkml; $(MAKE) VPATH=$(VPATH)/wkml check) ; \
	    touch wkml_lib_check ; \
	    touch wkml_lib_check.out ; \
	fi


common_lib: objsdir fsys_lib utils_lib
	(cd common; $(MAKE) VPATH=$(VPATH)/common)
common_lib_clean:
	if test -d common; then (cd common; $(MAKE) VPATH=$(VPATH)/common clean) fi
common_lib_check: $(BUILD_TARGETS) common_lib utils_lib
	if test -d examples && ! grep DUMMYLIB arch.make > /dev/null ; \
	then \
	    rm -f common_lib_check ; \
	    rm -f common_lib_check.out ; \
	    (cd common; $(MAKE) VPATH=$(VPATH)/common check) ; \
	    touch common_lib_check ; \
	    touch common_lib_check.out ; \
	fi
#
utils_lib: objsdir fsys_lib
	(cd utils; $(MAKE) VPATH=$(VPATH)/utils)
utils_lib_clean:
	if test -d utils; then (cd utils; $(MAKE) VPATH=$(VPATH)/utils clean) fi
utils_lib_check: $(BUILD_TARGETS) common_lib utils_lib
	if test -d examples && ! grep DUMMYLIB arch.make > /dev/null ; \
	then \
	    rm -f utils_lib_check ; \
	    rm -f utils_lib_check.out ; \
	    (cd utils; $(MAKE) VPATH=$(VPATH)/utils check) ; \
	    touch utils_lib_check ; \
	    touch utils_lib_check.out ; \
	fi
#
fsys_lib: objsdir
	(cd fsys; $(MAKE) VPATH=$(VPATH)/fsys)
fsys_lib_clean:
	if test -d fsys; then (cd fsys; $(MAKE) VPATH=$(VPATH)/fsys clean) fi
#
check: $(foreach target,$(BUILD_TARGETS) common_lib utils_lib,$(target)_check)
	@if ! test -d examples; then echo "You need to download the full version of FoX to run the testsuite"; \
	elif grep DUMMYLIB arch.make > /dev/null; then echo "You cannot run the testsuite on the dummy library"; \
	else \
	rm -f check.out; \
	touch check.out; \
	for i in $(BUILD_TARGETS) common_lib utils_lib; do cat $$i''_check.out >> check.out ; rm -f $$i''_check.out ; done; \
	grep RESULT check.out; \
	fi
#
DoX:
	(cd DoX; $(MAKE) VPATH=$(VPATH)/DoX)
#
cutdown:
	rm -rf .gitignore DoX/ config/aclocal.m4 config/autom4te.cache config/configure.ac config/m4/ config/makefile examples/ m4/ */test/ */*.m4 Changelog RELEASE release.sh
	rm -rf cmake/
	rm -rf CMakeLists.txt */CMakeLists.txt
	rm -rf FoX.vfproj Fox.vfproj.README
	for i in */makefile ; \
	do sed -e /m4/d $$i'' > $$i''.tmp ; \
	mv $$i''.tmp $$i'' ; \
	done

cutdown-wxml: cutdown
	rm -rf wcml/ wkml/ sax/ dom/

cutdown-wcml: cutdown
	rm -rf wkml/ sax/ dom/

cutdown-wkml: cutdown
	rm -rf sax/ dom/ wcml/

cutdown-sax: cutdown
	rm -rf wxml/ wcml/ wkml/ dom/

cutdown-dom: cutdown
	rm -rf wkml/ wcml/



#clean: wxml_lib_clean wcml_lib_clean common_lib_clean fsys_lib_clean sax_lib_clean dom_lib_clean utils_lib_clean
clean: wkml_lib_clean wxml_lib_clean wcml_lib_clean common_lib_clean fsys_lib_clean sax_lib_clean dom_lib_clean utils_lib_clean
	if test -d examples; then (cd examples; $(MAKE) VPATH=$(VPATH)/examples clean) fi
	rm -rf objs .FoX check.out *_check *_check.out
#
distclean: clean
	rm -f FoX-config arch.make config.log config.status .config check.out
	touch arch.make
