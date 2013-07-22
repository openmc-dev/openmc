xml_lib = lib/libxml.a

#===============================================================================
# Targets
#===============================================================================

$(xml_lib):
	mkdir -p include
	mkdir -p lib
	cd fsys; make MACHINE=$(MACHINE) F90=$(F90) F90FLAGS="$(F90FLAGS)"
	cd utils; make MACHINE=$(MACHINE) F90=$(F90) F90FLAGS="$(F90FLAGS)"
	cd common; make MACHINE=$(MACHINE) F90=$(F90) F90FLAGS="$(F90FLAGS)"
	cd wxml; make MACHINE=$(MACHINE) F90=$(F90) F90FLAGS="$(F90FLAGS)"
	cd sax; make MACHINE=$(MACHINE) F90=$(F90) F90FLAGS="$(F90FLAGS)"
	cd dom; make MACHINE=$(MACHINE) F90=$(F90) F90FLAGS="$(F90FLAGS)"
	cd lib; ar rcs libxml.a *.o; rm -f *.o

clean:
	cd fsys; make clean
	cd utils; make clean
	cd common; make clean
	cd wxml; make clean
	cd sax; make clean
	cd dom; make clean
	@rm -r -f include
	@rm -r -f lib

#===============================================================================
# Rules
#===============================================================================

.PHONY: clean
