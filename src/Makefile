all:
	mkdir -p build
	cmake -H. -Bbuild
	make -s -C build
clean:
	make -s -C build clean
distclean:
	rm -fr build

.PHONY: all clean distclean
