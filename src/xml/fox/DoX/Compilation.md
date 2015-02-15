#Configuration and compilation

You will have received the FoX source code as a tar.gz file.

Unpack it as normal, and change directory into the top-level directory, FoX-$VERSION.

### Requirements for use

FoX requires a Fortran 95 compiler - not just Fortran 90. All currently available versions of Fortran compilers claim to support F95. If your favoured compiler is not listed as working, I recommend the use of [g95](www.g95.org), which is free to download and use. And in such a case, please send a bug report to your compiler vendor.

In the event that you need to write a code targetted at multiple compilers, including some which have bugs preventing FoX compilation, please note the possibility of producing a [dummy library](#dummy_library).

##Configuration

* In order to generate the Makefile, make sure that you have a Fortran compiler in your `PATH`, and do:

    `./configure`

This should suffice for most installations. However:

1. You may not be interested in all of the modules that FoX supplies. For example, you may only be interested in output, not input. If so, you can select which modules you want using `--enable-MODULENAME` where MODULENAME is one of `wxml`, `wcml`, `wkml`, `sax`, `dom`. If none are explicitly enabled, then all will be built. (Alternatively, you can exclude modules one at a time with `--disable-MODULENAME`) Thus, for example, if you only care about CML output, and not anything else: `./configure --enable-wcml`

2. If you have more than one Fortran compiler available, or it is not on your `PATH`, you can force the choice by doing:

   `./configure FC=/path/to/compiler/of/choice`

3. It is possible that the configuration fails. In this case
	* please tell me about it so I can fix it
  	* all relevant compiler details are placed in the file arch.make; you may be able to edit that file to allow compilation. Again, if so, please let me know what you need to do.

4. By default the resultant files are installed under the objs directory. If you wish them to be installed elsewhere, you may do

    `./configure --prefix=/path/to/installation`

Note that the configure process encodes the current directory location in several
places.  If you move the FoX directory later on, you will need to re-run configure.

* You may be interested in [dummy compilation](#dummy_library). This is activated with the `--enable-dummy` switch (but only works for wxml/wcml currently).

    `./configure --enable-wcml --enable-dummy`

##Compilation

In order to compile the full library, now simply do:

    make

This will build all the requested FoX modules, and the relevant examples

##Testing

In the full version of the FoX library, there are several testsuites included.

To run them all, simply run `make check` from the top-level directory. This will run the individual testsuites, and collate their results.

If any failures occur (unrelated to known compiler issues, see the [up-to-date list](http://github.com/andreww/fox/issues)), please send a message to the mailing list (<fox-discuss@googlegroups.com>) with details of compiler, hardware platform, and the nature of the failure.

The testsuites for the SAX and DOM libraries are very extensive, and are somewhat fragile, so are not distributed with FoX. Please contact the author for details.

##Linking to an existing program

* The files all having been compiled and installed, you need to link them into your program.

A script is provided which will provide the appropriate compiler and linker flags for you; this will be created after configuration, in the top-level directory, and is called `FoX-config`. It may be taken from there and placed anywhere.

FoX-config takes the following arguments:

* `--fcflags`: return flags for compilation
* `--libs`: return flags for linking
* `--wxml`: return flags for compiling/linking against wxml
* `--wcml`: return flags for compiling/linking against wcml
* `--sax`: return flags for compiling/linking against sax

If it is called with no arguments, it will expand to compile & link flags, thusly:

 	f95 -o program program.f90 `FoX-config`

For compiling only against FoX, do the following:

 	f95 -c `FoX-config --fcflags` sourcefile.f90

For linking only to the FoX library, do:

  	f95 -o program `FoX-config --libs` *.o

or similar, according to your compilation scheme. 

Note that by default, `FoX-config` assumes you are using all modules of the library. If you are only using part, then this can be specified by also passing the name of each module required, like so:

	FoX-config --fcflags --wcml

## Compiling a dummy library

<a name="dummy_library"/>

Because of the shortcomings in some compilers, it is not possible to compile FoX everywhere. Equally, sometimes it is useful to be able to compile a code both with and without support for FoX (perhaps to reduce executable size). Especially where FoX is being used only for additional output, it is useful to be able to run the code and perform computations even without the possibility of XML output.

For this reason, it is possible to compile a dummy version of FoX. This includes all public interfaces, so that your code will compile and link correctly - however none of the subroutines do anything, so you can retain the same version of your code without having to comment out all FoX calls.

Because this dummy version of FoX contains nothing except empty subroutines, it compiles and links with all known Fortran 95 compilers, regardless of compiler bugs.

To compile the dummy code, use the `--enable-dummy` switch. Note that currently the dummy mode is not yet available for the DOM module.

## Alternative build methods

The "-full" versions of FoX are also shipped with files
to help compile the code on using other systems using CMake
or from within Microsoft Visual Studio. Brief instructions for
using these files are below.

### CMake

CMake does not build software itself but generates makefiles or
projectfiles (depending on the platform), that are then used to 
compile the software, it should thus be a cross platform 
method for building FoX (in theory at least).

Files needed for building FoX with CMake are included in the 
"-full" distribution. These can:

* do similar checks as current build tools to check for example how
ABORT and FLUSH work,
wether or not certain bugs in compilers are present,

* compile FoX into static libraries.

* generate Fortran files from m4 sources if m4 is present.

* do out-of-source build, does not interfere with current build system

* do a parallel build of fox (make -j)

However, CMake cannot, at present, build the run the test suite or
the packaging scripts used for release. To build FoX with CMake 
the following is needed:

* CMake version >= 2.6.0

* usual build tools: fortran compiler, make, ld, ...

* m4 to generate fortran files from m4 sources (optional).

**CMake Build instructions (linux):** 
Once you installed cmake, go to the main directory of fox
and create a build directory, and from there, execute cmake thus:

	cd fox/

	mkdir build/ && cd build/

	cmake ../

	make -j

Libaries and module files can then be found in the subdirectories of build.

### Windows

It is also possible to build FoX from within Microsoft Visual Studio
and the file FoX.vfproj contains a Visual Studio project for Intel Fortran
to simplify this process. At time of writing, it is compatible with 
Visual Studio 2011 and Intel Visual Fortran Composer XE 2011.

The project will build FoX in one of the four
 configurations: Win32/x64 and debug/release.
When building FoX for a specific configuration, 
an output library file Fox_debug.lib or Fox.lib 
and associated modules are created in a folder in 
a relative path  ../lib or ../libx64 respectively.

For a given configuration in in your application project
you will then need to:

1. In "Fortran" "General" "Additional Include Directories" add 
the respective modules folder (generated above)

2. In "Linker" "General" "Additional library directories" add the 
path to the respective lib or libx64 folder.

3. In "Linker" "Input" "Additional dependencies" add Fox_debug.lib 
or FoX.lib respectively. 

Your application should now be able to build and link with FoX.

