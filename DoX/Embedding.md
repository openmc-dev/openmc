#Using FoX in your own project.

The recommended way to use FoX is to embed the full source code as a subdirectory, into an existing project.

In order to do this, you need to do something like the following:

1. Put the full source code as a top-level subdirectory of the tree, called FoX.
2. Incorporate calls to FoX into the program.
3. Incorporate building FoX into your build process.

##To incorporate into the program

It is probably best to isolate use of XML facilities to a small part of the program. This is easily accomplished for XML input,
which will generally happen in only one or two places.

For XML output, this can be more complex. The easiest, and least intrusive way is probably to create a F90 module for your program, looking something like `example_xml_module.f90`

Then you must somewhere (probably in your main program), use this module, and call `initialize_xml_output()` at the start; and then `end_xml_output()` at the end of the program.

In any of the subroutines where you want to output data to the xml file, you should then insert `use example_xml_module` at the beginning of the subroutine. You can then use any of the xml output routines with no further worries, as shown in the examples.

It is easy to make the use of FoX optional, by the use of preprocessor defines. This can be done simply by wrapping each call to your XML wrapper routines in `#ifdef XML`, or similar. Alternatively, the use of the dummy FoX interfaces allows you to switch FoX on and off at compile time - see [Compilation](|Compilation|).

##To incorporate into the build process:

###Configuration

First, FoX must be configured, to ensure that it is set up correctly for your compiler.
(See [Compilation](|Compilation|))
If your main code has a `configure` step, then run FoX's `configure` as part of it.

If your code doesn't have its own configure step, then the first thing that "make" does
should be to configure FoX, if it's not already configured. But that should only happen
once; every time you make your code thereafter, you don't need to re-configure FoX,
because nothing has changed. To do that, put a target like the following in your 
Makefile.

    FoX/.config:
        	(cd FoX; ./configure FC=$(FC))

(Assuming that your `Makefile` already has a variable `FC` which sets the Fortran compiler)

When FoX configure completes, it "touch"es a file called `FoX/.config`. That means that
whenever you re-run your own make, it checks to see if `FoX/.config` exists - if it does,
then it knows FoX doesn't need to be re-configured, so it doesn't bother.

###Compilation of FoX

Then, FoX needs to be compiled before your code (because your modules will depend
on FoX's modules.) But again, it only needs to be compiled once. You won't be changing
FoX, you'll only be changing your own code, so recompiling your code doesn't require
recompiling FoX.

So, add another target like the following;

    FoX/.FoX: FoX/.config
        	(cd FoX; $(MAKE))

This has a dependency on the `configure` script as I showed above, but it will only run it 
if the `configure` script hasn't already been run.

When FoX is successfully compiled, the last thing its `Makefile` does is "touch" the file called
`FoX/.FoX`. So the above target checks to see if that file exists; and if it does, then it doesn't
bother recompiling FoX, because it's already compiled. On the very first time you compile
your code, it will `cd` into the FoX directory and compile it - but then never again.

You then need to have that rule be a dependency of your main target; like so:

      MyExecutable: FoX/.FoX

(or whatever your default `Makefile` rule is).

which will ensure that before `MyExecutable` is compiled, `make` will check to see that FoX
has been compiled (which most of the time it will be, so nothing further will happen).
But the first time you compile your code, it will call the FoX target, and FoX will be
configured & compiled.

###Compiling/linking your code

You should add this to your `FFLAGS` (or equivalent - the variable that holds
flags for compile-time use.

    FFLAGS=-g -O2 -whatever-else $$(FoX/FoX-config --fcflags)

to make sure that you get the path to your modules. (Different compilers have different flags for specifying module
paths; some use `-I`, some use `-M`, _etc_, if you use the above
construction it will pick the right one automatically for your compiler.)

Similarly, for linking, add the following to your `LDFLAGS` (or equivalent - the variable
that holds flags for link-time use.)

    LDFLAGS=-lwhatever $$(FoX/FoX-config --libs)

(For full details of the `FoX-config` script, see [Compilation](|Compilation|))

###Cleaning up

Finally - you probably have a `clean` target in your makefile. Don't tie FoX into this
target - most of the time when you `make clean`, you don't want to `make clean` with 
FoX as well, because there's no need - FoX won't have changed and
it'll take a couple of minutes to recompile.

However, you can add a `distclean` (or something) target, which you use before
moving your code to another machine, that looks like:

    distclean: clean
        	(cd FoX; $(MAKE) distclean)

and that will ensure that when you do `make distclean`, even FoX's object files are
cleaned up. But of course that will mean that you have to reconfigure & recompile
FoX next time you compile your code
