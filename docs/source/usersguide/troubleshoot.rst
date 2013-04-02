.. _usersguide_troubleshoot:

===============
Troubleshooting
===============

-------------------------
Problems with Compilation
-------------------------

If you are experiencing problems trying to compile OpenMC, first check if the
error you are receiving is among the following options.

undefined reference to `_vtab$...
*********************************

If you see this message when trying to compile, the most likely cause is that
you are using a compiler that does not support type-bound procedures from
Fortran 2003. This affects any version of gfortran prior to 4.6. Downloading and
installing the latest gfortran_ compiler should resolve this problem.

Fatal Error: Wrong module version '4' (expected '9') for file 'xml_data_cmfd_t.mod' opened at (1)
*************************************************************************************************

The `.mod` modules files that are created by gfortran are versioned and
sometimes are usually not backwards compatible. If gfortran is upgraded and the
modules files for xml-fortran source files are not deleted, this error may
occur. To fix this, clear out all module and object files with :program:`make
distclean` and then recompiling.

Fatal Error: File 'xml_data_cmfd_t.mod' opened at (1) is not a GFORTRAN module file
***********************************************************************************

When OpenMC compiles, the first thing it needs to do is compile source in the
xml-fortran subdirectory. If you compiled everything with a compiler other than
gfortran, performed a :program:`make clean`, and then tried to :program:`make`
with gfortran, the xml-fortran modules would have been compiled with a different
compiler. To fix this, try clearing out all module and object files with
:program:`make distclean` and then recompiling.

gfortran: unrecognized option '-cpp'
************************************

You are probably using a version of the gfortran compiler that is too
old. Download and install the latest version of gfortran_.

f951: error: unrecognized command line option "-fbacktrace"
***********************************************************

You are probably using a version of the gfortran compiler that is too
old. Download and install the latest version of gfortran_.


make[1]: ifort: Command not found
*********************************

You tried compiling with the Intel Fortran compiler and it was not found on your
:envvar:`PATH`. If you have the Intel compiler installed, make sure the shell
can locate it (this can be tested with :program:`which ifort`).

make[1]: pgf90: Command not found
*********************************

You tried compiling with the PGI Fortran compiler and it was not found on your
:envvar:`PATH`. If you have the PGI compiler installed, make sure the shell can
locate it (this can be tested with :program:`which pgf90`).

-------------------------
Problems with Simulations
-------------------------

Segmentation Fault
******************

A segmentation fault occurs when the program tries to access a variable in
memory that was outside the memory allocated for the program. The best way to
debug a segmentation fault is to re-compile OpenMC with debug options turned
on. First go to your ``openmc/src`` directory where OpenMC was compiled and type
the following commands:

.. code-block:: sh

    make distclean
    make DEBUG=yes

Now when you re-run your problem, it should report exactly where the program
failed. If after reading the debug output, you are still unsure why the program
failed, send an email to the OpenMC User's Group `mailing list`_.

ERROR: No cross_sections.xml file was specified in settings.xml or in the CROSS_SECTIONS environment variable.
**************************************************************************************************************

OpenMC needs to know where to find cross section data for each
nuclide. Information on what data is available and in what files is summarized
in a cross_sections.xml file. You need to tell OpenMC where to find the
cross_sections.xml file either with the :ref:`cross_sections` in settings.xml or
with the :envvar:`CROSS_SECTIONS` environment variable. It is recommended to add
a line in your ``.profile`` or ``.bash_profile`` setting the
:envvar:`CROSS_SECTIONS` environment variable.

ERROR: After particle __ crossed surface __ it could not be located in any cell and it did not leak.
****************************************************************************************************

This error can arise either if a problem is specified with no boundary
conditions or if there is an error in the geometry itself. First check to ensure
that all of the outer surfaces of your geometry have been given vacuum or
reflective boundary conditions. If proper boundary conditions have been applied
and you still receive this error, it means that a surface/cell/lattice in your
geometry has been specified incorrectly or is missing.

The best way to debug this error is to turn on a trace for the particle getting
lost. After the error message, the code will display what batch, generation, and
particle number caused the error. In your settings.xml, add a :ref:`trace`
followed by the batch, generation, and particle number. This will give you
detailed output every time that particle enters a cell, crosses a boundary, or
has a collision. For example, if you received this error at cycle 5, generation
1, particle 4032, you would enter:

.. code-block:: xml

    <trace>5 1 4032</trace>

.. _gfortran: http://gcc.gnu.org/wiki/GFortran
.. _mailing list: https://groups.google.com/forum/?fromgroups=#!forum/openmc-users
