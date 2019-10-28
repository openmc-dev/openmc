.. _usersguide_troubleshoot:

===============
Troubleshooting
===============

-------------------------
Problems with Compilation
-------------------------

If you are experiencing problems trying to compile OpenMC, first check if the
error you are receiving is among the following options.

-------------------------
Problems with Simulations
-------------------------

Segmentation Fault
******************

A segmentation fault occurs when the program tries to access a variable in
memory that was outside the memory allocated for the program. The best way to
debug a segmentation fault is to re-compile OpenMC with debug options turned
on. Create a new build directory and type the following commands:

.. code-block:: sh

    mkdir build-debug && cd build-debug
    cmake -Ddebug=on /path/to/openmc
    make

Now when you re-run your problem, it should report exactly where the program
failed. If after reading the debug output, you are still unsure why the program
failed, send an email to the OpenMC User's Group `mailing list`_.

ERROR: No cross_sections.xml file was specified in settings.xml or in the OPENMC_CROSS_SECTIONS environment variable.
*********************************************************************************************************************

OpenMC needs to know where to find cross section data for each nuclide.
Information on what data is available and in what files is summarized in a
cross_sections.xml file. You need to tell OpenMC where to find the
cross_sections.xml file either with the :ref:`cross_sections` in settings.xml or
with the :envvar:`OPENMC_CROSS_SECTIONS` environment variable. It is recommended
to add a line in your ``.profile`` or ``.bash_profile`` setting the
:envvar:`OPENMC_CROSS_SECTIONS` environment variable.

Geometry Debugging
******************

Overlapping Cells
^^^^^^^^^^^^^^^^^

For fast run times, normal simulations do not check if the geometry is
incorrectly defined to have overlapping cells.  This can lead to incorrect
results that may or may not be obvious when there are errors in the geometry
input file.  The built-in 2D and 3D plotters will check for cell overlaps at
the center of every pixel or voxel position they process, however this might
not be a sufficient check to ensure correctly defined geometry.  For instance,
if an overlap is of small aspect ratio, the plotting resolution might not be
high enough to produce any pixels in the overlapping area.

To reliably validate a geometry input, it is best to run the problem in
geometry debugging mode with the ``-g``, ``-geometry-debug``, or
``--geometry-debug`` command-line options.  This will enable checks for
overlapping cells at every move of esch simulated particle.  Depending on the
complexity of the geometry input file, this could add considerable overhead to
the run (these runs can still be done in parallel).  As a result, for this run
mode the user will probably want to run fewer particles than a normal
simulation run.  In this case it is important to be aware of how much coverage
each area of the geometry is getting.  For instance, if certain regions do not
have many particles travelling through them there will not be many locations
where overlaps are checked for in that region.  The user should refer to the
output after a geometry debug run to see how many checks were performed in each
cell, and then adjust the number of starting particles or starting source
distributions accordingly to achieve good coverage.

ERROR: After particle __ crossed surface __ it could not be located in any cell and it did not leak.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

For large runs it is often advantageous to run only the offending particle by
using particle restart mode with the ``-s``, ``-particle``, or ``--particle``
command-line options in conjunction with the particle restart files that are
created when particles are lost with this error.

.. _mailing list: https://groups.google.com/forum/?fromgroups=#!forum/openmc-users
