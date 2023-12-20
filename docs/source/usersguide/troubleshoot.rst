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

RuntimeError: OpenMC aborted unexpectedly.
******************************************

This error usually indicates that OpenMC experienced a segmentation fault. A
segmentation fault occurs when the program tries to access a variable in memory
that was outside the memory allocated for the program. The best way to debug a
segmentation fault is to :ref:`compile OpenMC from source <install_source>` with
debug options turned on. Create a new build directory and type the following
commands:

.. code-block:: sh

    mkdir build-debug && cd build-debug
    cmake -DCMAKE_BUILD_TYPE=Debug /path/to/openmc
    make

Now when you re-run your problem, it should report exactly where the program
failed. If after reading the debug output, you are still unsure why the program
failed, post a message on the `OpenMC Discourse Forum
<https://openmc.discourse.group/>`_.

.. _troubleshoot_lost_particles:

WARNING: After particle __ crossed surface __ it could not be located in any cell and it did not leak.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

During a simulation, particles can become "lost" if they reach a surface and
there is no cell defined on the other side of the surface. It is important to
ensure that 1) proper boundary conditions have been applied to the outer
surfaces of your model and 2) all space in your model is filled with a cell,
even regions that are void and have no material assigned to them.

Please see the instructions in :ref:`troubleshoot_geometry` on how to resolve
issues with lost particles.

ERROR: Maximum number of lost particles has been reached.
*********************************************************

See the above description regarding :ref:`lost particles
<troubleshoot_lost_particles>`. When too many particles are lost, the simulation
will abort altogether. Again, please see the instructions in
:ref:`troubleshoot_geometry` on how to resolve issues with lost particles.

ERROR: No cross_sections.xml file was specified in settings.xml or in the OPENMC_CROSS_SECTIONS environment variable.
*********************************************************************************************************************

OpenMC needs to know where to find cross section data for each nuclide.
Information on what data is available and in what files is summarized in a
cross_sections.xml file. You need to tell OpenMC where to find the
cross_sections.xml file either with the :ref:`cross_sections` in settings.xml or
with the :envvar:`OPENMC_CROSS_SECTIONS` environment variable. It is recommended
to add a line in your ``.profile`` or ``.bash_profile`` setting the
:envvar:`OPENMC_CROSS_SECTIONS` environment variable.

RuntimeError: Failed to open HDF5 file with mode 'w': summary.h5
****************************************************************

This often occurs when working with the Python API and executing multiple OpenMC
runs in a script. If an :class:`openmc.StatePoint` is open in the Python interpreter,
the file handle of the statepoint file as well as the linked `summary.h5` file will
be unavailable for writing, causing this error to appear. To avoid this situation,
it is recommended that data be extracted from statepoint files in a context manager:

.. code-block:: python

    with openmc.StatePoint('statepoint.10.h5') as sp:
        k_eff = sp.keff

or that the :meth:`StatePoint.close` method is called before executing a subsequent
OpenMC run.

.. _troubleshoot_geometry:

Geometry Debugging
******************

To identify issues in your geometry, it is highly recommended to use the `OpenMC
Plot Explorer <https://github.com/openmc-dev/plotter/>`_ GUI application. This
application enables you to interactively explore a model, identify regions that
may be missing a cell definition, and identify overlapping cells.

If you are having issues with lost particles, the following procedure may be
helpful. If OpenMC reports, for example, that a particle reaching surface 50
could not be located, look at your geometry.xml to see which cells have a region
definition that includes surface 50, e.g.:

.. code-block:: xml

    <cell id="10" material="5" region="-49 50 -51 52" />

This may indicate that you need to define a cell on the other side of cell 10.
At this point, using the OpenMC Plot Explorer to locate cell 10 may provide a
visual clue as to whether there is a missing or overlapping cell near cell 10.
Working with the unique integer IDs of cells may be cumbersome; if you provide
names to your cells, these names will show up in the Plot Explorer, which will
aid geometry debugging.

Another method to check for overlapping cells in a geometry is to run the problem in
geometry debugging mode with the ``-g``, ``-geometry-debug``, or
``--geometry-debug`` command-line options.  This will enable checks for
overlapping cells at every move of each simulated particle.  Depending on the
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

Depletion
*********

If you are running a depletion simulation and are experiencing random hangs or
crashes, you may need to set::

    openmc.deplete.pool.USE_MULTIPROCESSING = False

in your Python file before making any calls to the integrator. This can be
caused by an MPI implementation that is not compatible with forking (e.g., see
the `OpenMPI FAQ entry about forking
<https://www.open-mpi.org/faq/?category=tuning#fork-warning>`_).
