.. _usersguide_scripts:

=======================
Executables and Scripts
=======================

.. _scripts_openmc:

----------
``openmc``
----------

Once you have a model built (see :ref:`usersguide_basics`), you can either run
the openmc executable directly from the directory containing your XML input
files, or you can specify as a command-line argument the directory containing
the XML input files. For example, if your XML input files are in the directory
``/home/username/somemodel/``, one way to run the simulation would be:

.. code-block:: sh

    cd /home/username/somemodel
    openmc

Alternatively, you could run from any directory:

.. code-block:: sh

    openmc /home/username/somemodel

Note that in the latter case, any output files will be placed in the present
working directory which may be different from ``/home/username/somemodel``. If
you're using the Python API, :func:`openmc.run` is equivalent to running
``openmc`` from the command line. ``openmc`` accepts the following command line
flags:

-c, --volume           Run in stochastic volume calculation mode
-g, --geometry-debug   Run in geometry debugging mode, where cell overlaps are
                       checked for after each move of a particle
-n, --particles N      Use *N* particles per generation or batch
-p, --plot             Run in plotting mode
-r, --restart file     Restart a previous run from a state point or a particle
                       restart file
-s, --threads N        Run with *N* OpenMP threads
-t, --track            Write tracks for all particles
-v, --version          Show version information
-h, --help             Show help message

----------------------
``openmc-ace-to-hdf5``
----------------------

------------------------------
``openmc-convert-mcnp70-data``
------------------------------

------------------------------
``openmc-convert-mcnp71-data``
------------------------------

------------------------
``openmc-get-jeff-data``
------------------------

-----------------------------
``openmc-get-multipole-data``
-----------------------------

------------------------
``openmc-get-nndc-data``
------------------------

--------------------------
``openmc-plot-mesh-tally``
--------------------------

-----------------------
``openmc-track-to-vtk``
-----------------------

------------------------
``openmc-update-inputs``
------------------------

----------------------
``openmc-update-mgxs``
----------------------

.. _scripts_validate:

-----------------------
``openmc-validate-xml``
-----------------------

Input files can be checked before executing OpenMC using the
``openmc-validate-xml`` script which is installed alongside the Python API. Two
command line arguments can be set when running ``openmc-validate-xml``:

-i, --input-path      Location of OpenMC input files.
-r, --relaxng-path    Location of OpenMC RelaxNG files

If the RelaxNG path is not set, the script will search for these files because
it expects that the user is either running the script located in the install
directory ``bin`` folder or in ``src/utils``. Once executed, it will match
OpenMC XML files with their RelaxNG schema and check if they are valid.  Below
is a table of the messages that will be printed after each file is checked.

========================  ===================================
Message                   Description
========================  ===================================
[XML ERROR]               Cannot parse XML file.
[NO RELAXNG FOUND]        No RelaxNG file found for XML file.
[NOT VALID]               XML file does not match RelaxNG.
[VALID]                   XML file matches RelaxNG.
========================  ===================================

.. _scripts_voxel:

---------------------------
``openmc-voxel-to-silovtk``
---------------------------

When OpenMC generates :ref:`voxel plots <usersguide_voxel>`, they are in an
:ref:`HDF5 format <io_voxel>` that is not terribly useful by itself. The
``openmc-voxel-to-silovtk`` script converts a voxel HDF5 file to `VTK
<http://www.vtk.org/>`_ or `SILO
<https://wci.llnl.gov/simulation/computer-codes/silo>`_ file. For VTK, you need
to have the VTK Python bindings installed. For SILO, you need to have `silomesh
<https://github.com/nhorelik/silomesh>`_ installed. To convert a voxel file,
simply provide the path to the file:

.. code-block:: sh

   openmc-voxel-to-silovtk voxel_1.h5

The ``openmc-voxel-to-silovtk`` script also takes the following optional
command-line arguments:

-o, --output   Path to output VTK or SILO file
-s, --silo     Flag to convert to SILO instead of VTK
