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
working directory which may be different from
``/home/username/somemodel``. ``openmc`` accepts the following command line
flags:

-c, --volume           Run in stochastic volume calculation mode
-e, --event            Run using event-based parallelism
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

.. note:: If you're using the Python API, :func:`openmc.run` is equivalent to
          running ``openmc`` from the command line.

.. _scripts_ace:

----------------------
``openmc-ace-to-hdf5``
----------------------

This script can be used to create HDF5 nuclear data libraries used by OpenMC if
you have existing ACE files. There are four different ways you can specify ACE
libraries that are to be converted:

1. List each ACE library as a positional argument. This is very useful in
   conjunction with the usual shell utilities (``ls``, ``find``, etc.).
2. Use the ``--xml`` option to specify a pre-v0.9 cross_sections.xml file.
3. Use the ``--xsdir`` option to specify a MCNP xsdir file.
4. Use the ``--xsdata`` option to specify a Serpent xsdata file.

The script does not use any extra information from cross_sections.xml/ xsdir/
xsdata files to determine whether the nuclide is metastable. Instead, the
``--metastable`` argument can be used to specify whether the ZAID naming convention
follows the NNDC data convention (1000*Z + A + 300 + 100*m), or the MCNP data
convention (essentially the same as NNDC, except that the first metastable state
of Am242 is 95242 and the ground state is 95642).

The optional ``--fission_energy_release`` argument will accept an HDF5 file
containing a library of fission energy release (ENDF MF=1 MT=458) data. A
library built from ENDF/B-VII.1 data is released with OpenMC and can be found at
openmc/data/fission_Q_data_endb71.h5. This data is necessary for
'fission-q-prompt' and 'fission-q-recoverable' tallies, but is not needed
otherwise.

-h, --help            show help message and exit

-d DESTINATION, --destination DESTINATION
                      Directory to create new library in

-m META, --metastable META
                      How to interpret ZAIDs for metastable nuclides. META
                      can be either 'nndc' or 'mcnp'. (default: nndc)

--xml XML             Old-style cross_sections.xml that lists ACE libraries

--xsdir XSDIR         MCNP xsdir file that lists ACE libraries

--xsdata XSDATA       Serpent xsdata file that lists ACE libraries

--fission_energy_release FISSION_ENERGY_RELEASE
                      HDF5 file containing fission energy release data

.. _scripts_compton:

-----------------------
``openmc-make-compton``
-----------------------

This script generates an HDF5 file called ``compton_profiles.h5`` that contains
Compton profile data using an existing data library from `Geant4
<http://geant4.cern.ch/>`_. Note that OpenMC includes this data file by default
so it should not be necessary in practice to generate it yourself.


.. _scripts_depletion_chain:

-------------------------------
``openmc-make-depletion-chain``
-------------------------------

This script generates a depletion chain file called ``chain_endfb71.xml``
using ENDF/B-VII.1 nuclear data. If the :envvar:`OPENMC_ENDF_DATA` variable
is not set, and ``"neutron"``, ``"decay"``, ``"nfy"`` directories
do not exist, then ENDF/B-VII.1 data will be downloaded.

.. _scripts_depletion_chain_casl:

------------------------------------
``openmc-make-depletion-chain-casl``
------------------------------------

This script generates a depletion chain called ``chain_casl.xml``
using ENDF/B-VII.1 nuclear data for a simplified chain.
The nuclides were chosen by CASL-ORIGEN, which can be found in
Appendix A of Kang Seog Kim, `"Specification for the VERA Depletion
Benchmark Suite" <https://doi.org/10.2172/1256820>`_,
CASL-U-2015-1014-000, Rev. 0, ORNL/TM-2016/53, 2016.
``Te129`` has been added into this chain due to its link to
``I129`` production.

If the :envvar:`OPENMC_ENDF_DATA` variable is not set,
and ``"neutron"``, ``"decay"``, ``"nfy"`` directories
to not exist, then ENDF/B-VII.1 data will be downloaded.

.. _scripts_stopping:

-------------------------------
``openmc-make-stopping-powers``
-------------------------------

This script generates an HDF5 file called ``stopping_power.h5`` that contains
radiative and collision stopping powers and mean excitation energy pulled from
the `NIST ESTAR database
<https://physics.nist.gov/PhysRefData/Star/Text/ESTAR.html>`_. Note that OpenMC
includes this data file by default so it should not be necessary in practice to
generate it yourself.

.. _scripts_plot:

--------------------------
``openmc-plot-mesh-tally``
--------------------------

``openmc-plot-mesh-tally`` provides a graphical user interface for plotting mesh
tallies. The path to the statepoint file can be provided as an optional arugment
(if omitted, a file dialog will be presented).

.. _scripts_track:

-----------------------
``openmc-track-to-vtk``
-----------------------

This script converts HDF5 :ref:`particle track files <usersguide_track>` to VTK
poly data that can be viewed with ParaView or VisIt. The filenames of the
particle track files should be given as posititional arguments. The output
filename can also be changed with the ``-o`` flag:

-o OUT, --out OUT    Output VTK poly filename

------------------------
``openmc-update-inputs``
------------------------

If you have existing XML files that worked in a previous version of OpenMC that
no longer work with the current version, you can try to update these files using
``openmc-update-inputs``. If any of the given files do not match the most
up-to-date formatting, then they will be automatically rewritten.  The old
out-of-date files will not be deleted; they will be moved to a new file with
'.original' appended to their name.

Formatting changes that will be made:

geometry.xml
  Lattices containing 'outside' attributes/tags will be replaced with lattices
  containing 'outer' attributes, and the appropriate cells/universes will be
  added. Any 'surfaces' attributes/elements on a cell will be renamed 'region'.

materials.xml
  Nuclide names will be changed from ACE aliases (e.g., Am-242m) to HDF5/GND
  names (e.g., Am242_m1). Thermal scattering table names will be changed from
  ACE aliases (e.g., HH2O) to HDF5/GND names (e.g., c_H_in_H2O).

----------------------
``openmc-update-mgxs``
----------------------

This script updates OpenMC's deprecated multi-group cross section XML files to
the latest HDF5-based format.

-i IN, --input IN    Input XML file
-o OUT, --output OUT  Output file in HDF5 format

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
``openmc-voxel-to-vtk``
---------------------------

When OpenMC generates :ref:`voxel plots <usersguide_voxel>`, they are in an
:ref:`HDF5 format <io_voxel>` that is not terribly useful by itself. The
``openmc-voxel-to-vtk`` script converts a voxel HDF5 file to a `VTK
<http://www.vtk.org/>`_ file. To run this script, you will need to have the VTK
Python bindings installed. To convert a voxel file, simply provide the path to
the file:

.. code-block:: sh

   openmc-voxel-to-vtk voxel_1.h5

The ``openmc-voxel-to-vtk`` script also takes the following optional
command-line arguments:

-o, --output   Path to output VTK file
