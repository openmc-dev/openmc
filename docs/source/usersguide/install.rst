.. _usersguide_install:

==============================
Installation and Configuration
==============================

.. currentmodule:: openmc

----------------------------------------
Installing on Linux/Mac with conda-forge
----------------------------------------

`Conda <http://conda.pydata.org/docs/>`_ is an open source package management
system and environment management system for installing multiple versions of
software packages and their dependencies and switching easily between
them. `conda-forge <https://conda-forge.github.io/>`_ is a community-led conda
channel of installable packages. For instructions on installing conda, please
consult their `documentation
<http://conda.pydata.org/docs/install/quick.html>`_.

Once you have `conda` installed on your system, add the `conda-forge` channel to
your configuration with:

.. code-block:: sh

    conda config --add channels conda-forge

Once the `conda-forge` channel has been enabled, OpenMC can then be installed
with:

.. code-block:: sh

    conda install openmc

It is possible to list all of the versions of OpenMC available on your platform with:

.. code-block:: sh

    conda search openmc --channel conda-forge

-----------------------------
Installing on Ubuntu with PPA
-----------------------------

For users with Ubuntu 15.04 or later, a binary package for OpenMC is available
through a `Personal Package Archive`_ (PPA) and can be installed through the
`APT package manager`_. First, add the following PPA to the repository sources:

.. code-block:: sh

    sudo apt-add-repository ppa:paulromano/staging

Next, resynchronize the package index files:

.. code-block:: sh

    sudo apt update

Now OpenMC should be recognized within the repository and can be installed:

.. code-block:: sh

    sudo apt install openmc

Binary packages from this PPA may exist for earlier versions of Ubuntu, but they
are no longer supported.

.. _Personal Package Archive: https://launchpad.net/~paulromano/+archive/staging
.. _APT package manager: https://help.ubuntu.com/community/AptGet/Howto

--------------------
Building from Source
--------------------

.. _prerequisites:

Prerequisites
-------------

.. admonition:: Required
   :class: error

    * A Fortran compiler such as gfortran_

      In order to compile OpenMC, you will need to have a Fortran compiler
      installed on your machine. Since a number of Fortran 2003/2008 features
      are used in the code, it is recommended that you use the latest version of
      whatever compiler you choose. For gfortran_, it is necessary to use
      version 4.8.0 or above.

      If you are using Debian or a Debian derivative such as Ubuntu, you can
      install the gfortran compiler using the following command::

          sudo apt install gfortran

    * A C/C++ compiler such as gcc_

      OpenMC includes two libraries written in C and C++, respectively. These
      libraries have been tested to work with a wide variety of compilers. If
      you are using a Debian-based distribution, you can install the g++
      compiler using the following command::

          sudo apt install g++

    * CMake_ cross-platform build system

      The compiling and linking of source files is handled by CMake in a
      platform-independent manner. If you are using Debian or a Debian
      derivative such as Ubuntu, you can install CMake using the following
      command::

          sudo apt install cmake

    * HDF5_ Library for portable binary output format

      OpenMC uses HDF5 for binary output files. As such, you will need to have
      HDF5 installed on your computer. The installed version will need to have
      been compiled with the same compiler you intend to compile OpenMC with. If
      you are using HDF5 in conjunction with MPI, we recommend that your HDF5
      installation be built with parallel I/O features. An example of
      configuring HDF5_ is listed below::

           FC=/opt/mpich/3.1/bin/mpif90 CC=/opt/mpich/3.1/bin/mpicc \
           ./configure --prefix=/opt/hdf5/1.8.12 --enable-fortran \
                       --enable-fortran2003 --enable-parallel

      You may omit ``--enable-parallel`` if you want to compile HDF5_ in serial.

      .. important::

          OpenMC uses various parts of the HDF5 Fortran 2003 API; as such you
          must include ``--enable-fortran2003`` or else OpenMC will not be able
          to compile.

      On Debian derivatives, HDF5 and/or parallel HDF5 can be installed through
      the APT package manager:

      .. code-block:: sh

          sudo apt install libhdf5-dev hdf5-helpers

      Note that the exact package names may vary depending on your particular
      distribution and version.

.. admonition:: Optional
   :class: note

    * An MPI implementation for distributed-memory parallel runs

      To compile with support for parallel runs on a distributed-memory
      architecture, you will need to have a valid implementation of MPI
      installed on your machine. The code has been tested and is known to work
      with the latest versions of both OpenMPI_ and MPICH_. OpenMPI and/or MPICH
      can be installed on Debian derivatives with::

          sudo apt install mpich libmpich-dev
          sudo apt install openmpi-bin libopenmpi-dev

    * git_ version control software for obtaining source code

.. _gfortran: http://gcc.gnu.org/wiki/GFortran
.. _gcc: https://gcc.gnu.org/
.. _CMake: http://www.cmake.org
.. _OpenMPI: http://www.open-mpi.org
.. _MPICH: http://www.mpich.org
.. _HDF5: http://www.hdfgroup.org/HDF5/

Obtaining the Source
--------------------

All OpenMC source code is hosted on GitHub_. You can download the source code
directly from GitHub or, if you have the git_ version control software installed
on your computer, you can use git to obtain the source code. The latter method
has the benefit that it is easy to receive updates directly from the GitHub
repository. GitHub has a good set of `instructions
<http://help.github.com/set-up-git-redirect>`_ for how to set up git to work
with GitHub since this involves setting up ssh_ keys. With git installed and
setup, the following command will download the full source code from the GitHub
repository::

    git clone https://github.com/mit-crpg/openmc.git

By default, the cloned repository will be set to the development branch. To
switch to the source of the latest stable release, run the following commands::

    cd openmc
    git checkout master

.. _GitHub: https://github.com/mit-crpg/openmc
.. _git: http://git-scm.com
.. _ssh: http://en.wikipedia.org/wiki/Secure_Shell

.. _usersguide_build:

Build Configuration
-------------------

Compiling OpenMC with CMake is carried out in two steps. First, ``cmake`` is run
to determine the compiler, whether optional packages (MPI, HDF5) are available,
to generate a list of dependencies between source files so that they may be
compiled in the correct order, and to generate a normal Makefile. The Makefile
is then used by ``make`` to actually carry out the compile and linking
commands. A typical out-of-source build would thus look something like the
following

.. code-block:: sh

    mkdir build && cd build
    cmake ..
    make

Note that first a build directory is created as a subdirectory of the source
directory. The Makefile in the top-level directory will automatically perform an
out-of-source build with default options.

CMakeLists.txt Options
++++++++++++++++++++++

The following options are available in the CMakeLists.txt file:

debug
  Enables debugging when compiling. The flags added are dependent on which
  compiler is used.

profile
  Enables profiling using the GNU profiler, gprof.

optimize
  Enables high-optimization using compiler-dependent flags. For gfortran and
  Intel Fortran, this compiles with -O3.

openmp
  Enables shared-memory parallelism using the OpenMP API. The Fortran compiler
  being used must support OpenMP. (Default: on)

coverage
  Compile and link code instrumented for coverage analysis. This is typically
  used in conjunction with gcov_.

maxcoord
  Maximum number of nested coordinate levels in geometry. Defaults to 10.

To set any of these options (e.g. turning on debug mode), the following form
should be used:

.. code-block:: sh

    cmake -Ddebug=on /path/to/openmc

.. _gcov: https://gcc.gnu.org/onlinedocs/gcc/Gcov.html

.. _usersguide_compile_mpi:

Compiling with MPI
++++++++++++++++++

To compile with MPI, set the :envvar:`FC` and :envvar:`CC` environment variables
to the path to the MPI Fortran and C wrappers, respectively. For example, in a
bash shell:

.. code-block:: sh

    export FC=mpif90
    export CC=mpicc
    cmake /path/to/openmc

Note that in many shells, environment variables can be set for a single command,
i.e.

.. code-block:: sh

    FC=mpif90 CC=mpicc cmake /path/to/openmc

Selecting HDF5 Installation
+++++++++++++++++++++++++++

CMakeLists.txt searches for the ``h5fc`` or ``h5pfc`` HDF5 Fortran wrapper on
your PATH environment variable and subsequently uses it to determine library
locations and compile flags. If you have multiple installations of HDF5 or one
that does not appear on your PATH, you can set the HDF5_ROOT environment
variable to the root directory of the HDF5 installation, e.g.

.. code-block:: sh

    export HDF5_ROOT=/opt/hdf5/1.8.15
    cmake /path/to/openmc

This will cause CMake to search first in /opt/hdf5/1.8.15/bin for ``h5fc`` /
``h5pfc`` before it searches elsewhere. As noted above, an environment variable
can typically be set for a single command, i.e.

.. code-block:: sh

    HDF5_ROOT=/opt/hdf5/1.8.15 cmake /path/to/openmc

.. _compile_linux:

Compiling on Linux and Mac OS X
-------------------------------

To compile OpenMC on Linux or Max OS X, run the following commands from within
the root directory of the source code:

.. code-block:: sh

    mkdir build && cd build
    cmake ..
    make
    make install

This will build an executable named ``openmc`` and install it (by default in
/usr/local/bin). If you do not have administrative privileges, you can install
OpenMC locally by specifying an install prefix when running cmake:

.. code-block:: sh

    cmake -DCMAKE_INSTALL_PREFIX=$HOME/.local ..

The ``CMAKE_INSTALL_PREFIX`` variable can be changed to any path for which you
have write-access.

Compiling on Windows 10
-----------------------

Recent versions of Windows 10 include a subsystem for Linux that allows one to
run Bash within Ubuntu running in Windows. First, follow the installation guide
`here <https://msdn.microsoft.com/en-us/commandline/wsl/install_guide>`_ to get
Bash on Ubuntu on Windows setup. Once you are within bash, obtain the necessary
:ref:`prerequisites <prerequisites>` via ``apt``. Finally, follow the
:ref:`instructions for compiling on linux <compile_linux>`.

Compiling for the Intel Xeon Phi
--------------------------------

For the second generation Knights Landing architecture, nothing special is
required to compile OpenMC. You may wish to experiment with compiler flags that
control generation of vector instructions to see what configuration gives
optimal performance for your target problem.

For the first generation Knights Corner architecture, it is necessary to
cross-compile OpenMC. If you are using the Intel Fortran compiler, it is
necessary to specify that all objects be compiled with the ``-mmic`` flag as
follows:

.. code-block:: sh

    mkdir build && cd build
    FC=ifort CC=icc FFLAGS=-mmic cmake -Dopenmp=on ..
    make

Note that unless an HDF5 build for the Intel Xeon Phi (Knights Corner) is
already on your target machine, you will need to cross-compile HDF5 for the Xeon
Phi. An `example script`_ to build zlib and HDF5 provides several necessary
workarounds.

.. _example script: https://github.com/paulromano/install-scripts/blob/master/install-hdf5-mic

Testing Build
-------------

If you have ENDF/B-VII.1 cross sections from NNDC_ you can test your build.
Make sure the **OPENMC_CROSS_SECTIONS** environmental variable is set to the
*cross_sections.xml* file in the *data/nndc* directory.
There are two ways to run tests. The first is to use the Makefile present in
the source directory and run the following:

.. code-block:: sh

    make test

If you want more options for testing you can use ctest_ command. For example,
if we wanted to run only the plot tests with 4 processors, we run:

.. code-block:: sh

    cd build
    ctest -j 4 -R plot

If you want to run the full test suite with different build options please
refer to our :ref:`test suite` documentation.

--------------------
Python Prerequisites
--------------------

OpenMC's :ref:`Python API <pythonapi>` works with either Python 2.7 or Python
3.2+. In addition to Python itself, the API relies on a number of third-party
packages. All prerequisites can be installed using `conda
<http://conda.pydata.org/docs/>`_ (recommended), `pip
<https://pip.pypa.io/en/stable/>`_, or through the package manager in most Linux
distributions.

.. admonition:: Required
   :class: error

   `six <https://pythonhosted.org/six/>`_
      The Python API works with both Python 2.7+ and 3.2+. To do so, the six
      compatibility library is used.

   `NumPy <http://www.numpy.org/>`_
      NumPy is used extensively within the Python API for its powerful
      N-dimensional array.

   `SciPy <https://www.scipy.org/>`_
      SciPy's special functions, sparse matrices, and spatial data structures
      are used for several optional features in the API.

   `pandas <http://pandas.pydata.org/>`_
      Pandas is used to generate tally DataFrames as demonstrated in
      :ref:`examples_pandas` example notebook.

   `h5py <http://www.h5py.org/>`_
      h5py provides Python bindings to the HDF5 library. Since OpenMC outputs
      various HDF5 files, h5py is needed to provide access to data within these
      files from Python.

.. admonition:: Optional
   :class: note

   `Matplotlib <http://matplotlib.org/>`_
      Matplotlib is used to providing plotting functionality in the API like the
      :meth:`Universe.plot` method and the :func:`openmc.plot_xs` function.

   `uncertainties <https://pythonhosted.org/uncertainties/>`_
      Uncertainties are optionally used for decay data in the :mod:`openmc.data`
      module.

   `Cython <http://cython.org/>`_
      Cython is used for resonance reconstruction for ENDF data converted to
      :class:`openmc.data.IncidentNeutron`.

   `vtk <http://www.vtk.org/>`_
      The Python VTK bindings are needed to convert voxel and track files to VTK
      format.

   `silomesh <https://github.com/nhorelik/silomesh>`_
      The silomesh package is needed to convert voxel and track files to SILO
      format.

   `lxml <http://lxml.de/>`_
      lxml is used for the :ref:`scripts_validate` script.

.. _usersguide_nxml:

-----------------------------------------------------
Configuring Input Validation with GNU Emacs nXML mode
-----------------------------------------------------

The `GNU Emacs`_ text editor has a built-in mode that extends functionality for
editing XML files. One of the features in nXML mode is the ability to perform
real-time `validation`_ of XML files against a `RELAX NG`_ schema. The OpenMC
source contains RELAX NG schemas for each type of user input file. In order for
nXML mode to know about these schemas, you need to tell emacs where to find a
"locating files" description. Adding the following lines to your ``~/.emacs``
file will enable real-time validation of XML input files:

.. code-block:: common-lisp

    (require 'rng-loc)
    (add-to-list 'rng-schema-locating-files "~/openmc/schemas.xml")

Make sure to replace the last string on the second line with the path to the
schemas.xml file in your own OpenMC source directory.

.. _GNU Emacs: http://www.gnu.org/software/emacs/
.. _validation: http://en.wikipedia.org/wiki/XML_validation
.. _RELAX NG: http://relaxng.org/
.. _NNDC: http://www.nndc.bnl.gov/endf/b7.1/acefiles.html
.. _ctest: http://www.cmake.org/cmake/help/v2.8.12/ctest.html
