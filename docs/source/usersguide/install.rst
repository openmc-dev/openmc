.. _usersguide_install:

==============================
Installation and Configuration
==============================

.. currentmodule:: openmc

.. _install_conda:

--------------------------------------------------
Installing on Linux with Mamba and conda-forge
--------------------------------------------------

`Conda <https://conda.io/en/latest/>`_ is an open source package management
systems and environments management system for installing multiple versions of
software packages and their dependencies and switching easily between them.
`Mamba <https://mamba.readthedocs.io/en/latest/>`_ is a cross-platform package
manager and is compatible with `conda` packages.
OpenMC can be installed in a `conda` environment with `mamba`.
First, `conda` should be installed with one of the following installers:
`Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_,
`Anaconda <https://www.anaconda.com/>`_, or `Miniforge <https://github.com/conda-forge/miniforge>`_.
Once you have `conda` installed on your system, OpenMC can be installed via the
`conda-forge` channel with `mamba`.

Note: Conda will not find OpenMC. If you want to use OpenMC on Mac, 
we recommend installing from source.

First, add the `conda-forge` channel with:

.. code-block:: sh

    conda config --add channels conda-forge

Then create and activate a new conda enviroment called `openmc-env` in
which to install OpenMC.

.. code-block:: sh

    conda create -n openmc-env
    conda activate openmc-env

Then install `mamba`, which will be used to install OpenMC.

.. code-block:: sh

    conda install mamba

To list the versions of OpenMC that are available on the `conda-forge` channel,
in your terminal window or an Anaconda Prompt run:

.. code-block:: sh

    mamba search openmc

OpenMC can then be installed with:

.. code-block:: sh

    mamba install openmc

You are now in a conda environment called `openmc-env` that has OpenMC installed.

-------------------------------------------
Installing on Linux/Windows with Docker
-------------------------------------------

OpenMC can be easily deployed using `Docker <https://www.docker.com/>`_ on any
Windows, Mac, or Linux system. With Docker running, execute the following
command in the shell to download and run a `Docker image`_ with the most recent
release of OpenMC from `DockerHub <https://hub.docker.com/>`_:

.. code-block:: sh

    docker run openmc/openmc:latest

This will take several minutes to run depending on your internet download speed.
The command will place you in an interactive shell running in a `Docker
container`_ with OpenMC installed.

.. note:: The ``docker run`` command supports many `options`_ for spawning
          containers including `mounting volumes`_ from the host filesystem,
          which many users will find useful.

.. _Docker image: https://docs.docker.com/engine/reference/commandline/images/
.. _Docker container: https://www.docker.com/resources/what-container
.. _options: https://docs.docker.com/engine/reference/commandline/run/
.. _mounting volumes: https://docs.docker.com/storage/volumes/

.. _install-spack:

OpenMC does not currently support a Mac Docker image.

----------------------------------
Installing from Source using Spack (recommended for Mac)
----------------------------------

Though Spack might take longer to install OpenMC than Mamba and Code-Forge, 
it is currently the recommended way to install OpenMC on Mac. You can also install
OpenMC from Source as outlined in the next section.

Spack is a package management tool designed to support multiple versions and
configurations of software on a wide variety of platforms and environments.
Please follow Spack's `setup guide`_ to configure the Spack system.

The OpenMC Spack recipe has been configured with variants that match most
options provided in the CMakeLists.txt file. To see a list of these variants and
other information use:

.. code-block:: sh

    spack info openmc

.. note::

    It should be noted that by default OpenMC is built with
    `-DCMAKE_BUILD_TYPE=RelwithDebInfo`. In addition, MPI is OFF while OpenMP is
    ON.

It is recommended to install OpenMC with the Python API. Information about this
Spack recipe can be found with the following command:

.. code-block:: sh

    spack info py-openmc

.. note::

   The only variant for the Python API is ``mpi``.

The most basic installation of OpenMC can be accomplished by entering the
following command:

.. code-block::

    spack install py-openmc

.. caution::

    When installing any Spack package, dependencies are assumed to be at
    configured defaults unless otherwise specfied in the specification on the
    command line. In the above example, assuming the default options weren't
    changed in Spack's package configuration, py-openmc will link against a
    non-MPI non-release build of openmc. Even if a release build of openmc was
    built separately, it will rebuild openmc with the default build type. Thus,
    if you are trying to link against dependencies that were configured
    different than defaults, ``^openmc[variants]`` will have to be present in
    the command.

For a release build of OpenMC with MPI support on (provided by OpenMPI), the
following command can be used:

.. code-block:: sh

    spack install py-openmc +mpi ^openmpi ^openmc build_type=Release

.. note::

   ``+mpi`` is automatically forwarded to OpenMC.

.. tip::

    When installing py-openmc, it will use Spack's preferred Python. For
    example, assuming Spack's preferred Python is 3.8.7, to build py-openmc
    against the latest Python 3.7 instead, ``^python@3.7.0:3.7.99`` should be
    added to the specification on the command line. Additionally, a compiler
    type and version can be specified at the end of the command using
    ``%gcc@<version>``, ``%intel@<version>``, etc.

A useful tool in Spack is to look at the dependency tree before installation.
This can be observed using Spack's ``spec`` tool:

.. code-block::

    spack spec py-openmc +mpi ^openmc build_type=Release

Once installed, environment/lmod modules can be generated or Spack's ``load``
feature can be used to access the installed packages.

.. _Spack: https://spack.readthedocs.io/en/latest/
.. _setup guide: https://spack.readthedocs.io/en/latest/getting_started.html


.. _install_source:

----------------------
Installing from Source
----------------------

.. _prerequisites:

Prerequisites
-------------

.. admonition:: Required
   :class: error

    * A C/C++ compiler such as gcc_

      OpenMC's core codebase is written in C++. The source files have been
      tested to work with a wide variety of compilers. If you are using a
      Debian-based distribution, you can install the g++ compiler using the
      following command::

          sudo apt install g++

    * CMake_ cross-platform build system

      The compiling and linking of source files is handled by CMake in a
      platform-independent manner. If you are using Debian or a Debian
      derivative such as Ubuntu, you can install CMake using the following
      command::

          sudo apt install cmake

    * HDF5_ Library for portable binary output format

      OpenMC uses HDF5 for many input/output files. As such, you will need to
      have HDF5 installed on your computer. The installed version will need to
      have been compiled with the same compiler you intend to compile OpenMC
      with. If compiling with gcc from the APT repositories, users of Debian
      derivatives can install HDF5 and/or parallel HDF5 through the package
      manager::

          sudo apt install libhdf5-dev

      Parallel versions of the HDF5 library called `libhdf5-mpich-dev` and
      `libhdf5-openmpi-dev` exist which are built against MPICH and OpenMPI,
      respectively. To link against a parallel HDF5 library, make sure to set
      the HDF5_PREFER_PARALLEL CMake option, e.g.::

          cmake -DHDF5_PREFER_PARALLEL=on -DOPENMC_USE_MPI=on ..

      Note that the exact package names may vary depending on your particular
      distribution and version.

      If you are using building HDF5 from source in conjunction with MPI, we
      recommend that your HDF5 installation be built with parallel I/O
      features. An example of configuring HDF5_ is listed below::

          CC=mpicc ./configure --enable-parallel

      You may omit ``--enable-parallel`` if you want to compile HDF5_ in serial.

.. admonition:: Optional
   :class: note

    * libpng_ official reference PNG library

      OpenMC's built-in plotting capabilities use the libpng library to produce
      compressed PNG files. In the absence of this library, OpenMC will fallback
      to writing PPM files, which are uncompressed and only supported by select
      image viewers. libpng can be installed on Ddebian derivates with::

          sudo apt install libpng-dev

    * An MPI implementation for distributed-memory parallel runs

      To compile with support for parallel runs on a distributed-memory
      architecture, you will need to have a valid implementation of MPI
      installed on your machine. The code has been tested and is known to work
      with the latest versions of both OpenMPI_ and MPICH_. OpenMPI and/or MPICH
      can be installed on Debian derivatives with::

          sudo apt install mpich libmpich-dev
          sudo apt install openmpi-bin libopenmpi-dev

    * git_ version control software for obtaining source code

    * DAGMC_ toolkit for simulation using CAD-based geometries

      OpenMC supports particle tracking in CAD-based geometries via the Direct
      Accelerated Geometry Monte Carlo (DAGMC) toolkit (`installation
      instructions <https://svalinn.github.io/DAGMC/install/openmc.html>`_). For
      use in OpenMC, only the ``MOAB_DIR`` and ``BUILD_TALLY`` variables need to
      be specified in the CMake configuration step when building DAGMC. This
      option also allows unstructured mesh tallies on tetrahedral MOAB meshes.
      In addition to turning this option on, the path to the DAGMC installation
      should be specified as part of the ``CMAKE_PREFIX_PATH`` variable::

          cmake -DOPENMC_USE_DAGMC=on -DCMAKE_PREFIX_PATH=/path/to/dagmc/installation ..

    * MCPL_ library for reading and writing .mcpl files

      This option allows OpenMC to read and write MCPL (Monte Carlo Particle
      Lists) files instead of .h5 files for sources (external source
      distribution, k-eigenvalue source distribution, and surface sources). To
      turn this option on in the CMake configuration step, add the following
      option::

          cmake -DOPENMC_USE_MCPL=on ..

    * NCrystal_ library for defining materials with enhanced thermal neutron transport

      OpenMC supports the creation of materials from NCrystal, which replaces
      the scattering kernel treatment of ACE files with a modular, on-the-fly
      approach. OpenMC does not need any particular build option to use this,
      but NCrystal must be installed on the system. Refer to `NCrystal
      documentation
      <https://github.com/mctools/ncrystal/wiki/Get-NCrystal>`_ for how this is
      achieved.

    * libMesh_ mesh library framework for numerical simulations of partial differential equations

      This optional dependency enables support for unstructured mesh tally
      filters using libMesh meshes. Any 3D element type supported by libMesh can
      be used, but the implementation is currently restricted to collision
      estimators. In addition to turning this option on, the path to the libMesh
      installation should be specified as part of the ``CMAKE_PREFIX_PATH``
      variable::

          cmake -DOPENMC_USE_LIBMESH=on -DOPENMC_USE_MPI=on -DCMAKE_PREFIX_PATH=/path/to/libmesh/installation ..

      Note that libMesh is most commonly compiled with MPI support. If that
      is the case, then OpenMC should be compiled with MPI support as well.

.. _gcc: https://gcc.gnu.org/
.. _CMake: https://cmake.org
.. _OpenMPI: https://www.open-mpi.org
.. _MPICH: https://www.mpich.org
.. _HDF5: https://www.hdfgroup.org/solutions/hdf5/
.. _DAGMC: https://svalinn.github.io/DAGMC/index.html
.. _MOAB: https://bitbucket.org/fathomteam/moab
.. _libMesh: https://libmesh.github.io/
.. _libpng: http://www.libpng.org/pub/png/libpng.html
.. _MCPL: https://github.com/mctools/mcpl
.. _NCrystal: https://github.com/mctools/ncrystal

Obtaining the Source
--------------------

All OpenMC source code is hosted on GitHub_. You can download the source code
directly from GitHub or, if you have the git_ version control software installed
on your computer, you can use git to obtain the source code. The latter method
has the benefit that it is easy to receive updates directly from the GitHub
repository. GitHub has a good set of `instructions
<https://docs.github.com/en/github/getting-started-with-github/set-up-git>`_ for
how to set up git to work with GitHub since this involves setting up ssh_ keys.
With git installed and setup, the following command will download the full
source code from the GitHub repository::

    git clone --recurse-submodules https://github.com/openmc-dev/openmc.git

By default, the cloned repository will be set to the development branch. To
switch to the source of the latest stable release, run the following commands::

    cd openmc
    git checkout master

.. _GitHub: https://github.com/openmc-dev/openmc
.. _git: https://git-scm.com
.. _ssh: https://en.wikipedia.org/wiki/Secure_Shell

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

OPENMC_ENABLE_COVERAGE
  Compile and link code instrumented for coverage analysis. This is typically
  used in conjunction with gcov_. (Default: off)

OPENMC_ENABLE_PROFILE
  Enables profiling using the GNU profiler, gprof. (Default: off)

OPENMC_USE_OPENMP
  Enables shared-memory parallelism using the OpenMP API. The C++ compiler
  being used must support OpenMP. (Default: on)

OPENMC_USE_DAGMC
  Enables use of CAD-based DAGMC_ geometries and MOAB_ unstructured mesh
  tallies. Please see the note about DAGMC in the optional dependencies list
  for more information on this feature. The installation directory for DAGMC
  should also be defined as `DAGMC_ROOT` in the CMake configuration command.
  (Default: off)

OPENMC_USE_MCPL
  Turns on support for reading MCPL_ source files and writing MCPL source points
  and surface sources. (Default: off)

OPENMC_USE_LIBMESH
  Enables the use of unstructured mesh tallies with libMesh_. (Default: off)

OPENMC_USE_MPI
  Turns on compiling with MPI (Default: off). For further information on MPI
  options, please see the `FindMPI.cmake documentation
  <https://cmake.org/cmake/help/latest/module/FindMPI.html>`_.

To set any of these options (e.g., turning on profiling), the following form
should be used:

.. code-block:: sh

    cmake -DOPENMC_ENABLE_PROFILE=on /path/to/openmc

.. _gcov: https://gcc.gnu.org/onlinedocs/gcc/Gcov.html

.. _usersguide_compile_mpi:

Specifying the Build Type
+++++++++++++++++++++++++

OpenMC can be configured for debug, release, or release with debug info by setting
the `CMAKE_BUILD_TYPE` option.

Debug
  Enable debug compiler flags with no optimization. On most platforms/compilers,
  this is equivalent to `-O0 -g`.

Release
  Disable debug and enable optimization. On most platforms/compilers, this is
  equivalent to `-O3 -DNDEBUG`.

RelWithDebInfo
  (Default if no type is specified.) Enable optimization and debug. On most
  platforms/compilers, this is equivalent to `-O2 -g`.

Example of configuring for Debug mode:

.. code-block:: sh

    cmake -DCMAKE_BUILD_TYPE=Debug /path/to/openmc

Selecting HDF5 Installation
+++++++++++++++++++++++++++

CMakeLists.txt searches for the ``h5cc`` or ``h5pcc`` HDF5 C wrapper on
your PATH environment variable and subsequently uses it to determine library
locations and compile flags. If you have multiple installations of HDF5 or one
that does not appear on your PATH, you can set the HDF5_ROOT environment
variable to the root directory of the HDF5 installation, e.g.

.. code-block:: sh

    export HDF5_ROOT=/opt/hdf5/1.8.15
    cmake /path/to/openmc

This will cause CMake to search first in /opt/hdf5/1.8.15/bin for ``h5cc`` /
``h5pcc`` before it searches elsewhere. As noted above, an environment variable
can typically be set for a single command, i.e.

.. code-block:: sh

    HDF5_ROOT=/opt/hdf5/1.8.15 cmake /path/to/openmc

.. _compile_linux:

Compiling on Linux and macOS
----------------------------

To compile OpenMC on Linux or macOS, run the following commands from within the
root directory of the source code:

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

Compiling on Windows
--------------------

Recent versions of Windows include a subsystem for Linux that allows one to run
Bash within Ubuntu running in Windows. First, follow the installation guide
`here <https://learn.microsoft.com/en-us/windows/wsl/install>`_ to get Bash on
Ubuntu on Windows set up. Once you are within bash, obtain the necessary
:ref:`prerequisites <prerequisites>` via ``apt``. Finally, follow the
:ref:`instructions for compiling on linux <compile_linux>`.

Testing Build
-------------

To run the test suite, you will first need to download a pre-generated cross
section library along with windowed multipole data. Please refer to our
:ref:`devguide_tests` documentation for further details.

---------------------
Installing Python API
---------------------

If you installed OpenMC using :ref:`Conda <install_conda>`, no further steps are
necessary in order to use OpenMC's :ref:`Python API <pythonapi>`. However, if
you are :ref:`installing from source <install_source>`, the Python API is not
installed by default when ``make install`` is run because in many situations it
doesn't make sense to install a Python package in the same location as the
``openmc`` executable (for example, if you are installing the package into a
`virtual environment <https://docs.python.org/3/tutorial/venv.html>`_). The
easiest way to install the :mod:`openmc` Python package is to use pip_, which is
included by default in Python 3.4+. From the root directory of the OpenMC
distribution/repository, run:

.. code-block:: sh

    python -m pip install .

pip will first check that all :ref:`required third-party packages
<usersguide_python_prereqs>` have been installed, and if they are not present,
they will be installed by downloading the appropriate packages from the Python
Package Index (`PyPI <https://pypi.org/>`_).

Installing in "Development" Mode
--------------------------------

If you are primarily doing development with OpenMC, it is strongly recommended
to install the Python package in :ref:`"editable" mode <devguide_editable>`.

.. _usersguide_python_prereqs:

Prerequisites
-------------

The Python API works with Python 3.8+. In addition to Python itself, the API
relies on a number of third-party packages. All prerequisites can be installed
using Conda_ (recommended), pip_, or through the package manager in most Linux
distributions.

.. admonition:: Required
   :class: error

   `NumPy <https://numpy.org/>`_
      NumPy is used extensively within the Python API for its powerful
      N-dimensional array.

   `SciPy <https://www.scipy.org/>`_
      SciPy's special functions, sparse matrices, and spatial data structures
      are used for several optional features in the API.

   `pandas <https://pandas.pydata.org/>`_
      Pandas is used to generate tally DataFrames as demonstrated in an `example
      notebook
      <https://nbviewer.jupyter.org/github/openmc-dev/openmc-notebooks/blob/main/pandas-dataframes.ipynb>`_.

   `h5py <http://www.h5py.org/>`_
      h5py provides Python bindings to the HDF5 library. Since OpenMC outputs
      various HDF5 files, h5py is needed to provide access to data within these
      files from Python.

   `Matplotlib <https://matplotlib.org/>`_
      Matplotlib is used to providing plotting functionality in the API like the
      :meth:`Universe.plot` method and the :func:`openmc.plot_xs` function.

   `uncertainties <https://pythonhosted.org/uncertainties/>`_
      Uncertainties are used for decay data in the :mod:`openmc.data` module.

   `lxml <https://lxml.de/>`_
      lxml is used for various parts of the Python API.

.. admonition:: Optional
   :class: note

   `mpi4py <https://mpi4py.readthedocs.io/en/stable/>`_
      mpi4py provides Python bindings to MPI for running distributed-memory
      parallel runs. This package is needed if you plan on running depletion
      simulations in parallel using MPI.

   `vtk <https://vtk.org/>`_
      The Python VTK bindings are needed to convert voxel and track files to VTK
      format.

   `pytest <https://docs.pytest.org>`_
      The pytest framework is used for unit testing the Python API.

If you are running simulations that require OpenMC's Python bindings to the C
API (including depletion and CMFD), it is recommended to build ``h5py`` (and
``mpi4py``, if you are using MPI) using the same compilers and HDF5 version as
for OpenMC. Thus, the install process would proceed as follows:

.. code-block:: sh

    mkdir build && cd build
    HDF5_ROOT=<path to HDF5> CXX=<path to mpicxx> cmake ..
    make
    make install

    cd ..
    MPICC=<path to mpicc> python -m pip install mpi4py
    HDF5_DIR=<path to HDF5> python -m pip install --no-binary=h5py h5py

If you are using parallel HDF5, you'll also need to make sure the right MPI
wrapper is used when installing h5py:

.. code-block:: sh

    CC=<path to mpicc> HDF5_MPI=ON HDF5_DIR=<path to HDF5> python -m pip install --no-binary=h5py h5py

.. _Conda: https://conda.io/en/latest/
.. _pip: https://pip.pypa.io/en/stable/
