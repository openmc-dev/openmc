.. _usersguide_install:

==============================
Installation and Configuration
==============================

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

    sudo apt-get update

Now OpenMC should be recognized within the repository and can be installed:

.. code-block:: sh

    sudo apt-get install openmc

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

    * A Fortran compiler such as gfortran_

      In order to compile OpenMC, you will need to have a Fortran compiler
      installed on your machine. Since a number of Fortran 2003/2008 features
      are used in the code, it is recommended that you use the latest version of
      whatever compiler you choose. For gfortran_, it is necessary to use
      version 4.6.0 or above.

      If you are using Debian or a Debian derivative such as Ubuntu, you can
      install the gfortran compiler using the following command::

          sudo apt-get install gfortran

    * CMake_ cross-platform build system

      The compiling and linking of source files is handled by CMake in a
      platform-independent manner. If you are using Debian or a Debian
      derivative such as Ubuntu, you can install CMake using the following
      command::

          sudo apt-get install cmake

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

          sudo apt-get install libhdf5-8 libhdf5-dev hdf5-helpers

      Note that the exact package names may vary depending on your particular
      distribution and version.

.. admonition:: Optional

    * An MPI implementation for distributed-memory parallel runs

      To compile with support for parallel runs on a distributed-memory
      architecture, you will need to have a valid implementation of MPI
      installed on your machine. The code has been tested and is known to work
      with the latest versions of both OpenMPI_ and MPICH_. OpenMPI and/or MPICH
      can be installed on Debian derivatives with::

          sudo apt-get install mpich libmpich-dev
          sudo apt-get install openmpi-bin libopenmpi1.6 libopenmpi-dev

    * git_ version control software for obtaining source code

.. _gfortran: http://gcc.gnu.org/wiki/GFortran
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
  being used must support OpenMP.

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

Compiling with MPI
++++++++++++++++++

To compile with MPI, set the :envvar:`FC` environment variable to the path to
the MPI Fortran wrapper. For example, in a bash shell:

.. code-block:: sh

    export FC=mpif90
    cmake /path/to/openmc

Note that in many shells, an environment variable can be set for a single
command, i.e.

.. code-block:: sh

    FC=mpif90 cmake /path/to/openmc

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

Compiling on Windows
--------------------

Using Cygwin
++++++++++++

One option for compiling OpenMC on a Windows operating system is to use Cygwin_,
a Linux-like environment for Windows. You will need to first `install
Cygwin`_. When you are asked to select packages, make sure the following are
selected:

* Devel: gcc-core
* Devel: gcc-fortran
* Devel: make
* Devel: cmake

If you plan on obtaining the source code directly using git, select the
following packages:

* Devel: git
* Devel: git-completion (Optional)
* Devel: gitk (Optional)

In order to use the Python scripts provided with OpenMC, you will also need to
install Python. This can be done within Cygwin or directly in Windows. To
install within Cygwin, select the following packages:

* Python: python (Version > 2.7 recommended)

Once you have obtained the source code, run the following commands from within
the source code root directory:

.. code-block:: sh

    mkdir build && cd build
    cmake ..
    make

This will build an executable named ``openmc``.

.. _Cygwin: http://cygwin.com/
.. _install Cygwin: http://cygwin.com/setup.exe

Using MinGW
+++++++++++

An alternate option for installing OpenMC on Windows is using MinGW_, which
stands for Minimalist GNU for Windows. An executable for installing the MinGW
distribution is available on SourceForge_. When installing MinGW, make sure the
following components are selected:

* MinGW Compiler Suite: Fortran Compiler
* MSYS Basic System

Once MinGW is installed, copy the OpenMC source distribution to your MinGW home
directory (usually C:\\MinGW\\msys\\1.0\\home\\YourUsername). Once you have
the source code in place, run the following commands from within the MinGW shell
in the root directory of the OpenMC distribution:

.. code-block:: sh

    make

This will build an executable named ``openmc``.

.. _MinGW: http://www.mingw.org
.. _SourceForge: http://sourceforge.net/projects/mingw

Compiling for the Intel Xeon Phi
--------------------------------

In order to build OpenMC for the Intel Xeon Phi using the Intel Fortran
compiler, it is necessary to specify that all objects be compiled with the
``-mmic`` flag as follows:

.. code-block:: sh

    mkdir build && cd build
    FC=ifort FFLAGS=-mmic cmake -Dopenmp=on ..
    make

Note that unless an HDF5 build for the Intel Xeon Phi is already on your target
machine, you will need to cross-compile HDF5 for the Xeon Phi. An `example
script`_ to build zlib and HDF5 provides several necessary workarounds.

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

---------------------------
Cross Section Configuration
---------------------------

In order to run a simulation with OpenMC, you will need cross section data for
each nuclide or material in your problem. OpenMC can be run in
continuous-energy or multi-group mode.

In continuous-energy mode OpenMC uses ACE format cross sections; in this case
you can use nuclear data that was processed with NJOY_, such as that
distributed with MCNP_ or Serpent_.  Several sources provide free processed
ACE data as described below. The TALYS-based evaluated nuclear data library,
TENDL_, is also openly available in ACE format.

In multi-group mode, OpenMC utilizes an XML-based library format which can be
used to describe nuclide- or material-specific quantities.

Using ENDF/B-VII.1 Cross Sections from NNDC
-------------------------------------------

The NNDC_ provides ACE data from the ENDF/B-VII.1 neutron and thermal scattering
sublibraries at four temperatures processed using NJOY_. To use this data with
OpenMC, a script is provided with OpenMC that will automatically download,
extract, and set up a confiuration file:

.. code-block:: sh

    cd openmc/data
    python get_nndc_data.py

At this point, you should set the :envvar:`OPENMC_CROSS_SECTIONS` environment
variable to the absolute path of the file
``openmc/data/nndc/cross_sections.xml``. This cross section set is used by the
test suite.

Using JEFF Cross Sections from OECD/NEA
---------------------------------------

The NEA_ provides processed ACE data from the JEFF_ nuclear library upon
request. A DVD of the data can be requested here_. To use this data with OpenMC,
the following steps must be taken:

1. Copy and unzip the data on the DVD to a directory on your computer.
2. In the root directory, a file named ``xsdir``, or some variant thereof,
   should be present. This file contains a listing of all the cross sections and
   is used by MCNP. This file should be converted to a ``cross_sections.xml``
   file for use with OpenMC. A utility is provided in the OpenMC distribution
   for this purpose:

   .. code-block:: sh

       openmc/scripts/openmc-xsdir-to-xml xsdir31 cross_sections.xml

3. In the converted ``cross_sections.xml`` file, change the contents of the
   <directory> element to the absolute path of the directory containing the
   actual ACE files.
4. Additionally, you may need to change any occurrences of upper-case "ACE"
   within the ``cross_sections.xml`` file to lower-case.
5. Either set the :ref:`cross_sections` in a settings.xml file or the
   :envvar:`OPENMC_CROSS_SECTIONS` environment variable to the absolute path of
   the ``cross_sections.xml`` file.

Using Cross Sections from MCNP
------------------------------

To use cross sections distributed with MCNP, change the <directory> element in
the ``cross_sections.xml`` file in the root directory of the OpenMC distribution
to the location of the MCNP cross sections. Then, either set the
:ref:`cross_sections` in a settings.xml file or the
:envvar:`OPENMC_CROSS_SECTIONS` environment variable to the absolute path of
the ``cross_sections.xml`` file.

Using Cross Sections from Serpent
---------------------------------

To use cross sections distributed with Serpent, change the <directory> element
in the ``cross_sections_serpent.xml`` file in the root directory of the OpenMC
distribution to the location of the Serpent cross sections. Then, either set the
:ref:`cross_sections` in a settings.xml file or the
:envvar:`OPENMC_CROSS_SECTIONS` environment variable to the absolute path of
the ``cross_sections_serpent.xml``
file.

Using Multi-Group Cross Sections
--------------------------------

Multi-group cross section libraries are generally tailored to the specific
calculation to be performed.  Therefore, at this point in time, OpenMC is not
distributed with any pre-existing multi-group cross section libraries.
However, if the user has obtained or generated their own library, the user
should set the :envvar:`OPENMC_MG_CROSS_SECTIONS` environment variable
to the absolute path of the file library expected to used most frequently.

.. _NJOY: http://t2.lanl.gov/nis/codes.shtml
.. _NNDC: http://www.nndc.bnl.gov/endf/b7.1/acefiles.html
.. _NEA: http://www.oecd-nea.org
.. _JEFF: http://www.oecd-nea.org/dbdata/jeff/
.. _here: http://www.oecd-nea.org/dbdata/pubs/jeff312-cd.html
.. _MCNP: http://mcnp.lanl.gov
.. _Serpent: http://montecarlo.vtt.fi
.. _TENDL: ftp://ftp.nrg.eu/pub/www/talys/tendl2012/tendl2012.html

--------------
Running OpenMC
--------------

Once you have a model built (see :ref:`usersguide_input`), you can either run
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
working directory which may be different from ``/home/username/somemodel``.

Command-Line Flags
------------------

OpenMC accepts the following command line flags:

-g, --geometry-debug   Run in geometry debugging mode, where cell overlaps are
                       checked for after each move of a particle
-n, --particles N      Use *N* particles per generation or batch
-p, --plot             Run in plotting mode
-r, --restart file     Restart a previous run from a state point or a particle
                       restart file
-s, --threads N        Run with *N* OpenMP threads
-t, --track            Write tracks for all particles
-v, --version          Show version information

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
.. _ctest: http://www.cmake.org/cmake/help/v2.8.12/ctest.html
