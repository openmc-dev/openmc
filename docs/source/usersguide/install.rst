.. _usersguide_install:

==============================
Installation and Configuration
==============================

-----------------------------
Installing on Ubuntu with PPA
-----------------------------

For users with Ubuntu 11.10 or later, a binary package for OpenMC is available
through a Personal Package Archive (PPA) and can be installed through the APT
package manager. First, add the following PPA to the repository sources:

.. code-block:: sh

    sudo apt-add-repository ppa:paulromano/staging

Next, resynchronize the package index files:

.. code-block:: sh

    sudo apt-get update

Now OpenMC should be recognized within the repository and can be installed:

.. code-block:: sh

    sudo apt-get install openmc

--------------------
Building from Source
--------------------

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

.. admonition:: Optional

    * An MPI implementation for distributed-memory parallel runs

      To compile with support for parallel runs on a distributed-memory
      architecture, you will need to have a valid implementation of MPI
      installed on your machine. The code has been tested and is known to work
      with the latest versions of both OpenMPI_ and MPICH_. Note that if using
      OpenMPI, make sure that --with-mpi-f90-size is not set to medium or large
      since this may prevent MPI calls from completing successfully in
      OpenMC. OpenMPI and/or MPICH can be installed on Debian derivatives
      with::

          sudo apt-get install mpich2 libmpich2-dev
          sudo apt-get install openmpi1.6-bin libopenmpi1.6-dev

    * HDF5_ Library for portable binary output format

      To compile with support for HDF5_ output (highly recommended), you will
      need to have HDF5 installed on your computer. The installed version will
      need to have been compiled with the same compiler you intend to compile
      OpenMC with. HDF5_ must be built with parallel I/O features if you intend
      to use HDF5_ with MPI. An example of configuring HDF5_ is listed below::

           FC=/opt/mpich/3.0.4-gnu/bin/mpif90 CC=/opt/mpich/3.0.4-gnu/bin/mpicc \
           ./configure --prefix=/opt/hdf5/1.8.11-gnu --enable-fortran \
                       --enable-fortran2003 --enable-parallel

      You may omit ``--enable-parallel`` if you want to compile HDF5_ in serial.

    * PETSc_ for CMFD acceleration

      To enable CMFD acceleration, you will need to have PETSc_ (3.4.2 or higher)
      installed on your computer. The installed version will need to have been
      compiled with the same compiler you intend to compile OpenMC with. OpenMC
      requires PETSc_ to be configured with Fortran datatypes. An example of
      configuring PETSc_ is listed below::

           ./configure --prefix=/opt/petsc/3.4.2-gnu --download-f-blas-lapack \
                       --with-mpi-dir=/opt/mpich/3.0.4-gnu/ --with-shared-libraries=0 \
                       --with-fortran-datatypes

      The BLAS/LAPACK library is not required to be downloaded and can be linked
      explicitly (e.g., Intel MKL library).

    * git_ version control software for obtaining source code

.. _gfortran: http://gcc.gnu.org/wiki/GFortran
.. _OpenMPI: http://www.open-mpi.org
.. _MPICH: http://www.mpich.org
.. _HDF5: http://www.hdfgroup.org/HDF5/
.. _PETSc: http://www.mcs.anl.gov/petsc/

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

    git clone git://github.com/mit-crpg/openmc.git

By default, the cloned repository will be set to the development branch. To
switch to the source of the latest stable release, run the following commands::

    cd openmc/src
    git checkout master

.. _GitHub: https://github.com/mit-crpg/openmc
.. _git: http://git-scm.com
.. _ssh: http://en.wikipedia.org/wiki/Secure_Shell

Build Configuration
-------------------

All configuration for OpenMC is done within the Makefile located in
``src/Makefile``. In the Makefile, you will see that there are a number of User
Options which can be changed. It is recommended that you do not change anything
else in the Makefile unless you are experienced with compiling and building
software using Makefiles. The following parameters can be set from the User
Options sections in the Makefile:

COMPILER
  This variable tells the Makefile which compiler to use. Valid options are
  gnu, intel, pgi, ibm, and cray. The default is gnu (gfortran).

DEBUG
  Enables debugging when compiling. The flags added are dependent on which
  compiler is used.

PROFILE
  Enables profiling using the GNU profiler, gprof.

OPTIMIZE
  Enables high-optimization using compiler-dependent flags. For gfortran and
  Intel Fortran, this compiles with -O3.

MPI
  Enables parallel runs using the Message Passing Interface. The MPI_DIR
  variable should be set to the base directory of the MPI implementation.

OPENMP
  Enables shared-memory parallelism using the OpenMP API. The Fortran compiler
  being used must support OpenMP.

HDF5
  Enables HDF5 output in addition to normal screen and text file output. The
  HDF5_DIR variable should be set to the base directory of the HDF5
  installation.

PETSC
  Enables PETSc for use in CMFD acceleration. The PETSC_DIR variable should be
  set to the base directory of the PETSc installation.

It is also possible to change these options from the command line itself. For
example, if you want to compile with DEBUG turned on without actually change the
Makefile, you can enter the following from a terminal::

    make DEBUG=yes

Compiling on Linux and Mac OS X
-------------------------------

To compile OpenMC on Linux or Max OS X, run the following commands from within
the root directory of the source code:

.. code-block:: sh

    cd src
    make
    sudo make install

This will build an executable named ``openmc`` and install it (by default in
/usr/local/bin). If you do not have administrative privileges, you can install
OpenMC locally by replacing the last command with:

.. code-block:: sh

    make install -e prefix=$HOME/.local

The ``prefix`` variable can be changed to any path for which you have
write-access.

Compiling on Windows
--------------------

Using Cygwin
++++++++++++

One option for compiling OpenMC on a Windows operating system is to use Cygwin_,
a Linux-like environment for Windows. You will need to first `install
Cygwin`_. When you are asked to select packages, make sure the following are
selected:

* Devel: gcc4-core
* Devel: gcc4-fortran
* Devel: make

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

    cd src
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

    cd src
    make

This will build an executable named ``openmc``.

.. _MinGW: http://www.mingw.org
.. _SourceForge: http://sourceforge.net/projects/mingw

---------------------------
Cross Section Configuration
---------------------------

In order to run a simulation with OpenMC, you will need cross section data for
each nuclide in your problem. Since OpenMC uses ACE format cross sections, you
can use nuclear data that was processed with NJOY_, such as that distributed
with MCNP_ or Serpent_. Several sources provide free processed ACE data as
described below. The TALYS-based evaluated nuclear data library, TENDL_, is also
openly available in ACE format.

Using ENDF/B-VII.1 Cross Sections from NNDC
-------------------------------------------

The NNDC_ provides ACE data from the ENDF/B-VII.1 neutron and thermal scattering
sublibraries at four temperatures processed using NJOY_. To use this data with
OpenMC, a script is provided with OpenMC that will automatically download,
extract, and set up a confiuration file:

.. code-block:: sh

    cd openmc/data
    python get_nndc_data.py

At this point, you should set the :envvar:`CROSS_SECTIONS` environment variable
to the absolute path of the file ``openmc/data/nndc/cross_sections.xml``.

Using JEFF Cross Sections from OECD/NEA
---------------------------------------

The NEA_ provides processed ACE data from the JEFF_ nuclear library upon
request. A DVD of the data can be requested here_. To use this data with OpenMC,
the following steps must be taken:

1. Copy and unzip the data on the DVD to a directory on your computer.
2. In the root directory, a file named ``xsdir``, or some variant thereof,
   should be present. This file contains a listing of all the cross sections and
   is used by MCNP. This file should be converted to a ``cross_sections.xml``
   file for use with OpenMC. A Python script is provided in the OpenMC
   distribution for this purpose:

   .. code-block:: sh

       openmc/src/utils/convert_xsdir.py xsdir31 cross_sections.xml

3. In the converted ``cross_sections.xml`` file, change the contents of the
   <directory> element to the absolute path of the directory containing the
   actual ACE files.
4. Additionally, you may need to change any occurrences of upper-case "ACE"
   within the ``cross_sections.xml`` file to lower-case.
5. Either set the :ref:`cross_sections` in a settings.xml file or the
   :envvar:`CROSS_SECTIONS` environment variable to the absolute path of the
   ``cross_sections.xml`` file.

Using Cross Sections from MCNP
------------------------------

To use cross sections distributed with MCNP, change the <directory> element in
the ``cross_sections.xml`` file in the root directory of the OpenMC distribution
to the location of the MCNP cross sections. Then, either set the
:ref:`cross_sections` in a settings.xml file or the :envvar:`CROSS_SECTIONS`
environment variable to the absolute path of the ``cross_sections.xml`` file.

Using Cross Sections from Serpent
---------------------------------

To use cross sections distributed with Serpent, change the <directory> element
in the ``cross_sections_serpent.xml`` file in the root directory of the OpenMC
distribution to the location of the Serpent cross sections. Then, either set the
:ref:`cross_sections` in a settings.xml file or the :envvar:`CROSS_SECTIONS`
environment variable to the absolute path of the ``cross_sections_serpent.xml``
file.

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
the XML input files. For example, if the path of your OpenMC executable is
``/home/username/openmc/src/openmc`` and your XML input files are in the
directory ``/home/username/somemodel/``, one way to run the simulation would be:

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
