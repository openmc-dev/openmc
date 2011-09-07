.. _install:

============
Installation
============

-------------
Prerequisites
-------------

In order to compile OpenMC, you will need to have a Fortran compiler installed
on your machine. Since a number of Fortran 2003 features are used in the code,
it is recommended that you use the latest version of whatever compiler you
choose. For gfortran_, it is recommended that you use version 4.5.0 or above.

To compile with support for parallel runs on a distributed-memory architecture,
you will need to have a valid implementation of MPI installed on your
machine. The code has been tested and is known to work with the latest versions
of both OpenMPI_ and MPICH2_. You may use older versions of MPI implementations at
your own risk.

.. _gfortran: http://gcc.gnu.org/wiki/GFortran
.. _OpenMPI: http://www.open-mpi.org
.. _MPICH2: http://www.mcs.anl.gov/mpi/mpich/

-------------
Configuration
-------------

All configuration for OpenMC is done within the Makefile located
``src/Makefile``. In the Makefile, you will see that there are a number of User
Options which can be changed. It is recommended that you do not change anything
else in the Makefile unless you are experienced with compiling and building
software using Makefiles. The following parameters can be set from the User
Options sections in the Makefile:

COMPILER
  This variable tells the Makefile which compiler to use. Valid options are
  gfortran, intel, and pgi.

DEBUG
  Enables debugging when compiling. The flags added are dependent on which
  compiler is used.

PROFILE
  Enables profiling using the GNU profiler, gprof.

OPTIMIZE
  Enables high-optimization using compiler-dependent flags. For gfortran,
  this compiles with -O3. For Intel Fortran, this compiles with -O3 as well as
  interprocedural optimization.

USE_MPI
  Enables parallel runs using the Message Passing Interface. Users should
  also set the MPI root directory by setting the MPI variable further down in
  the Makefile.

USE_MPI
  Enables parallel runs on shared-memory architecture using OpenMP threading.

USE_COARRAY
  Enables parallel runs using Fortran 2008 coarrays.

---------
Compiling
---------

To compile the code, run the following commands from within the root directory
for OpenMC::

    cd openmc/src
    make -f Makefile

This will build an executable named ``openmc``.
