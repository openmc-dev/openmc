==========================================
OpenMC Monte Carlo Particle Transport Code
==========================================

The OpenMC project aims to provide a fully-featured Monte Carlo particle
transport code based on modern methods. The project started under the
Computational Reactor Physics Group at MIT.

If you are interested in the project or would like to help and contribute,
please contact `Paul Romano`_.

---------------
Building OpenMC
---------------

Building the executable for OpenMC is as easy as going to the src directory and
running make. By default, the Makefile is set to use the GNU Fortran compiler,
gfortran. This can be changed by the F90 variable in the Makefile.

Compiling for multiple processors can be controlled with the USE_MPI and MPI
variables. The MPI variable should be set to the base directory of the MPI
implementation installed on your computer.

OpenMC has been tested on Linux, Mac, and Windows platforms with Intel, GNU,
IBM, Cray, and PGI compilers. It is recommended to use the latest version of
whatever compiler you should choose to use as there are a number of Fortran 2003
and 2008 intrinsics that are used in OpenMC.

--------------
Running OpenMC
--------------

To run OpenMC after building, you'll need to write an input file according to
the specifications in this documentation. With your input file complete, simply
run::

	openmc inputFile

where *inputFile* is the name of your input file. To run in parallel, use the
mpiexec or mpirun script provided by your MPI implementation. Assuming this is
on your PATH, you may run::

   mpiexec -n numProcs openmc inputFile

where *numProcs* is the number of processors you desire to run on.

.. _Paul Romano: mailto:paul.k.romano@gmail.com
