.. _quickinstall:

===================
Quick Install Guide
===================

This quick install guide outlines the basic steps needed to install OpenMC on
your computer. For more detailed instructions on configuring and installing
OpenMC, see :ref:`usersguide_install` in the User's Manual.

-------------
Prerequisites
-------------

In order to compile and run OpenMC, a number of prerequisite software packages
and libraries may be needed. These include:

- A Fortran compiler such as gfortran_
- git_ version control software for obtaining source code (*optional*)
- An MPI implementation for parallel runs (*optional*)
- HDF5_ Library for improved output format (*optional*)

.. _gfortran: http://gcc.gnu.org/wiki/GFortran
.. _git: http://git-scm.com
.. _HDF5: http://www.hdfgroup.org/HDF5/

--------------------
Obtaining the Source
--------------------

All OpenMC source code is hosted on GitHub_. The latest version of the OpenMC
source code can be downloaded from GitHub with the following command::

    git clone git@github.com:mit-crpg/openmc.git

If you don't have git installed, you can download the source as a zip or tar
file directly from GitHub.

.. _GitHub: https://github.com/mit-crpg/openmc

-------------------------------
Compiling on Linux and Mac OS X
-------------------------------

To compile OpenMC on Linux or Max OS X, run the following commands from within
the root directory of the source code:

.. code-block:: sh

    cd src
    make

This will build an executable named ``openmc``.

--------------------
Compiling on Windows
--------------------

To compile OpenMC on a Windows operating system, you will need to first install
Cygwin_, a Linux-like environment for Windows. When configuring Cygwin, make
sure you install both the gfortran compiler as well as git. Once you have
obtained the source code, run the following commands from within the source code
root directory:

.. code-block:: sh

    cd src
    make

This will build an executable named ``openmc``.

.. _Cygwin: http://www.cygwin.com/
