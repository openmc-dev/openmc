.. _quickinstall:

===================
Quick Install Guide
===================

This quick install guide outlines the basic steps needed to install OpenMC on
your computer. For more detailed instructions on configuring and installing
OpenMC, see :ref:`usersguide_install` in the User's Manual.

--------------------------------
Installing on Ubuntu through PPA
--------------------------------

For users with Ubuntu 11.10 or later, a binary package for OpenMC is available
through a `Personal Package Archive`_ (PPA) and can be installed through the `APT
package manager`_. Simply enter the following commands into the terminal:

.. code-block:: sh

    sudo apt-add-repository ppa:paulromano/staging
    sudo apt-get update
    sudo apt-get install openmc

Currently, the binary package does not allow for parallel simulations, HDF5_, or
CMFD acceleration through PETSc_. Users who need such capabilities should build
OpenMC from source as is described in :ref:`usersguide_install`.

.. _Personal Package Archive: https://launchpad.net/~paulromano/+archive/staging
.. _APT package manager: https://help.ubuntu.com/community/AptGet/Howto
.. _HDF5: http://www.hdfgroup.org/HDF5/
.. _PETSc: http://www.mcs.anl.gov/petsc/

-------------------------------------------
Installing from Source on Linux or Mac OS X
-------------------------------------------

All OpenMC source code is hosted on GitHub_. If you have git_ and the gfortran_
compiler installed, you can download and install OpenMC be entering the
following commands in a terminal:

.. code-block:: sh

    git clone git://github.com/mit-crpg/openmc.git
    cd openmc/src
    git checkout master
    make
    sudo make install

This will build an executable named ``openmc`` and install it (by default in
/usr/local/bin). If you do not have administrator privileges, the last command
can be replaced with a local install, e.g.

.. code-block:: sh

    make install -e prefix=$HOME/.local

.. _GitHub: https://github.com/mit-crpg/openmc
.. _git: http://git-scm.com
.. _gfortran: http://gcc.gnu.org/wiki/GFortran
