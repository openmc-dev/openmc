.. _quickinstall:

===================
Quick Install Guide
===================

This quick install guide outlines the basic steps needed to install OpenMC on
your computer. For more detailed instructions on configuring and installing
OpenMC, see :ref:`usersguide_install` in the User's Manual.

----------------------------------------
Installing on Linux/Mac with conda-forge
----------------------------------------

`Conda <http://conda.pydata.org/docs/>`_ is an open source package management
system and environment management system for installing multiple versions of
software packages and their dependencies and switching easily between them. If
you have `conda` installed on your system, OpenMC can be installed via the
`conda-forge` channel. First, add the `conda-forge` channel with:

.. code-block:: sh

    conda config --add channels conda-forge

OpenMC can then be installed with:

.. code-block:: sh

    conda install openmc

--------------------------------
Installing on Ubuntu through PPA
--------------------------------

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

-------------------------------------------
Installing on Linux/Mac/Windows with Docker
-------------------------------------------

OpenMC can be easily deployed using `Docker <https://www.docker.com/>`_ on any
Windows, Mac or Linux system. With Docker running, execute the following
command in the shell to download and run a `Docker image`_ with the most recent release of OpenMC from `DockerHub <https://hub.docker.com/>`_ called ``openmc/openmc:v0.10.0``:

.. code-block:: sh

    docker run openmc/openmc:v0.10.0

This will take several minutes to run depending on your internet download speed. The command will place you in an interactive shell running in a `Docker container`_ with OpenMC installed.

.. note:: The ``docker run`` command supports many `options`_ for spawning
          containers -- including `mounting volumes`_ from the host
	  filesystem -- which many users will find useful.

.. _Docker image: https://docs.docker.com/engine/reference/commandline/images/
.. _Docker container: https://www.docker.com/resources/what-container
.. _options: https://docs.docker.com/engine/reference/commandline/run/
.. _mounting volumes: https://docs.docker.com/storage/volumes/

---------------------------------------
Installing from Source on Ubuntu 15.04+
---------------------------------------

To build OpenMC from source, several :ref:`prerequisites <prerequisites>` are
needed. If you are using Ubuntu 15.04 or higher, all prerequisites can be
installed directly from the package manager.

.. code-block:: sh

    sudo apt install g++ cmake libhdf5-dev

After the packages have been installed, follow the instructions below for
building and installing OpenMC from source.

-------------------------------------------
Installing from Source on Linux or Mac OS X
-------------------------------------------

All OpenMC source code is hosted on `GitHub
<https://github.com/openmc-dev/openmc>`_. If you have `git
<https://git-scm.com>`_, the `gcc <https://gcc.gnu.org/>`_ compiler suite,
`CMake <http://www.cmake.org>`_, and `HDF5 <https://www.hdfgroup.org/HDF5/>`_
installed, you can download and install OpenMC be entering the following
commands in a terminal:

.. code-block:: sh

    git clone --recurse-submodules https://github.com/openmc-dev/openmc.git
    cd openmc
    mkdir build && cd build
    cmake ..
    make
    sudo make install

This will build an executable named ``openmc`` and install it (by default in
/usr/local/bin). If you do not have administrator privileges, the cmake command
should specify an installation directory where you have write access, e.g.

.. code-block:: sh

    cmake -DCMAKE_INSTALL_PREFIX=$HOME/.local ..

The :mod:`openmc` Python package must be installed separately. The easiest way
to install it is using `pip <https://pip.pypa.io/en/stable/>`_, which is
included by default in Python 2.7 and Python 3.4+. From the root directory of
the OpenMC distribution/repository, run:

.. code-block:: sh

    pip install .

If you want to build a parallel version of OpenMC (using OpenMP or MPI),
directions can be found in the :ref:`detailed installation instructions
<usersguide_build>`.
