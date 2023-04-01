.. _quickinstall:

===================
Quick Install Guide
===================

This quick install guide outlines the basic steps needed to install OpenMC on
your computer. For more detailed instructions on configuring and installing
OpenMC, see :ref:`usersguide_install` in the User's Manual.

--------------------------------------------------
Installing on Linux/Mac with Mamba and conda-forge
--------------------------------------------------

`Conda <https://conda.io/en/latest/>`_ is an open source package management
system and environments management system for installing multiple versions of
software packages and their dependencies and switching easily between them.
`Mamba <https://mamba.readthedocs.io/en/latest/>`_ is a cross-platform package
manager and is compatible with `conda` packages.
OpenMC can be installed in a `conda` environment with `mamba`.
First, `conda` should be installed with one of the following installers:
`Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_,
`Anaconda <https://www.anaconda.com/>`_, or `Miniforge <https://github.com/conda-forge/miniforge>`_.
Once you have `conda` installed on your system, OpenMC can be installed via the
`conda-forge` channel with `mamba`.

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
Installing on Linux/Mac/Windows with Docker
-------------------------------------------

OpenMC can be easily deployed using `Docker <https://www.docker.com/>`_ on any
Windows, Mac, or Linux system. With Docker running, execute the following command
in the shell to download and run a `Docker image`_ with the most recent release
of OpenMC from `DockerHub <https://hub.docker.com/>`_:

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

----------------------------------
Installing from Source using Spack
----------------------------------

Spack_ is a package management tool designed to support multiple versions and
configurations of software on a wide variety of platforms and environments.
Please follow Spack's `setup guide`_ to configure the Spack system.

To install the latest OpenMC with the Python API, use the following command:

.. code-block:: sh

    spack install py-openmc

For more information about customizations including MPI, see the
:ref:`detailed installation instructions using Spack <install-spack>`.
Once installed, environment/lmod modules can be generated or Spack's `load` feature
can be used to access the installed packages.

.. _Spack: https://spack.readthedocs.io/en/latest/
.. _setup guide: https://spack.readthedocs.io/en/latest/getting_started.html

--------------------------------
Installing from Source on Ubuntu
--------------------------------

To build OpenMC from source, several :ref:`prerequisites <prerequisites>` are
needed. If you are using Ubuntu or higher, all prerequisites can be installed
directly from the package manager:

.. code-block:: sh

    sudo apt install g++ cmake libhdf5-dev libpng-dev

After the packages have been installed, follow the instructions below for
building and installing OpenMC from source.

-------------------------------------------
Installing from Source on Linux or Mac OS X
-------------------------------------------

All OpenMC source code is hosted on `GitHub
<https://github.com/openmc-dev/openmc>`_. If you have `git
<https://git-scm.com>`_, the `gcc <https://gcc.gnu.org/>`_ compiler suite,
`CMake <https://cmake.org>`_, and `HDF5
<https://www.hdfgroup.org/solutions/hdf5/>`_ installed, you can download and
install OpenMC be entering the following commands in a terminal:

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
included by default in Python 3.4+. From the root directory of the OpenMC
distribution/repository, run:

.. code-block:: sh

    pip install .

If you want to build a parallel version of OpenMC (using OpenMP or MPI),
directions can be found in the :ref:`detailed installation instructions
<usersguide_build>`.
