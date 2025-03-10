.. _quickinstall:

===================
Quick Install Guide
===================

This quick install guide outlines the basic steps needed to install OpenMC on
your computer. For more detailed instructions on configuring and installing
OpenMC, see :ref:`usersguide_install` in the User's Manual.

----------------------------------
Installing on Linux/Mac with Conda
----------------------------------

`Conda <https://docs.conda.io/en/latest/>`_ is an open source package management
system and environments management system for installing multiple versions of
software packages and their dependencies and switching easily between them.
OpenMC can be installed in a `conda` environment. First, `conda` should be
`installed <https://www.anaconda.com/docs/getting-started/getting-started>`_
with either Anaconda Distribution or Miniconda. Once you have `conda` installed
on your system, OpenMC can be installed via the `conda-forge` channel.

First, add the `conda-forge` channel with:

.. code-block:: sh

    conda config --add channels conda-forge
    conda config --set channel_priority strict

Then create and activate a new conda enviroment called `openmc-env` (or whatever
you wish) with OpenMC installed.

.. code-block:: sh

    conda create --name openmc-env openmc
    conda activate openmc-env

You are now in a conda environment called `openmc-env` that has OpenMC
installed.

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

-------------------------------
Manually Installing from Source
-------------------------------

Obtaining prerequisites on Ubuntu
---------------------------------

When building OpenMC from source, all :ref:`prerequisites <prerequisites>` can
be installed using the package manager:

.. code-block:: sh

    sudo apt install g++ cmake libhdf5-dev libpng-dev

After the packages have been installed, follow the instructions to build from
source below.

Obtaining prerequisites on macOS
--------------------------------

For an OpenMC build with multithreading enabled, a package manager like
`Homebrew <https://brew.sh>`_ should first be installed. Then, the following
packages should be installed, for example in Homebrew via:

.. code-block:: sh

   brew install llvm cmake xtensor hdf5 python libomp libpng

The compiler provided by the above LLVM package should be used in place of the
one provisioned by XCode, which does not support the multithreading library used
by OpenMC. To ensure CMake picks up the correct compiler, make sure that either
the :envvar:`CXX` environment variable is set to the brew-installed ``clang++``
or that the directory containing it is on your :envvar:`PATH` environment
variable. Common locations for the brew-installed compiler are
``/opt/homebrew/opt/llvm/bin`` and ``/usr/local/opt/llvm/bin``.

After the packages have been installed, follow the instructions to build from
source below.

Building Source on Linux or macOS
---------------------------------

All OpenMC source code is hosted on `GitHub
<https://github.com/openmc-dev/openmc>`_. If you have `git
<https://git-scm.com>`_, a modern C++ compiler, `CMake <https://cmake.org>`_,
and `HDF5 <https://www.hdfgroup.org/solutions/hdf5/>`_ installed, you can
download and install OpenMC by entering the following commands in a terminal:

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
to install it is using `pip <https://pip.pypa.io/en/stable/>`_.
From the root directory of the OpenMC repository, run:

.. code-block:: sh

    python -m pip install .

By default, OpenMC will be built with multithreading support. To build
distributed-memory parallel versions of OpenMC using MPI or to configure other
options, directions can be found in the :ref:`detailed installation instructions
<usersguide_build>`.
