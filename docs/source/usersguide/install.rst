.. _usersguide_install:

==============================
Installation and Configuration
==============================

.. currentmodule:: openmc

.. _install_conda:

----------------------------------------
Installing on Linux/Mac with conda-forge
----------------------------------------

Conda_ is an open source package management system and environment management
system for installing multiple versions of software packages and their
dependencies and switching easily between them. If you have `conda` installed on
your system, OpenMC can be installed via the `conda-forge` channel. First, add
the `conda-forge` channel with:

.. code-block:: sh

    conda config --add channels conda-forge

To list the versions of OpenMC that are available on the `conda-forge` channel,
in your terminal window or an Anaconda Prompt run:

.. code-block:: sh

    conda search openmc

OpenMC can then be installed with:

.. code-block:: sh

    conda create -n openmc-env openmc

This will install OpenMC in a conda environment called `openmc-env`. To activate
the environment, run:

.. code-block:: sh

    conda activate openmc-env

.. _install_ppa:

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

    sudo apt update

Now OpenMC should be recognized within the repository and can be installed:

.. code-block:: sh

    sudo apt install openmc

Binary packages from this PPA may exist for earlier versions of Ubuntu, but they
are no longer supported.

.. _Personal Package Archive: https://launchpad.net/~paulromano/+archive/staging
.. _APT package manager: https://help.ubuntu.com/community/AptGet/Howto

.. _install-spack:

----------------------------------
Installing from Source using Spack
----------------------------------

Spack_ is a package management tool designed to support multiple versions and
configurations of software on a wide variety of platforms and environments.
Please following Spack's `setup guide`_ to configure the Spack system.

The OpenMC Spack recipe has been configured with variants that match most options
provided in the CMakeLists.txt file. To see a list of these variants and other
information:

.. code-block:: sh

    spack info openmc

    CMakePackage:   openmc

    Description:
        The OpenMC project aims to provide a fully-featured Monte Carlo particle
        transport code based on modern methods. It is a constructive solid
        geometry, continuous-energy transport code that uses ACE format cross
        sections. The project started under the Computational Reactor Physics
        Group at MIT.

    Homepage: http://openmc.readthedocs.io/

    Tags: 
        None

    Preferred version:  
        0.12.0     [git] https://github.com/openmc-dev/openmc.git at tag v0.12.0

    Safe versions:  
        develop    [git] https://github.com/openmc-dev/openmc.git on branch develop
        master     [git] https://github.com/openmc-dev/openmc.git on branch master
        0.12.0     [git] https://github.com/openmc-dev/openmc.git at tag v0.12.0
        0.11.0     https://github.com/openmc-dev/openmc/tarball/v0.11.0
        0.10.0     https://github.com/openmc-dev/openmc/tarball/v0.10.0

    Variants:
        Name [Default]                 Allowed values          Description
        ===========================    ====================    ==================================

        build_type [RelWithDebInfo]    Debug, Release,         CMake build type
                                       RelWithDebInfo,         
                                       MinSizeRel              
        debug [off]                    on, off                 Enable debug flags
        ipo [off]                      on, off                 CMake interprocedural optimization
        mpi [off]                      on, off                 Enable MPI support
        openmp [on]                    on, off                 Enable OpenMP support
        optimize [off]                 on, off                 Enable optimization flags

    Installation Phases:
        cmake    build    install

    Build Dependencies:
        cmake  git  hdf5  mpi

    Link Dependencies:
        hdf5  mpi

    Run Dependencies:
        None

    Virtual Packages: 
        None

.. note::

    It should be noted that by default OpenMC builds with ``-O2 -g`` flags which
    are equivalent to a CMake build type of `RelwithDebInfo`. In addition, MPI
    is OFF while OpenMP is ON.

Assuming one had configured Spack with a GNU 9.3.0 compiler, to build OpenMC
with optimization, OpenMP support and OpenMPI, use:

.. code-block:: sh

    spack install openmc+mpi+optimize+openmp ^openmpi %gcc@9.3.0

Although the OpenMP variant is by default ON, one can explicitly write it in
the Spack specification. If a user wanted OpenMC without any MPI or OpenMP support, 
the variants can be deactivated:

.. code-block:: sh

    spack install openmc~mpi+optimize~openmp

The Python API for OpenMC can be installed in a similar fashion. Looking at the
information page presented with Spack:

.. code-block:: sh

    spack info py-openmc

    PythonPackage:   py-openmc

    Description:
        The OpenMC project aims to provide a fully-featured Monte Carlo particle
        transport code based on modern methods. It is a constructive solid
        geometry, continuous-energy transport code that uses ACE format cross
        sections. The project started under the Computational Reactor Physics
        Group at MIT.

    Homepage: http://openmc.readthedocs.io/

    Tags: 
        None

    Preferred version:  
        0.12.0     [git] https://github.com/openmc-dev/openmc.git at tag v0.12.0

    Safe versions:  
        develop    [git] https://github.com/openmc-dev/openmc.git on branch develop
        master     [git] https://github.com/openmc-dev/openmc.git on branch master
        0.12.0     [git] https://github.com/openmc-dev/openmc.git at tag v0.12.0
        0.11.0     https://github.com/openmc-dev/openmc/tarball/v0.11.0

    Variants:
        Name [Default]    Allowed values    Description
        ==============    ==============    ==================

        mpi [off]         on, off           Enable MPI support

    Installation Phases:
        build    install

    Build Dependencies:
        git  openmc  py-cython  py-h5py  py-ipython  py-lxml  py-matplotlib  py-mpi4py  py-numpy  py-pandas  py-scipy  py-setuptools  py-uncertainties  python

    Link Dependencies:
        python

    Run Dependencies:
        openmc  py-h5py  py-ipython  py-lxml  py-matplotlib  py-mpi4py  py-numpy  py-pandas  py-scipy  py-uncertainties  python

    Virtual Packages: 
        None

The only variant for the Python API is MPI. To configure a serial version of the Python API against a serial version of OpenMC:

.. code-block::

    spack install py-openmc

.. caution::

    When installing any Spack package, dependencies are assumed to be at configured defaults unless otherwise specfied in the
    specification on the command line. In the above example, assuming the default options weren't changed in Spack's package
    configuration, py-openmc will want to link against a non-optimized openmc. Even if you have a built an optimized openmc,
    it will rebuild openmc with optimization OFF. Thus, if you are trying to link against dependencies that were configured
    different than defaults, ^openmc[variants] will have to be present on command line.

To build a parallel version of py-openmc that links against a parallel version of openmc that was used in the previous example:

.. code-block::

    spack install py-openmc+mpi ^openmc+optimize ^openmpi

.. note::

    If py-openmc is given the +mpi variant, it is automatically passed to openmc so there is no need to specify it explicitly.

.. tip::

    When installing py-openmc it will use Spack's preferred Python. For example, assuming Spack's preferred Python
    is 3.8.7, to build py-openmc against the latest Python 3.7 instead, ``^python@3.7.0:3.7.99`` should be added to the
    specification on the command line.

A useful tool in Spack is to look at the dependency tree before installation. This can be observed using
Spack's `spec` tool:

.. code-block::

    spack spec py-openmc+mpi ^openmc+optimize %gcc@9.3.0

    Input spec
    --------------------------------
    py-openmc+mpi
        ^openmc%gcc@9.3.0+optimize

    Concretized
    --------------------------------
    py-openmc@0.12.0%gcc@9.3.0+mpi arch=linux-centos7-x86_64
        ^git@2.29.0%gcc@9.3.0+svn~tcltk arch=linux-centos7-x86_64
            ^autoconf@2.69%gcc@9.3.0 arch=linux-centos7-x86_64
                ^m4@1.4.18%gcc@9.3.0+sigsegv patches=3877ab548f88597ab2327a2230ee048d2d07ace1062efe81fc92e91b7f39cd00,fc9b61654a3ba1a8d6cd78ce087e7c96366c290bc8d2c299f09828d793b853c8 arch=linux-centos7-x86_64
                    ^libsigsegv@2.12%gcc@9.3.0 arch=linux-centos7-x86_64
                ^perl@5.32.0%gcc@9.3.0+cpanm+shared+threads arch=linux-centos7-x86_64
                    ^berkeley-db@18.1.40%gcc@9.3.0 arch=linux-centos7-x86_64
                    ^gdbm@1.18.1%gcc@9.3.0 arch=linux-centos7-x86_64
                        ^readline@8.0%gcc@9.3.0 arch=linux-centos7-x86_64
                            ^ncurses@6.2%gcc@9.3.0~symlinks+termlib arch=linux-centos7-x86_64
                                ^pkgconf@1.7.3%gcc@9.3.0 arch=linux-centos7-x86_64
            ^automake@1.16.2%gcc@9.3.0 arch=linux-centos7-x86_64
            ^curl@7.72.0%gcc@9.3.0~darwinssl~gssapi~libssh~libssh2~nghttp2 arch=linux-centos7-x86_64
                ^libidn2@2.3.0%gcc@9.3.0 arch=linux-centos7-x86_64
                    ^libunistring@0.9.10%gcc@9.3.0 arch=linux-centos7-x86_64
                        ^libiconv@1.16%gcc@9.3.0 arch=linux-centos7-x86_64
                ^openssl@1.1.1h%gcc@9.3.0+systemcerts arch=linux-centos7-x86_64
                    ^zlib@1.2.11%gcc@9.3.0+optimize+pic+shared arch=linux-centos7-x86_64
            ^expat@2.2.10%gcc@9.3.0+libbsd arch=linux-centos7-x86_64
                ^libbsd@0.10.0%gcc@9.3.0 arch=linux-centos7-x86_64
            ^gettext@0.21%gcc@9.3.0+bzip2+curses+git~libunistring+libxml2+tar+xz arch=linux-centos7-x86_64
                ^bzip2@1.0.8%gcc@9.3.0+shared arch=linux-centos7-x86_64
                    ^diffutils@3.7%gcc@9.3.0 arch=linux-centos7-x86_64
                ^libxml2@2.9.10%gcc@9.3.0~python arch=linux-centos7-x86_64
                    ^xz@5.2.5%gcc@9.3.0~pic arch=linux-centos7-x86_64
                ^tar@1.32%gcc@9.3.0 arch=linux-centos7-x86_64
            ^libtool@2.4.6%gcc@9.3.0 arch=linux-centos7-x86_64
            ^openssh@8.4p1%gcc@9.3.0 arch=linux-centos7-x86_64
                ^libedit@3.1-20191231%gcc@9.3.0 arch=linux-centos7-x86_64
            ^pcre2@10.35%gcc@9.3.0~jit+multibyte arch=linux-centos7-x86_64
            ^perl-alien-svn@1.8.11.0%gcc@9.3.0 arch=linux-centos7-x86_64
                ^apr@1.6.2%gcc@9.3.0 arch=linux-centos7-x86_64
                ^apr-util@1.6.1%gcc@9.3.0+crypto~gdbm~odbc~pgsql~sqlite arch=linux-centos7-x86_64
                ^perl-module-build@0.4224%gcc@9.3.0 arch=linux-centos7-x86_64
                ^sqlite@3.33.0%gcc@9.3.0+column_metadata+fts~functions~rtree arch=linux-centos7-x86_64
        ^openmc@0.12.0%gcc@9.3.0~debug~ipo+mpi+openmp+optimize build_type=Release arch=linux-centos7-x86_64
            ^cmake@3.18.4%gcc@9.3.0~doc+ncurses+openssl+ownlibs~qt patches=bf695e3febb222da2ed94b3beea600650e4318975da90e4a71d6f31a6d5d8c3d arch=linux-centos7-x86_64
            ^hdf5@1.10.7%gcc@9.3.0+cxx~debug+fortran+hl+java+mpi+pic+shared+szip+threadsafe api=none arch=linux-centos7-x86_64
                ^libszip@2.1.1%gcc@9.3.0 arch=linux-centos7-x86_64
                ^numactl@2.0.14%gcc@9.3.0 patches=4e1d78cbbb85de625bad28705e748856033eaafab92a66dffd383a3d7e00cc94 arch=linux-centos7-x86_64
                ^openjdk@11.0.2%gcc@9.3.0 arch=linux-centos7-x86_64
                ^openmpi@3.1.6%gcc@9.3.0~atomics~cuda~cxx~cxx_exceptions+gpfs~java~legacylaunchers~lustre~memchecker~pmi~singularity~sqlite3+static+thread_multiple+vt+wrapper-rpath fabrics=ofi,verbs schedulers=none arch=linux-centos7-x86_64
                    ^hwloc@1.11.11%gcc@9.3.0~cairo~cuda~gl~libudev+libxml2~netloc~nvml+pci+shared arch=linux-centos7-x86_64
                        ^libpciaccess@0.16%gcc@9.3.0 arch=linux-centos7-x86_64
                            ^util-macros@1.19.1%gcc@9.3.0 arch=linux-centos7-x86_64
                    ^libfabric@1.11.0%gcc@9.3.0~kdreg fabrics=shm,tcp,udp,verbs arch=linux-centos7-x86_64
                        ^rdma-core@20%gcc@9.3.0~ipo build_type=Release arch=linux-centos7-x86_64
                            ^libnl@3.3.0%gcc@9.3.0 arch=linux-centos7-x86_64
                                ^bison@3.6.4%gcc@9.3.0 arch=linux-centos7-x86_64
                                    ^help2man@1.47.11%gcc@9.3.0 arch=linux-centos7-x86_64
                                ^flex@2.6.4%gcc@9.3.0+lex patches=09c22e5c6fef327d3e48eb23f0d610dcd3a35ab9207f12e0f875701c677978d3 arch=linux-centos7-x86_64
                                    ^findutils@4.6.0%gcc@9.3.0 patches=84b916c0bf8c51b7e7b28417692f0ad3e7030d1f3c248ba77c42ede5c1c5d11e,bd9e4e5cc280f9753ae14956c4e4aa17fe7a210f55dd6c84aa60b12d106d47a2 arch=linux-centos7-x86_64
                                        ^texinfo@6.5%gcc@9.3.0 patches=12f6edb0c6b270b8c8dba2ce17998c580db01182d871ee32b7b6e4129bd1d23a,1732115f651cff98989cb0215d8f64da5e0f7911ebf0c13b064920f088f2ffe1 arch=linux-centos7-x86_64
                            ^py-docutils@0.15.2%gcc@9.3.0 arch=linux-centos7-x86_64
                                ^py-setuptools@50.3.2%gcc@9.3.0 arch=linux-centos7-x86_64
                                    ^python@3.8.7%gcc@9.3.0+bz2+ctypes+dbm~debug+libxml2+lzma~nis~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tix~tkinter~ucs4+uuid+zlib patches=0d98e93189bc278fbc37a50ed7f183bd8aaf249a8e1670a465f0db6bb4f8cf87 arch=linux-centos7-x86_64
                                        ^libffi@3.3%gcc@9.3.0 patches=26f26c6f29a7ce9bf370ad3ab2610f99365b4bdd7b82e7c31df41a3370d685c0 arch=linux-centos7-x86_64
                                        ^uuid@1.6.2%gcc@9.3.0 arch=linux-centos7-x86_64
        ^py-cython@0.29.21%gcc@9.3.0 arch=linux-centos7-x86_64
        ^py-h5py@2.10.0%gcc@9.3.0+mpi arch=linux-centos7-x86_64
            ^py-mpi4py@3.0.3%gcc@9.3.0 arch=linux-centos7-x86_64
            ^py-numpy@1.19.4%gcc@9.3.0+blas+lapack arch=linux-centos7-x86_64
                ^openblas@0.3.12%gcc@9.3.0~consistent_fpcsr~ilp64+pic+shared threads=none arch=linux-centos7-x86_64
            ^py-pkgconfig@1.5.1%gcc@9.3.0 arch=linux-centos7-x86_64
                ^py-nose@1.3.7%gcc@9.3.0 arch=linux-centos7-x86_64
            ^py-six@1.14.0%gcc@9.3.0 arch=linux-centos7-x86_64
        ^py-ipython@7.18.1%gcc@9.3.0 arch=linux-centos7-x86_64
            ^py-backcall@0.1.0%gcc@9.3.0 arch=linux-centos7-x86_64
            ^py-decorator@4.4.2%gcc@9.3.0 arch=linux-centos7-x86_64
            ^py-jedi@0.13.3%gcc@9.3.0 arch=linux-centos7-x86_64
                ^py-parso@0.6.1%gcc@9.3.0 arch=linux-centos7-x86_64
            ^py-pexpect@4.7.0%gcc@9.3.0 arch=linux-centos7-x86_64
                ^py-ptyprocess@0.6.0%gcc@9.3.0 arch=linux-centos7-x86_64
            ^py-pickleshare@0.7.5%gcc@9.3.0 arch=linux-centos7-x86_64
            ^py-prompt-toolkit@2.0.9%gcc@9.3.0 arch=linux-centos7-x86_64
                ^py-wcwidth@0.1.7%gcc@9.3.0 arch=linux-centos7-x86_64
            ^py-pygments@2.6.1%gcc@9.3.0 arch=linux-centos7-x86_64
            ^py-traitlets@5.0.4%gcc@9.3.0 arch=linux-centos7-x86_64
                ^py-ipython-genutils@0.2.0%gcc@9.3.0 arch=linux-centos7-x86_64
        ^py-lxml@4.5.2%gcc@9.3.0~cssselect~html5~htmlsoup arch=linux-centos7-x86_64
            ^libxslt@1.1.33%gcc@9.3.0+crypto~python arch=linux-centos7-x86_64
                ^libgcrypt@1.8.5%gcc@9.3.0 arch=linux-centos7-x86_64
                    ^libgpg-error@1.37%gcc@9.3.0 arch=linux-centos7-x86_64
                        ^gawk@5.0.1%gcc@9.3.0 arch=linux-centos7-x86_64
                            ^gmp@6.1.2%gcc@9.3.0 arch=linux-centos7-x86_64
                            ^mpfr@4.0.2%gcc@9.3.0 patches=3f80b836948aa96f8d1cb9cc7f3f55973f19285482a96f9a4e1623d460bcccf0 arch=linux-centos7-x86_64
                                ^autoconf-archive@2019.01.06%gcc@9.3.0 arch=linux-centos7-x86_64
        ^py-matplotlib@3.3.3%gcc@9.3.0~animation~fonts+image~latex~movies backend=agg arch=linux-centos7-x86_64
            ^freetype@2.10.1%gcc@9.3.0 arch=linux-centos7-x86_64
                ^libpng@1.6.37%gcc@9.3.0 arch=linux-centos7-x86_64
            ^py-certifi@2020.6.20%gcc@9.3.0 arch=linux-centos7-x86_64
            ^py-cycler@0.10.0%gcc@9.3.0 arch=linux-centos7-x86_64
            ^py-kiwisolver@1.1.0%gcc@9.3.0 arch=linux-centos7-x86_64
            ^py-pillow@7.2.0%gcc@9.3.0~freetype~imagequant+jpeg~jpeg2000~lcms~tiff~webp~webpmux~xcb+zlib arch=linux-centos7-x86_64
                ^libjpeg-turbo@2.0.4%gcc@9.3.0 arch=linux-centos7-x86_64
                    ^nasm@2.15.05%gcc@9.3.0 arch=linux-centos7-x86_64
            ^py-pyparsing@2.4.7%gcc@9.3.0 arch=linux-centos7-x86_64
            ^py-python-dateutil@2.8.0%gcc@9.3.0 arch=linux-centos7-x86_64
                ^py-setuptools-scm@4.1.2%gcc@9.3.0+toml arch=linux-centos7-x86_64
                    ^py-toml@0.10.2%gcc@9.3.0 arch=linux-centos7-x86_64
            ^qhull@2020.1%gcc@9.3.0~ipo build_type=Release arch=linux-centos7-x86_64
        ^py-pandas@1.1.4%gcc@9.3.0 arch=linux-centos7-x86_64
            ^py-bottleneck@1.2.1%gcc@9.3.0 arch=linux-centos7-x86_64
            ^py-numexpr@2.7.0%gcc@9.3.0 arch=linux-centos7-x86_64
            ^py-pytz@2020.1%gcc@9.3.0 arch=linux-centos7-x86_64
        ^py-scipy@1.5.4%gcc@9.3.0 arch=linux-centos7-x86_64
            ^py-pybind11@2.5.0%gcc@9.3.0~ipo build_type=Release arch=linux-centos7-x86_64
        ^py-uncertainties@3.1.4%gcc@9.3.0~docs~optional arch=linux-centos7-x86_64
            ^py-future@0.18.2%gcc@9.3.0 arch=linux-centos7-x86_64

Once installed, environment/lmod modules can be generated or Spack's `load` feature
can be used to access the installed packages. 

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

          CXX=mpicxx.mpich cmake -DHDF5_PREFER_PARALLEL=on ..

      Note that the exact package names may vary depending on your particular
      distribution and version.

      If you are using building HDF5 from source in conjunction with MPI, we
      recommend that your HDF5 installation be built with parallel I/O
      features. An example of configuring HDF5_ is listed below::

          CC=mpicc ./configure --enable-parallel

      You may omit ``--enable-parallel`` if you want to compile HDF5_ in serial.

.. admonition:: Optional
   :class: note

    * An MPI implementation for distributed-memory parallel runs

      To compile with support for parallel runs on a distributed-memory
      architecture, you will need to have a valid implementation of MPI
      installed on your machine. The code has been tested and is known to work
      with the latest versions of both OpenMPI_ and MPICH_. OpenMPI and/or MPICH
      can be installed on Debian derivatives with::

          sudo apt install mpich libmpich-dev
          sudo apt install openmpi-bin libopenmpi-dev

    * DAGMC_ toolkit for simulation using CAD-based geometries

      OpenMC supports particle tracking in CAD-based geometries via the Direct
      Accelerated Geometry Monte Carlo (DAGMC) toolkit (`installation
      instructions
      <https://svalinn.github.io/DAGMC/install/dag_multiple.html>`_). For use in
      OpenMC, only the ``MOAB_DIR`` and ``BUILD_TALLY`` variables need to be
      specified in the CMake configuration step.

    * git_ version control software for obtaining source code


.. _gcc: https://gcc.gnu.org/
.. _CMake: http://www.cmake.org
.. _OpenMPI: http://www.open-mpi.org
.. _MPICH: http://www.mpich.org
.. _HDF5: https://www.hdfgroup.org/solutions/hdf5/
.. _DAGMC: https://svalinn.github.io/DAGMC/index.html

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

debug
  Enables debugging when compiling. The flags added are dependent on which
  compiler is used.

profile
  Enables profiling using the GNU profiler, gprof.

optimize
  Enables high-optimization using compiler-dependent flags. For gcc and
  Intel C++, this compiles with -O3.

openmp
  Enables shared-memory parallelism using the OpenMP API. The C++ compiler
  being used must support OpenMP. (Default: on)

dagmc
  Enables use of CAD-based DAGMC_ geometries. Please see the note about DAGMC in
  the optional dependencies list for more information on this feature. The
  installation directory for DAGMC should also be defined as `DAGMC_ROOT` in the
  CMake configuration command. (Default: off)

coverage
  Compile and link code instrumented for coverage analysis. This is typically
  used in conjunction with gcov_.

To set any of these options (e.g. turning on debug mode), the following form
should be used:

.. code-block:: sh

    cmake -Ddebug=on /path/to/openmc

.. _gcov: https://gcc.gnu.org/onlinedocs/gcc/Gcov.html

.. _usersguide_compile_mpi:

Compiling with MPI
++++++++++++++++++

To compile with MPI, set the :envvar:`CXX` environment variable to the path to
the MPI C++ wrapper. For example, in a bash shell:

.. code-block:: sh

    export CXX=mpicxx
    cmake /path/to/openmc

Note that in many shells, environment variables can be set for a single command,
i.e.

.. code-block:: sh

    CXX=mpicxx cmake /path/to/openmc

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

Compiling on Windows 10
-----------------------

Recent versions of Windows 10 include a subsystem for Linux that allows one to
run Bash within Ubuntu running in Windows. First, follow the installation guide
`here <https://msdn.microsoft.com/en-us/commandline/wsl/install_guide>`_ to get
Bash on Ubuntu on Windows setup. Once you are within bash, obtain the necessary
:ref:`prerequisites <prerequisites>` via ``apt``. Finally, follow the
:ref:`instructions for compiling on linux <compile_linux>`.

Compiling for the Intel Xeon Phi
--------------------------------

For the second generation Knights Landing architecture, nothing special is
required to compile OpenMC. You may wish to experiment with compiler flags that
control generation of vector instructions to see what configuration gives
optimal performance for your target problem.

For the first generation Knights Corner architecture, it is necessary to
cross-compile OpenMC. If you are using the Intel compiler, it is necessary to
specify that all objects be compiled with the ``-mmic`` flag as follows:

.. code-block:: sh

    mkdir build && cd build
    CXX=icpc CXXFLAGS=-mmic cmake -Dopenmp=on ..
    make

Note that unless an HDF5 build for the Intel Xeon Phi (Knights Corner) is
already on your target machine, you will need to cross-compile HDF5 for the Xeon
Phi. An `example script`_ to build zlib and HDF5 provides several necessary
workarounds.

.. _example script: https://github.com/paulromano/install-scripts/blob/master/install-hdf5-mic

Testing Build
-------------

To run the test suite, you will first need to download a pre-generated cross
section library along with windowed multipole data. Please refer to our
:ref:`devguide_tests` documentation for further details.

---------------------
Installing Python API
---------------------

If you installed OpenMC using :ref:`Conda <install_conda>` or :ref:`PPA
<install_ppa>`, no further steps are necessary in order to use OpenMC's
:ref:`Python API <pythonapi>`. However, if you are :ref:`installing from source
<install_source>`, the Python API is not installed by default when ``make
install`` is run because in many situations it doesn't make sense to install a
Python package in the same location as the ``openmc`` executable (for example,
if you are installing the package into a `virtual environment
<https://docs.python.org/3/tutorial/venv.html>`_). The easiest way to install
the :mod:`openmc` Python package is to use pip_, which is included by default in
Python 3.4+. From the root directory of the OpenMC distribution/repository, run:

.. code-block:: sh

    pip install .

pip will first check that all :ref:`required third-party packages
<usersguide_python_prereqs>` have been installed, and if they are not present,
they will be installed by downloading the appropriate packages from the Python
Package Index (`PyPI <https://pypi.org/>`_). However, do note that since pip
runs the ``setup.py`` script which requires NumPy, you will have to first
install NumPy:

.. code-block:: sh

    pip install numpy

Installing in "Development" Mode
--------------------------------

If you are primarily doing development with OpenMC, it is strongly recommended
to install the Python package in :ref:`"editable" mode <devguide_editable>`.

.. _usersguide_python_prereqs:

Prerequisites
-------------

The Python API works with Python 3.5+. In addition to Python itself, the API
relies on a number of third-party packages. All prerequisites can be installed
using Conda_ (recommended), pip_, or through the package manager in most Linux
distributions.

.. admonition:: Required
   :class: error

   `NumPy <http://www.numpy.org/>`_
      NumPy is used extensively within the Python API for its powerful
      N-dimensional array.

   `SciPy <https://www.scipy.org/>`_
      SciPy's special functions, sparse matrices, and spatial data structures
      are used for several optional features in the API.

   `pandas <http://pandas.pydata.org/>`_
      Pandas is used to generate tally DataFrames as demonstrated in
      an `example notebook <../examples/pandas-dataframes.ipynb>`_.

   `h5py <http://www.h5py.org/>`_
      h5py provides Python bindings to the HDF5 library. Since OpenMC outputs
      various HDF5 files, h5py is needed to provide access to data within these
      files from Python.

   `Matplotlib <http://matplotlib.org/>`_
      Matplotlib is used to providing plotting functionality in the API like the
      :meth:`Universe.plot` method and the :func:`openmc.plot_xs` function.

   `uncertainties <https://pythonhosted.org/uncertainties/>`_
      Uncertainties are used for decay data in the :mod:`openmc.data` module.

   `lxml <http://lxml.de/>`_
      lxml is used for the :ref:`scripts_validate` script and various other
      parts of the Python API.

.. admonition:: Optional
   :class: note

   `mpi4py <https://mpi4py.readthedocs.io/en/stable/>`_
      mpi4py provides Python bindings to MPI for running distributed-memory
      parallel runs. This package is needed if you plan on running depletion
      simulations in parallel using MPI.

   `Cython <http://cython.org/>`_
      Cython is used for resonance reconstruction for ENDF data converted to
      :class:`openmc.data.IncidentNeutron`.

   `vtk <http://www.vtk.org/>`_
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
    MPICC=<path to mpicc> pip install mpi4py
    HDF5_DIR=<path to HDF5> pip install --no-binary=h5py h5py

If you are using parallel HDF5, you'll also need to make sure the right MPI
wrapper is used when installing h5py:

.. code-block:: sh

    CC=<path to mpicc> HDF5_MPI=ON HDF5_DIR=<path to HDF5> pip install --no-binary=h5py h5py

.. _usersguide_nxml:

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
.. _validation: https://en.wikipedia.org/wiki/XML_validation
.. _RELAX NG: http://relaxng.org/
.. _ctest: https://cmake.org/cmake/help/latest/manual/ctest.1.html
.. _Conda: https://conda.io/en/latest/
.. _pip: https://pip.pypa.io/en/stable/
