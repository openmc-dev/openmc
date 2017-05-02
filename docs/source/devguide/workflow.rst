.. _devguide_workflow:

====================
Development Workflow
====================

Anyone wishing to make contributions to OpenMC should be fully acquianted and
comfortable working with git_ and GitHub_. We assume here that you have git
installed on your system, have a GitHub account, and have setup SSH keys to be
able to create/push to repositories on GitHub.

Overview
--------

Development of OpenMC relies heavily on branching; specifically, we use a
branching model sometimes referred to as `git flow`_. If you plan to contribute
to OpenMC development, we highly recommend that you read the linked blog post to
get a sense of how the branching model works. There are two main branches that
always exist: *master* and *develop*. The *master* branch is a stable branch
that contains the latest release of the code. The *develop* branch is where any
ongoing development takes place prior to a release and is not guaranteed to be
stable. When the development team decides that a release should occur, the
*develop* branch is merged into *master*.

Trivial changes to the code may be committed directly to the *develop* branch by
a trusted developer. However, most new features should be developed on a branch
that branches off of *develop*. When the feature is completed, a `pull request`_
is initiated on GitHub that is then reviewed by a trusted developer. If the pull
request is satisfactory, it is then merged into *develop*. Note that a trusted
developer may not review their own pull request (i.e., an independent code
review is required).

Code Review Criteria
--------------------

In order to be considered suitable for inclusion in the *develop* branch, the
following criteria must be satisfied for all proposed changes:

- Changes have a clear purpose and are useful.
- Compiles and passes the regression suite with all configurations (This is
  checked by Travis CI).
- If appropriate, test cases are added to regression suite.
- No memory leaks (checked with valgrind_).
- Conforms to the OpenMC `style guide`_.
- No degradation of performance or greatly increased memory usage. This is not a
  hard rule -- in certain circumstances, a performance loss might be acceptable
  if there are compelling reasons.
- New features/input are documented.
- No unnecessary external software dependencies are introduced.

Contributing
------------

Now that you understand the basic development workflow, let's discuss how an
individual to contribute to development. Note that this would apply to both new
features and bug fixes. The general steps for contributing are as follows:

1. Fork the main openmc repository from `mit-crpg/openmc`_. This will create a
   repository with the same name under your personal account. As such, you can
   commit to it as you please without disrupting other developers.

   .. image:: ../_images/fork.png

2. Clone your fork of OpenMC and create a branch that branches off of *develop*:

   .. code-block:: sh

       git clone git@github.com:yourusername/openmc.git
       cd openmc
       git checkout -b newbranch develop

3. Make your changes on the new branch that you intend to have included in
   *develop*. If you have made other changes that should not be merged back,
   ensure that those changes are made on a different branch.

4. Issue a pull request from GitHub and select the *develop* branch of
   mit-crpg/openmc as the target.

   .. image:: ../_images/pullrequest.png

   At a minimum, you should describe what the changes you've made are and why
   you are making them. If the changes are related to an oustanding issue, make
   sure it is cross-referenced.

5. A trusted developer will review your pull request based on the criteria
   above. Any issues with the pull request can be discussed directly on the pull
   request page itself.

6. After the pull request has been thoroughly vetted, it is merged back into the
   *develop* branch of mit-crpg/openmc.

.. _test suite:

OpenMC Test Suite
-----------------

The purpose of this test suite is to ensure that OpenMC compiles using various
combinations of compiler flags and options, and that all user input options can
be used successfully without breaking the code. The test suite is comprised of
regression tests where different types of input files are configured and the
full OpenMC code is executed. Results from simulations are compared with
expected results. The test suite is comprised of many build configurations
(e.g. debug, mpi, hdf5) and the actual tests which reside in sub-directories
in the tests directory. We recommend to developers to test their branches
before submitting a formal pull request using gfortran and Intel compilers
if available.

The test suite is designed to integrate with cmake using ctest_.  It is
configured to run with cross sections from NNDC_ augmented with 0 K elastic
scattering data for select nuclides as well as multipole data. To download the
proper data, run the following commands:

.. code-block:: sh

    wget -O nndc_hdf5.tar.xz $(cat <openmc_root>/.travis.yml | grep anl.box | awk '{print $2}')
    tar xJvf nndc_hdf5.tar.xz
    export OPENMC_CROSS_SECTIONS=$(pwd)/nndc_hdf5/cross_sections.xml

    git clone --branch=master git://github.com/smharper/windowed_multipole_library.git wmp_lib
    tar xzvf wmp_lib/multipole_lib.tar.gz
    export OPENMC_MULTIPOLE_LIBRARY=$(pwd)/multipole_lib

The test suite can be run on an already existing build using:

.. code-block:: sh

    cd build
    make test

or

.. code-block:: sh

    cd build
    ctest

There are numerous ctest_ command line options that can be set to have
more control over which tests are executed.

Before running the test suite python script, the following environmental
variables should be set if the default paths are incorrect:

    * **FC** - The command for a Fortran compiler (e.g. gfotran, ifort).

        * Default - *gfortran*

    * **CC** - The command for a C compiler (e.g. gcc, icc).

        * Default - *gcc*

    * **CXX** - The command for a C++ compiler (e.g. g++, icpc).

        * Default - *g++*

    * **MPI_DIR** - The path to the MPI directory.

        * Default - */opt/mpich/3.2-gnu*

    * **HDF5_DIR** - The path to the HDF5 directory.

        * Default - */opt/hdf5/1.8.16-gnu*

    * **PHDF5_DIR** - The path to the parallel HDF5 directory.

        * Default - */opt/phdf5/1.8.16-gnu*

To run the full test suite, the following command can be executed in the
tests directory:

.. code-block:: sh

    python run_tests.py

A subset of build configurations and/or tests can be run. To see how to use
the script run:

.. code-block:: sh

    python run_tests.py --help

As an example, say we want to run all tests with debug flags only on tests
that have cone and plot in their name. Also, we would like to run this on
4 processors. We can run:

.. code-block:: sh

    python run_tests.py -j 4 -C debug -R "cone|plot"

Note that standard regular expression syntax is used for selecting build
configurations and tests. To print out a list of build configurations, we
can run:

.. code-block:: sh

    python run_tests.py -p

Adding tests to test suite
++++++++++++++++++++++++++

To add a new test to the test suite, create a sub-directory in the tests
directory that conforms to the regular expression *test_*. To configure
a test you need to add the following files to your new test directory,
*test_name* for example:

    * OpenMC input XML files
    * **test_name.py** - Python test driver script, please refer to other
      tests to see how to construct. Any output files that are generated
      during testing must be removed at the end of this script.
    * **inputs_true.dat** - ASCII file that contains Python API-generated XML
      files concatenated together. When the test is run, inputs that are
      generated are compared to this file.
    * **results_true.dat** - ASCII file that contains the expected results
      from the test. The file *results_test.dat* is compared to this file
      during the execution of the python test driver script. When the
      above files have been created, generate a *results_test.dat* file and
      copy it to this name and commit. It should be noted that this file
      should be generated with basic compiler options during openmc
      configuration and build (e.g., no MPI/HDF5, no debug/optimization).

In addition to this description, please see the various types of tests that
are already included in the test suite to see how to create them. If all is
implemented correctly, the new test directory will automatically be added
to the CTest framework.

Private Development
-------------------

While the process above depends on the fork of the OpenMC repository being
publicly available on GitHub, you may also wish to do development on a private
repository for research or commercial purposes. The proper way to do this is to
create a complete copy of the OpenMC repository (not a fork from GitHub). The
private repository can then either be stored just locally or in conjunction with
a private repository on Github (this requires a `paid plan`_). Alternatively,
`Bitbucket`_ offers private repositories for free. If you want to merge some
changes you've made in your private repository back to mit-crpg/openmc
repository, simply follow the steps above with an extra step of pulling a branch
from your private repository into a public fork.

.. _git: http://git-scm.com/
.. _GitHub: https://github.com/
.. _git flow: http://nvie.com/git-model
.. _valgrind: http://valgrind.org/
.. _style guide: http://mit-crpg.github.io/openmc/devguide/styleguide.html
.. _pull request: https://help.github.com/articles/using-pull-requests
.. _mit-crpg/openmc: https://github.com/mit-crpg/openmc
.. _paid plan: https://github.com/plans
.. _Bitbucket: https://bitbucket.org
.. _ctest: http://www.cmake.org/cmake/help/v2.8.12/ctest.html
.. _NNDC:  http://www.nndc.bnl.gov/endf/b7.1/acefiles.html
