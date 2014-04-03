=================
OpenMC Test Suite
=================

The purpose of this test suite is to ensure that OpenMC compiles using various
combinations of compiler flags and options and that all user input options can
be used successfully without breaking the code. The test suite is based on
regression or integrated testing where different types of input files are 
configured and the full OpenMC code is executed. Results from simulations
are compared with expected results. The test suite is comprised of many
build configurations (e.g. debug, mpi, hdf5) and the actual tests which
reside in sub-directories in the tests directory.

The test suite is designed to integrate with cmake using ctest_. To run the
full test suite run:

.. sh::
    python run_tests.py

A subset of build configurations and/or tests can be run. To see how to use
the script run:

.. sh::
    python run_tests.py --help

As an example, say we want to run all tests with debug flags only on tests
that have cone and plot in their name. Also, we would like to run this on 
4 processors. We can run:

.. sh::
    python run_tests.py -j 4 -C debug -R "cone|plot"

Note that standard regular expression syntax is used for selecting build
configurations and tests. To print out a list of build configurations, we
can run:

.. sh::
    python run_tests.py -p

.. _ctest: http://www.cmake.org/cmake/help/v2.8.8/ctest.html
