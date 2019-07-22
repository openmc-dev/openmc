.. _devguide_tests:

==========
Test Suite
==========

The OpenMC test suite consists of two parts, a regression test suite and a unit
test suite. The regression test suite is based on regression or integrated
testing where different types of input files are configured and the full OpenMC
code is executed. Results from simulations are compared with expected
results. The unit tests are primarily intended to test individual
functions/classes in the OpenMC Python API.

Prerequisites
-------------

- The test suite relies on the third-party `pytest <https://pytest.org>`_
  package. To run either or both the regression and unit test suites, it is
  assumed that you have OpenMC fully installed, i.e., the :ref:`scripts_openmc`
  executable is available on your :envvar:`PATH` and the :mod:`openmc` Python
  module is importable. In development where it would be onerous to continually
  install OpenMC every time a small change is made, it is recommended to install
  OpenMC in development/editable mode. With setuptools, this is accomplished by
  running::

      python setup.py develop

  or using pip (recommended)::

      pip install -e .[test]

- The test suite requires a specific set of cross section data in order for
  tests to pass. A download URL for the data that OpenMC expects can be found
  within ``tools/ci/download-xs.sh``.
- In addition to the HDF5 data, some tests rely on ENDF files. A download URL
  for those can also be found in ``tools/ci/download-xs.sh``.
- Some tests require `NJOY <https://www.njoy21.io/NJOY2016>`_ to preprocess
  cross section data. The test suite assumes that you have an ``njoy``
  executable available on your :envvar:`PATH`.

Running Tests
-------------

To execute the test suite, go to the ``tests/`` directory and run::

    pytest

If you want to collect information about source line coverage in the Python API,
you must have the `pytest-cov <https://pypi.python.org/pypi/pytest-cov>`_ plugin
installed and run::

    pytest --cov=../openmc --cov-report=html

Generating XML Inputs
---------------------

Many of the regression tests rely on the Python API to build an appropriate
model. However, it can sometimes be desirable to work directly with the XML
input files rather than having to run a script in order to run the problem/test.
To build the input files for a test without actually running the test, you can
run::

    pytest --build-inputs <name-of-test>

Adding Tests to the Regression Suite
------------------------------------

To add a new test to the regression test suite, create a sub-directory in the
``tests/regression_tests/`` directory. To configure a test you need to add the
following files to your new test directory:

    * OpenMC input XML files, if they are not generated through the Python API
    * **test.py** - Python test driver script; please refer to other tests to
      see how to construct. Any output files that are generated during testing
      must be removed at the end of this script.
    * **inputs_true.dat** - ASCII file that contains Python API-generated XML
      files concatenated together. When the test is run, inputs that are
      generated are compared to this file.
    * **results_true.dat** - ASCII file that contains the expected results from
      the test. The file *results_test.dat* is compared to this file during the
      execution of the python test driver script. When the above files have been
      created, generate a *results_test.dat* file and copy it to this name and
      commit. It should be noted that this file should be generated with basic
      compiler options during openmc configuration and build (e.g., no MPI, no
      debug/optimization).

In addition to this description, please see the various types of tests that are
already included in the test suite to see how to create them. If all is
implemented correctly, the new test will automatically be discovered by pytest.
