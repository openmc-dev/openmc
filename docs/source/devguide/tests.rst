.. _devguide_tests:

==========
Test Suite
==========

Running Tests
-------------

The OpenMC test suite consists of two parts, a regression test suite and a unit
test suite. The regression test suite is based on regression or integrated
testing where different types of input files are configured and the full OpenMC
code is executed. Results from simulations are compared with expected
results. The unit tests are primarily intended to test individual
functions/classes in the OpenMC Python API.

The test suite relies on the third-party `pytest <https://pytest.org>`_
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

It is also assumed that you have cross section data available that is pointed to
by the :envvar:`OPENMC_CROSS_SECTIONS` and :envvar:`OPENMC_MULTIPOLE_LIBRARY`
environment variables. Furthermore, to run unit tests for the :mod:`openmc.data`
module, it is necessary to have ENDF/B-VII.1 data available and pointed to by
the :envvar:`OPENMC_ENDF_DATA` environment variable. All data sources can be
obtained using the ``tools/ci/travis-before-script.sh`` script.

To execute the test suite, go to the ``tests/`` directory and run::

    pytest

If you want to collect information about source line coverage in the Python API,
you must have the `pytest-cov <https://pypi.python.org/pypi/pytest-cov>`_ plugin
installed and run::

    pytest --cov=../openmc --cov-report=html

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
