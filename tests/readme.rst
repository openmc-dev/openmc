=================
OpenMC Test Suite
=================

The purpose of this test suite is to ensure that OpenMC compiles using various
combinations of compiler flags and options and that all user input options can
be used successfully without breaking the code. The test suite is by no means
complete and should not be viewed as a comprehensive unit test suite will full
coverage. Until more effort can be put into actual unit testing, this suite is a
simple means of making sure that new features added into the code don't break
existing features.

The test suite is designed to run with the third-party Python package
nose_. Running the test suite is as simple as going to the tests/ directory and
running:

.. sh::
    nosetests

However, usually testing is split into two parts: compilation and running. To
run the compilation tests, use:

.. sh::
    nosetests test_compile

Then, to run all the normal tests (which require that an OpenMC executable is
already built):

.. sh::
    nosetests --exclude-dir=test_compile

.. _nose: https://nose.readthedocs.org
