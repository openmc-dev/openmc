#!/bin/sh

set -ev

# Run all debug tests
./check_source.py
if [ "$TRAVIS_PULL_REQUEST" != "false" ]; then
  ./run_tests.py -C "^hdf5-debug$|^omp-hdf5-debug|^mpi-hdf5-debug|^phdf5-debug$|^phdf5-omp-debug$" -j 2 -s
else
  ./run_tests.py -C "^hdf5-debug$" -j 2
fi
