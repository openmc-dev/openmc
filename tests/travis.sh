#!/bin/sh

set -ev

# Run all debug tests
if [ "$TRAVIS_PULL_REQUEST" != "false" ]; then
  ./run_tests.py -C "^basic-debug$|^hdf5-debug$|^mpi-omp-debug$|^phdf5-omp-debug$" -j 2 -s
else
  ./run_tests.py -C "^basic-debug$" -j 2
fi
