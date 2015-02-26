#!/bin/sh

set -ev

# Run all debug tests
if [ "$TRAVIS_PULL_REQUEST" != "false" ]; then
  ./run_tests.py -C "^basic-debug$|^hdf5-debug$|^mpi-omp-debug$|^phdf5-omp-debug$|^omp-phdf5-petsc-debug$" -j 4 -s
else
  ./run_tests.py -C "^basic-debug$" -j 4
fi
