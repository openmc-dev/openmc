#!/bin/bash
set -ex

# Run source check
cd tests
if [[ $TRAVIS_PYTHON_VERSION == "3.4" && $OMP == 'n' && $MPI == 'n' ]]; then
    ./check_source.py
fi

# Run regression and unit tests
if [[ $MPI == 'y' ]]; then
    pytest --cov=../openmc -v --mpi
else
    pytest --cov=../openmc -v
fi
