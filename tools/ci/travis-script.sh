#!/bin/bash
set -ex

# Run source check
if [[ $TRAVIS_PYTHON_VERSION == "3.4" && $OMP == 'n' && $MPI == 'n' ]]; then
    pushd tests && python check_source.py && popd
fi

# Run regression and unit tests
if [[ $MPI == 'y' ]]; then
    pytest --cov=openmc -v --mpi tests
else
    pytest --cov=openmc -v tests
fi
