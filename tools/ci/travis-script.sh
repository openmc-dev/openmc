#!/bin/bash
set -ex

# Run regression and unit tests
if [[ $MPI == 'y' ]]; then
    pytest --cov=openmc -v --mpi tests
else
    pytest --cov=openmc -v tests
fi
