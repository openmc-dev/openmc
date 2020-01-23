#!/bin/bash
set -ex

# Run regression and unit tests
if [[ $MPI == 'y' ]]; then
  if [[ $EVENT == 'y' ]]; then
    pytest --cov=openmc -v --mpi --event tests
  else
    pytest --cov=openmc -v --mpi tests
  fi
else
  if [[ $EVENT == 'y' ]]; then
    pytest --cov=openmc -v --event tests
  else
    pytest --cov=openmc -v tests
  fi
fi
