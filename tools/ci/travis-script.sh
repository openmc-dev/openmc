#!/bin/bash
set -ex

# Compile argument list
args=" "

if [[ $MPI == 'y' ]]; then
  args="${args} --mpi "
fi
  
if [[ $EVENT == 'y' ]]; then
  args="${args} --event "
fi

# Run regression and unit tests
pytest --cov=openmc -v $args tests
