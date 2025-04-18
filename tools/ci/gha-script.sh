#!/bin/bash
set -ex

# Argument List
args=" "

# Check for MPI
if [[ $MPI == 'y' ]]; then
  args="${args} --mpi "
fi

# Check for event-based
if [[ $EVENT == 'y' ]]; then
  args="${args} --event "
fi

# Skip source directory warning
export OPENMC_DEV_MODE=1

# Run regression and unit tests
pytest --cov=openmc -v $args tests

