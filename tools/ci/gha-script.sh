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

# Check NCrystal installation
if [[ $NCRYSTAL = 'y' ]]; then
  nctool --test
fi

# Run regression and unit tests
pytest --cov=openmc -v $args tests
