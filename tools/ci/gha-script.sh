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

# Rename openmc to openmc-test
mv openmc openmc-test

# Run regression and unit tests
pytest --cov=openmc -v $args tests

# Rename to openmc back
mv openmc-test openmc
