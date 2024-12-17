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
  # Change environmental variables
  eval $( "${HOME}/ncrystal_inst/bin/ncrystal-config" --setup )
  nctool --test
fi

# Move openmc to python/openmc
mv openmc python/openmc

# Run regression and unit tests
pytest --cov=openmc -v $args tests

# Move openmc back
mv python/openmc openmc
