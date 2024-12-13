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

# Remove openmc directory to avoid warnings
rm -rf openmc/*

# Run regression and unit tests
pytest --cov=openmc -v $args tests
