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

# Run unit tests and then regression tests
pytest -v $args \
  tests/test_matplotlib_import.py \
  tests/unit_tests \
  tests/regression_tests
