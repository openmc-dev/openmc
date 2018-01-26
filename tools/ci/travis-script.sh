#!/bin/bash
set -ex

# Run regression test suite
cd build
ctest

# Run source check
cd ../tests
if [[ $TRAVIS_PYTHON_VERSION == "3.4" && $OMP == 'n' && $MPI == 'n' ]]; then
    ./check_source.py
fi

# Run unit tests
pytest --cov=../openmc -v unit_tests/
