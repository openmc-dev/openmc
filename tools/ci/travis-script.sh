#!/bin/bash
set -ex
cd tests
if [[ $TRAVIS_PYTHON_VERSION == "3.4" && $OPENMC_CONFIG == '^hdf5-debug$' ]]; then
    ./check_source.py
fi
./run_tests.py -C $OPENMC_CONFIG -j 2
pytest --cov=../openmc -v unit_tests/
