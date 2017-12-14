#!/bin/bash
set -ex
cd tests
if [[ $OPENMC_CONFIG == "check_source" ]]; then
    ./check_source.py;
else
    ./run_tests.py -C $OPENMC_CONFIG -j 2
    pytest --cov=../openmc -v unit_tests/
fi
