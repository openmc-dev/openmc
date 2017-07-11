#!/bin/sh

set -ev

# Run all debug tests
./check_source.py
./run_tests.py -C $OPENMC_CONFIG -j 2
