#!/bin/bash
set -ex

# Install NJOY 2016
./tools/ci/travis-install-njoy.sh

# Running OpenMC's setup.py requires numpy/cython already
pip install numpy cython

# pytest installed by default -- make sure we get latest
pip install --upgrade pytest

# Pandas stopped supporting Python 3.4 with version 0.21
if [[ "$TRAVIS_PYTHON_VERSION" == "3.4" ]]; then
    pip install pandas==0.20.3
fi

# Install OpenMC in editable mode
pip install -e .[test]
