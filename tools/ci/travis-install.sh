#!/bin/bash
set -ex

# Install NJOY 2016
./tools/ci/travis-install-njoy.sh

# Running OpenMC's setup.py requires numpy/cython already
pip install numpy cython

# pytest installed by default -- make sure we get latest
pip install --upgrade pytest

# IPython stopped supporting Python 2.7 with version 6
if [[ "$TRAVIS_PYTHON_VERSION" == "2.7" ]]; then
    pip install "ipython<6"
fi

# Pandas stopped supporting Python 3.4 with version 0.21
if [[ "$TRAVIS_PYTHON_VERSION" == "3.4" ]]; then
    pip install pandas==0.20.3
fi

# Build and install OpenMC executable
python tools/ci/travis-install.py

# Install Python API in editable mode
pip install -e .[test]

# For uploading to coveralls
pip install python-coveralls
