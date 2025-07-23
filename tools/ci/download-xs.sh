#!/bin/bash
# touching this script invalidate the test data xs cache
set -ex

# Download HDF5 data
if [[ ! -e $HOME/nndc_hdf5/cross_sections.xml ]]; then
    wget -q -O - https://github.com/GuySten/data/releases/download/v1.0.0/nndc_hdf5_test.tar.xz | tar -C $HOME -xJ
fi

# Download ENDF/B-VII.1 distribution
ENDF=$HOME/endf-b-vii.1
if [[ ! -d $ENDF/neutrons || ! -d $ENDF/photoat || ! -d $ENDF/atomic_relax || ! -d $ENDF/gammas ]]; then
    wget -q -O - https://anl.box.com/shared/static/4kd2gxnf4gtk4w1c8eua5fsua22kvgjb.xz | tar -C $HOME -xJ
fi
