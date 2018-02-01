#!/bin/bash
set -ex

# Allow tests that require GUI as described at:
# https://docs.travis-ci.com/user/gui-and-headless-browsers/#Using-xvfb-to-Run-Tests-That-Require-a-GUI
sh -e /etc/init.d/xvfb start

# Download NNDC HDF5 data
if [[ ! -e $HOME/nndc_hdf5/cross_sections.xml ]]; then
    wget https://anl.box.com/shared/static/a6sw2cep34wlz6b9i9jwiotaqoayxcxt.xz -O - | tar -C $HOME -xvJ
fi

# Download ENDF/B-VII.1 distribution
if [[ ! -d $HOME/endf-b-vii.1/neutrons ]]; then
    wget https://anl.box.com/shared/static/4kd2gxnf4gtk4w1c8eua5fsua22kvgjb.xz -O - | tar -C $HOME -xvJ
fi

# Download multipole library
git clone --branch=master git://github.com/smharper/windowed_multipole_library.git wmp_lib
tar -C $HOME -xzvf wmp_lib/multipole_lib.tar.gz
