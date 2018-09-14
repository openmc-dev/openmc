#!/bin/bash
set -ex

# Allow tests that require GUI as described at:
# https://docs.travis-ci.com/user/gui-and-headless-browsers/#Using-xvfb-to-Run-Tests-That-Require-a-GUI
sh -e /etc/init.d/xvfb start

# Download NNDC HDF5 data
if [[ ! -e $HOME/nndc_hdf5/cross_sections.xml ]]; then
    wget https://anl.box.com/shared/static/a0eflty17atnpd0pp7460exagr3nuhm7.xz -O - | tar -C $HOME -xvJ
fi

# Download ENDF/B-VII.1 distribution
ENDF=$HOME/endf-b-vii.1/
if [[ ! -d $ENDF/neutrons || ! -d $ENDF/photoat || ! -d $ENDF/atomic_relax ]]; then
    wget https://anl.box.com/shared/static/4kd2gxnf4gtk4w1c8eua5fsua22kvgjb.xz -O - | tar -C $HOME -xvJ
fi

# Download multipole library
if [[ ! -e $HOME/WMP_Library/092235.h5 ]]; then
    wget https://github.com/mit-crpg/WMP_Library/releases/download/v1.0/WMP_Library_v1.0.tar.gz
    tar -C $HOME -xzvf WMP_Library_v1.0.tar.gz
fi
