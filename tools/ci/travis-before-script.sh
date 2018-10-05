#!/bin/bash
set -ex

# Allow tests that require GUI as described at:
# https://docs.travis-ci.com/user/gui-and-headless-browsers/#Using-xvfb-to-Run-Tests-That-Require-a-GUI
sh -e /etc/init.d/xvfb start

# Download NNDC HDF5 data, ENDF/B-VII.1 distribution, multipole library
sh ./tools/ci/download-xs.sh
