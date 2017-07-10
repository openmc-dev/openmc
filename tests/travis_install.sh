#!/bin/bash

set -ev

# Build PHDF5
if [[ ! -e $HOME/phdf5_install/bin/h5pfc ]]; then
    wget -q https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.19/src/hdf5-1.8.19.tar.bz2
    tar -xjf hdf5-1.8.19.tar.bz2
    mv hdf5-1.8.19 phdf5-1.8.19
    cd phdf5-1.8.19
    CC=/usr/bin/mpicc FC=/usr/bin/mpif90 \
      ./configure --prefix=$HOME/phdf5_install -q \
      --enable-fortran --enable-fortran2003 --enable-parallel
    make -j 2 &> /dev/null
    make install &> /dev/null
    cd ..
fi

# Build HDF5
if [[ ! -e $HOME/hdf5_install/bin/h5fc ]]; then
    tar -xjf hdf5-1.8.19.tar.bz2
    cd hdf5-1.8.19
    CC=gcc FC=gfortran ./configure --prefix=$HOME/hdf5_install -q \
      --enable-fortran --enable-fortran2003
    make -j 2 &> /dev/null
    make install &> /dev/null
    cd ..
fi
