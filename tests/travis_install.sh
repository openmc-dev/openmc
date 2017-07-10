#!/bin/bash

set -ev

# Build PHDF5
if [[ ! -e $HOME/phdf5_install/bin/h5pfc ]]; then
    wget -q https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.1/src/hdf5-1.10.1.tar.bz2
    tar -xjf hdf5-1.10.1.tar.bz2
    mv hdf5-1.10.1 phdf5-1.10.1
    cd phdf5-1.10.1
    CC=/usr/bin/mpicc FC=/usr/bin/mpif90 \
      ./configure --prefix=$HOME/phdf5_install -q \
      --enable-fortran --enable-parallel
    make -j 2 &> /dev/null
    make install &> /dev/null
    cd ..
fi

# Build HDF5
if [[ ! -e $HOME/hdf5_install/bin/h5fc ]]; then
    tar -xjf hdf5-1.10.1.tar.bz2
    cd hdf5-1.10.1
    CC=gcc FC=gfortran ./configure --prefix=$HOME/hdf5_install -q \
                                   --enable-fortran
    make -j 2 &> /dev/null
    make install &> /dev/null
    cd ..
fi
