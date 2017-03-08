#!/bin/bash

set -ev

# Build MPICH
if [[ ! -e $HOME/mpich_install/bin/mpiexec ]]; then
    wget -q http://www.mpich.org/static/downloads/3.1.3/mpich-3.1.3.tar.gz
    tar -xzvf mpich-3.1.3.tar.gz >/dev/null 2>&1
    cd mpich-3.1.3
    ./configure --prefix=$HOME/mpich_install -q
    make -j 2 >/dev/null 2>&1
    make install >/dev/null 2>&1
    cd ..
fi

# Build PHDF5
if [[ ! -e $HOME/phdf5_install/bin/h5pfc ]]; then
    wget -q http://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.15/src/hdf5-1.8.15.tar.gz
    tar -xzvf hdf5-1.8.15.tar.gz >/dev/null 2>&1
    mv hdf5-1.8.15 phdf5-1.8.15; cd phdf5-1.8.15
    CC=$HOME/mpich_install/bin/mpicc FC=$HOME/mpich_install/bin/mpif90 \
       ./configure \
       --prefix=$HOME/phdf5_install -q --enable-fortran \
       --enable-fortran2003 --enable-parallel
    make -j 2 >/dev/null 2>&1
    make install >/dev/null 2>&1
    cd ..
fi

# Build HDF5
if [[ ! -e $HOME/hdf5_install/bin/h5fc ]]; then
    tar -xzvf hdf5-1.8.15.tar.gz >/dev/null 2>&1
    cd hdf5-1.8.15
    CC=gcc FC=gfortran ./configure --prefix=$HOME/hdf5_install -q \
                                   --enable-fortran --enable-fortran2003
    make -j 2 >/dev/null 2>&1
    make install >/dev/null 2>&1
    cd ..
fi
