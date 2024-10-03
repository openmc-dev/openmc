# ----------------------------------------------------------------------------
# Dockerfile for building the OpenMC project with support for various dependencies
# and configurations. This Dockerfile allows you to build OpenMC with support
# for GCC or OpenMPI compilers, along with several libraries like NJOY2016, HDF5,
# NetCDF, MOAB, EMBREE, Double Down, DAGMC, NCrystal, PyBind11, Xtensor, Vectfit,
# libMesh, and MCPL. Each of these dependencies is installed from their respective
# repositories and tags.

# The build process is split into stages:
# 1. Base Stage: Sets up a Manylinux base image and installs necessary dependencies.
# 2. Compiler Configuration: Defines the compilers (GCC or OpenMPI) to be used.
# 3. Dependencies Stage: Downloads and builds all external dependencies.
# 4. Python Dependencies Stage: Downloads and installs Python dependencies.
# 5. OpenMC Stage: Copies OpenMC source code, build wheel with specific
#    flags and installs it in the container. Then, runs tests.

# Arguments and environment variables can be customized for different compiler and
# dependency versions.
#
# To build the Docker image, use the following command from the repository's root:
# docker build -t openmc -f tools/ci/manylinux.dockerfile .
#
# For more information about each step, refer to the inline comments.
# ----------------------------------------------------------------------------

# Configure base image
ARG MANYLINUX_IMAGE=manylinux_2_28_x86_64

# Configure Compiler to use (gcc or openmpi)
ARG COMPILER="gcc"

# Configure Python ABI to use
ARG Python_ABI="cp312-cp312"

# OpenMC options
ARG OPENMC_USE_OPENMP="ON"
ARG OPENMC_BUILD_TESTS="ON"
ARG OPENMC_ENABLE_PROFILE="OFF"
ARG OPENMC_ENABLE_COVERAGE="OFF"
ARG OPENMC_USE_DAGMC="ON"
ARG OPENMC_USE_LIBMESH="ON"
ARG OPENMC_USE_MCPL="ON"
ARG OPENMC_USE_NCRYSTAL="ON"
ARG OPENMC_USE_UWUW="OFF"

# Configure dependencies tags
ARG NJOY2016_TAG="2016.76"
ARG HDF5_TAG="hdf5_1.14.4.3"
ARG NETCDF_TAG="v4.9.2"
ARG MOAB_TAG="5.5.1"
ARG EMBREE_TAG="v4.3.3"
ARG DD_TAG="v1.1.0"
ARG DAGMC_TAG="v3.2.3"
ARG NCrystal_TAG="v3.9.7"
ARG PYBIND_TAG="v2.13.6"
ARG XTL_TAG="0.7.7"
ARG XTENSOR_TAG="0.25.0"
ARG XTENSOR_PYTHON_TAG="0.27.0"
ARG XTENSOR_BLAS_TAG="0.21.0"
ARG VECTFIT_TAG="master"
ARG LIBMESH_TAG="v1.7.2"
ARG MCPL_TAG="v1.6.2"


# Base stage
FROM quay.io/pypa/${MANYLINUX_IMAGE} AS base

# Set timezone
ENV TZ=America/Chicago
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Set Home directory
ENV HOME /root
WORKDIR $HOME

# Setup Epel repository and install build dependencies
RUN yum install -y epel-release && \
    yum config-manager --enable epel && \
    yum install -y \ 
        wget \
        git \
        gcc \
        gcc-c++ \
        gcc-gfortran \
        make \
        python3.12-devel \
        zlib-devel \
        curl-devel \
        eigen3-devel \
        lapack-devel \
        libpng-devel && \
    yum clean all

# Set up environment variables for shared libraries
ENV LD_LIBRARY_PATH=/usr/local/lib64:$LD_LIBRARY_PATH


# Compiler configuration stage: gcc
FROM base AS compiler-gcc

ENV CC=gcc
ENV CXX=g++
ENV FC=gfortran
ENV F77=gfortran


# Compiler configuration stage: openmpi
FROM base AS compiler-openmpi

# Install OpenMPI
RUN yum install -y \
        openmpi-devel && \
    yum clean all

ENV CC=mpicc
ENV CXX=mpicxx
ENV FC=mpif90
ENV F77=mpif77

# Set up OpenMPI environment variables
ENV PATH=/usr/lib64/openmpi/bin:$PATH
ENV LD_LIBRARY_PATH=/usr/lib64/openmpi/lib:$LD_LIBRARY_PATH


# Dependencies stage
FROM compiler-${COMPILER} AS dependencies

ARG COMPILER

# Build and install NJOY2016
ARG NJOY2016_TAG
RUN git clone --depth 1 -b ${NJOY2016_TAG} https://github.com/njoy/njoy2016.git njoy && \
    cd njoy && \
    mkdir build && cd build && \
    cmake .. \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -Dstatic=ON && \
    make -j$(nproc) && make install && \
    cd ../.. && \
    rm -rf njoy

# Build and install HDF5
ARG HDF5_TAG
RUN git clone --depth 1 -b ${HDF5_TAG} https://github.com/HDFGroup/hdf5.git hdf5 && \
    cd hdf5 && \
    mkdir build && cd build && \
    cmake .. \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DHDF5_ENABLE_PARALLEL=$([ ${COMPILER} == "openmpi" ] && echo "ON" || echo "OFF") \
        -DHDF5_BUILD_HL_LIB=ON \
        -DBUILD_SHARED_LIBS=ON && \
    make -j$(nproc) && make install && \
    cd ../.. && \
    rm -rf hdf5

# Build and install NetCDF
ARG NETCDF_TAG
RUN git clone --depth 1 -b ${NETCDF_TAG} https://github.com/Unidata/netcdf-c.git netcdf && \
    cd netcdf && \
    mkdir build && cd build && \
    cmake .. \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DBUILD_SHARED_LIBS=ON \
        -DENABLE_DAP=ON \
        -DENABLE_TESTS=OFF && \
    make -j$(nproc) && make install && \
    cd ../.. && \
    rm -rf netcdf

# Build and install MOAB
ARG MOAB_TAG
RUN git clone --depth 1 -b ${MOAB_TAG} https://bitbucket.org/fathomteam/moab.git moab && \
    cd moab && \
    mkdir build && cd build && \
    cmake .. \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DENABLE_MPI=$([ ${COMPILER} == "openmpi" ] && echo "ON" || echo "OFF") \
        -DENABLE_HDF5=ON \
        -DHDF5_ROOT=/usr/local \
        -DENABLE_NETCDF=ON \
        -DNETCDF_ROOT=/usr/local \
        -DBUILD_SHARED_LIBS=ON \
        -DENABLE_BLASLAPACK=OFF \
        -DENABLE_PYMOAB=OFF && \
    make -j$(nproc) && make install && \
    cd ../.. && \
    rm -rf moab

# Build and install Embree
ARG EMBREE_TAG
RUN git clone --depth 1 -b ${EMBREE_TAG} https://github.com/embree/embree.git embree && \
    cd embree && \
    mkdir build && cd build && \
    cmake .. \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DEMBREE_TASKING_SYSTEM=INTERNAL \
        -DEMBREE_ISPC_SUPPORT=OFF \
        -DEMBREE_TUTORIALS=OFF && \
    make -j$(nproc) && make install && \
    cd ../.. && \
    rm -rf embree

# Build and install Double Down
ARG DD_TAG
RUN git clone --depth 1 -b ${DD_TAG} https://github.com/pshriwise/double-down.git dd && \
    cd dd && \
    mkdir build && cd build && \
    cmake .. \
        -DCMAKE_INSTALL_PREFIX=/usr/local && \
    make -j$(nproc) && make install && \
    cd ../.. && \
    rm -rf dd

# Build and install DAGMC
ARG DAGMC_TAG
RUN git clone --depth 1 -b ${DAGMC_TAG} https://github.com/svalinn/DAGMC.git dagmc && \
    cd dagmc && \
    mkdir build && cd build && \
    cmake .. \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DMOAB_DIR=/usr/local \
        -Ddd_ROOT=/usr/local \
        -DBUILD_TALLY=ON \
        -DBUILD_UWUW=ON \
        -DDOUBLE_DOWN=ON \
        -DBUILD_STATIC_LIBS=OFF \
        -DBUILD_RPATH=OFF && \
    make -j$(nproc) && make install && \
    cd ../.. && \
    rm -rf dagmc

# Build and install libMesh
ARG LIBMESH_TAG
RUN git clone --depth 1 -b ${LIBMESH_TAG} https://github.com/libMesh/libmesh.git libmesh && \
    cd libmesh && \
    git submodule update --init --recursive && \
    mkdir build && cd build && \
    export METHODS="opt" && \
    ../configure \
        $([ ${COMPILER} = 'openmpi' ] && echo '--enable-mpi' || echo '--disable-mpi') \
        --prefix=/usr/local \
        --enable-exodus \
        --disable-netcdf-4 \
        --disable-eigen \
        --disable-lapack && \
    make -j$(nproc) && make install && \
    cd ../.. && \
    rm -rf libmesh

# Build and install MCPL
ARG MCPL_TAG
RUN git clone --depth 1 -b ${MCPL_TAG} https://github.com/mctools/mcpl.git mcpl && \
    cd mcpl && \
    mkdir build && cd build && \
    cmake .. \
        -DCMAKE_INSTALL_PREFIX=/usr/local && \
    make -j$(nproc) && make install && \
    cd ../.. && \
    rm -rf mcpl

# Download and extract HDF5 data
RUN wget -q -O - https://anl.box.com/shared/static/teaup95cqv8s9nn56hfn7ku8mmelr95p.xz | tar -C $HOME -xJ
ENV OPENMC_CROSS_SECTIONS=$HOME/nndc_hdf5/cross_sections.xml

# Download and extract ENDF/B-VII.1 distribution
RUN wget -q -O - https://anl.box.com/shared/static/4kd2gxnf4gtk4w1c8eua5fsua22kvgjb.xz | tar -C $HOME -xJ
ENV OPENMC_ENDF_DATA=$HOME/endf-b-vii.1


# Python dependencies stage
FROM dependencies AS python-dependencies

ARG Python_ABI

# Use Python from manylinux as the default Python
ENV PATH="/opt/python/${Python_ABI}/bin:${PATH}"
RUN ln -sf /opt/python/${Python_ABI}/bin/python3 /usr/bin/python

# Build and install NCrystal
ARG NCrystal_TAG
RUN git clone --depth 1 -b ${NCrystal_TAG} https://github.com/mctools/ncrystal.git ncrystal && \
    cd ncrystal && \
    mkdir build && cd build && \
    cmake .. \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DBUILD_SHARED_LIBS=ON \
        -DNCRYSTAL_NOTOUCH_CMAKE_BUILD_TYPE=ON \
        -DNCRYSTAL_MODIFY_RPATH=OFF \
        -DCMAKE_BUILD_TYPE=Release \
        -DNCRYSTAL_ENABLE_EXAMPLES=OFF \
        -DNCRYSTAL_ENABLE_SETUPSH=OFF \
        -DNCRYSTAL_ENABLE_DATA=EMBED \
        -DPython3_EXECUTABLE=$(which python) && \
    make -j$(nproc) && make install && \
    cd ../.. && \
    ncrystal-config --setup && \
    rm -rf ncrystal

# Build and install pybind
ARG PYBIND_TAG
RUN git clone --depth 1 -b ${PYBIND_TAG} https://github.com/pybind/pybind11.git pybind11 && \
    cd pybind11 && \
    mkdir build && cd build && \
    cmake .. && \
    make -j$(nproc) && make install && \
    cd .. && \
    python -m pip install . && \
    cd .. && \
    rm -rf pybind11 

# Build and install xtl
ARG XTL_TAG
RUN git clone --depth 1 -b ${XTL_TAG} https://github.com/xtensor-stack/xtl.git xtl && \
    cd xtl && \
    mkdir build && cd build && \
    cmake .. && \
    make -j$(nproc) && make install && \
    cd ../.. && \
    rm -rf xtl

# Build and install xtensor
ARG XTENSOR_TAG
RUN git clone --depth 1 -b ${XTENSOR_TAG} https://github.com/xtensor-stack/xtensor.git xtensor && \
    cd xtensor && \
    mkdir build && cd build && \
    cmake .. && \
    make -j$(nproc) && make install && \
    cd ../.. && \
    rm -rf xtensor

# Build and install xtensor-python
ARG XTENSOR_PYTHON_TAG
RUN git clone --depth 1 -b ${XTENSOR_PYTHON_TAG} https://github.com/xtensor-stack/xtensor-python.git xtensor-python && \
    cd xtensor-python && \
    mkdir build && cd build && \
    python -m pip install numpy && \
    cmake .. \
    -DNUMPY_INCLUDE_DIRS=$(python -c "import numpy; print(numpy.get_include())") && \
    make -j$(nproc) && make install && \
    cd ../.. && \
    rm -rf xtensor-python

# Build and install xtensor-blas
ARG XTENSOR_BLAS_TAG
RUN git clone --depth 1 -b ${XTENSOR_BLAS_TAG} https://github.com/xtensor-stack/xtensor-blas.git xtensor-blas && \
    cd xtensor-blas && \
    mkdir build && cd build && \
    cmake .. && \
    make -j$(nproc) && make install && \
    cd ../.. && \
    rm -rf xtensor-blas

# Build and install vectfit
ARG VECTFIT_TAG
RUN git clone --depth 1 -b ${VECTFIT_TAG} https://github.com/liangjg/vectfit.git vectfit && \
    cd vectfit && \
    python -m pip install . && \
    cd .. && \
    rm -rf vectfit


# OpenMC stage
FROM python-dependencies AS openmc

ARG COMPILER
ARG Python_ABI
ARG OPENMC_USE_OPENMP
ARG OPENMC_BUILD_TESTS
ARG OPENMC_ENABLE_PROFILE
ARG OPENMC_ENABLE_COVERAGE
ARG OPENMC_USE_DAGMC
ARG OPENMC_USE_LIBMESH
ARG OPENMC_USE_MCPL
ARG OPENMC_USE_NCRYSTAL
ARG OPENMC_USE_UWUW

# Copy OpenMC source to docker image
COPY . $HOME/openmc

# Configure SKBUILD CMake arguments
RUN export SKBUILD_CMAKE_ARGS="-DOPENMC_USE_MPI=$([ ${COMPILER} == 'openmpi' ] && echo 'ON' || echo 'OFF'); \
                        -DOPENMC_USE_OPENMP=${OPENMC_USE_OPENMP}; \
                        -DOPENMC_BUILD_TESTS=${OPENMC_BUILD_TESTS}; \
                        -DOPENMC_ENABLE_PROFILE=${OPENMC_ENABLE_PROFILE}; \
                        -DOPENMC_ENABLE_COVERAGE=${OPENMC_ENABLE_COVERAGE}; \
                        -DOPENMC_USE_DAGMC=${OPENMC_USE_DAGMC}; \
                        -DOPENMC_USE_LIBMESH=${OPENMC_USE_LIBMESH}; \
                        -DOPENMC_USE_MCPL=${OPENMC_USE_MCPL}; \
                        -DOPENMC_USE_NCRYSTAL=${OPENMC_USE_NCRYSTAL}; \
                        -DOPENMC_USE_UWUW=${OPENMC_USE_UWUW}" && \
    cd $HOME/openmc && \
    python -m build . -w

# Repair wheel
RUN auditwheel repair $HOME/openmc/dist/openmc-*.whl -w $HOME/openmc/dist/

# Install OpenMC wheel
RUN python -m pip install \
        "$(echo $HOME/openmc/dist/*manylinux**.whl)[$([ ${COMPILER} == 'openmpi' ] && echo 'depletion-mpi,')test,ci,vtk]"

# Test OpenMC
RUN cd $HOME/openmc && \
    nctool --test && \
    pytest --cov=openmc -v $([ ${COMPILER} == 'openmpi' ] && echo '--mpi') --event tests
