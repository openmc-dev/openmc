# Configure base image
ARG MANYLINUX_IMAGE=manylinux_2_28_x86_64
ARG Python_ABI="cp312-cp312"

# Configure Compiler to use (gcc or openmpi)
ARG COMPILER="gcc"

# Configure dependencies tags
ARG NJOY2016_TAG="2016.76"
ARG HDF5_TAG="hdf5_1.14.4.3"
ARG NETCDF_TAG="v4.9.2"
ARG MOAB_TAG="master"
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


# Build base stage
FROM quay.io/pypa/${MANYLINUX_IMAGE} AS base

ARG Python_ABI

# Set timezone
ENV TZ=America/Chicago
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Set Home directory
ENV HOME /root
WORKDIR $HOME

# Setup Epel
RUN yum install -y epel-release && \
    yum config-manager --enable epel

# Install basic dependencies
RUN yum install -y \ 
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
        libpng-devel \
        openmpi-devel && \
    yum clean all

# Use Python from manylinux as the default Python
ENV PATH="/opt/python/${Python_ABI}/bin:${PATH}"
RUN ln -sf /opt/python/${Python_ABI}/bin/python3 /usr/bin/python

# Set up general library environment variables
ENV LD_LIBRARY_PATH=/usr/lib64:$LD_LIBRARY_PATH

# Install necessary Python packages
RUN python -m pip install --upgrade \
        scikit-build-core \
        setuptools \  
        numpy \
        cmake \
        ninja \
        h5py \
        scipy \
        ipython \
        matplotlib \
        pandas \
        lxml \
        uncertainties \
        endf \
        vtk \
        mpi4py \
        packaging \
        pytest \
        pytest-cov \
        colorama \
        openpyxl

FROM base AS compiler-gcc

ENV CC=gcc
ENV CXX=g++
ENV FC=gfortran
ENV F77=gfortran

FROM base AS compiler-openmpi

ENV CC=mpicc
ENV CXX=mpicxx
ENV FC=mpif90
ENV F77=mpif77

# Set up OpenMPI environment variables
ENV PATH=/usr/lib64/openmpi/bin:$PATH
ENV LD_LIBRARY_PATH=/usr/lib64/openmpi/lib:$LD_LIBRARY_PATH


FROM compiler-${COMPILER} AS dependencies

ARG COMPILER

# Set up NJOY2016
ARG NJOY2016_TAG

# Build and install NJOY2016
RUN git clone --depth 1 -b ${NJOY2016_TAG} https://github.com/njoy/njoy2016.git njoy && \
    cd njoy && \
    mkdir build && cd build && \
    cmake .. \
        -Dstatic=ON && \
    make -j$(nproc) && make install && \
    cd ../.. && \
    rm -rf njoy


# Set up HDF5
ARG HDF5_TAG

# Build and install HDF5
RUN git clone --depth 1 -b ${HDF5_TAG} https://github.com/HDFGroup/hdf5.git hdf5 && \
    cd hdf5 && \
    mkdir build && cd build && \
    cmake .. \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DHDF5_ENABLE_PARALLEL=$([ ${COMPILER} == "openmpi" ] && echo "ON" || echo "OFF") \
        -DHDF5_BUILD_HL_LIB=ON \
        -DBUILD_SHARED_LIBS=ON \
        -DBUILD_EXAMPLES=OFF && \
    make -j$(nproc) && make install && \
    cd ../.. && \
    rm -rf hdf5


# Set up NetCDF
ARG NETCDF_TAG

# Build and install NetCDF
RUN git clone --depth 1 -b ${NETCDF_TAG} https://github.com/Unidata/netcdf-c.git netcdf && \
    cd netcdf && \
    mkdir build && cd build && \
    cmake .. \
        -DBUILD_SHARED_LIBS=ON \
        -DENABLE_DAP=ON \
        -DENABLE_TESTS=OFF && \
    make -j$(nproc) && make install && \
    cd ../.. && \
    rm -rf netcdf


# Set up MOAB
ARG MOAB_TAG

# Build and install MOAB
RUN git clone --depth 1 -b ${MOAB_TAG} https://bitbucket.org/fathomteam/moab.git moab && \
    cd moab && \
    mkdir build && cd build && \
    cmake .. \
        -DENABLE_MPI=$([ ${COMPILER} == "openmpi" ] && echo "ON" || echo "OFF") \
        -DENABLE_HDF5=ON \
        -DHDF5_ROOT=/usr/local \
        -DENABLE_NETCDF=ON \
        -DNETCDF_ROOT=/usr/local \
        -DBUILD_SHARED_LIBS=ON \
        -DENABLE_BLASLAPACK=OFF && \
    make -j$(nproc) && make install && \
    cd ../.. && \
    rm -rf moab

# Set up EMBREE
ARG EMBREE_TAG

# Build and install EMBREE
RUN git clone --depth 1 -b ${EMBREE_TAG} https://github.com/embree/embree.git embree && \
    cd embree && \
    mkdir build && cd build && \
    cmake .. \
        -DEMBREE_TASKING_SYSTEM=INTERNAL \
        -DEMBREE_ISPC_SUPPORT=OFF \
        -DEMBREE_TUTORIALS=OFF && \
    make -j$(nproc) && make install && \
    cd ../.. && \
    rm -rf embree

# Set up Double Down
ARG DD_TAG

# Build and install Double Down
RUN git clone --depth 1 -b ${DD_TAG} https://github.com/pshriwise/double-down.git dd && \
    cd dd && \
    mkdir build && cd build && \
    cmake .. && \
    make -j$(nproc) && make install && \
    cd ../.. && \
    rm -rf dd

# Set up DAGMC
ARG DAGMC_TAG

# Build and install DAGMC
RUN git clone --depth 1 -b ${DAGMC_TAG} https://github.com/svalinn/DAGMC.git dagmc && \
    cd dagmc && \
    mkdir build && cd build && \
    cmake .. \
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

# Set up NCrystal
ARG NCrystal_TAG

# Build and install NCrystal
RUN git clone --depth 1 -b ${NCrystal_TAG} https://github.com/mctools/ncrystal.git ncrystal && \
    cd ncrystal && \
    mkdir build && cd build && \
    cmake .. \
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

# Set up pybind
ARG PYBIND_TAG

# Build and install pybind
RUN git clone --depth 1 -b ${PYBIND_TAG} https://github.com/pybind/pybind11.git pybind11 && \
    cd pybind11 && \
    mkdir build && cd build && \
    cmake .. && \
    make -j$(nproc) && make install && \
    cd .. && \
    python -m pip install . && \
    cd .. && \
    rm -rf pybind11 

# Set up xtl
ARG XTL_TAG

# Build and install xtl
RUN git clone --depth 1 -b ${XTL_TAG} https://github.com/xtensor-stack/xtl.git xtl && \
    cd xtl && \
    mkdir build && cd build && \
    cmake .. && \
    make -j$(nproc) && make install && \
    cd ../.. && \
    rm -rf xtl

# Set up xtensor
ARG XTENSOR_TAG

# Build and install xtensor
RUN git clone --depth 1 -b ${XTENSOR_TAG} https://github.com/xtensor-stack/xtensor.git xtensor && \
    cd xtensor && \
    mkdir build && cd build && \
    cmake .. && \
    make -j$(nproc) && make install && \
    cd ../.. && \
    rm -rf xtensor

# Set up xtensor-python
ARG XTENSOR_PYTHON_TAG

# Build and install xtensor-python
RUN git clone --depth 1 -b ${XTENSOR_PYTHON_TAG} https://github.com/xtensor-stack/xtensor-python.git xtensor-python && \
    cd xtensor-python && \
    mkdir build && cd build && \
    cmake .. \
    -DNUMPY_INCLUDE_DIRS=$(python -c "import numpy; print(numpy.get_include())") && \
    make -j$(nproc) && make install && \
    cd ../.. && \
    rm -rf xtensor-python

# Set up xtensor-blas
ARG XTENSOR_BLAS_TAG

# Build and install xtensor-blas
RUN git clone --depth 1 -b ${XTENSOR_BLAS_TAG} https://github.com/xtensor-stack/xtensor-blas.git xtensor-blas && \
    cd xtensor-blas && \
    mkdir build && cd build && \
    cmake .. && \
    make -j$(nproc) && make install && \
    cd ../.. && \
    rm -rf xtensor-blas

# Set up vectfit
ARG VECTFIT_TAG

# Build and install vectfit
RUN git clone --depth 1 -b ${VECTFIT_TAG} https://github.com/liangjg/vectfit.git vectfit && \
    cd vectfit && \
    python -m pip install . && \
    cd .. && \
    rm -rf vectfit

# Set up libMesh
ARG LIBMESH_TAG

# Build and install libMesh
RUN git clone --depth 1 -b ${LIBMESH_TAG} https://github.com/libMesh/libmesh.git libmesh && \
    cd libmesh && \
    git submodule update --init --recursive && \
    mkdir build && cd build && \
    export METHODS="opt" && \
    ../configure \
        --enable-exodus \
        --disable-netcdf-4 \
        --disable-eigen \
        --disable-lapack && \
    make -j$(nproc) && make install && \
    cd .. && \
    rm -rf libmesh

# Set up MCPL
ARG MCPL_TAG

# Build and install MCPL
RUN git clone --depth 1 --single-branch -b ${MCPL_TAG} https://github.com/mctools/mcpl.git mcpl && \
    cd mcpl && \
    mkdir build && cd build && \
    cmake .. && \
    make -j$(nproc) && make install && \
    cd .. && \
    rm -rf mcpl


# Download and extract HDF5 data
RUN wget -q -O - https://anl.box.com/shared/static/teaup95cqv8s9nn56hfn7ku8mmelr95p.xz | tar -C $HOME -xJ

# Download and extract ENDF/B-VII.1 distribution
RUN wget -q -O - https://anl.box.com/shared/static/4kd2gxnf4gtk4w1c8eua5fsua22kvgjb.xz | tar -C $HOME -xJ


FROM dependencies AS openmc

ARG COMPILER

# Copy OpenMC source
COPY . $HOME/openmc

# Configure SKBUILD CMake arguments
ENV SKBUILD_CMAKE_ARGS "-DOPENMC_USE_OPENMP=ON; \
                        -DOPENMC_USE_DAGMC=ON; \
                        -DOPENMC_USE_LIBMESH=ON; \
                        -DOPENMC_USE_MPI=$([ ${COMPILER} = 'openmpi' ] && echo 'ON' || echo 'OFF'); \
                        -DOPENMC_USE_MCPL=ON; \
                        -DOPENMC_USE_NCRYSTAL=ON;"

# Build OpenMC wheel
RUN cd $HOME/openmc && python -m build . -w

# Repair wheel
RUN auditwheel repair $HOME/openmc/dist/openmc-*.whl -w $HOME/openmc/dist

# Install OpenMC wheel
RUN python -m pip install $HOME/openmc/dist/*manylinux**.whl


FROM openmc AS openmc-test

ARG COMPILER

# Test OpenMC
RUN cd $HOME/openmc && \
    nctool --test && \
    pytest --cov=openmc -v $([ ${COMPILER} = 'openmpi' ] && echo '--mpi') --event tests
