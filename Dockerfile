# To build with OpenMC
# docker build -t openmc .

# To build with OpenMC and DAGMC enabled
# docker build -t openmc_dagmc --build-arg include_dagmc=true .

FROM ubuntu:latest

# By default this Dockerfile builds OpenMC without dagmc
ARG include_dagmc=false

# Setup environment variables for Docker image
ENV CC=/usr/bin/mpicc CXX=/usr/bin/mpicxx \
    PATH=/opt/openmc/bin:/opt/NJOY2016/build:$PATH \
    LD_LIBRARY_PATH=/opt/openmc/lib:$LD_LIBRARY_PATH \
    OPENMC_CROSS_SECTIONS=/root/nndc_hdf5/cross_sections.xml \
    OPENMC_ENDF_DATA=/root/endf-b-vii.1 \
    DEBIAN_FRONTEND=noninteractive

# Install dependencies from Debian package manager
RUN apt-get update -y && \
    apt-get upgrade -y && \
    apt-get install -y \
        python3-pip python-is-python3 wget git gfortran g++ cmake \
        mpich libmpich-dev libhdf5-serial-dev libhdf5-mpich-dev \
        imagemagick && \
    apt-get autoremove

# Update system-provided pip
RUN pip3 install --upgrade pip

# Clone and install NJOY2016
RUN git clone https://github.com/njoy/NJOY2016 /opt/NJOY2016 && \
    cd /opt/NJOY2016 && \
    mkdir build && cd build && \
    cmake -Dstatic=on .. && make 2>/dev/null && make install

# Clone and install OpenMC without DAGMC
RUN if [ "$include_dagmc" = "false" ] ; \
    then git clone --recurse-submodules https://github.com/openmc-dev/openmc.git \
    /opt/openmc ; \
    cd /opt/openmc ; \
    mkdir -p build ; \
    cd build ; \
    cmake -Doptimize=on -DHDF5_PREFER_PARALLEL=on .. ; \
    make  ; \
    make install ; \
    cd ..  ; \
    pip install -e .[test] ; \
    fi

# install addition packages required for DAGMC
RUN if [ "$include_dagmc" = "true" ] ; \
    then apt-get --yes install libeigen3-dev ; \
    apt-get --yes install libblas-dev ; \
    apt-get --yes install liblapack-dev ; \
    apt-get --yes install libnetcdf-dev ; \
    #apt-get --yes install libnetcdf13 ; \
    fi

ENV LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial:$LD_LIBRARY_PATH

# Clone and install MOAB
RUN if [ "$include_dagmc" = "true" ] ; \
    then mkdir -p MOAB/bld ; \
    cd MOAB ; \
    git clone --depth 1 https://bitbucket.org/fathomteam/moab -b Version5.1.0 ; \
    cd bld ; \
    cmake ../moab -DENABLE_HDF5=ON ; \
                -DENABLE_NETCDF=ON ; \
                -DBUILD_SHARED_LIBS=OFF ; \
                -DENABLE_FORTRAN=OFF ; \
                -DCMAKE_INSTALL_PREFIX=$HOME/MOAB ; \
    make -j4; \
    make install ; \
    rm -rf * ; \
    cmake ../moab -DBUILD_SHARED_LIBS=ON ; \
                -DENABLE_HDF5=ON ; \
                -DENABLE_PYMOAB=ON ; \
                -DENABLE_BLASLAPACK=OFF ; \
                -DENABLE_FORTRAN=OFF ; \
                -DCMAKE_INSTALL_PREFIX=$HOME/MOAB ; \
    make -j4; \
    make install ; \
    fi

# Clone and install DAGMC
RUN if [ "$include_dagmc" = "true" ] ; \
    then mkdir DAGMC ; \
    cd DAGMC ; \
    git clone -b develop https://github.com/svalinn/dagmc ; \
    mkdir build ; \
    cd build ; \
    cmake ../dagmc -DBUILD_TALLY=ON ; \
                -DCMAKE_INSTALL_PREFIX=$HOME/DAGMC/ ; \
                -DMOAB_DIR=$HOME/MOAB  ; \
                -DBUILD_STATIC_LIBS=OFF ; \
                -DBUILD_STATIC_EXE=OFF ; \
    make -j install ; \
    rm -rf $HOME/DAGMC/dagmc $HOME/DAGMC/build ; \
    fi

# Clone and install OpenMC with DAGMC
RUN if [ "$include_dagmc" = "true" ] ; \
    then git clone --recurse-submodules https://github.com/openmc-dev/openmc.git \
    /opt/openmc ; \
    cd /opt/openmc ; \
    mkdir -p build ; \
    cd build ; \
    cmake -Doptimize=on -Ddagmc=ON ; \
                    -DDAGMC_DIR=$HOME/DAGMC/ ; \
                    -DHDF5_PREFER_PARALLEL=on ..  ; \
    make  ; \
    make install ; \
    cd ..  ; \
    pip install -e .[test] ; \
    fi

# Download cross sections (NNDC and WMP) and ENDF data needed by test suite
RUN /opt/openmc/tools/ci/download-xs.sh
