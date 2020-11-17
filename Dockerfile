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
    make ; \
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
    apt-get --yes install libtbb-dev ; \
    apt-get --yes install libglfw3-dev ; \
    fi

# perhaps not needed#
# ENV LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial:$LD_LIBRARY_PATH

RUN git clone https://github.com/embree/embree
RUN git clone https://github.com/pshriwise/double-down
RUN git clone -b develop https://github.com/svalinn/dagmc
RUN git clone --recurse-submodules https://github.com/openmc-dev/openmc.git
RUN cd $HOME && \
    mkdir MOAB && \
    cd MOAB && \
    git clone  --single-branch --branch develop https://bitbucket.org/fathomteam/moab/


RUN apt-get --yes install  


# Clone and install Embree
RUN if [ "$include_dagmc" = "true" ] ; \
    then echo installing embree ; \
    # git clone https://github.com/embree/embree ; \
    cd embree ; \
    mkdir build ; \
    cd build ; \
    cmake .. -DCMAKE_INSTALL_PREFIX=.. \
        -DEMBREE_ISPC_SUPPORT=OFF ; \
    make ; \
    make install ; \
    fi


# Clone and install MOAB
RUN if [ "$include_dagmc" = "true" ] ; \
    then mkdir MOAB ; \
    mkdir build ; \
    cd MOAB ; \
    # git clone --depth 1 https://bitbucket.org/fathomteam/moab -b develop ; \
    cd build ; \
    cmake ../moab -DENABLE_HDF5=ON ; \
                -DENABLE_NETCDF=ON ; \
                -DBUILD_SHARED_LIBS=OFF ; \
                -DENABLE_FORTRAN=OFF ; \
                -DCMAKE_INSTALL_PREFIX=/MOAB ; \
    make -j4; \
    make install ; \
    rm -rf * ; \
    cmake ../moab -DBUILD_SHARED_LIBS=ON ; \
                -DENABLE_HDF5=ON ; \
                -DENABLE_PYMOAB=ON ; \
                -DENABLE_BLASLAPACK=OFF ; \
                -DENABLE_FORTRAN=OFF ; \
                -DCMAKE_INSTALL_PREFIX=/MOAB ; \
    make -j4; \
    make install ; \
    fi


# Clone and install double-down
RUN if [ "$include_dagmc" = "true" ] ; \
    then echo installing double-down ; \
    # git clone https://github.com/pshriwise/double-down ; \
    cd double-down ; \
    mkdir build ; \
    cd build ; \
    cmake .. -DCMAKE_INSTALL_PREFIX=.. ; \
        -DMOAB_DIR=/MOAB ; \
        -DEMBREE_DIR=/embree ; \
        -DEMBREE_ROOT=/embree ; \
    make ; \
    make install ; \
    fi

# Clone and install DAGMC
RUN if [ "$include_dagmc" = "true" ] ; \
    then mkdir DAGMC ; \
    cd DAGMC ; \
    # git clone -b develop https://github.com/svalinn/dagmc ; \
    mkdir build ; \
    cd build ; \
    cmake ../dagmc -DBUILD_TALLY=ON ; \
        -DCMAKE_INSTALL_PREFIX=/DAGMC/ ; \
        -DMOAB_DIR=/MOAB  ; \
        -DBUILD_STATIC_LIBS=OFF ; \
        -DBUILD_STATIC_EXE=OFF ; \
    make -j install ; \
    rm -rf /DAGMC/dagmc /DAGMC/build ; \
    fi


# Clone and install OpenMC with DAGMC
RUN if [ "$include_dagmc" = "true" ] ; \
    then echo installing openmc with dagmc ; \
    # then git clone --recurse-submodules https://github.com/openmc-dev/openmc.git ; \
    /opt/openmc ; \
    cd /opt/openmc ; \
    mkdir -p build ; \
    cd build ; \
    cmake -Doptimize=on -Ddagmc=ON ; \
        -DDAGMC_DIR=/DAGMC/ ; \
        -DHDF5_PREFER_PARALLEL=on ..  ; \
    make  ; \
    make install ; \
    cd ..  ; \
    pip install -e .[test] ; \
    fi

# Download cross sections (NNDC and WMP) and ENDF data needed by test suite
RUN /opt/openmc/tools/ci/download-xs.sh
