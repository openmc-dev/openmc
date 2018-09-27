FROM ubuntu:latest

# Setup environment variables for Docker image
ENV FC=/usr/bin/mpif90 CC=/usr/bin/mpicc CXX=/usr/bin/mpicxx \
    PATH=/opt/openmc/bin:/opt/conda/bin:/opt/NJOY2016/build:$PATH \
    LD_LIBRARY_PATH=/opt/openmc/lib:$LD_LIBRARY_PATH \
    OPENMC_CROSS_SECTIONS=/root/nndc_hdf5/cross_sections.xml \
    OPENMC_MULTIPOLE_LIBRARY=/root/WMP_Library \
    OPENMC_ENDF_DATA=/root/endf-b-vii.1

# Install dependencies from Debian package manager
RUN apt-get update -y && \
    apt-get upgrade -y && \
    apt-get install -y python3-pip && \
    apt-get install -y wget git emacs && \
    apt-get install -y gfortran g++ cmake && \
    apt-get install -y mpich libmpich-dev && \
    apt-get install -y libhdf5-serial-dev libhdf5-mpich-dev && \
    apt-get install -y imagemagick && \
    apt-get autoremove

# Download Miniconda3 and install Python dependencies
#RUN wget --quiet https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
#    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
#    rm -rf ~/miniconda.sh
#RUN pip install --upgrade pip

# Clone and install NJOY2016
RUN git clone https://github.com/njoy/NJOY2016 /opt/NJOY2016 && \
    cd /opt/NJOY2016 && \
    mkdir build && cd build && \
    cmake -Dstatic=on .. && make 2>/dev/null && make install

# Clone and install OpenMC
RUN git clone https://github.com/wbinventor/openmc.git /opt/openmc && \
    cd /opt/openmc && git checkout docker && mkdir -p build && cd build && \
    cmake  -Doptimize=on -DHDF5_PREFER_PARALLEL=on .. && \
    make && make install && \
    cd .. && pip3 install -e .[test]

# Download cross sections (NNDC and WMP) and ENDF data needed by test suite
#RUN ./opt/openmc/docker/docker-download-xs.sh
RUN ./opt/openmc/tools/ci/download-xs.sh