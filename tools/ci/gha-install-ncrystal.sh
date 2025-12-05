#!/bin/bash
set -ex
cd $HOME

#Use the NCrystal develop branch (in the near future we can move this to master):
git clone https://github.com/mctools/ncrystal --branch develop --single-branch --depth 1 ncrystal_src

SRC_DIR="$PWD/ncrystal_src"
BLD_DIR="$PWD/ncrystal_bld"
INST_DIR="$PWD/ncrystal_inst"
PYTHON=$(which python3)

CPU_COUNT=1

mkdir "$BLD_DIR"
cd ncrystal_bld

cmake \
    "${SRC_DIR}" \
    -DBUILD_SHARED_LIBS=ON \
    -DNCRYSTAL_NOTOUCH_CMAKE_BUILD_TYPE=ON \
    -DNCRYSTAL_MODIFY_RPATH=OFF \
    -DCMAKE_BUILD_TYPE=Release \
    -DNCRYSTAL_ENABLE_EXAMPLES=OFF \
    -DNCRYSTAL_ENABLE_SETUPSH=OFF \
    -DNCRYSTAL_ENABLE_DATA=EMBED \
    -DCMAKE_INSTALL_PREFIX="${INST_DIR}" \
    -DPython3_EXECUTABLE="$PYTHON"

make -j${CPU_COUNT:-1}
make install

#Note: There is no "make test" or "make ctest" functionality for NCrystal
#      yet. If it appears in the future, we should add it here.

# Output the configuration to the log
"${INST_DIR}/bin/ncrystal-config" --setup

# Change environmental variables

eval $( "${INST_DIR}/bin/ncrystal-config" --setup )

# Check installation worked

nctool --test
