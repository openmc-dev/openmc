#!/bin/bash
set -ex

XDG_BRANCH='main'
XDG_REPO='https://github.com/xdg-org/xdg.git'
XDG_INSTALL_DIR=$HOME/XDG/

pushd $HOME

mkdir -p XDG
cd XDG
git clone -b $XDG_BRANCH $XDG_REPO xdg
mkdir build
cd build

cmake_args=(
    ../xdg
    -DCMAKE_INSTALL_PREFIX=$XDG_INSTALL_DIR
    -DXDG_BUILD_TESTS=OFF
    -DXDG_BUILD_TOOLS=OFF
    -DXDG_ENABLE_MOAB=ON
    -DMOAB_DIR=$HOME/MOAB
)

if [[ $LIBMESH == 'y' || $XDG == 'y' ]]; then
    cmake_args+=(
        -DXDG_ENABLE_LIBMESH=ON
        -DCMAKE_PREFIX_PATH=$HOME/LIBMESH
    )
else
    cmake_args+=(-DXDG_ENABLE_LIBMESH=OFF)
fi

if [[ $MPI == 'y' ]]; then
    cmake_args+=(-DXDG_LINK_MPI=ON)
fi

cmake "${cmake_args[@]}"
make -j4
make -j4 install
rm -rf $HOME/XDG/xdg $HOME/XDG/build

popd
